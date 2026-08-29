# Manually sourcing compound data for `tpt-thermo-data`

`tpt-thermo-data`'s seed dataset (`crates/tpt-thermo-data/data/seed.toml`) ships
2719 `[[components]]` entries, but only **69** carry real, cited critical
constants (`source = "NIST / Poling et al. 2001"` or similar). The other 2650
are `source = "estimated"` — and on inspection these aren't rough estimates of
real chemicals, they're a mechanically generated homologous-series filler
(`"tetratriacontylbenzene"`, `"butyl argenticyanide"`, `"C350-ane"`, ...) that
pads the count, not compounds worth curating. There's also a known data bug: a
bogus `"butane"` entry (Tc=290K, Pc=2.52MPa — physically wrong; real Tc≈425K)
sitting alongside the correctly-curated `"n-butane"` — worth deleting whenever
the seed set is next touched.

This guide is for whoever (you, a contributor, a future session) wants to grow
the curated set by hand. It can't be automated: the values themselves are
public science, but there's no free bulk source. NIST's site disallows large
automated scraping, DDBST's bulk data is a paid product, and an LLM
recalling these numbers from memory risks quietly wrong values — which is
worse than the current honest placeholders, because this crate is built
around every value being traceable to a real citation.

## What a compound needs

Every `[[components]]` entry is a `ComponentRecord`
(`crates/tpt-thermo-data/src/record.rs`):

| Field | Unit | Required? | Notes |
|---|---|---|---|
| `name` | — | yes | unique, lowercase, e.g. `"ethylene glycol"` |
| `formula` | — | no | e.g. `"C2H6O2"` |
| `cas` | — | no | CAS registry number |
| `critical_temperature_k` | K | yes | validated to (0, 2000] |
| `critical_pressure_pa` | Pa | yes | validated to (0, 1.0e9] |
| `acentric_factor` | — | yes | Pitzer ω, validated to [-0.6, 1.6] |
| `molar_mass_kg_per_mol` | kg/mol | yes | validated to [1e-4, 1.0] |
| `normal_boiling_point_k` | K | no | must be ≤ `critical_temperature_k` |
| `source` | — | **yes** | provenance string — never blank, never guessed |

`ComponentRecord::validate()` enforces the ranges above at load time, and the
loader rejects duplicate names — a malformed entry fails loudly rather than
silently corrupting the dataset.

## Where to find real values

1. **NIST Chemistry WebBook** (webbook.nist.gov) — the primary free source.
   Search by name, CAS number, or formula; the compound's page lists critical
   temperature/pressure directly when available (look under "Phase change
   data" / "Fluid properties"). Acentric factor usually isn't listed directly
   — derive it from the Pitzer definition using the vapor pressure at
   Tr = 0.7 if NIST gives a vapor-pressure correlation, or take it from Poling
   et al. (below) if already tabulated there.
   **Do not automate this** — NIST's site policy doesn't permit bulk automated
   scraping; this is a one-compound-at-a-time manual lookup.
2. **Poling, Prausnitz & O'Connell**, *The Properties of Gases and Liquids*
   (5th ed., 2001) — the source already cited for all 69 existing curated
   entries. Appendix A tabulates Tc/Pc/ω/Tb for several hundred common
   compounds directly — often faster than deriving ω by hand from NIST alone.
3. **DDBST / Dortmund Data Bank** (ddbst.com) — commercial; the free tier is a
   single-compound lookup form, not bulk export. Useful for spot-checking or
   filling a specific gap, but respect their terms — don't attempt bulk
   pulls.
4. **CoolProp** (MIT-licensed, github.com/CoolProp/CoolProp) — the one source
   that's legitimately bulk-pullable: `pip install CoolProp` and read
   `CoolProp.CoolProp.PropsSI("Tcrit"/"pcrit"/"acentric", fluid)`, or read its
   fluid JSON files directly on GitHub. Covers ~120-130 industrial fluids with
   citations already attached per fluid. Good for a quick batch of common
   fluids; carry CoolProp's copyright notice in the `source` field or a
   nearby attribution comment when you use it (MIT requires the notice be
   preserved, not a full relicense).

## Priority order

Given there's no realistic path to 2000+ compounds by hand, work top-down by
industrial relevance rather than alphabetically. The 69 already curated cover
light gases, C1-C10 n-alkanes, a few branched alkanes, water/ammonia/CO2/H2S,
BTX minus two isomers, a handful of common solvents (acetone, esters, ethers,
chlorinated solvents), a few common refrigerants, and some acids/amines
(see `crates/tpt-thermo-data/data/seed.toml`'s non-`estimated` entries for the
exact list). The tiers below are what's still missing, roughly in the order a
process-simulation databank would prioritize:

**Tier 1 — high-frequency industrial chemicals**
- Glycols: ethylene glycol, propylene glycol, diethylene glycol — dominant in
  gas dehydration and antifreeze systems
- Alkanolamines: MEA, DEA, MDEA — dominant amines in CO2/H2S gas sweetening
- Polymer monomers: styrene, vinyl chloride, vinyl acetate, propylene oxide,
  ethylene oxide, acrylonitrile — huge production-volume feedstocks
- Remaining BTX-family isomers: o-xylene, m-xylene, cumene (isopropylbenzene)

**Tier 2 — common process/refrigeration chemicals**
- Refrigerants: R22, R134a, R32, R125 (R245fa is already present)
- Glycol ethers: 2-methoxyethanol, 2-ethoxyethanol
- C4-C8 alcohols/acids/esters not yet covered: n-butanol, isobutanol,
  propionic acid, butyric acid, methyl acetate, propyl acetate
- Additional amines: diethylamine, triethylamine, morpholine

**Tier 3 — broader solvent/specialty coverage**
- More aromatics/naphthenes: styrene oxide, indene, methylcyclohexane,
  decalin
- Halogenated solvents: 1,2-dichloroethane, trichloroethylene,
  perchloroethylene
- Nitriles/sulfur compounds beyond what's curated: propionitrile,
  dimethyl sulfide, dimethyl sulfoxide

**Tier 4 — everything else**
- Lower-priority/niche compounds, worked through as time allows. Do **not**
  try to work through the 2650 auto-generated placeholder names — they're
  filler, not a priority queue.

## Adding data to the dataset

1. Copy `crates/tpt-thermo-data/data/component_import_template.toml`'s
   `[[components]]` block once per compound and fill it in with your looked-up
   values, including a real `source` string (e.g. `"NIST WebBook"` or
   `"Poling et al. 2001"`).
2. Choose how to load it:
   - **Permanent addition (typical case)**: paste the filled block(s) into
     `crates/tpt-thermo-data/data/seed.toml`. It's compiled into the crate via
     `include_str!` (`crates/tpt-thermo-data/src/seed.rs`), so
     `cargo build -p tpt-thermo-data` / `cargo test -p tpt-thermo-data`
     picks it up automatically — no code changes needed.
   - **Staging/standalone file**: keep new entries in a separate `.toml` and
     load them at runtime with
     `SeedComponentDatabase::from_toml_str(&contents)` (or `from_json_str`
     for JSON) — useful for validating a batch before deciding to merge it
     into the shipped seed.
3. If you also have a literature/fitted binary interaction parameter for a
   pair already in the database, use the template's commented-out
   `[[binary_interactions]]` block (`crates/tpt-thermo-data/src/bip.rs`) the
   same way.
4. Run `cargo test -p tpt-thermo-data` — `ComponentRecord::validate()` and the
   duplicate-name check will catch malformed entries immediately.

## Provenance rule

Every value added this way must carry a real `source` string naming exactly
where it came from. Never leave `source` blank, and never fill it with a
guess or an LLM-recalled number presented as if it were looked up — the
dataset's usefulness depends on every value being traceable to a citation.
