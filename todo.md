# tpt-thermodynamics — Build Todo

> Tracks bootstrap + full 12-crate build-out for the tpt-thermodynamics repo, per
> `spec.txt`. This is a **multi-session** effort (~16-18 sessions estimated). Crates
> are consumed/published like every other `tpt-*` repo: upstream deps
> (`tpt-math-*`, `tpt-eng-props`) are crates.io version strings, not path/git deps.
> **Per explicit decision, this repo does not publish to crates.io itself** during
> this build-out — that is left to the user to do separately later. License for
> every crate: `MIT OR Apache-2.0`. Author: TPT Solutions. Full phase rationale
> lives in the approved plan file (see project notes); this file is the actual
> per-session progress tracker.

## Phase 0 — Repo Bootstrap

- [x] Root `Cargo.toml` (`[workspace]` resolver "2", `[workspace.package]`,
      `[workspace.dependencies]` seeded with upstream version-string deps,
      `[workspace.lints]`, `[profile.release]`)
- [x] `rust-toolchain.toml` (stable + rustfmt/clippy + `wasm32-unknown-unknown`)
- [x] `rustfmt.toml`
- [x] `deny.toml`
- [x] `.cargo/config.toml` (`xtask` alias)
- [x] `.gitignore`
- [x] `xtask/` crate (`Cargo.toml` + `src/main.rs`: `fmt`, `clippy`, `test`, `deny`,
       `wasm` build check, `check`, `new-crate`, `all` subcommands)
- [x] `justfile` mirroring `cargo xtask` commands
- [x] `LICENSE-MIT`, `LICENSE-APACHE` (copied verbatim from sibling repo)
- [x] Root `README.md` (purpose, crate inventory table, build-order note, links to
       `spec.txt`/`todo.md`)
- [x] `CONTRIBUTING.md`, `SECURITY.md`, `CHANGELOG.md`
- [x] `.github/workflows/ci.yml` (self-contained: fmt / clippy / nextest test+doc /
       wasm / cargo-deny / msrv — matches actual sibling-repo practice, not shared
       reusable workflows)
- [x] `.github/dependabot.yml`
- [x] `examples/` crate scaffold (empty `Cargo.toml` + `src/`, populated per-phase)
- [x] `git init`; initial commit
- [x] This `todo.md`
- [x] Sanity check: `cargo build` succeeds on the empty (xtask + examples only)
       workspace; `cargo xtask check` passes

**Done when:** `cargo xtask check` is green on a workspace containing only
`xtask`/`examples`, all standard docs exist, and this file reflects Phases 0-13 as
tracked TODOs.

---

## Phase 1 — Upstream prerequisite work in `tpt-math` (sibling repo)

> **This phase's work happens in `d:/Programming/1PRODUCTION/Open Source/tpt-math`,
> not this repo.** `tpt-thermo-core`'s `EquationOfState` trait and every
> `Pressure`/`Temperature`/composition type depend on these two crates directly, and
> both are currently thin, already-published 0.1.0 wrappers missing the
> thermodynamics-specific surface (confirmed by reading their source — see below).

### 1.1 — `tpt-math-numeric`: root-finding, nonlinear-system, and ODE solvers

Currently (78 lines): only `num_traits` re-exports + a `Scalar` supertrait. **Missing:**
no root finder, no nonlinear-system solver, no ODE/IVP solver.

- [ ] `src/root.rs`: bisection, Newton-Raphson, **Brent's method** (spec's explicit
      Newton fallback for flash/phase solvers)
- [ ] `src/nonlinear.rs`: thin wrapper over `tpt-math-optimize-general`'s existing
      `minimize_newton`/`minimize_newton_with` for n-dim `F(x)=0` systems, damped-step
      (ω relaxation) support
- [ ] `src/ode.rs`: RK4 + adaptive RKF45/RK45 IVP solver (needed for arc-length
      continuation in Phases 8-9, HKF path integration in Phase 11)
- [ ] Wire deps: `tpt-math-optimize-general`, `nalgebra`/`tpt-math-linalg` as needed
      for vector state in ODE solves (verify against tpt-math's actual linalg backend
      — its README states Apache-2.0-only deps incl. `nalgebra` were replaced with
      in-house `tpt-math-linalg*`; confirm which crate `tpt-math-optimize-general`
      actually uses before assuming `nalgebra` is available)
- [ ] Define a local, domain-agnostic `Convergence`-shaped result type (do NOT
      depend on `tpt-thermo-core` — keep this crate generic/reusable)
- [ ] Unit tests: root finders vs. known analytical roots; nonlinear solve vs. a
      small 2-3 equation system; ODE solvers vs. a known analytical IVP
- [ ] Doctests + rustdoc for every new public item
- [ ] `cargo fmt` / `clippy` / `cargo deny check` clean
- [ ] no_std verify (if kept no_std)
- [ ] Bump `0.1.0 -> 0.2.0`; `CHANGELOG.md` entry
- [ ] **Do not publish to crates.io** (explicit user decision — local commit only)

### 1.2 — `tpt-math-units`: thermodynamic quantity aliases

Currently (76 lines): aliases `Area`, `Length`, `Mass`, `Ratio`,
`ThermodynamicTemperature`, `Time`, `Velocity`, `Volume` from `uom`. **Missing:**
`Pressure`, `MolarVolume`, `EnergyPerMol`/`MolarEnergy`, `AmountOfSubstance`,
`MolarMass`, `MolarHeatCapacity`, `DynamicViscosity`, `ThermalConductivity`,
`DiffusionCoefficient` — all already exist as `uom` ISQ quantities, this is aliasing
only, no new dimensional math.

- [ ] Extend `pub mod q` with the aliases listed above
- [ ] Alias `uom`'s `ThermodynamicTemperature` as `Temperature` deliberately,
      documented clearly re: the absolute-vs-interval-temperature footgun
- [ ] `EnergyPerMol` type alias (`= MolarEnergy`) matching spec's exact
      `EquationOfState` trait signature name
- [ ] Extend `prelude` to bring the new aliases into scope
- [ ] Verify no gap exists vs. `uom`'s ISQ before assuming any custom
      `system!`/`quantity!` extension is needed
- [ ] Unit tests: conversion round-trip per new alias (incl. degC<->K respecting the
      interval/absolute distinction)
- [ ] Doctests + rustdoc; `cargo fmt`/`clippy`/`deny` clean; no_std verify
- [ ] Bump version; `CHANGELOG.md` entry
- [ ] **Do not publish to crates.io** (explicit user decision — local commit only)

### 1.3 — Land

- [ ] Run tpt-math's full `cargo xtask check` + workspace test suite (regression
      check against ~20 other `tpt-math-*` crates depending on these two)
- [ ] Commit in `tpt-math` (separate commits per crate, per its existing convention)
- [ ] **Open item:** since these stay unpublished, `tpt-thermodynamics`'s
      `Cargo.toml` version-string deps for `tpt-math-numeric`/`tpt-math-units` will
      not resolve until the user publishes them (or a temporary, clearly-labeled
      `[patch.crates-io]` path override is added and removed before any real
      release) — confirm approach with the user at the start of Phase 2.

**Done when:** both crates compile/test/lint clean in `tpt-math` and expose the
surface above, committed locally.

---

## Per-Crate Checklist Template (Phases 2-13)

**Standard crate:**
1. Scaffold `crates/<name>/` (`Cargo.toml` inheriting workspace fields, `lib.rs` stub)
2. Wire dependencies (in-workspace path deps + upstream version-string deps)
3. Implement scope (see phase-specific breakdown below)
4. Unit tests + doctests
5. Rustdoc (crate-level + public API)
6. `cargo fmt --check` / `cargo clippy --all-targets --all-features -- -D warnings` clean
7. `cargo deny check` clean
8. no_std verification (only if the crate is `no_std`)
9. Validation targets met (spec sec6, against curated seed dataset — see Deferred Scope)
10. Changelog entry

**Umbrella crate (`tpt-thermo`, Phase 13):** scaffold with feature-gated optional
deps -> wire optional deps + matching flat feature flags -> re-export each
constituent's public API behind its feature -> rustdoc documenting the feature
matrix -> fmt/clippy/deny clean across representative feature combinations
(default, `--all-features`, each Tier-2 consumption profile from spec sec7).

---

## Phase 2 — `tpt-thermo-core`

*Foundation layer. Build order 1/12. `no_std` + `alloc`. Depends on:
`tpt-math-units`, `tpt-math-numeric`, `tpt-math-linalg` (all upstream).*

- [x] Scaffold `crates/tpt-thermo-core/`
- [x] Composition types (`src/composition.rs`): `MoleFraction`, `MassFraction`,
      `Molality` newtypes + normalization + conversion utilities
- [x] `ConvergenceStatus` / `DivergenceReason` / `NumericalIssueReason` enums
      (`src/convergence.rs`, per spec sec4's exact shape)
- [x] `EquationOfState` trait (`src/eos.rs`): pressure, fugacity coefficient,
      enthalpy, entropy, heat capacity, speed of sound, compressibility /
      thermal-expansion / molar-volume-solve (Brent); numerical-default
      implementations where the spec allows; ideal-gas reference impl
- [x] Mixing-rule surface (`src/mixing.rs`): `MixingRule` trait + **forward-declared**
      `ExcessGibbsModel` and `StabilityTest` traits. Concrete vdW1f / Huron-Vidal /
      Wong-Sandler combiners land in Phase 4 (cubic) and Phase 5 (activity).
- [x] `ComponentDatabase` trait (`src/component.rs`) with unit-safe accessors
- [x] `BipParameter` / `ParameterSource` / provenance structs (`src/provenance.rs`)
      — `chrono::NaiveDate` replaced by a `no_std`-friendly `SourceDate` (year/month/day)
- [x] Forward-declared `ExcessGibbsModel` trait (Phase 4 defines usage, Phase 5
      implements it — avoids cyclic cubic<->activity crate dependency)
- [x] Forward-declared `StabilityTest` trait (Phase 7 defines usage, Phase 8
      implements it — avoids cyclic flash<->phase crate dependency)
- [x] `no_std`/`alloc` split: `Vec`-returning methods gated behind `alloc`
      (verified `cargo build --no-default-features --features alloc`)
- [x] Toy ideal-gas `EquationOfState` impl (tests + living documentation)
- [x] Unit tests, doctests, rustdoc, fmt/clippy/deny clean, no_std verify
- [x] `examples/` entry: composition conversion + ideal-gas EoS toy

**Done when:** compiles standalone + no_std-verified, full trait/enum/struct
surface above exposed, ideal-gas reference impl passing.

---

## Phase 3 — `tpt-thermo-data`

*Build order 2/12. Depends on: `tpt-thermo-core` (path), `tpt-eng-props` (upstream).*

- [x] Scaffold `crates/tpt-thermo-data/`
- [x] 3a: `ComponentRecord` schema + TOML/JSON (de)serialization (serde) + physical-
      constraint validation
- [x] 3b: `ComponentDatabase` impl backed by a curated **seed set** (24 compounds:
      water, CO2, methane, ethane, propane, n-alkanes C4–C8, N2, O2, H2, Ar, He,
      benzene, toluene, ethanol, methanol, NH3, H2S, ethylene, propylene, HCl).
      Expanding to ~50-100 is straightforward and tracked as Deferred Scope.
- [x] 3c: BIP tables — `BipTable` structure + loader shipped; fitted values seeded
      alongside Phase 4/5 (every pair defaults to 0.0 until then)
- [ ] 3d: Parameter estimation utilities — deferred to Phase 4+ (per spec)
- [x] Simple schema-version field for data versioning (not a full audit-log system)
- [x] Unit tests: schema validation edge cases, TOML/JSON round-trip, seed-dataset
      sanity checks vs. literature values
- [x] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** `ComponentDatabase` implemented for the seed dataset sufficient to
validate Phases 4-9, user-supplied TOML/JSON components load/validate, provenance
attached to every seeded value.

---

## Phase 4 — `tpt-thermo-eos-cubic`

*Build order 3/12. Depends on: `tpt-thermo-core`, `tpt-thermo-data`.*

- [ ] Scaffold `crates/tpt-thermo-eos-cubic/`
- [ ] PR (+ modified variants) `src/pr.rs`, SRK `src/srk.rs`, volume-translated PR
      (Peneloux) `src/volume_translation.rs` — all implementing `EquationOfState`
- [ ] Alpha functions (`src/alpha.rs`): Soave, Twu, Mathias-Copeman via
      `AlphaFunction` trait
- [ ] van der Waals 1-fluid mixing with T-dependent BIPs (`k_ij = a + b/T + c*ln(T)`)
- [ ] Huron-Vidal (MHV1, MHV2, PSRK), generic over `tpt-thermo-core`'s
      `ExcessGibbsModel` trait
- [ ] Wong-Sandler mixing rules (`src/wong_sandler.rs`)
- [ ] Cardano's method cubic root solver + physically-meaningful-root selection via
      stability criteria (`src/cubic_solver.rs`)
- [ ] Critical point detection, spinodal curve, mechanical stability (`src/critical.rs`)
- [ ] Validation: pure-component density <1%, enthalpy <2%, vapor pressure <3% vs.
      seed compounds (spec sec6)
- [ ] Doctests, rustdoc, fmt/clippy/deny clean, `examples/` entry (PR P-V-T calc)

**Done when:** PR/SRK/vPR pass pure-component validation targets for the seed set,
Cardano root selection robust across 2/3-real-root cases, full `EquationOfState`
trait implemented.

---

## Phase 5 — `tpt-thermo-eos-activity`

*Build order 4/12. Depends on: `tpt-thermo-core` (implements `ExcessGibbsModel`),
`tpt-thermo-data`.*

- [ ] Scaffold `crates/tpt-thermo-eos-activity/`
- [ ] NRTL (`src/nrtl.rs`), UNIQUAC (`src/uniquac.rs`)
- [ ] UNIFAC original + Dortmund modified (`src/unifac/`) — seed group-parameter
      table only; full group coverage tracked as Deferred Scope
- [ ] Wilson (`src/wilson.rs`)
- [ ] eNRTL/Pitzer electrolyte extensions **explicitly deferred to Phase 11**
- [ ] Temperature-dependent parameter helper (`A + B/T + C*ln(T)`), infinite-
      dilution limiting-law tests
- [ ] Validation: pressure <5%, temperature <2K, vapor composition <0.02 mole
      fraction vs. 10-20 seed binary pairs (spec sec6, full 100+ tracked as
      Deferred Scope)
- [ ] Integration test: Huron-Vidal (Phase 4) consuming this crate's models via
      `ExcessGibbsModel`
- [ ] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** NRTL/UNIQUAC/Wilson pass infinite-dilution + VLE validation for the
seed set, UNIFAC predicts without fitting, Huron-Vidal cross-crate coupling tested.

---

## Phase 6 — `tpt-thermo-eos-saft`

*Build order 5/12. Expect its own multi-session sub-effort. Depends on:
`tpt-thermo-core`, `tpt-thermo-data`, `tpt-math-numeric`,
`tpt-math-optimize-general`/`-convex`.*

- [ ] Scaffold `crates/tpt-thermo-eos-saft/`
- [ ] PC-SAFT (`src/pc_saft.rs`): hard-chain reference, dispersion term
- [ ] Association term (`src/association.rs`): 1/2/3/4-site schemes,
      cross-association (water-alcohol style), solvation; Newton-Raphson w/
      analytical Jacobian, returns `ConvergenceStatus`
- [ ] SAFT-VR Mie (`src/saft_vr_mie.rs`)
- [ ] eSAFT electrolyte extension: basic ion-ion/ion-solvent/ion-segment term if
      schedule allows, else explicit Deferred Scope item (not silently dropped)
- [ ] Full derivative set (analytical per Gross & Sadowski 2001 where practical,
      numerical-default fallback from `tpt-thermo-core`)
- [ ] Parameter estimation utilities (fit to pure-component/binary data)
- [ ] Validation: density/enthalpy vs. REFPROP-style data for associating fluids
      (water, alcohols) in the seed set
- [ ] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** PC-SAFT reproduces literature density/vapor-pressure within spec
sec6 tolerances for seed compounds, association solver converges for
cross-associating mixtures (e.g. water-ethanol) with `ConvergenceStatus` reporting.

---

## Phase 7 — `tpt-thermo-flash`

*Build order 6/12. Depends on: `tpt-thermo-core`, `tpt-thermo-eos-cubic`,
`tpt-thermo-eos-activity`, `tpt-thermo-eos-saft`, `tpt-thermo-data`.*

- [ ] Scaffold `crates/tpt-thermo-flash/`
- [ ] Rachford-Rice successive substitution (`src/rachford_rice.rs`) + Wilson/UNIFAC
      K-value initialization (`src/initialization.rs`)
- [ ] Newton-Raphson w/ full Jacobian (`src/newton_flash.rs`)
- [ ] PT, PH, TV, TS, PU, PV flash variants (`src/{pt,ph,tv,ts,pu,pv}.rs`) — PT first
- [ ] Near-critical density-based fallback (`src/density_based.rs`)
- [ ] Trace-component `ln K_i` handling (threaded through RR/Newton, not separate)
- [ ] LLE isoactivity flash
- [ ] **VLL nested-loop flash deferred to Phase 8** (needs `tpt-thermo-phase`'s
      `StabilityTest` — hard sequencing dependency, tracked explicitly, closed out
      in Phase 8, not dropped)
- [ ] Convergence acceleration (`src/acceleration.rs`): dominant eigenvalue, volume
      substitution, GDEM
- [ ] `flash_pt_batch` — straightforward per-composition loop first; explicit SIMD
      tracked as Deferred Scope follow-up
- [ ] Validation: phase fraction <1%, composition <0.01 mole fraction vs. 3-5 seed
      multi-component systems (spec sec6, full 20+ tracked as Deferred Scope)
- [ ] Criterion benchmark harness (`benches/flash_pt.rs`) targeting <1ms/10-component
- [ ] Doctests, rustdoc, fmt/clippy/deny clean, `examples/` entry (PT flash, PR,
      natural-gas-like mixture)

**Done when:** PT/PH/TV/TS/PU/PV all converge (`ConvergenceStatus::Converged`) on
seed systems within tolerance, LLE works for a known partially-miscible pair,
criterion benchmark exists.

---

## Phase 8 — `tpt-thermo-phase`

*Build order 7/12. Depends on: `tpt-thermo-core`, `tpt-thermo-flash`. Implements
`StabilityTest`.*

- [ ] Scaffold `crates/tpt-thermo-phase/`
- [ ] TPD minimization (`src/tpd.rs`): Michelsen method, successive substitution ->
      Newton-Raphson refinement
- [ ] Multiple trial-composition initialization strategies (`src/trial_compositions.rs`)
- [ ] `StabilityResult` struct (phase count, compositions, status)
- [ ] Multiphase equilibrium V-L-L / V-L-L-L / L-L-L (`src/multiphase.rs`)
- [ ] SLE with T-dependent solubility (`src/sle.rs`)
- [ ] Mixture critical point calculation (Heidemann-Rahal) + continuation
      (`src/critical_locus.rs`)
- [ ] Phase boundary arc-length continuation (`src/continuation.rs`)
- [ ] **Close out Phase 7's deferred VLL flash** as a cross-crate integration test
- [ ] Validation: TPD correctly classifies stable/unstable for seed
      azeotrope/miscibility-gap/near-critical systems (spec sec6, full 30+ tracked
      as Deferred Scope)
- [ ] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** TPD classifies correctly for the seed set, VLL flash works
end-to-end via this crate, mixture critical point located for >=1 known binary.

---

## Phase 9 — `tpt-thermo-bubble-dew`

*Build order 8/12. Depends on: `tpt-thermo-core`, `tpt-thermo-flash`,
`tpt-thermo-phase`.*

- [ ] Scaffold `crates/tpt-thermo-bubble-dew/`
- [ ] Bubble point (`src/bubble.rs`): Newton on Σ K_i x_i = 1, both "find T" / "find P"
- [ ] Dew point (`src/dew.rs`): Newton on Σ x_i = Σ y_i/K_i = 1
- [ ] Phase envelope continuation (`src/envelope.rs`): P-T, P-x-y, T-x-y, reusing
      Phase 8's continuation machinery
- [ ] Azeotrope detection (`src/azeotrope.rs`)
- [ ] Cricondenbar/cricondentherm detection (`src/cricondentherm.rs`)
- [ ] Reactive distillation: implement only if trivially composable with existing
      flash machinery (reaction kinetics out-of-scope per spec sec2), else Deferred
      Scope item
- [ ] Validation: pressure <5%, temperature <2K, vapor composition <0.02 vs. seed set
- [ ] Doctests, rustdoc, fmt/clippy/deny clean, `examples/` entry (phase envelope,
      e.g. ethanol-water)

**Done when:** bubble/dew solvers converge for the seed binary set, >=1 full phase
envelope traced end-to-end, azeotrope detection flags a known pair.

---

## Phase 10 — `tpt-thermo-transport`

*Build order 9/12. Depends on: `tpt-thermo-core`, `tpt-thermo-data`. Lower-coupling
— can parallelize against Phases 8-9 with multiple contributors.*

- [ ] Scaffold `crates/tpt-thermo-transport/`
- [ ] Viscosity (`src/viscosity.rs`): Chung et al., Lucas, corresponding-states
- [ ] Thermal conductivity (`src/conductivity.rs`): Chung et al., Ely-Hanley,
      corresponding-states
- [ ] Diffusivity (`src/diffusivity.rs`): Fuller-Schettler-Giddings, Darken, Vignes
- [ ] Mixture averaging (`src/mixing_rules.rs`): Wilke, Mason-Saxena, Filippov, Darken
- [ ] Residual entropy scaling (`src/residual_entropy_scaling.rs`)
- [ ] Unit-safe throughout (`DynamicViscosity`, `ThermalConductivity`,
      `DiffusionCoefficient` aliases from Phase 1.2)
- [ ] Validation: viscosity <10%, conductivity <15%, diffusivity <20% vs. seed set
- [ ] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** each correlation family passes its spec sec6 target for the seed
set, ideal-gas/incompressible-liquid limiting behavior verified.

---

## Phase 11 — `tpt-thermo-electrolyte`

*Build order 10/12. Depends on: `tpt-thermo-core`, `tpt-thermo-eos-activity`
(eNRTL reuses Phase 5's NRTL), `tpt-thermo-data`.*

- [ ] Scaffold `crates/tpt-thermo-electrolyte/`
- [ ] Pitzer (`src/pitzer.rs`): β⁰/β¹/β²/C_φ, mixed-electrolyte θ/ψ, to 6 molal
- [ ] eNRTL (`src/enrtl.rs`): long-range Pitzer-Debye-Hückel + short-range NRTL
      (composes Phase 5's NRTL — this is where Phase 5's deferred electrolyte
      extension actually lands)
- [ ] HKF (`src/hkf.rs`): high-T/P (to 1000°C, 5 kbar) standard partial molal
      properties, leans on Phase 1.1 ODE solvers for path integration
- [ ] Ion association (`src/ion_association.rs`): Bjerrum criterion, mass-action
- [ ] Gas solubility (`src/solubility.rs`): Setschenow equation
- [ ] Debye-Hückel limiting-law infinite-dilution tests across all three models
- [ ] Validation: activity coefficient <5%, osmotic coefficient <2%, to 6 molal vs.
      seed single-electrolyte set (e.g. NaCl-H2O) + a handful of mixed systems
- [ ] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** Pitzer reproduces literature values for seed set, eNRTL composes
correctly with Phase 5's NRTL, HKF produces sane properties across a documented
T/P path.

---

## Phase 12 — `tpt-thermo-polymer`

*Build order 11/12. Depends on: `tpt-thermo-core`, `tpt-thermo-eos-saft` (chain
m->infinity limit), `tpt-thermo-phase` (cloud point).*

- [ ] Scaffold `crates/tpt-thermo-polymer/`
- [ ] Flory-Huggins (`src/flory_huggins.rs`): combinatorial entropy + χ parameter
- [ ] Sanchez-Lacombe lattice-fluid EoS (`src/sanchez_lacombe.rs`), implements
      `EquationOfState`
- [ ] PC-SAFT-for-polymers (`src/pc_saft_polymer.rs`): thin specialization of Phase
      6, regression-tested to reduce correctly to Phase 6's PC-SAFT
- [ ] Cloud point (`src/cloud_point.rs`): UCST/LCST via Phase 8's stability machinery
- [ ] Molecular weight distribution (`src/mwd.rs`): Schulz-Zimm, most-probable
- [ ] Polymer-solvent χ parameter estimation from VLE/LLE/osmotic-pressure data
      (`src/parameter_estimation.rs`)
- [ ] Crystallization / melting-point depression via Flory equation
      (`src/crystallization.rs`)
- [ ] Validation: cloud-point prediction vs. >=1 literature UCST/LCST system
- [ ] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** Flory-Huggins + Sanchez-Lacombe implement the workspace's
EoS/activity interfaces consistently, cloud-point locates a known UCST/LCST,
PC-SAFT-for-polymers regression-tested against Phase 6's limit.

---

## Phase 13 — `tpt-thermo` (umbrella)

*Build order 12/12. Mirrors `tpt-eng-props`'s umbrella pattern
(`tpt-engineering/crates/tpt-eng-props/`).*

- [ ] Scaffold `crates/tpt-thermo/` — optional deps on all 11 constituent crates
- [ ] Flat feature tree (`cubic`, `activity`, `saft`, `flash`, `phase`,
      `bubble-dew`, `transport`, `electrolyte`, `polymer`, `data`) — no
      auto-implied features, per spec sec3's explicit flat-tree mandate
- [ ] `default = []`; `tpt-thermo-core`/`-data` always re-exported (non-optional)
- [ ] `src/lib.rs`: `pub use tpt_thermo_core as core;` / `pub use tpt_thermo_data as
      data;` always; feature-gated `pub use` for every other constituent crate
- [ ] Unified `ThermoError` (`src/error.rs`) with feature-gated variants
- [ ] Composition conversion utilities re-exposed (implemented in Phase 2, not
      reimplemented here)
- [ ] High-level convenience API (`src/api.rs`): `FlashCalculator` builder,
      `bubble_point(...)`, `flash_pt(...)` top-level functions
- [ ] Doctests for every public convenience function; rustdoc documenting the full
      feature matrix
- [ ] fmt/clippy/deny clean across: default, `--all-features`, each Tier-2
      consumption profile (spec sec7)
- [ ] `examples/` entries per Tier-2 profile (`tpt-process`/`tpt-materials`/
      `tpt-earth`-shaped minimal builds)

**Done when:** `cargo build --no-default-features` yields only core+data
re-exports, each feature flag builds standalone, `--all-features` builds, each
Tier-2 example compiles with only its listed feature subset.

---

## Known Deferred Scope

Explicit tracking so intentionally-reduced scope isn't silently lost:

- [ ] `tpt-thermo-data`: full 2000+ compound coverage (Phase 3 ships ~50-100 seed
      compounds only)
- [ ] `tpt-thermo-eos-activity`: full UNIFAC group table (Phase 5 ships seed groups
      only)
- [ ] `tpt-thermo-eos-saft`: full eSAFT electrolyte extension, if not completed
      alongside Phase 6/11
- [ ] `tpt-thermo-flash`: explicit SIMD vectorization for `flash_pt_batch`
      (loop-level ships first)
- [ ] `tpt-thermo-bubble-dew`: reactive distillation (likely skipped — needs
      out-of-scope reaction-equilibrium machinery)
- [ ] `tpt-thermo-data`: parameter estimation utilities beyond the Phase 4+ minimal
      set
- [ ] Full spec sec6 validation breadth (100+ binary VLE pairs, 20+ multicomponent
      flash systems, 30+ stability systems) — every phase validates against a
      curated seed set first; expanding to full spec sec6 counts is a tracked
      per-crate follow-up
- [ ] Publish `tpt-thermodynamics` (and the bumped `tpt-math-numeric`/
      `tpt-math-units`) to crates.io — intentionally left to the user, not done as
      part of this build-out
- [ ] Phase 1 (`tpt-math` sibling repo) was **not** performed: the root-finding /
      nonlinear / ODE solvers and the extended units aliases were instead
      implemented locally in `tpt-thermo-core` (`src/numerics.rs` + the
      `Temperature` / `EnergyPerMol` / `MolarEntropy` aliases in `src/quantities.rs`)
      because the published `tpt-math-numeric` / `tpt-math-units` 0.1.0 crates are
      thin wrappers. If Phase 1 is later done upstream, `tpt-thermo-core` should be
      re-pointed at the published surface and its local copies removed.

## Suggested session granularity

Roughly one phase per session (Phase 1 as two sub-sessions, one per `tpt-math`
crate), with Phase 6 (`tpt-thermo-eos-saft`) and Phase 11
(`tpt-thermo-electrolyte`) likely needing 2+ sessions each. ~16-18 sessions
minimum for the full build-out plus upstream prerequisite work.
