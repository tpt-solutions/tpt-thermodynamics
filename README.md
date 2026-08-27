# tpt-thermodynamics

Multi-component mixture properties, phase equilibria, and equations of state
for the TPT Solutions engineering / materials / earth / energy stack.

This workspace is the thermodynamics layer referenced by downstream `tpt-*`
crates (process modelling, materials, earth sciences, energy systems). It is
built bottom-up: a `no_std`-capable core trait surface, then data, then the
equation-of-state families, then flash / phase-stability / bubble-dew solvers,
transport properties, electrolytes, and polymers, finally re-exported through a
flat-feature umbrella crate.

## Crate inventory

Build order follows `todo.md` (Phases 2-13). Each crate depends only on earlier
layers plus the upstream `tpt-math-*` / `tpt-eng-props` substrate.

| Order | Crate | Purpose |
|------:|-------|---------|
| 1 | `tpt-thermo-core` | Foundation: `EquationOfState` trait, composition types, mixing-rule traits, convergence enums. `no_std` + `alloc`. |
| 2 | `tpt-thermo-data` | Component database, TOML/JSON (de)serialization, curated seed compound set, BIP tables, provenance. |
| 3 | `tpt-thermo-eos-cubic` | Cubic EoS: PR / SRK / volume-translated PR, alpha functions, vdW1f / Huron-Vidal / Wong-Sandler mixing. |
| 4 | `tpt-thermo-eos-activity` | Activity models: NRTL, UNIQUAC, Wilson, UNIFAC (+ Dortmund). Implements `ExcessGibbsModel`. |
| 5 | `tpt-thermo-eos-saft` | PC-SAFT, association term, SAFT-VR Mie. |
| 6 | `tpt-thermo-flash` | Rachford-Rice / Newton flash: PT, PH, TV, TS, PU, PV, LLE. |
| 7 | `tpt-thermo-phase` | Phase stability (TPD), multiphase equilibrium, SLE, critical locus, continuation. Implements `StabilityTest`. |
| 8 | `tpt-thermo-bubble-dew` | Bubble / dew point solvers, phase envelopes, azeotrope + cricondenbar/therm detection. |
| 9 | `tpt-thermo-transport` | Viscosity, thermal conductivity, diffusivity, mixture averaging, residual-entropy scaling. |
| 10 | `tpt-thermo-electrolyte` | Pitzer, eNRTL, HKF, ion association, gas solubility. |
| 11 | `tpt-thermo-polymer` | Flory-Huggins, Sanchez-Lacombe, PC-SAFT-for-polymers, cloud point, MWD. |
| 12 | `tpt-thermo` | Umbrella crate: flat feature tree re-exporting every constituent. |

## Build order & dependencies

The crates are consumed/published like every other `tpt-*` repo: upstream
dependencies (`tpt-math-*`, `tpt-eng-props`) are crates.io **version strings**,
never path/git deps across repo boundaries. In-workspace path deps are added to
`[workspace.dependencies]` as each phase lands.

> **Note:** this repo does **not** publish to crates.io as part of the
> build-out (left to the user). See `todo.md` "Known Deferred Scope".

A `no_std` + `alloc` split lives only in `tpt-thermo-core`; every other crate is
`std`-based.

## Quick start

```bash
cargo xtask check      # fast compile check across the workspace
cargo xtask test       # tests + doctests
cargo xtask all        # fmt, clippy, test, deny, wasm, check
```

`just` recipes mirror the `cargo xtask` commands above.

## Tracking

- `spec.txt` — full design specification (phases, validation targets, feature
  matrix).
- `todo.md` — per-session progress tracker for the bootstrap + 12-crate
  build-out.

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at
your option. © TPT Solutions.
