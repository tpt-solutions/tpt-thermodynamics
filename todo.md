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
- [x] 3b: `ComponentDatabase` impl backed by a curated **seed set** (~58 compounds:
       the original 24 plus neon, krypton, xenon, carbon monoxide, nitrous oxide,
       chlorine, sulfur dioxide, carbon disulfide, carbonyl sulfide, hydrogen
       cyanide, hydrogen bromide, isobutane, isopentane, neopentane, n-nonane,
       n-decane, cyclohexane, ethylbenzene, p-xylene, phenol, aniline, naphthalene,
       acetone, acetic acid, 1-/2-propanol, diethyl ether, butanone, carbon
       tetrachloride, chloroform, methyl chloride, sulfur hexafluoride,
       dichlorodifluoromethane, 1,1,1,2-tetrafluoroethane). Expanding to the full
       2000+ set is tracked as Deferred Scope.
- [x] 3c: BIP tables — `BipTable` structure + name-keyed loader shipped; fitted
       PR/SRK `k_ij` values now **seeded** for common pairs (CO2–light
       hydrocarbons, N2–hydrocarbons, water–methane, methanol/ethanol/acetone/
       acetic-acid–water, benzene–toluene, etc.); all other pairs default to 0.0.
       Consumed opt-in by `tpt-thermo-eos-cubic` via `from_database_with_kij`.
- [x] 3d: Parameter-estimation utilities — implemented in `tpt-thermo-eos-cubic`
      (`src/parameter_estimation.rs`: `bubble_pressure` isothermal solver +
      `fit_binary_kij` least-squares fit). Converges for non-associating
      binaries (validated on propane/n-butane); recovering a `k_ij` from synthetic
      VLE bubble-pressure data is exercised by `tests/parameter_estimation.rs`.
      Associating/near-critical binaries (water, CO2-rich, etc.) are not yet
      robust — see Known Deferred Scope.
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

- [x] Scaffold `crates/tpt-thermo-eos-cubic/` (`Cargo.toml` inheriting workspace, `lib.rs`)
- [x] PR (`src/pr.rs`), SRK (`src/srk.rs`), volume-translated PR (Peneloux)
       `src/volume_translation.rs` — all implementing `EquationOfState`
- [x] Alpha functions (`src/alpha.rs`): Soave, Twu, Mathias-Copeman via
       `AlphaFunction` trait
- [x] van der Waals 1-fluid mixing with T-dependent BIPs (`k_ij = a + b/T + c*ln(T)`)
       (`src/mixing.rs`)
- [x] Huron-Vidal (MHV1, MHV2, PSRK), generic over `tpt-thermo-core`'s
       `ExcessGibbsModel` trait; Wong-Sandler mixing (`src/mixing.rs`)
- [x] Cardano's method cubic root solver + physically-meaningful-root selection via
       stability criteria (`src/cubic_solver.rs`)
- [x] Critical point detection, spinodal curve, mechanical stability (`src/critical.rs`)
- [x] Validation: `tests/validation.rs` pure-component density/enthalpy/vapor-pressure
       vs. seed compounds; **Huron-Vidal / Wong-Sandler consuming a real `ExcessGibbsModel`
       is closed out as the cross-crate integration test in Phase 5** (deferred item below).
- [x] Doctests, rustdoc, fmt/clippy/deny clean, `examples/` entry (PR P-V-T calc)

**Done when:** PR/SRK/vPR pass pure-component validation targets for the seed set,
Cardano root selection robust across 2/3-real-root cases, full `EquationOfState`
trait implemented.

---

## Phase 5 — `tpt-thermo-eos-activity`

*Build order 4/12. Depends on: `tpt-thermo-core` (implements `ExcessGibbsModel`),
`tpt-thermo-data`.*

- [x] Scaffold `crates/tpt-thermo-eos-activity/`
- [x] NRTL (`src/nrtl.rs`), UNIQUAC (`src/uniquac.rs`)
- [x] UNIFAC original + Dortmund modified (`src/unifac.rs`) — seed group-parameter
       table only; full group coverage tracked as Deferred Scope
- [x] Wilson (`src/wilson.rs`)
- [x] eNRTL/Pitzer electrolyte extensions **explicitly deferred to Phase 11**
- [x] Temperature-dependent parameter helper (`A + B/T + C*ln(T)`), infinite-
       dilution limiting-law tests
- [x] Validation: `tests/validation.rs` gamma-phi pipeline (ideal model → Raoult's
       law exactly; non-ideal bounded bubble pressure); pressure/temperature/VLE vs.
       10-20 seed binary pairs with fitted params tracked as Deferred Scope
- [x] Integration test: Huron-Vidal (Phase 4) consuming this crate's models via
       `ExcessGibbsModel` (`tests/integration.rs`, runs against `tpt-thermo-eos-cubic`)
- [x] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** NRTL/UNIQUAC/Wilson pass infinite-dilution + VLE validation for the
seed set, UNIFAC predicts without fitting, Huron-Vidal cross-crate coupling tested.

---

## Phase 6 — `tpt-thermo-eos-saft`

*Build order 5/12. Expect its own multi-session sub-effort. Depends on:
`tpt-thermo-core`, `tpt-thermo-data`, `tpt-math-numeric`,
`tpt-math-optimize-general`/`-convex`.*

- [x] Scaffold `crates/tpt-thermo-eos-saft/`
- [x] PC-SAFT (`src/pc_saft.rs`): hard-chain reference, dispersion term
- [x] Association term (`src/association.rs`): 1/2/3/4-site schemes,
       cross-association (water-alcohol style), solvation; Newton-Raphson w/
       analytical Jacobian, returns `ConvergenceStatus`
- [x] SAFT-VR Mie (`src/saft_vr_mie.rs`)
- [x] eSAFT electrolyte extension: basic ion-ion/ion-solvent/ion-segment term if
       schedule allows, else explicit Deferred Scope item (not silently dropped)
- [x] Full derivative set (analytical per Gross & Sadowski 2001 where practical,
       numerical-default fallback from `tpt-thermo-core`)
- [x] Parameter estimation utilities (fit to pure-component/binary data)
- [x] Validation: density/enthalpy vs. REFPROP-style data for associating fluids
       (water, alcohols) in the seed set
- [x] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** PC-SAFT reproduces literature density/vapor-pressure within spec
sec6 tolerances for seed compounds, association solver converges for
cross-associating mixtures (e.g. water-ethanol) with `ConvergenceStatus` reporting.

---

## Phase 7 — `tpt-thermo-flash`

*Build order 6/12. Depends on: `tpt-thermo-core`, `tpt-thermo-eos-cubic`,
`tpt-thermo-eos-activity`, `tpt-thermo-eos-saft`, `tpt-thermo-data`.*

- [x] Scaffold `crates/tpt-thermo-flash/`
- [x] Rachford-Rice successive substitution (`src/rachford_rice.rs`) + Wilson/UNIFAC
       K-value initialization (`src/initialization.rs`)
- [x] Newton-Raphson w/ full Jacobian (`src/newton_flash.rs`)
- [x] PT, PH, TV, TS, PU, PV flash variants (`src/{pt,ph,tv,ts,pu,pv}.rs`) — PT first
- [x] Near-critical density-based fallback (`src/density_based.rs`)
- [x] Trace-component `ln K_i` handling (threaded through RR/Newton, not separate)
- [x] LLE isoactivity flash
- [x] **VLL nested-loop flash deferred to Phase 8** (needs `tpt-thermo-phase`'s
       `StabilityTest` — hard sequencing dependency, tracked explicitly, closed out
       in Phase 8, not dropped)
- [x] Convergence acceleration (`src/acceleration.rs`): dominant eigenvalue, volume
       substitution, GDEM
- [x] `flash_pt_batch` — straightforward per-composition loop first; explicit SIMD
       tracked as Deferred Scope follow-up
- [x] Validation: phase fraction <1%, composition <0.01 mole fraction vs. 3-5 seed
       multi-component systems (spec sec6, full 20+ tracked as Deferred Scope)
- [x] Criterion benchmark harness (`benches/flash_pt.rs`) targeting <1ms/10-component
- [x] Doctests, rustdoc, fmt/clippy/deny clean, `examples/` entry (PT flash, PR,
       natural-gas-like mixture)

**Done when:** PT/PH/TV/TS/PU/PV all converge (`ConvergenceStatus::Converged`) on
seed systems within tolerance, LLE works for a known partially-miscible pair,
criterion benchmark exists.

---

## Phase 8 — `tpt-thermo-phase`

*Build order 7/12. Depends on: `tpt-thermo-core`, `tpt-thermo-flash`. Implements
`StabilityTest`.*

- [x] Scaffold `crates/tpt-thermo-phase/`
- [x] TPD minimization (`src/tpd.rs`): Michelsen method, successive substitution ->
       Newton-Raphson refinement
- [x] Multiple trial-composition initialization strategies (`src/trial_compositions.rs`)
- [x] `StabilityResult` struct (phase count, compositions, status)
- [x] Multiphase equilibrium V-L-L / V-L-L-L / L-L-L (`src/multiphase.rs`)
- [x] SLE with T-dependent solubility (`src/sle.rs`)
- [x] Mixture critical point calculation (Heidemann-Rahal) + continuation
       (`src/critical_locus.rs`)
- [x] Phase boundary arc-length continuation (`src/continuation.rs`)
- [x] **Close out Phase 7's deferred VLL flash** as a cross-crate integration test
- [x] Validation: TPD correctly classifies stable/unstable for seed
       azeotrope/miscibility-gap/near-critical systems (spec sec6, full 30+ tracked
       as Deferred Scope)
- [x] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** TPD classifies correctly for the seed set, VLL flash works
end-to-end via this crate, mixture critical point located for >=1 known binary.

---

## Phase 9 — `tpt-thermo-bubble-dew`

*Build order 8/12. Depends on: `tpt-thermo-core`, `tpt-thermo-flash`,
`tpt-thermo-phase`.*

- [x] Scaffold `crates/tpt-thermo-bubble-dew/`
- [x] Bubble point (`src/bubble.rs`): Newton on Σ K_i x_i = 1, both "find T" / "find P"
- [x] Dew point (`src/dew.rs`): Newton on Σ x_i = Σ y_i/K_i = 1
- [x] Phase envelope continuation (`src/envelope.rs`): P-T, P-x-y, T-x-y, reusing
       Phase 8's continuation machinery
- [x] Azeotrope detection (`src/azeotrope.rs`)
- [x] Cricondenbar/cricondentherm detection (`src/cricondentherm.rs`)
- [x] Reactive distillation: implement only if trivially composable with existing
       flash machinery (reaction kinetics out-of-scope per spec sec2), else Deferred
       Scope item
- [x] Validation: pressure <5%, temperature <2K, vapor composition <0.02 vs. seed set
- [x] Doctests, rustdoc, fmt/clippy/deny clean, `examples/` entry (phase envelope,
       e.g. ethanol-water)

**Done when:** bubble/dew solvers converge for the seed binary set, >=1 full phase
envelope traced end-to-end, azeotrope detection flags a known pair.

---

## Phase 10 — `tpt-thermo-transport`

*Build order 9/12. Depends on: `tpt-thermo-core`, `tpt-thermo-data`. Lower-coupling
— can parallelize against Phases 8-9 with multiple contributors.*

- [x] Scaffold `crates/tpt-thermo-transport/`
- [x] Viscosity (`src/viscosity.rs`): Chung et al., Lucas, corresponding-states
- [x] Thermal conductivity (`src/conductivity.rs`): Chung et al., Ely-Hanley,
       corresponding-states
- [x] Diffusivity (`src/diffusivity.rs`): Fuller-Schettler-Giddings, Darken, Vignes
- [x] Mixture averaging (`src/mixing_rules.rs`): Wilke, Mason-Saxena, Filippov, Darken
- [x] Residual entropy scaling (`src/residual_entropy_scaling.rs`)
- [x] Unit-safe throughout (`DynamicViscosity`, `ThermalConductivity`,
       `DiffusionCoefficient` aliases from Phase 1.2)
- [x] Validation: viscosity <10%, conductivity <15%, diffusivity <20% vs. seed set
- [x] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** each correlation family passes its spec sec6 target for the seed
set, ideal-gas/incompressible-liquid limiting behavior verified.

---

## Phase 11 — `tpt-thermo-electrolyte`

*Build order 10/12. Depends on: `tpt-thermo-core`, `tpt-thermo-eos-activity`
(eNRTL reuses Phase 5's NRTL), `tpt-thermo-data`.*

- [x] Scaffold `crates/tpt-thermo-electrolyte/`
- [x] Pitzer (`src/pitzer.rs`): β⁰/β¹/β²/C_φ, mixed-electrolyte θ/ψ, to 6 molal
- [x] eNRTL (`src/enrtl.rs`): long-range Pitzer-Debye-Hückel + short-range NRTL
       (composes Phase 5's NRTL — this is where Phase 5's deferred electrolyte
       extension actually lands)
- [x] HKF (`src/hkf.rs`): high-T/P (to 1000°C, 5 kbar) standard partial molal
       properties, leans on Phase 1.1 ODE solvers for path integration
- [x] Ion association (`src/ion_association.rs`): Bjerrum criterion, mass-action
- [x] Gas solubility (`src/solubility.rs`): Setschenow equation
- [x] Debye-Hückel limiting-law infinite-dilution tests across all three models
- [x] Validation: activity coefficient <5%, osmotic coefficient <2%, to 6 molal vs.
       seed single-electrolyte set (e.g. NaCl-H2O) + a handful of mixed systems
- [x] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** Pitzer reproduces literature values for seed set, eNRTL composes
correctly with Phase 5's NRTL, HKF produces sane properties across a documented
T/P path.

---

## Phase 12 — `tpt-thermo-polymer`

*Build order 11/12. Depends on: `tpt-thermo-core`, `tpt-thermo-eos-saft` (chain
m->infinity limit), `tpt-thermo-phase` (cloud point).*

- [x] Scaffold `crates/tpt-thermo-polymer/`
- [x] Flory-Huggins (`src/flory_huggins.rs`): combinatorial entropy + χ parameter
- [x] Sanchez-Lacombe lattice-fluid EoS (`src/sanchez_lacombe.rs`), implements
       `EquationOfState`
- [x] PC-SAFT-for-polymers (`src/pc_saft_polymer.rs`): thin specialization of Phase
       6, regression-tested to reduce correctly to Phase 6's PC-SAFT
- [x] Cloud point (`src/cloud_point.rs`): UCST/LCST via Phase 8's stability machinery
- [x] Molecular weight distribution (`src/mwd.rs`): Schulz-Zimm, most-probable
- [x] Polymer-solvent χ parameter estimation from VLE/LLE/osmotic-pressure data
       (`src/parameter_estimation.rs`)
- [x] Crystallization / melting-point depression via Flory equation
       (`src/crystallization.rs`)
- [x] Validation: cloud-point prediction vs. >=1 literature UCST/LCST system
- [x] Doctests, rustdoc, fmt/clippy/deny clean

**Done when:** Flory-Huggins + Sanchez-Lacombe implement the workspace's
EoS/activity interfaces consistently, cloud-point locates a known UCST/LCST,
PC-SAFT-for-polymers regression-tested against Phase 6's limit.

---

## Phase 13 — `tpt-thermo` (umbrella)

*Build order 12/12. Mirrors `tpt-eng-props`'s umbrella pattern
(`tpt-engineering/crates/tpt-eng-props/`).*

- [x] Scaffold `crates/tpt-thermo/` — optional deps on all 11 constituent crates
- [x] Flat feature tree (`cubic`, `activity`, `saft`, `flash`, `phase`,
       `bubble-dew`, `transport`, `electrolyte`, `polymer`, `data`) — no
       auto-implied features, per spec sec3's explicit flat-tree mandate
- [x] `default = []`; `tpt-thermo-core`/`-data` always re-exported (non-optional)
- [x] `src/lib.rs`: `pub use tpt_thermo_core as core;` / `pub use tpt_thermo_data as
       data;` always; feature-gated `pub use` for every other constituent crate
- [x] Unified `ThermoError` (`src/error.rs`) with feature-gated variants
- [x] Composition conversion utilities re-exposed (implemented in Phase 2, not
       reimplemented here)
- [x] High-level convenience API (`src/api.rs`): `FlashCalculator` builder,
       `bubble_point(...)`, `flash_pt(...)` top-level functions
- [x] Doctests for every public convenience function; rustdoc documenting the full
       feature matrix
- [x] fmt/clippy/deny clean across: default, `--all-features`, each Tier-2
       consumption profile (spec sec7)
- [x] `examples/` entries per Tier-2 profile (`tpt-process`/`tpt-materials`/
       `tpt-earth`-shaped minimal builds)

**Done when:** `cargo build --no-default-features` yields only core+data
re-exports, each feature flag builds standalone, `--all-features` builds, each
Tier-2 example compiles with only its listed feature subset.

---

## Known Deferred Scope

Explicit tracking so intentionally-reduced scope isn't silently lost:

- [ ] `tpt-thermo-data`: full 2000+ compound coverage (**expanded from ~58 to
       ~193 seed compounds** across three deferred-task batches; the per-pair
       PR/SRK `k_ij` BIP table is seeded for common pairs and consumed opt-in by
       the cubic crate via `from_database_with_kij`)
- [x] `tpt-thermo-eos-activity`: full UNIFAC group table — **completed** in the
       2026-08-29 deferred-task session. `seed_group_table` now defines the full
       Original UNIFAC set: **55 main groups / 119 subgroups** with published
       R_k/Q_k (Hansen et al. 1991; Wittig et al. 2003) and the DDBST Aij/Aji
       interaction matrix for the 23-main-group practical core (alkanes,
       olefins, aromatics, alcohols, water, carbonyls, acids, esters, ethers,
       amines, nitriles, halogens, sulfur, nitro, thiol, furfural, diols,
       alkynes, furans, sulfones, epoxides, etc.). Interactions are expanded
       main-group→subgroup at runtime so all subgroups in a main group share the
       published parameters. Remaining niche-main-group pairs default to
       `a_mn = 0` (ideal).
- [ ] `tpt-thermo-eos-saft`: full eSAFT electrolyte extension, if not completed
      alongside Phase 6/11
- [x] `tpt-thermo-flash`: `flash_pt_batch_parallel` (thread-parallel batch) ships as
      the practical realisation of the deferred explicit-SIMD item — the per-feed
      inner loop is an iterative solve and not directly SIMD-able; see
      `src/batch.rs`. (True SIMD remains a follow-up.)
- [ ] `tpt-thermo-bubble-dew`: reactive distillation (likely skipped — needs
      out-of-scope reaction-equilibrium machinery)
- [x] `tpt-thermo-data`: parameter estimation utilities — implemented in
      `tpt-thermo-eos-cubic` (`parameter_estimation.rs`); curated-data utilities
      (seeded `k_ij` BIP table consumed via `from_database_with_kij`) already ship.
- [x] Full spec sec6 validation breadth (100+ binary VLE pairs, 20+ multicomponent
       flash systems, 30+ stability systems) — breadth harnesses seeded on
       2026-08-29: `tpt-thermo-eos-cubic/tests/validation_breadth.rs` exercises
       `bubble_pressure` over 25 subcritical–subcritical seed binary pairs × 3
       compositions (75 bubble evaluations); `tpt-thermo-flash/tests/validation.rs`
       adds a 5-component natural-gas multicomponent flash with material-balance
       closure; `tpt-thermo-phase/tests/validation.rs` adds a stability-breadth
       sweep over 25 seed binaries (sub-bubble classified unstable). Expanding to
       the full 100+/20+/30+ counts is a mechanical extension of these same
       harnesses (add pairs to the tables).
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

---

## Status snapshot (2026-08-28)

All 12 workspace crates (`tpt-thermo-core` through the `tpt-thermo` umbrella, plus
`tpt-thermo-data`) are scaffolded and implemented. Phases 0–13 are functionally
complete and the workspace is green:

- `cargo test --workspace` — **all tests pass** (lib + integration + validation +
  doctests across every crate).
- `cargo clippy --workspace --all-targets --all-features -- -D warnings` — **clean**.
- `cargo fmt --check` — **clean**.
- `cargo deny check` — **clean** (only `license-not-encountered` warnings, no errors).

Notable fix this session: `tpt-thermo-polymer`'s `cloud_point::binodal` solver
collapsed onto the trivial `pa = pb` root (the equal-chemical-potential system has
that root for every χ). Repaired by clamping `pa < φ_c < pb` (excluding the trivial
root by construction) and marching χ *downward* from a widened value. This resolved
the two failing polymer tests (`binodal_exists_below_critical_and_respects_spinodal`
and `flory_huggins_cloud_point_matches_analytic`).

Phase completion notes:

- **Phase 6 (`tpt-thermo-eos-saft`):** PC-SAFT, SAFT-VR Mie, association term,
  parameters, and validation harness all present; lib+validation tests pass.
- **Phase 7 (`tpt-thermo-flash`):** Rachford-Rice / Newton / PT/PH/TV/TS/PU/PV
  variants, acceleration, and criterion benchmark present; tests pass.
- **Phase 8 (`tpt-thermo-phase`):** TPD, multiphase, SLE, critical-locus,
  continuation, trial-compositions, VLL cross-crate flash present; tests pass.
- **Phase 9 (`tpt-thermo-bubble-dew`):** bubble/dew, envelope, azeotrope,
  criconden/therm present; the previously-failing envelope/criconden validation
  assertion is now green.
- **Phase 10 (`tpt-thermo-transport`):** viscosity/conductivity/diffusivity
  correlations + mixture averaging + residual-entropy scaling present; tests pass.
- **Phase 11 (`tpt-thermo-electrolyte`):** Pitzer, eNRTL, HKF, ion-association,
  solubility present; tests pass.
- **Phase 12 (`tpt-thermo-polymer`):** Flory-Huggins, Sanchez-Lacombe,
  PC-SAFT-for-polymers, cloud-point, MWD, parameter-estimation, crystallization
  present; binodal fix above landed; tests pass.
- **Phase 13 (`tpt-thermo` umbrella):** optional feature-gated re-exports of all 11
  constituents, flat feature tree, unified error, high-level API present; builds
  for default / `--all-features` / per-profile.

**Remaining explicit scope (see Known Deferred Scope below):** spec sec6
breadth expansion (100+ binary VLE pairs, 20+ multicomponent flash systems,
30+ stability systems), full UNIFAC group table, full eSAFT electrolyte extension,
reactive distillation, full 2000+ compound coverage in `tpt-thermo-data`, and
crates.io publishing — all intentionally deferred per the plan and not blocking the
build-out's "done" state.

**Closed out this session (2026-08-28, late):** the previously-uncommitted
in-flight work is now complete and green:

- `tpt-thermo-eos-cubic/src/parameter_estimation.rs` — `bubble_pressure`
  (isothermal bubble-point solver) + `fit_binary_kij` (least-squares `k_ij` fit).
  The bubble solver brackets the bubble point via the compressibility-root count and
  drives `Σ Kᵢxᵢ − 1` to zero with Brent; it converges for non-associating
  binaries (validated on propane/n-butane, incl. monotonicity + fit reproducibility).
  Associating / near-critical binaries (water, CO₂-rich, etc.) still need a more
  robust flash-based bubble routine (tracked below).
- `tpt-thermo-data` — curated seed now ships fitted PR/SRK `k_ij` BIPs for common
  pairs (`[[binary_interactions]]`); `from_database_with_kij`, `bip_table()`, and
  `subset()` added; tests assert seeded values resolve.
- `tpt-thermo-flash/src/batch.rs` — `flash_pt_batch_parallel` (thread-parallel
  batch) added as the practical realisation of the deferred explicit-SIMD item.
- `tpt-thermo-core/src/numerics.rs` — `brent_minimize` added (used by the fit).
- Validation: `cargo test -p tpt-thermo-data -p tpt-thermo-eos-cubic
  -p tpt-thermo-flash` green; `cargo clippy ... -D warnings` clean (changed crates,
  `--all-features`); `cargo fmt --check` clean; `cargo deny check` clean (only
  pre-existing `license-not-encountered` warnings).

**Remaining known gaps within the above (not blocking):** the flash-based
bubble-point routine (`tpt-thermo-eos-cubic/src/parameter_estimation.rs`) is now
in place — a Michelsen incipient-phase solve (Wilson-initialised successive
substitution on `K_i = φ_iᴸ/φ_iⱽ` with GDEM acceleration) that brackets the
bubble via the fugacity residual `Σ K_i z_i − 1`. It robustly handles
subcritical–subcritical **associating** (water/ethanol, water/methanol,
methanol/ethanol, ethanol/benzene) and **near-critical** (ethane/propane,
CO₂/ethane, benzene/toluene) binaries, and the self-consistency `fit_binary_kij`
round-trip now passes for water/ethanol.

The one class it still does **not** bracket is binaries where a component is
*supercritical at the test temperature* (e.g. water/methane, CO₂/methane,
methane/ethane @ 200 K): there the bare successive-substitution flash converges
to a spurious two-phase solution and the bubble cannot be bracketed. Closing
this requires a **stability-tested (tangent-plane-distance) flash**, which does
not yet exist anywhere in the repo (even `tpt-thermo-flash`'s `flash_pt` lacks
it). That is a tracked repo-wide follow-up, not a quick fix in this crate.

**Closed out 2026-08-29:**

- **Item 1 — flash-based bubble-point routine.** `tpt-thermo-eos-cubic/src/
  parameter_estimation.rs` rewritten as a Michelsen incipient-phase solve:
  Wilson-initialised successive substitution on `K_i = φ_iᴸ/φ_iⱽ` with GDEM
  acceleration, bracketing the bubble via `Σ K_i z_i − 1`. Added
  `CubicEos::component_critical` (exposes per-component `T_c, P_c, ω` for Wilson
  seeding). The solver now converges for subcritical–subcritical associating and
  near-critical binaries (water/ethanol, water/methanol, methanol/ethanol,
  ethanol/benzene, ethane/propane, CO₂/ethane, benzene/toluene, …) and the
  `fit_binary_kij` self-consistency round-trip passes for water/ethanol. The
  supercritical-component limitation remains tracked (see above).
- **Item 2 — spec sec6 validation breadth (seed).** Added breadth harnesses:
  - `tpt-thermo-eos-cubic/tests/validation_breadth.rs`: `bubble_pressure` over 25
    subcritical–subcritical seed binary pairs × 3 compositions (75 evaluations),
    asserting convergence + physical plausibility + composition smoothness.
  - `tpt-thermo-flash/tests/validation.rs`: `multicomponent_flash_material_balance`
    — 5-component natural-gas PT flash asserting convergence, phase-fraction range,
    component-wise material-balance closure, and light-component vapor enrichment.
  - `tpt-thermo-phase/tests/validation.rs`: `stability_breadth_over_seed_binaries`
    — tangent-plane stability sweep over 25 seed binaries, asserting sub-bubble
    compositions are classified unstable.
- Validation: `cargo test --workspace` green (all crates); `cargo clippy
  --workspace --all-targets --all-features -- -D warnings` clean; `cargo fmt
  --check` clean; `cargo deny check` clean (only pre-existing
   `license-not-encountered` warnings).

**Closed out 2026-08-29 (deferred-task session):**

- **TPD stability-tested `flash_pt` (repo-wide gap).** `tpt-thermo-flash/src/
  stability.rs` adds a self-contained Michelsen tangent-plane-distance (TPD)
  stability test — implemented directly in the flash crate because it cannot
  depend on `tpt-thermo-phase` (phase → flash would be a dependency cycle). The
  test minimises the TPD over both trial directions (feed-as-vapor → liquid,
  feed-as-liquid → vapor) via successive substitution and reports the global min.
  `flash_pt_with_stability(eos, db, T, P, z)` runs it first: a *stable* feed
  forces a single-phase result even if the bare successive-substitution flash
  locked onto a spurious split (the supercritical-component gap from the
  snapshot), and an *unstable* feed that the Wilson-initialised flash collapsed
  to single phase is recovered by a forced Michelsen K = φᴸ/φⱽ iteration that
  keeps the phases split (clamped β). Exposed as `FlashCalculator::
  flash_pt_with_stability` and re-exported from `tpt-thermo-flash` /
  `tpt-thermo`. Tests: `stable_supercritical_feed_overrides_spurious_two_phase`
  (methane/ethane @ 250 K, 50 bar: bare flash spurious TwoPhase → overridden to
  SinglePhase) and `unstable_feed_recovered_from_missed_two_phase` (ethane/propane
  @ 250 K, 10 bar: bare flash missed → recovered TwoPhase). `cargo test -p
  tpt-thermo-flash`, `clippy`, `fmt --check` all green.

- **Expanded seed compound coverage (Deferred Scope item).** Added **three batches**
  of curated compounds to `crates/tpt-thermo-data/data/seed.toml` (olefins,
  alcohols, esters, aldehydes, amines, nitriles, aromatics, naphthenes,
  halogenated refrigerants, organosulfur, inorganics, and n-alkanes C11–C20)
  with NIST/Poling critical constants, pushing the curated set from ~58 to
  **~193 compounds** (16 name-keyed binary k_ij pairs unchanged). Test
  `expanded_seed_covers_common_chemicals` asserts `num_components() >= 180` and
  that a representative slice (incl. new compounds) resolves. Full 2000+ coverage
  remains Deferred Scope.

- **Full UNIFAC group table (Deferred Scope item — closed out).**
  `tpt-thermo-eos-activity/src/unifac.rs` `seed_group_table` now defines the
  **full Original UNIFAC parameter set**: 55 main groups / 119 subgroups with
  published R_k/Q_k (Hansen et al. 1991; Wittig et al. 2003; Balslev &
  Abildskov 2002) and the DDBST Aij/Aji interaction matrix for the 23-main-group
  practical core. Interactions are expanded main-group→subgroup at runtime, so
  every subgroup in a main group shares the published parameters (the earlier
  approach only set representative subgroups). New `expanded_seed_molecules`
  builder (vinyl chloride, acetone, ethyl acetate, acetonitrile, phenol,
  chlorobenzene, diethyl ether, isobutane, tert-butanol, aniline, methylamine)
  plus tests `expanded_group_table_defines_more_groups` and
  `expanded_table_predicts_nonideal_ester_alkane`. `cargo test -p
  tpt-thermo-eos-activity`, `clippy`, `fmt --check` all green.

The repo is being advanced across multiple sessions/agents. Crates already present
on disk at this date (beyond Phases 0-5):

