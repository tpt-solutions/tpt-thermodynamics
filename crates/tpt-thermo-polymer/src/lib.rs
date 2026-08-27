//! `tpt-thermo-polymer` — polymer thermodynamics for `tpt-thermodynamics`.
//!
//! Covers the models a process/材料 stack needs for polymer systems:
//!
//! * [`flory_huggins`] — Flory–Huggins activity model (combinatorial entropy plus a
//!   `χ` interaction term), binary and multicomponent.
//! * [`mwd`] — molecular-weight distributions (Schulz–Zimm, most-probable) and their
//!   number/weight/moment averages.
//! * [`crystallization`] — Flory melting-point depression.
//! * [`sanchez_lacombe`] — Sanchez–Lacombe lattice-fluid equation of state
//!   ([`tpt_thermo_core::EquationOfState`] implementation).
//! * [`pc_saft_polymer`] — a thin PC-SAFT-for-polymers specialisation reusing
//!   [`tpt_thermo_eos_saft::PcSaft`], regression-tested to reduce to plain PC-SAFT.
//! * [`cloud_point`] — L–L phase-split (cloud-point) detection via the
//!   [`tpt_thermo_phase`] tangent-plane-distance machinery.
//! * [`parameter_estimation`] — estimate a Flory–Huggins `χ` from a single LLE
//!   tie-line.
//!
//! Quantitative validation against full literature datasets (UCST/LCST cloud points,
//! polymer-solvent VLE/LLE, osmotic-pressure `χ`) is tracked as Deferred Scope
//! (consistent with the rest of this build-out); the implementations below are
//! internally consistent and checked against analytical limiting cases.

pub mod cloud_point;
pub mod crystallization;
pub mod flory_huggins;
pub mod mwd;
pub mod parameter_estimation;
pub mod pc_saft_polymer;
pub mod sanchez_lacombe;
