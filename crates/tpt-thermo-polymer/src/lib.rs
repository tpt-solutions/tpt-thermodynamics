//! `tpt-thermo-polymer` — polymer thermodynamics for `tpt-thermodynamics`.
//!
//! Phase 12 of the build-out. This crate collects the polymer-relevant models and
//! provides them through the same workspace interfaces used everywhere else:
//!
//! * [`FloryHuggins`] — the classic Flory-Huggins combinatorial + `χ` activity
//!   model, implementing [`tpt_thermo_core::mixing::ExcessGibbsModel`].
//! * [`SanchezLacombe`] — the Sanchez-Lacombe lattice-fluid equation of state,
//!   implementing [`tpt_thermo_core::EquationOfState`].
//! * [`pc_saft_polymer::PolymerPcSaft`] — a thin specialization of Phase 6's
//!   PC-SAFT for polymer chains (large segment count `m`), re-using the existing
//!   [`tpt_thermo_eos_saft::PcSaft`] engine.
//! * [`cloud_point`] — UCST/LCST cloud-point (binodal/spinodal) computation for
//!   binary polymer solutions.
//! * [`mwd`] — molecular-weight distributions (Schulz-Zimm, most-probable) and
//!   their moments.
//! * [`parameter_estimation`] — fitting `χ` from osmotic-pressure / activity data.
//! * [`crystallization`] — Flory melting-point depression.
//!
//! All models are `no_std`/`alloc` and unit-safe (uom-backed).

#![cfg_attr(not(feature = "std"), no_std)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::many_single_char_names)]

extern crate alloc;

pub mod cloud_point;
pub mod crystallization;
pub mod flory_huggins;
pub mod mwd;
pub mod parameter_estimation;
pub mod pc_saft_polymer;
pub mod sanchez_lacombe;

pub use cloud_point::{binodal, critical_point, ChiTemperature, CloudPointKind, CloudPointResult};
pub use crystallization::flory_melting_depression;
pub use flory_huggins::FloryHuggins;
pub use mwd::{most_probable, schulz_zimm, MolecularWeightDistribution};
pub use parameter_estimation::chi_from_osmotic_pressure;
pub use pc_saft_polymer::PolymerPcSaft;
pub use sanchez_lacombe::SanchezLacombe;
