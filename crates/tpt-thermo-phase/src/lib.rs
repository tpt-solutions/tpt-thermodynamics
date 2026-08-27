//! `tpt-thermo-phase` — phase stability, multiphase equilibrium, solid–liquid
//! equilibrium, and mixture critical-point / arc-length-continuation tools.
//!
//! Phase 8 of the `tpt-thermodynamics` build-out. This crate implements the
//! core [`tpt_thermo_core::StabilityTest`] trait (forward-declared in Phase 2)
//! via Michelsen tangent-plane-distance (TPD) minimisation and builds
//! higher-level multiphase / SLE / critical-locus machinery on top of it.
//!
//! # Scope note
//!
//! The VLL nested-loop *flash* that the Phase-7 todo defers to this crate
//! (`tpt-thermo-phase`) requires the `tpt-thermo-flash` crate (Phase 7), which is
//! not yet present in this workspace. That cross-crate integration test is
//! tracked as deferred here (see [`multiphase`]): the stability / TPD core this
//! crate provides is the component it needs.

#![cfg_attr(not(feature = "std"), no_std)]

// Numerical code necessarily uses index-addressed loops over small matrices and
// single-letter quantities (a, b, v, t, …); these are stylistic, not bugs.
#![allow(clippy::needless_range_loop)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::similar_names)]
#![allow(clippy::too_many_arguments)]

extern crate alloc;

pub mod continuation;
pub mod critical_locus;
pub mod linalg;
pub mod multiphase;
pub mod phase_volume;
pub mod sle;
pub mod tpd;
pub mod trial_compositions;

pub use continuation::arc_length_step;
pub use critical_locus::{critical_locus_binary, mixture_critical_point, CriticalGuess};
pub use multiphase::{detect_phases, MultiphaseResult};
pub use phase_volume::{BrentPhaseVolume, PhaseVolume};
pub use sle::solid_liquid_solubility;
pub use tpd::{StabilityAnalyzer, TangentPlaneDistance, TpdSolution};
pub use trial_compositions::{pure_component_trials, regular_grid_trials, wilson_k_values};
