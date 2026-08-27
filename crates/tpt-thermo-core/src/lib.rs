//! `tpt-thermo-core` — foundation layer for `tpt-thermodynamics`.
//!
//! This crate defines the domain surface every other `tpt-thermo-*` crate is
//! built on, with no dependence on any other in-repo crate:
//!
//! * [`quantities`] — compile-time-typed thermodynamic quantity aliases over
//!   `uom` (the thin `tpt-math-units` 0.1.0 surface is extended here with the
//!   aliases the spec requires: [`quantities::Temperature`],
//!   [`quantities::EnergyPerMol`], [`quantities::MolarEntropy`]).
//! * [`composition`] — mole/mass fraction and molality newtypes plus a
//!   [`composition::Composition`] helper with conversions.
//! * [`convergence`] — [`convergence::ConvergenceStatus`] and its reason enums.
//! * [`numerics`] — root finders (bisection, Newton, Brent) generic over the
//!   `tpt-math-numeric` [`Scalar`](tpt_math_numeric::Scalar) trait. The spec
//!   expected these in `tpt-math-numeric`; they are provided here because the
//!   published 0.1.0 is a thin wrapper (see `todo.md` Phase 1 open item).
//! * [`eos`] — the [`eos::EquationOfState`] trait, a [`eos::State`] value, and a
//!   fully-working ideal-gas reference implementation.
//! * [`mixing`] — mixing-rule and excess-Gibbs/phase-stability traits (the
//!   latter two forward-declared so cubic and activity crates can couple
//!   without a cyclic dependency).
//! * [`component`] — the [`component::ComponentDatabase`] trait.
//! * [`provenance`] — parameter/value provenance structs.
//!
//! # `no_std`
//!
//! The crate is `#![no_std]` and works with only `alloc` available; `Vec`- and
//! `String`-returning APIs are gated behind the `alloc` feature (enabled by
//! `std`).

#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

pub mod component;
pub mod composition;
pub mod convergence;
pub mod eos;
pub mod error;
pub mod mixing;
pub mod numerics;
pub mod provenance;
pub mod quantities;

pub use component::ComponentDatabase;
pub use composition::{Composition, CompositionError, MassFraction, Molality, MoleFraction};
pub use convergence::{ConvergenceStatus, DivergenceReason, NumericalIssueReason};
pub use eos::{EquationOfState, IdealGas, State, StateBuilder};
pub use error::ThermoError;
pub use mixing::{ExcessGibbsModel, MixingRule, StabilityResult, StabilityTest};
pub use numerics::{bisection, brent, newton};
pub use provenance::{BipParameter, ParameterSource, Provenance, SourceDate};

/// Universal gas constant, `R`, in J·mol⁻¹·K⁻¹.
pub const R: f64 = 8.314_462_618_153_24;
