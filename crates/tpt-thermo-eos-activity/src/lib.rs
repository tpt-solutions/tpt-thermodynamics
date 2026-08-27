//! `tpt-thermo-eos-activity` — liquid-phase activity models for
//! `tpt-thermodynamics`.
//!
//! This crate implements the [`tpt_thermo_core::mixing::ExcessGibbsModel`] trait
//! (forward-declared in the core) for the classical local-composition and
//! group-contribution models, so the cubic EoS mixing rules (Huron-Vidal,
//! Wong-Sandler — Phase 4) can consume them directly:
//!
//! * [`Nrtl`] — Non-Random Two-Liquid (Renon & Prausnitz, 1968).
//! * [`Wilson`] — Wilson's equation (1964).
//! * [`Uniquac`] — UNIQUAC (Abrams & Prausnitz, 1975).
//! * [`unifac::UnifacModel`] — original and Dortmund-modified UNIFAC
//!   (Fredenslund et al., 1975), with a seed group table.
//!
//! All models are `no_std`/`alloc` and take their parameters explicitly (the
//! seed dataset does not yet carry activity-parameter tables — tracked as
//! Deferred Scope in `todo.md`).
//!
//! # Example
//!
//! ```
//! use tpt_thermo_core::mixing::ExcessGibbsModel;
//! use tpt_thermo_core::quantities::{Temperature, Pressure};
//! use tpt_thermo_eos_activity::{Nrtl, parameters::TdParam};
//! use uom::si::{thermodynamic_temperature::kelvin, pressure::pascal};
//!
//! let nrtl = Nrtl::binary(
//!     TdParam::new(0.5, 100.0, 0.0),
//!     TdParam::new(-0.2, 50.0, 0.0),
//!     0.3,
//! ).unwrap();
//! let t = Temperature::new::<kelvin>(333.15);
//! let p = Pressure::new::<pascal>(1.0e5);
//! let ge = nrtl.reduced_excess_gibbs(t, p, &[0.5, 0.5]).unwrap();
//! let lng1 = nrtl.ln_gamma(t, p, &[0.5, 0.5], 0).unwrap();
//! assert!(ge.is_finite());
//! assert!(lng1.is_finite());
//! ```

#![cfg_attr(not(feature = "std"), no_std)]
// Matrix/tensor algebra in this crate is clearest with explicit index loops
// over `i`/`j`; `needless_range_loop` is a false positive for double-indexed
// access (e.g. `a[i][j]` with `x[j]`), so it is allowed crate-wide.
#![allow(clippy::needless_range_loop)]

extern crate alloc;

pub mod nrtl;
pub mod parameters;
pub mod unifac;
pub mod uniquac;
pub mod wilson;

pub use nrtl::Nrtl;
pub use unifac::{UnifacModel, UnifacVariant};
pub use uniquac::{StructuralParams, Uniquac};
pub use wilson::Wilson;
