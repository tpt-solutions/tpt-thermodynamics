//! `tpt-thermo-electrolyte` — aqueous electrolyte thermodynamics for
//! `tpt-thermodynamics`.
//!
//! * [`pitzer`] — Pitzer (1973) molality-scale activity/osmotic coefficients for
//!   single and mixed electrolytes to ~6 molal.
//! * [`enrtl`] — eNRTL: long-range Pitzer–Debye–Hückel + short-range NRTL
//!   (composing `tpt-thermo-eos-activity`'s [`Nrtl`](tpt_thermo_eos_activity::Nrtl)).
//! * [`hkf`] — Helgeson–Kirkham–Flowers (1981) standard partial-molal properties
//!   as a function of `T, P`.
//! * [`ion_association`] — Bjerrum criterion and mass-action association constants.
//! * [`solubility`] — Setschenow gas solubility in electrolyte solutions.
//!
//! All routines target the Debye–Hückel limiting law at infinite dilution and
//! validate against it in the seed tests.

#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

pub mod enrtl;
pub mod hkf;
pub mod ion_association;
pub mod parameters;
pub mod pitzer;
pub mod solubility;

pub use parameters::{HkfParams, PitzerParams};
