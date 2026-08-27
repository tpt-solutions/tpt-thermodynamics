//! `tpt-thermo-eos-saft` — SAFT (statistical-associating-fluid-theory)
//! equations of state for `tpt-thermodynamics`.
//!
//! This crate implements the SAFT family on top of the
//! [`tpt_thermo_core::EquationOfState`] trait:
//!
//! * [`PcSaft`] — perturbed-chain SAFT (Gross & Sadowski 2001): hard-chain
//!   reference, dispersion perturbation, and the association (hydrogen-bonding)
//!   term via a Newton-Raphson site-fraction solver
//!   ([`association`](crate::association)).
//! * [`SaftVrMie`] — SAFT-VR Mie (Lafitte et al. 2013) sharing the same
//!   framework with per-component Mie repulsion/attraction ranges.
//!
//! Pure-component results are exact PC-SAFT; the mixture hard-chain uses the
//! Carnahan-Starling one-fluid approximation (documented refinement over full
//! bmcsL). Pressure is recovered from the packing-fraction derivative and
//! fugacity / enthalpy / entropy use the numerical-default composition and
//! temperature derivatives permitted by the build-out spec.
//!
//! # Example
//!
//! ```
//! use tpt_thermo_core::component::ComponentDatabase;
//! use tpt_thermo_core::EquationOfState;
//! use tpt_thermo_eos_saft::{PcSaft, parameters::SEED_SAFT_PARAMETERS};
//! use tpt_thermo_core::quantities::{Temperature, MolarVolume, Pressure};
//! use uom::si::{thermodynamic_temperature::kelvin, molar_volume::cubic_meter_per_mole, pressure::pascal};
//!
//! // Build PC-SAFT for methane from the curated parameter table.
//! let params = tpt_thermo_eos_saft::parameters::SaftParameters::new(
//!     SEED_SAFT_PARAMETERS.iter().filter(|c| c.name == "methane").copied().collect(),
//! );
//! let eos = PcSaft::new(params, vec![0.016_043]);
//! let t = Temperature::new::<kelvin>(300.0);
//! let v = MolarVolume::new::<cubic_meter_per_mole>(0.025);
//! let p = eos.pressure(t, v, &[1.0]).unwrap();
//! let phi = eos.ln_fugacity_coefficient(t, v, &[1.0], 0).unwrap();
//! // At 1 bar and 300 K methane is near-ideal: pressure ≈ R·T/v and ln φ ≈ 0.
//! assert!((p.value / 1.0e5 - 1.0).abs() < 0.1);
//! assert!(phi.abs() < 0.1);
//! ```

#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

pub mod association;
pub mod engine;
pub mod parameters;
pub mod pc_saft;
pub mod saft_vr_mie;

pub use association::AssociationResult;
pub use parameters::{
    AssociationParams, AssociationScheme, SaftComponent, SaftParameters, SEED_SAFT_PARAMETERS,
};
pub use pc_saft::PcSaft;
pub use saft_vr_mie::SaftVrMie;
