//! `tpt-thermo-flash` — multiphase flash calculators for `tpt-thermodynamics`.
//!
//! This crate implements the classical flash algorithms on top of any
//! [`tpt_thermo_core::EquationOfState`]:
//!
//! * [`pt::flash_pt`] — isothermal/isobaric (Rachford–Rice) flash with successive
//!   substitution and convergence acceleration ([`acceleration`]).
//! * [`variants::flash_ph`], [`variants::flash_tv`], [`variants::flash_ts`],
//!   [`variants::flash_pu`], [`variants::flash_pv`] — the other four spec flash
//!   specifications, each a robust outer loop around [`pt::flash_pt`].
//! * [`lle::lle_isoactivity`] — liquid–liquid (LLE) flash driven by an
//!   [`tpt_thermo_core::mixing::ExcessGibbsModel`] (e.g. the NRTL/Wilson/UNIQUAC
//!   models in `tpt-thermo-eos-activity`).
//! * [`batch::flash_pt_batch`] — a per-composition loop over a feed table.
//! * [`batch::flash_pt_batch_parallel`] — a `std`-feature, cross-thread variant of
//!   the batch (the practical realisation of the deferred explicit-SIMD batch).
//!
//! K-values are built from the equilibrium fugacity equality `K_i = φ_i^L/φ_i^V`,
//! where each phase's molar volume is recovered by root-solving the EoS pressure
//! ([`phase_volume::phase_volume`]) — so the solvers are model-agnostic and work
//! for any cubic, SAFT, or activity-coupled EoS that implements the core trait.
//!
//! # Example
//!
//! ```
//! use tpt_thermo_core::component::ComponentDatabase;
//! use tpt_thermo_core::quantities::{Pressure, Temperature};
//! use tpt_thermo_eos_cubic::PengRobinson;
//! use tpt_thermo_data::SeedComponentDatabase;
//! use tpt_thermo_flash::FlashCalculator;
//! use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};
//!
//! let db = SeedComponentDatabase::from_seed();
//! let eos = PengRobinson::from_database(&db).unwrap();
//! let calc = FlashCalculator::with_db(&eos, &db);
//! let methane = db.index_of("methane").unwrap();
//! let ethane = db.index_of("ethane").unwrap();
//! let mut z = vec![0.0; db.num_components()];
//! z[methane] = 0.7; z[ethane] = 0.3;
//! let t = Temperature::new::<kelvin>(280.0);
//! let p = Pressure::new::<pascal>(2.0e6);
//! let res = calc.flash_pt(t, p, &z).unwrap();
//! // A two-phase methane/ethane mixture at 280 K, 20 bar splits into vapour + liquid.
//! assert!(res.vapor_fraction >= 0.0 && res.vapor_fraction <= 1.0);
//! ```

#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

pub mod acceleration;
pub mod batch;
pub mod error;
pub mod initialization;
pub mod lle;
pub mod phase_volume;
pub mod pt;
pub mod rachford_rice;
pub mod variants;

pub use acceleration::AccelerationMemory;
pub use batch::flash_pt_batch;
#[cfg(feature = "std")]
pub use batch::flash_pt_batch_parallel;
pub use error::FlashError;
pub use initialization::wilson_k_values;
pub use lle::{lle_isoactivity, LleResult};
pub use phase_volume::{phase_volume, Phase};
pub use pt::{flash_pt, FlashCalculator, FlashResult};
pub use rachford_rice::{rachford_rice, RachfordRiceResult};
pub use variants::{flash_ph, flash_pu, flash_pv, flash_ts, flash_tv};
