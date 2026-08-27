//! `tpt-thermo-eos-cubic` — cubic equations of state for `tpt-thermodynamics`.
//!
//! This crate implements the classical cubic EoS family on top of the
//! [`tpt_thermo_core::EquationOfState`] trait:
//!
//! * [`PengRobinson`] (1976) and [`SoaveRedlichKwong`] (1972), with selectable
//!   alpha functions ([`SoaveAlpha`], [`TwuAlpha`], [`MathiasCopemanAlpha`]) and
//!   mixing rules.
//! * [`VolumeTranslated`] — the Peneloux volume-translation correction layered
//!   on a cubic EoS, improving liquid densities.
//! * Mixing rules: [`VdwMixing`] (van der Waals one-fluid, with optional
//!   `k_ij(T)`), and the excess-Gibbs-coupled [`HuronVidal`] (MHV1/MHV2/PSRK)
//!   and [`WongSandler`] combiners, generic over the core's
//!   [`ExcessGibbsModel`] trait (implemented by `tpt-thermo-eos-activity`,
//!   Phase 5).
//! * The shared [`CubicEos`] engine and the Cardano [`cubic_real_roots`]
//!   solver with [`Phase`]-aware [`select_root`].
//!
//! # Example
//!
//! ```
//! use tpt_thermo_core::component::ComponentDatabase;
//! use tpt_thermo_eos_cubic::{PengRobinson, cubic_solver::Phase};
//! use tpt_thermo_data::SeedComponentDatabase;
//! use tpt_thermo_core::quantities::{Temperature, Pressure, MolarVolume};
//! use uom::si::{thermodynamic_temperature::kelvin, pressure::pascal, molar_volume::cubic_meter_per_mole};
//!
//! let db = SeedComponentDatabase::from_seed();
//! let eos = PengRobinson::from_database(&db).unwrap();
//! let t = Temperature::new::<kelvin>(300.0);
//! let p = Pressure::new::<pascal>(1.0e6);
//! let v = eos.solve_phase(t, p, &[1.0], Phase::Vapor).unwrap();
//! assert!(v.value > 0.0);
//! ```

#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

pub mod alpha;
pub mod critical;
pub mod cubic_solver;
pub mod engine;
pub mod mixing;
pub mod pr;
pub mod srk;
pub mod volume_translation;

pub use alpha::{AlphaFunction, MathiasCopemanAlpha, SoaveAlpha, TwuAlpha};
pub use critical::{critical_point, mechanical_stability, spinodal_roots};
pub use cubic_solver::{compressibility_roots, cubic_real_roots, select_root, CubicModel, Phase};
pub use engine::CubicEos;
pub use mixing::{CubicMixing, HuronVidal, HvVariant, VdwMixing, WongSandler};
pub use pr::PengRobinson;
pub use srk::SoaveRedlichKwong;
pub use volume_translation::{peneloux_c, VolumeTranslated};
