//! `tpt-thermo-transport` — transport properties for `tpt-thermodynamics`.
//!
//! Low-pressure (dilute-gas) and liquid correlations:
//!
//! * [`viscosity`] — Chung et al. (1988) dilute-gas viscosity and the Lucas
//!   (1981) liquid-viscosity correlation, with Wilke/Mason–Saxena mixture rules.
//! * [`conductivity`] — Chung et al. (1988) gas thermal conductivity and a
//!   liquid corresponding-states estimate, with Filippov mixing.
//! * [`diffusivity`] — Fuller–Schettler–Giddings binary gas diffusion and Vignes
//!   liquid interdiffusion.
//! * [`mixing_rules`] — Wilke, Mason–Saxena, Filippov, and Darken mixture rules.
//! * [`residual_entropy_scaling`] — entropy-scaling collapse of transport data.
//!
//! All routines are unit-safe (returning `uom` quantities) and read critical
//! constants from a [`tpt_thermo_core::component::ComponentDatabase`].

#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

pub mod conductivity;
pub mod diffusivity;
pub mod mixing_rules;
pub mod parameters;
pub mod residual_entropy_scaling;
pub mod viscosity;

pub use parameters::{lj_params_for, LjParams};
