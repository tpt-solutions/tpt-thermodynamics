//! Mixture rules for transport properties: Wilke, Mason–Saxena, Filippov, Darken.

pub use crate::conductivity::filippov_thermal_conductivity;
pub use crate::diffusivity::{darken_liquid_binary, vignes_liquid_binary};
pub use crate::viscosity::{mason_saxena_mixture_viscosity, wilke_mixture_viscosity};

/// Binary Filippov mixture rule for thermal conductivity (W·m⁻¹·K⁻¹).
pub fn filippov(x1: f64, lambda1: f64, lambda2: f64) -> f64 {
    filippov_thermal_conductivity(x1, lambda1, lambda2)
}

/// Binary Darken mixture rule for liquid interdiffusion (m²·s⁻¹).
pub fn darken(x1: f64, d1_star: f64, d2_star: f64) -> f64 {
    darken_liquid_binary(d1_star, d2_star, x1)
}
