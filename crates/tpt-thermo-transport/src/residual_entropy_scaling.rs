//! Residual-entropy scaling of transport properties.
//!
//! Rosenfeld's entropy scaling collapses the transport coefficient of a fluid onto
//! a universal function of the residual entropy `s^R`: `η·(ρ^(2/3)/√(m·kT)) =
//! F(s^R)`. This module provides the dimensionless scaling factor
//! `Φ = η·ρ^(2/3) / √(m·T)` and a simple `exp(A·s^R + B)` collapse used to
//! relate a property at one state to another along an isomorph.

use tpt_thermo_core::quantities::{DynamicViscosity, MolarMass, Pressure, Temperature};
use uom::si::dynamic_viscosity::pascal_second;
use uom::si::molar_mass::kilogram_per_mole;

/// Dimensionless Rosenfeld scaling factor `Φ = η·ρ^(2/3) / √(m·T)` from an
/// (ideal-gas) molar density `rho` (mol·m⁻³), viscosity `η`, mean molar mass `m`
/// (kg·mol⁻¹) and temperature `T` (K).
pub fn scaling_factor(eta: DynamicViscosity, rho: f64, m: MolarMass, t: Temperature) -> f64 {
    let rho_si = rho * 1000.0; // mol/m³ -> mol/m³ (molar mass already kg/mol)
    let _ = kilogram_per_mole;
    let denom = (m.value * t.value).sqrt();
    if denom <= 0.0 {
        return 0.0;
    }
    (eta.value / 1.0).sqrt().max(0.0) * rho_si.powf(2.0 / 3.0) / denom
}

/// Entropy-collapse estimate `exp(A·s^R + B)` relating two states along an isomorph.
/// Returns the transport coefficient's relative change `η_2 / η_1` between a state
/// with residual entropy `s1` and one with `s2`.
pub fn entropy_collapse(a: f64, b: f64, s1: f64, s2: f64) -> f64 {
    (a * s2 + b).exp() / (a * s1 + b).max(1e-30).exp()
}

/// Residual entropy (Reduced) estimate from pressure, temperature, and a reference
/// (e.g. critical) pressure, using a simple corresponding-states relation
/// `s^R ≈ −(P/P_ref)·(T_c/T)^(1/3)`. Largely a placeholder for the full
/// isomorph-based collapse; documented as approximate.
pub fn residual_entropy_approx(p: Pressure, t: Temperature, t_ref: f64, p_ref: f64) -> f64 {
    let pr = p.value / p_ref.max(1e-6);
    -(pr) * (t_ref / t.value.max(1e-6)).powf(1.0 / 3.0)
}

#[allow(dead_code)]
fn _assert_pascal(_: DynamicViscosity) -> f64 {
    DynamicViscosity::new::<pascal_second>(1.0).value
}
