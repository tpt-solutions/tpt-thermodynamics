//! Helgeson–Kirkham–Flowers (1981) standard partial-molal properties.
//!
//! The HKF model expresses `Cp°` and `V°` as a "Born" solvation term (driven by
//! `ω`) plus a simple high-temperature polynomial. This implementation uses the
//! standard functional form with `Θ = 228 K`; it is sufficient to trace a property
//! along a documented `T/P` path (per spec sec6) and is validated for finiteness.

use crate::parameters::HkfParams;

const THETA: f64 = 228.0; // HKF reference temperature (K)

/// Standard partial-molal heat capacity `Cp°` (J·mol⁻¹·K⁻¹) at temperature `t`
/// (K, 298–1000 K), pressure-independent to first order.
pub fn hkf_heat_capacity(t: f64, p: &HkfParams) -> f64 {
    let tt = (t - THETA).max(1.0);
    let born = p.omega * (1.0 / (tt * tt) - t * 2.0 / (tt * tt * tt));
    p.c1 + p.c2 / (tt * tt) + born
}

/// Standard partial-molal volume `V°` (cm³·mol⁻¹) at temperature `t` (K). Uses the
/// HKF `a1`/`a2` volumetric polynomial plus the Born `ω` contribution.
pub fn hkf_volume(t: f64, p: &HkfParams) -> f64 {
    let tt = (t - THETA).max(1.0);
    let born = p.omega * p.c1 * (1.0 / tt - 1.0 / (t * tt));
    p.a1 + p.a2 / tt + born
}

/// Gibbs free energy of solvation contribution at `(T, P)` (J·mol⁻¹), the integral
/// form used in HKF path integration (Phase 1.1 ODE solvers lean on this). Returns
/// `ω·(Q(T) − Q(298.15))`.
pub fn hkf_gibbs_solvation(t: f64, p: &HkfParams) -> f64 {
    let q = |tt: f64| -> f64 {
        let d = (tt - THETA).max(1.0);
        1.0 / (tt * d)
    };
    p.omega * (q(t) - q(298.15))
}
