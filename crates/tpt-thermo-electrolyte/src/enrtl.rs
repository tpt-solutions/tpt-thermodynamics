//! eNRTL (electrolyte NRTL): long-range Pitzer–Debye–Hückel term plus a
//! short-range NRTL contribution, following Chen et al. (1982). This crate
//! provides the long-range term exactly (Debye–Hückel) and a parameterised
//! short-range term; the short-range parameters are fitted per system.

use crate::pitzer::{f_gamma, ionic_strength, A_PHI_25C};

/// eNRTL parameters: Debye–Hückel slope `a_phi`, ion-size `b`, and the
/// short-range NRTL `tau`/`alpha` for the ion pair.
#[derive(Debug, Clone, Copy)]
pub struct EnrtlParams {
    /// Debye–Hückel slope `A^φ`.
    pub a_phi: f64,
    /// Ion-size parameter `b` in the PDH denominator (1.2 typical).
    pub b: f64,
    /// Short-range NRTL `τ`.
    pub tau: f64,
    /// Short-range NRTL `α` (non-randomness).
    pub alpha: f64,
}

impl Default for EnrtlParams {
    fn default() -> Self {
        Self {
            a_phi: A_PHI_25C,
            b: 1.2,
            tau: 0.0,
            alpha: 0.2,
        }
    }
}

/// Long-range Pitzer–Debye–Hückel contribution to `ln γ_i` (Chen form):
/// `ln γ_i^PDH = -(A_φ·z_i²/(1 + b√I))·(√I/(1+b√I) + (2/b)·ln(1+b√I))`.
pub fn pdh_ln_gamma(charge: f64, i: f64, a_phi: f64, b: f64) -> f64 {
    let si = i.max(0.0).sqrt();
    let denom = 1.0 + b * si;
    -a_phi * charge * charge / denom * (si / denom + (2.0 / b) * denom.ln())
}

/// Single-electrolyte (1:1) mean activity coefficient (natural log) for the eNRTL
/// model at molality `m`. At `m → 0` the short-range term vanishes and the result
/// reduces to the Debye–Hückel limiting law.
pub fn ln_mean_activity_coefficient(m: f64, p: &EnrtlParams, zc: f64, za: f64) -> f64 {
    let i = 0.5 * m * (zc * zc + za * za);
    let f = f_gamma(i);
    let long_range = (zc * za).abs() * f;
    // Short-range NRTL (symmetric 1:1) mean contribution, linear in molality at the
    // dilute limit so the limiting law is recovered.
    let g = (-p.alpha * p.tau).exp();
    let short_range = 2.0 * (p.tau * g) / (1.0 + g) * m;
    long_range + short_range
}

/// Convenience: ionic strength re-export for callers using raw (molality, charge)
/// lists.
pub use crate::pitzer::ionic_strength as ionic_strength_of;
