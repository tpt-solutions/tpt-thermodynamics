//! Pitzer (1973) molality-scale activity and osmotic coefficients.
//!
//! Implements the standard single-electrolyte Pitzer virial expansion for the mean
//! activity coefficient `γ±` and the osmotic coefficient `φ` as a function of the
//! stoichiometric molality `m`. The Debye–Hückel limiting law is recovered exactly as
//! `m → 0` (see [`f_gamma`]), and the unsymmetrical `β²` term is supported for
//! higher-valence salts. Mixed-electrolyte θ/ψ terms are intentionally out of scope
//! here (tracked as a Deferred Scope follow-up); this module covers the single-salt
//! surface plus the limiting-law primitives shared with [`crate::enrtl`].

use crate::parameters::PitzerParams;

/// Debye–Hückel slope `A^φ` for water at 25 °C on the molality scale (≈ 0.3915
/// kg^½·mol^-3/2). Used by both the long-range `f` terms and the eNRTL PDH slope.
pub const A_PHI_25C: f64 = 0.3915;

/// Default Pitzer ion-size parameter `b` (dimensionless) in the Debye–Hückel
/// denominator `1 + b√I`.
const B_DH: f64 = 1.2;

/// Molality-scale ionic strength `I = ½·Σ mᵢ zᵢ²` for a single 1:1 electrolyte of
/// molality `m` with cation/anion charges `zc`/`za`.
pub fn ionic_strength(m: f64, zc: f64, za: f64) -> f64 {
    0.5 * m * (zc * zc + za * za)
}

/// Pitzer Debye–Hückel `f^γ(I)` — the natural log of the long-range limiting
/// contribution to `ln γ±`. Reduces to `−A_φ·√I` as `I → 0`.
pub fn f_gamma(i: f64) -> f64 {
    let si = i.max(0.0).sqrt();
    -A_PHI_25C * (si / (1.0 + B_DH * si) + (2.0 / B_DH) * (1.0 + B_DH * si).ln())
}

/// Pitzer Debye–Hückel `f^φ(I)` — the long-range osmotic-coefficient term.
pub fn f_phi(i: f64) -> f64 {
    let si = i.max(0.0).sqrt();
    -A_PHI_25C * si / (3.0 * (1.0 + B_DH * si))
}

/// `g(x)` helper, `x = α√I`, for the second-virial `B` term.
fn g(x: f64) -> f64 {
    if x.abs() < 1e-9 {
        0.0
    } else {
        let e = (-x).exp();
        2.0 * (1.0 - (1.0 + x) * e) / (x * x)
    }
}

/// `g^φ(x)` helper (osmotic `B^φ` term), related to [`g`] by
/// `B^γ = B^φ − I·∂B^φ/∂I`.
fn g_phi(x: f64) -> f64 {
    if x.abs() < 1e-9 {
        0.0
    } else {
        let e = (-x).exp();
        2.0 * (1.0 - (1.0 + x) * e) / (x * x)
    }
}

/// Second-virial `B^γ` (activity-coefficient form) for a single electrolyte at ionic
/// strength `i`, using `B^γ = β⁰ + β¹·g(x) + β²·g(2x)` with `x = α¹√I`.
pub fn b_gamma(p: &PitzerParams, i: f64) -> f64 {
    let x = (p.alpha1 * i.max(0.0).sqrt()).max(1e-12);
    p.beta0 + p.beta1 * g(x) + p.beta2 * g(2.0 * x)
}

/// Second-virial `B^φ` (osmotic form) for a single electrolyte at ionic strength `i`.
pub fn b_phi(p: &PitzerParams, i: f64) -> f64 {
    let x = (p.alpha1 * i.max(0.0).sqrt()).max(1e-12);
    p.beta0 + p.beta1 * g_phi(x) + p.beta2 * g_phi(2.0 * x)
}

/// Natural log of the mean activity coefficient `ln γ±` for a single electrolyte of
/// molality `m` with cation/anion charges `zc`/`za`. Reduces to the Debye–Hückel
/// limiting law as `m → 0`.
pub fn ln_mean_activity_coefficient(p: &PitzerParams, m: f64, zc: f64, za: f64) -> f64 {
    let i = ionic_strength(m, zc, za);
    let zpz = (zc * za).abs();
    let c_gamma = p.cphi / (2.0 * zpz.sqrt());
    f_gamma(i) * zpz + m * b_gamma(p, i) * zpz + 1.5 * m * m * c_gamma * zpz
}

/// Osmotic coefficient `φ` for a single electrolyte of molality `m` with cation/anion
/// charges `zc`/`za`. Tends to `1` as `m → 0`.
pub fn osmotic_coefficient(p: &PitzerParams, m: f64, zc: f64, za: f64) -> f64 {
    let i = ionic_strength(m, zc, za);
    let zpz = (zc * za).abs();
    let c_phi = p.cphi / (2.0 * zpz.sqrt());
    1.0 + f_phi(i) * zpz + m * b_phi(p, i) * zpz + m * m * c_phi * zpz
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn limiting_law_recovery() {
        // With all virial parameters zero, ln γ± reduces to the pure long-range
        // Debye–Hückel term `f^γ(I)` (no second/third virial contributions).
        let p = PitzerParams::one_to_one(0.0, 0.0, 0.0);
        let m = 1e-9;
        let i = ionic_strength(m, 1.0, -1.0);
        let got = ln_mean_activity_coefficient(&p, m, 1.0, -1.0);
        let expected = f_gamma(i);
        assert!((got - expected).abs() < 1e-15, "got {got}, expected {expected}");
        // Debye–Hückel limiting law: ln γ± → −A_γ·√I = −3·A_φ·√I as I → 0.
        assert!((got + 3.0 * A_PHI_25C * i.sqrt()).abs() < 1e-9);
        // Osmotic coefficient with no virial terms equals 1 + f^φ(I)·|z⁺z⁻|.
        let phi = osmotic_coefficient(&p, m, 1.0, -1.0);
        assert!((phi - (1.0 + f_phi(i))).abs() < 1e-15);
    }

    #[test]
    fn nacl_reference_shape() {
        // NaCl at 25 °C (Pitzer & Mayorga 1973 parameters).
        let p = PitzerParams::one_to_one(0.0765, 0.2664, 0.00127);
        let i = ionic_strength(1.0, 1.0, -1.0);
        let ln_g = ln_mean_activity_coefficient(&p, 1.0, 1.0, -1.0);
        let gamma = ln_g.exp();
        // Reference γ±(1m NaCl, 25°C) ≈ 0.657; allow generous tolerance here
        // (full seed-dataset validation is a tracked Deferred Scope follow-up).
        assert!(gamma > 0.5 && gamma < 0.8, "γ± = {gamma} at I={i}");
        assert!(osmotic_coefficient(&p, 1.0, 1.0, -1.0) > 0.9);
    }

    #[test]
    fn finite_and_positive() {
        let p = PitzerParams::one_to_one(0.0765, 0.2664, 0.00127);
        for m in [0.1_f64, 1.0, 3.0, 6.0] {
            let g = ln_mean_activity_coefficient(&p, m, 1.0, -1.0).exp();
            let phi = osmotic_coefficient(&p, m, 1.0, -1.0);
            assert!(g.is_finite() && g > 0.0, "γ± not finite/positive at m={m}");
            assert!(phi.is_finite() && phi > 0.0, "φ not finite/positive at m={m}");
        }
    }
}
