//! Crystallization and melting-point depression (Flory equation).
//!
//! For a semi-crystalline polymer diluted by a diluent, the equilibrium melting point
//! `T_m` is depressed relative to the pure-polymer melting point `T_m⁰` by the Flory
//! relation:
//!
//! ```text
//! 1/T_m − 1/T_m⁰ = −(R/ΔH_f)·(V_u/V_1)·(φ_1 − χ·φ_1²)
//! ```
//!
//! where `φ_1` is the volume fraction of diluent, `V_u`/`V_1` the repeat-unit/diluent
//! molar-volume ratio, `ΔH_f` the heat of fusion per repeat unit, and `χ` the
//! Flory–Huggins interaction parameter.

/// Equilibrium melting temperature of a diluted polymer.
///
/// * `t_m0` — pure-polymer melting point (K).
/// * `delta_h_f` — heat of fusion per repeat unit (J·mol⁻¹).
/// * `v_u_over_v1` — ratio of repeat-unit to diluent molar volume.
/// * `phi_1` — diluent volume fraction.
/// * `chi` — Flory–Huggins interaction parameter.
///
/// Returns `T_m` in kelvin.
pub fn melting_point_depression(
    t_m0: f64,
    delta_h_f: f64,
    v_u_over_v1: f64,
    phi_1: f64,
    chi: f64,
) -> f64 {
    assert!(t_m0 > 0.0 && delta_h_f > 0.0, "non-physical inputs");
    let r = 8.314462618_f64;
    let term = (r / delta_h_f) * v_u_over_v1 * (phi_1 - chi * phi_1 * phi_1);
    let inv_t = 1.0 / t_m0 - term;
    1.0 / inv_t
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_polymer_no_depression() {
        // φ_1 = 0 → T_m = T_m⁰ exactly.
        let t = melting_point_depression(400.0, 8000.0, 0.5, 0.0, 0.0);
        assert!((t - 400.0).abs() < 1e-9);
    }

    #[test]
    fn depression_is_monotonic_in_phi() {
        let t_low = melting_point_depression(400.0, 8000.0, 0.5, 0.1, 0.0);
        let t_high = melting_point_depression(400.0, 8000.0, 0.5, 0.4, 0.0);
        assert!(t_low < 400.0 && t_high < t_low, "more diluent → lower T_m");
    }
}
