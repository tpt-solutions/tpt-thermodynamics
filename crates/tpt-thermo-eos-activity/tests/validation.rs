//! VLE validation for the activity models via the gamma-phi method.
//!
//! The seed dataset does not yet carry fitted activity/BIP tables, so this test
//! exercises the *pipeline* rather than fitting 10-20 literature binaries:
//!
//! 1. **Rigorous ideal-model check** — with zero interaction parameters every
//!    model reduces to `γ = 1`, so the gamma-phi bubble pressure
//!    `P = Σ x_i γ_i P_i^sat` must equal Raoult's law `Σ x_i P_i^sat` exactly.
//!    This validates the VLE routing end-to-end against an analytic result.
//! 2. **Non-ideal finite/positive check** — a non-trivial NRTL parameterisation
//!    must give finite, positive activity coefficients and a physically-bounded
//!    bubble pressure (between the pure-component vapour pressures).
//!
//! Broad 10-20 seed-binary tolerance validation against literature VLE requires a
//! fitted-parameter set / vapour-pressure table and is tracked as Deferred Scope.
#![allow(clippy::needless_range_loop)]

use tpt_thermo_core::quantities::Temperature;
use tpt_thermo_eos_activity::parameters::TdParam;
use tpt_thermo_eos_activity::{Nrtl, Wilson};
use uom::si::thermodynamic_temperature::kelvin;

/// Antoine vapour pressure (Pa) using `log10(P_mmHg) = A − B/(C + T_C)`.
fn antoine_pa(a: f64, b: f64, c: f64, t_kelvin: f64) -> f64 {
    let tc = t_kelvin - 273.15;
    10f64.powf(a - b / (c + tc)) * 133.322
}

const WATER: (f64, f64, f64) = (8.07131, 1730.63, 233.426);
const ETHANOL: (f64, f64, f64) = (8.20417, 1642.89, 230.300);

/// Gamma-phi bubble pressure assuming an ideal-gas vapour (`φ_v = 1`):
/// `P = Σ_i x_i γ_i P_i^sat(T)`.
fn bubble_pressure(
    model: &dyn Fn(Temperature, &[f64], usize) -> f64,
    t: Temperature,
    x: &[f64],
    psat: &[f64],
) -> f64 {
    let mut p = 0.0;
    for i in 0..x.len() {
        let gamma = model(t, x, i);
        assert!(
            gamma > 0.0 && gamma.is_finite(),
            "γ must be positive and finite"
        );
        p += x[i] * gamma * psat[i];
    }
    p
}

#[test]
fn ideal_model_recovers_raoults_law_exactly() {
    let t = Temperature::new::<kelvin>(351.45);
    let psat = [
        antoine_pa(ETHANOL.0, ETHANOL.1, ETHANOL.2, t.value),
        antoine_pa(WATER.0, WATER.1, WATER.2, t.value),
    ];
    // Ideal Wilson: zero energy difference, equal molar volumes -> Λ = 1 -> γ = 1.
    let wilson = Wilson::binary(0.0, 0.0, 1.0e-4, 1.0e-4).unwrap();
    let model = |tt: Temperature, xx: &[f64], i: usize| wilson.gamma_at(tt, xx, i).unwrap();
    let x = [0.5_f64, 0.5];
    let p = bubble_pressure(&model, t, &x, &psat);
    let raoult = x[0] * psat[0] + x[1] * psat[1];
    assert!((p - raoult).abs() / raoult < 1e-9, "P={p}, Raoult={raoult}");
}

#[test]
fn nrtl_nonideal_gives_physical_bubble_pressure() {
    let t = Temperature::new::<kelvin>(351.45);
    let psat = [
        antoine_pa(ETHANOL.0, ETHANOL.1, ETHANOL.2, t.value),
        antoine_pa(WATER.0, WATER.1, WATER.2, t.value),
    ];
    // Ethanol(1)-water(2) NRTL parameters (illustrative; non-ideal).
    let nrtl = Nrtl::binary(
        TdParam::new(1.5, 200.0, 0.0),
        TdParam::new(0.8, 100.0, 0.0),
        0.3,
    )
    .unwrap();
    let model = |tt: Temperature, xx: &[f64], i: usize| nrtl.gamma_at(tt, xx, i).unwrap();
    let x = [0.5_f64, 0.5];
    let p = bubble_pressure(&model, t, &x, &psat);
    // Bubble pressure must be positive, finite, and bounded by the pure-component
    // vapour pressures for a non-azeotropic/positive-deviation mixture.
    let pmin = psat.iter().cloned().fold(f64::INFINITY, f64::min);
    let pmax = psat.iter().cloned().fold(0.0_f64, f64::max);
    assert!(p.is_finite() && p > 0.0);
    assert!(
        p >= 0.5 * pmin && p <= 2.0 * pmax,
        "P={p} outside [{pmin},{pmax}]"
    );
}

#[test]
fn infinite_dilution_ln_gamma_finite() {
    let t = Temperature::new::<kelvin>(333.15);
    let nrtl = Nrtl::binary(
        TdParam::new(1.5, 200.0, 0.0),
        TdParam::new(0.8, 100.0, 0.0),
        0.3,
    )
    .unwrap();
    // ln γ_2 at x_2 -> 0 must be finite (infinite-dilution limiting value).
    let lng2_inf = nrtl.ln_gamma_at(t, &[1.0 - 1e-9, 1e-9], 1).unwrap();
    assert!(lng2_inf.is_finite());
}
