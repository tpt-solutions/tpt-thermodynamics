//! Rachford–Rice flash equation: solve for the vapor fraction `β` and the phase
//! compositions given K-values and the overall feed `z`.

use alloc::vec::Vec;
use tpt_thermo_core::brent;
use tpt_thermo_core::convergence::ConvergenceStatus;
use tpt_thermo_core::error::ThermoError;

/// Result of a Rachford–Rice solve.
#[derive(Debug, Clone, PartialEq)]
pub struct RachfordRiceResult {
    /// Vapor mole fraction `β ∈ [0, 1]`.
    pub beta: f64,
    /// Liquid-phase composition `x_i = z_i / (1 + β(K_i − 1))`.
    pub x: Vec<f64>,
    /// Vapor-phase composition `y_i = K_i x_i`.
    pub y: Vec<f64>,
}

/// Evaluate the Rachford–Rice function `g(β) = Σ z_i (K_i − 1)/(1 + β(K_i − 1))`.
fn g(beta: f64, k: &[f64], z: &[f64]) -> f64 {
    let mut s = 0.0_f64;
    for i in 0..z.len() {
        let denom = 1.0 + beta * (k[i] - 1.0);
        if denom.abs() > 1e-15 {
            s += z[i] * (k[i] - 1.0) / denom;
        }
    }
    s
}

/// Solve the Rachford–Rice equation for `β, x, y` given K-values `k` and feed `z`.
///
/// If all `K_i ≤ 1` the mixture is entirely liquid (`β = 0`); if all `K_i ≥ 1` it
/// is entirely vapor (`β = 1`). Otherwise `β` is the unique root of `g(β) = 0` on
/// `[0, 1]`, found with Brent's method.
pub fn rachford_rice(k: &[f64], z: &[f64]) -> Result<RachfordRiceResult, ThermoError> {
    if k.len() != z.len() || z.is_empty() {
        return Err(ThermoError::InvalidInput("mismatched K / z lengths"));
    }
    let sum_z: f64 = z.iter().sum();
    if (sum_z - 1.0).abs() > 1e-6 {
        return Err(ThermoError::InvalidInput("feed does not sum to 1"));
    }

    let n = z.len();
    let all_le = k.iter().all(|&ki| ki <= 1.0 + 1e-12);
    let all_ge = k.iter().all(|&ki| ki >= 1.0 - 1e-12);

    let beta = if all_le {
        0.0
    } else if all_ge {
        1.0
    } else {
        // g(0) > 0 and g(1) < 0 for a two-phase mixture.
        let g0 = g(0.0, k, z);
        let g1 = g(1.0, k, z);
        if !(g0 > 0.0 && g1 < 0.0) {
            // Degenerate: fall back to the boundary nearest the sign change.
            if g0 <= 0.0 {
                0.0
            } else {
                1.0
            }
        } else {
            brent(
                |b| g(b, k, z),
                0.0,
                1.0,
                1e-13,
                200,
            )
            .map_err(ThermoError::Numerical)?
        }
    };

    let mut x = alloc::vec![0.0_f64; n];
    let mut y = alloc::vec![0.0_f64; n];
    for i in 0..n {
        let denom = 1.0 + beta * (k[i] - 1.0);
        if denom.abs() > 1e-15 {
            x[i] = z[i] / denom;
            y[i] = k[i] * x[i];
        } else {
            x[i] = 0.0;
            y[i] = 0.0;
        }
    }
    // Normalise to protect against tiny drift.
    let sx: f64 = x.iter().sum();
    let sy: f64 = y.iter().sum();
    if sx > 0.0 {
        for xi in x.iter_mut() {
            *xi /= sx;
        }
    }
    if sy > 0.0 {
        for yi in y.iter_mut() {
            *yi /= sy;
        }
    }
    if beta.is_nan() {
        return Err(ThermoError::Numerical(ConvergenceStatus::NotConverged));
    }
    Ok(RachfordRiceResult { beta, x, y })
}
