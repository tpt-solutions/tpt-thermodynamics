//! Liquid–liquid (LLE) flash via isoactivity, driven by an
//! [`ExcessGibbsModel`] (e.g. NRTL/Wilson/UNIQUAC from `tpt-thermo-eos-activity`).

use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::mixing::ExcessGibbsModel;
use tpt_thermo_core::quantities::{Pressure, Temperature};

use crate::rachford_rice::rachford_rice;

/// Result of an LLE flash: compositions of the two co-existing liquid phases and
/// the phase fraction of phase II (`β`).
#[derive(Debug, Clone, PartialEq)]
pub struct LleResult {
    /// Composition of liquid phase I.
    pub phase_i: Vec<f64>,
    /// Composition of liquid phase II (the incipient phase).
    pub phase_ii: Vec<f64>,
    /// Mole fraction of phase II.
    pub beta: f64,
    /// Iterations performed.
    pub iterations: usize,
    /// Whether the iteration met tolerance.
    pub converged: bool,
}

/// LLE flash for a partially-miscible system using an activity model.
///
/// The two liquid phases are in equilibrium when their activities are equal:
/// `γ_i^I·x_i^I = γ_i^II·x_i^II`, so the pseudo-K-values are `K_i = γ_i^II/γ_i^I`.
/// Successive substitution iterates these against the Rachford–Rice mass balance.
pub fn lle_isoactivity(
    model: &dyn ExcessGibbsModel,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> Result<LleResult, ThermoError> {
    let n = z.len();
    if n < 2 {
        return Err(ThermoError::InvalidInput("LLE requires >= 2 components"));
    }
    let sum_z: f64 = z.iter().sum();
    if (sum_z - 1.0).abs() > 1e-6 {
        return Err(ThermoError::InvalidInput("feed does not sum to 1"));
    }

    // Break symmetry with two distinct initial phase compositions.
    let mut x_i = z.to_vec();
    let mut x_ii = z.to_vec();
    x_i[0] *= 1.1;
    x_ii[0] *= 0.9;
    let norm = |v: &mut [f64]| {
        let s: f64 = v.iter().sum();
        if s > 0.0 {
            for v in v.iter_mut() {
                *v /= s;
            }
        }
    };
    norm(&mut x_i);
    norm(&mut x_ii);

    let tol = 1e-9;
    let max_iter = 200;
    let mut converged = false;
    let mut beta = 0.5;

    for it in 0..max_iter {
        // K_i = gamma_ii / gamma_i
        let mut k = alloc::vec![1.0_f64; n];
        for (i, ki) in k.iter_mut().enumerate() {
            let gln_i = model.ln_gamma(t, p, &x_i, i)?;
            let gln_ii = model.ln_gamma(t, p, &x_ii, i)?;
            *ki = (gln_ii - gln_i).exp();
        }
        let rr = rachford_rice(&k, z)?;
        beta = rr.beta.clamp(0.0, 1.0);
        let new_xi = rr.x.clone();
        let new_xii = rr.y.clone();

        let d = relative_change(&x_i, &new_xi).max(relative_change(&x_ii, &new_xii));
        x_i = new_xi;
        x_ii = new_xii;
        if d < tol {
            converged = true;
            return Ok(LleResult {
                phase_i: x_i,
                phase_ii: x_ii,
                beta,
                iterations: it + 1,
                converged,
            });
        }
    }

    Ok(LleResult {
        phase_i: x_i,
        phase_ii: x_ii,
        beta,
        iterations: max_iter,
        converged,
    })
}

fn relative_change(a: &[f64], b: &[f64]) -> f64 {
    let mut m = 0.0_f64;
    for i in 0..a.len().min(b.len()) {
        let denom = a[i].abs().max(1e-12);
        m = m.max((a[i] - b[i]).abs() / denom);
    }
    m
}
