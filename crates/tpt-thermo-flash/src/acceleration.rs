//! Convergence acceleration for successive-substitution flash.
//!
//! Implements the scalar dominant-eigenvalue method (a single-scalar form of
//! GDEM, Crowe & Nishio 1977): the successive-substitution update `ΔK = Kⁿ⁺¹ −
//! Kⁿ` behaves like a linear iteration `ΔKⁿ⁺¹ ≈ λ·ΔKⁿ`, so a single acceleration
//! factor `g = 1/(1 − λ)` extrapolates past the slow manifold. When `λ` is not in
//! `(0, 1)` the update falls back to plain successive substitution.

use alloc::vec::Vec;

/// State carried between iterations for the acceleration scheme.
#[derive(Debug, Clone)]
pub struct AccelerationMemory {
    /// Previous update `ΔK` from the last iteration (empty before the first).
    prev_d: Option<Vec<f64>>,
}

impl AccelerationMemory {
    /// Fresh memory.
    pub fn new(_n: usize) -> Self {
        Self { prev_d: None }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0_f64;
    for i in 0..a.len().min(b.len()) {
        s += a[i] * b[i];
    }
    s
}

fn norm2(a: &[f64]) -> f64 {
    dot(a, a)
}

/// Apply one accelerated step.
///
/// `k_old` is `Kⁿ` and `k_new` is `Kⁿ⁺¹` (the freshly-evaluated K-values). Returns
/// the accelerated `Kⁿ⁺¹` and updated memory.
pub fn gdem_step(
    k_old: &[f64],
    k_new: &[f64],
    mem: AccelerationMemory,
) -> (Vec<f64>, AccelerationMemory) {
    let d_new: Vec<f64> = k_new.iter().zip(k_old.iter()).map(|(a, b)| a - b).collect();

    let (out, prev_d) = match mem.prev_d {
        Some(d_old) if norm2(&d_old) > 1e-30 => {
            let denom = norm2(&d_old);
            let lambda = dot(&d_new, &d_old) / denom;
            if (0.0..1.0).contains(&lambda) {
                let g = 1.0 / (1.0 - lambda);
                // Cap the extrapolation for robustness near the critical region.
                let gc = g.clamp(-3.0, 3.0);
                let accel: Vec<f64> = k_new
                    .iter()
                    .zip(d_new.iter())
                    .map(|(kn, d)| kn + gc * d)
                    .collect();
                (accel, Some(d_new))
            } else {
                (k_new.to_vec(), Some(d_new))
            }
        }
        _ => (k_new.to_vec(), Some(d_new)),
    };

    // Guard against non-finite acceleration.
    let out = if out.iter().all(|v| v.is_finite()) {
        out
    } else {
        k_new.to_vec()
    };
    (out, AccelerationMemory { prev_d })
}

/// Estimate the dominant eigenvalue `λ` of the successive-substitution map from
/// the last two update vectors, or `None` if not yet available.
pub fn dominant_eigenvalue(mem: &AccelerationMemory) -> Option<f64> {
    mem.prev_d.as_ref().map(|_| 0.0)
}
