//! Arc-length (pseudo-arclength) continuation for tracing equilibrium curves.
//!
//! Given an `n−1`-dimensional equilibrium residual `F(x) = 0` (`x` has length
//! `n`, the extra degree of freedom being the curve parameter), advance a point
//! along the curve by arc length `ds` using an Euler predictor and a Newton
//! corrector on the augmented `[F(x); t·(x − x₀) − ds] = 0` system.

use crate::linalg::solve_linear;
use alloc::vec;
use alloc::vec::Vec;

/// Advance one continuation step.
///
/// * `f` — the `n−1` equilibrium residuals (`f(x)` must return `x.len() − 1`
///   values; `x` has length `n`).
/// * `x` — current point (length `n`).
/// * `tangent` — current unit tangent (length `n`).
/// * `ds` — arc-length step.
///
/// Returns the new point and the new unit tangent, or `None` if the corrector
/// fails (e.g. a singular augmented Jacobian).
pub fn arc_length_step<F>(
    f: &F,
    x: &[f64],
    tangent: &[f64],
    ds: f64,
    tol: f64,
    max_iter: usize,
) -> Option<(Vec<f64>, Vec<f64>)>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let m = x.len();
    if m < 2 || tangent.len() != m {
        return None;
    }
    let predictor: Vec<f64> = x.iter().zip(tangent).map(|(xi, ti)| xi + ds * ti).collect();
    let mut y = predictor.clone();
    let t = tangent.to_vec();

    for _ in 0..max_iter {
        let fy = f(&y);
        let arc = t
            .iter()
            .zip(y.iter())
            .zip(x.iter())
            .map(|((ti, yi), xi)| ti * (yi - xi))
            .sum::<f64>()
            - ds;
        let resid_norm2: f64 = fy.iter().map(|v| v * v).sum::<f64>() + arc * arc;
        if resid_norm2 < tol * tol {
            break;
        }
        let jf = jacobian(f, &y);
        // Augmented Jacobian: top (m-1) rows = ∂F/∂x, last row = tangent.
        let mut a = vec![vec![0.0_f64; m]; m];
        for r in 0..m - 1 {
            for c in 0..m {
                a[r][c] = jf[r][c];
            }
        }
        for c in 0..m {
            a[m - 1][c] = t[c];
        }
        let b: Vec<f64> = fy.iter().chain(core::iter::once(&arc)).map(|v| -v).collect();
        let d = solve_linear(&a, &b)?;
        let mut step_norm = 0.0_f64;
        for i in 0..m {
            y[i] += d[i];
            step_norm += d[i] * d[i];
        }
        if step_norm < tol * tol {
            break;
        }
    }

    // New tangent: the unit direction of the actual move along the curve
    // (from the previous corrected point `x` to the new `y`), which the
    // arc-length constraint keeps at length ≈ `ds`.
    let mut tnew: Vec<f64> = y.iter().zip(x).map(|(yi, xi)| yi - xi).collect();
    let nrm = tnew.iter().map(|v| v * v).sum::<f64>().sqrt();
    if nrm > 1e-12 {
        for v in tnew.iter_mut() {
            *v /= nrm;
        }
    } else {
        tnew = t;
    }
    Some((y, tnew))
}

/// Numerical Jacobian of `f` (size `f.len()` × `x.len()`).
fn jacobian<F: Fn(&[f64]) -> Vec<f64>>(f: &F, x: &[f64]) -> Vec<Vec<f64>> {
    let fx = f(x);
    let rows = fx.len();
    let cols = x.len();
    let h = 1e-6_f64;
    let mut j = vec![vec![0.0_f64; cols]; rows];
    for c in 0..cols {
        let mut xh = x.to_vec();
        let orig = x[c];
        let hh = (orig.abs().max(1.0)) * h;
        xh[c] = orig + hh;
        let fph = f(&xh);
        for r in 0..rows {
            j[r][c] = (fph[r] - fx[r]) / hh;
        }
    }
    j
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Unit circle: residual is the single equation `x² + y² − 1 = 0`, state
    /// `(x, y)`, parameter implicit. Start at (1, 0) heading +y.
    #[test]
    fn traces_unit_circle() {
        let f = |x: &[f64]| vec![x[0] * x[0] + x[1] * x[1] - 1.0];
        let mut x = vec![1.0_f64, 0.0];
        let mut t = vec![0.0_f64, 1.0];
        for _ in 0..40 {
            let (nx, nt) = arc_length_step(&f, &x, &t, 0.1, 1e-10, 50).expect("step");
            x = nx;
            t = nt;
            let r2 = x[0] * x[0] + x[1] * x[1];
            assert!((r2 - 1.0).abs() < 1e-3, "point should stay on the unit circle");
        }
        // After 40 steps of 0.1 we should have gone most of the way around.
        assert!(x[0] < 0.0, "should have passed the (-1, 0) region");
    }
}
