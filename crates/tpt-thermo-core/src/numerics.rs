//! One-dimensional root finders, generic over the
//! [`Scalar`](tpt_math_numeric::Scalar) trait.
//!
//! The build-out spec expected these in `tpt-math-numeric` (Phase 1), but the
//! published `tpt-math-numeric` 0.1.0 is a thin `num_traits` wrapper without
//! them (see `todo.md` Phase 1 open item). They are therefore provided here so
//! downstream `tpt-thermo-*` crates have a single, in-repo home for them.

use crate::convergence::{ConvergenceStatus, DivergenceReason, NumericalIssueReason};
use tpt_math_numeric::Scalar;

/// Default tolerance for the root finders below.
pub const ROOT_TOL: f64 = 1e-10;
/// Default iteration budget for the root finders below.
pub const ROOT_MAX_ITER: usize = 200;

/// Bisection root finder on a bracket `[a, b]` with `f(a)·f(b) ≤ 0`.
///
/// Returns the approximate root or a [`ConvergenceStatus`].
pub fn bisection<F, Func>(
    mut f: Func,
    mut a: F,
    mut b: F,
    tol: F,
    max_iter: usize,
) -> Result<F, ConvergenceStatus>
where
    F: Scalar,
    Func: FnMut(F) -> F,
{
    let two = F::one() + F::one();
    let mut fa = f(a);
    let fb = f(b);
    if fa * fb > F::zero() {
        return Err(ConvergenceStatus::NumericalIssue(
            NumericalIssueReason::OutOfDomain,
        ));
    }
    for _ in 0..max_iter {
        let m = (a + b) / two;
        let fm = f(m);
        if fm.abs() < tol || (b - a).abs() < tol {
            return Ok(m);
        }
        if fa * fm < F::zero() {
            b = m;
        } else {
            a = m;
            fa = fm;
        }
    }
    Err(ConvergenceStatus::Diverged(DivergenceReason::MaxIterations))
}

/// Newton-Raphson root finder given `f` and its derivative `df`.
pub fn newton<F, Func, Deriv>(
    mut f: Func,
    mut df: Deriv,
    mut x: F,
    tol: F,
    max_iter: usize,
) -> Result<F, ConvergenceStatus>
where
    F: Scalar,
    Func: FnMut(F) -> F,
    Deriv: FnMut(F) -> F,
{
    for _ in 0..max_iter {
        let fx = f(x);
        if fx.abs() < tol {
            return Ok(x);
        }
        let dfx = df(x);
        if dfx.abs() < F::from(1e-14).unwrap_or(F::zero()) {
            return Err(ConvergenceStatus::NumericalIssue(
                NumericalIssueReason::SingularJacobian,
            ));
        }
        let step = fx / dfx;
        x = x - step;
        if step.abs() < tol {
            return Ok(x);
        }
    }
    Err(ConvergenceStatus::Diverged(DivergenceReason::MaxIterations))
}

/// Brent's method: bisection safeguarded by inverse quadratic interpolation.
///
/// Combines the robustness of bisection with the speed of interpolation; the
/// spec names this as the explicit Newton fallback for flash/phase solvers.
pub fn brent<F, Func>(
    mut f: Func,
    mut a: F,
    mut b: F,
    tol: F,
    max_iter: usize,
) -> Result<F, ConvergenceStatus>
where
    F: Scalar,
    Func: FnMut(F) -> F,
{
    let eps = F::from(1e-14).unwrap_or(F::zero());
    let two = F::one() + F::one();
    let mut fa = f(a);
    let mut fb = f(b);
    if fa * fb > F::zero() {
        return Err(ConvergenceStatus::NumericalIssue(
            NumericalIssueReason::OutOfDomain,
        ));
    }
    // Ensure |f(a)| <= |f(b)| so `a` is the "better" endpoint.
    if fa.abs() < fb.abs() {
        core::mem::swap(&mut a, &mut b);
        core::mem::swap(&mut fa, &mut fb);
    }
    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;

    for _ in 0..max_iter {
        let tol_act = tol * b.abs() + eps;
        // Inverse quadratic interpolation when three distinct points are
        // available; otherwise secant. Result is `s`, with proposed step `d_new`.
        let (mut s, mut d_new) = if (fa - fb).abs() > eps && (fc - fb).abs() > eps {
            let num = (c - b) * (b - a) * (fc - fa) - (a - b) * (b - c) * (fa - fc);
            let den = (c - b) * (fa - fc) - (a - b) * (fc - fa);
            if den.abs() > eps {
                (num / den, c - b)
            } else {
                (a - fb / ((fb - fa) / (b - a)), b - a)
            }
        } else {
            (a - fb / ((fb - fa) / (b - a)), b - a)
        };

        // Accept the interpolation only if it sits inside a safe bisection step.
        let half = d.abs() / two;
        if s.abs() < two * tol_act || s.abs() >= two * half || (s - d).abs() * two >= d.abs() {
            s = (a + b) / two;
            d_new = b - a;
        }
        d = d_new;

        let fs = f(s);
        if fs.abs() < tol_act {
            return Ok(s);
        }
        if fa * fs < F::zero() {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            core::mem::swap(&mut a, &mut b);
            core::mem::swap(&mut fa, &mut fb);
        }
        c = a;
        fc = fa;
    }
    Err(ConvergenceStatus::Diverged(DivergenceReason::MaxIterations))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bisection_finds_root() {
        // f(x) = x^2 - 2, bracket [1, 2].
        let r = bisection(|x: f64| x * x - 2.0, 1.0, 2.0, 1e-12, 100).unwrap();
        assert!((r - core::f64::consts::SQRT_2).abs() < 1e-10);
    }

    #[test]
    fn newton_finds_root() {
        let r = newton(
            |x: f64| x - x.cbrt() - 1.0,
            |x: f64| 1.0 - 1.0 / (3.0 * x.powf(2.0 / 3.0)),
            2.0,
            1e-12,
            100,
        )
        .unwrap();
        assert!((r - 2.0_f64.cbrt() - 1.0).abs() < 1e-10 || (r - (r.cbrt() + 1.0)).abs() < 1e-9);
    }

    #[test]
    fn brent_finds_root() {
        // f(x) = cos(x) - x, bracket [0, 1].
        let r = brent(|x: f64| x.cos() - x, 0.0, 1.0, 1e-12, 100).unwrap();
        assert!((r.cos() - r).abs() < 1e-9);
    }

    #[test]
    fn unbracketed_returns_error() {
        let r = bisection(|x: f64| x * x + 1.0, 0.0, 2.0, 1e-12, 10);
        assert!(r.is_err());
    }
}
