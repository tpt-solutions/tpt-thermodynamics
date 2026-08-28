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

/// Brent's method for **minimisation** (golden-section search safeguarded by
/// inverse-parabolic interpolation; Brent 1973). Finds the minimum of a
/// unimodal `f` on the bracket `[a, b]`, returning `(x_min, f_min)`.
///
/// This is the companion to [`brent`] (root finding) and underpins the
/// parameter-estimation utilities (e.g. fitting a binary interaction
/// parameter to VLE data).
pub fn brent_minimize<F, Func>(mut f: Func, ax: F, bx: F, tol: F, max_iter: usize) -> (F, F)
where
    F: Scalar,
    Func: FnMut(F) -> F,
{
    let zero = F::zero();
    let two = F::one() + F::one();
    let half = F::from(0.5).unwrap_or(F::one());
    let gold = F::from(0.3819660112501051).unwrap_or(half);
    let tiny = F::from(1e-20).unwrap_or(F::zero());

    let mut a = if ax < bx { ax } else { bx };
    let mut b = if ax < bx { bx } else { ax };
    let mut v = a + gold * (b - a);
    let mut w = v;
    let mut x = v;
    let mut fx = f(x);
    let mut fv = fx;
    let mut fw = fx;
    let mut d = b - a;
    let mut e = d;

    for _ in 0..max_iter {
        let xm = half * (a + b);
        let tol1 = tol * x.abs() + tiny;
        let tol2 = two * tol1;
        if (x - xm).abs() <= tol2 - half * (b - a) {
            return (x, fx);
        }
        let (mut u, fu);
        if e.abs() > tol1 {
            // Inverse-parabolic interpolation through (v, fv), (w, fw), (x, fx).
            let r = (x - w) * (fx - fv);
            let q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            let mut qq = two * (q - r);
            if qq > zero {
                p = -p;
            }
            qq = qq.abs();
            let etemp = e;
            e = d;
            let par_ok = p.abs() < half * qq * etemp.abs() && p > qq * (a - x) && p < qq * (b - x);
            if par_ok {
                d = p / qq;
                u = x + d;
                if u - a < tol2 || b - u < tol2 {
                    let du = if xm >= x { tol1 } else { -tol1 };
                    u = x + du;
                }
                fu = f(u);
            } else {
                if x >= xm {
                    e = a - x;
                } else {
                    e = b - x;
                }
                d = gold * e;
                u = x + d;
                fu = f(u);
            }
        } else {
            if x >= xm {
                e = a - x;
            } else {
                e = b - x;
            }
            d = gold * e;
            u = x + d;
            fu = f(u);
        }
        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if w == x || fu <= fw {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if v == x || v == w || fu <= fv {
                v = u;
                fv = fu;
            }
        }
    }
    (x, fx)
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
    fn brent_minimize_finds_quadratic() {
        // f(x) = (x - 2)^2, minimum at x = 2, f = 0.
        let (x, fx) = brent_minimize(|x: f64| (x - 2.0).powi(2), 0.0, 5.0, 1e-10, 200);
        assert!((x - 2.0).abs() < 1e-6, "x = {x}");
        assert!(fx.abs() < 1e-9);
    }

    #[test]
    fn brent_minimize_finds_cosine() {
        // f(x) = cos(x), minimum at x = π, f = -1.
        let (x, fx) = brent_minimize(|x: f64| x.cos(), 0.0, 6.0, 1e-10, 200);
        assert!((x - core::f64::consts::PI).abs() < 1e-5, "x = {x}");
        assert!((fx + 1.0).abs() < 1e-6);
    }

    #[test]
    fn brent_minimize_recovers_shifted_parabola() {
        let (x, _fx) = brent_minimize(|x: f64| (x + 3.0).powi(2) + 1.0, -10.0, 0.0, 1e-9, 200);
        assert!((x + 3.0).abs() < 1e-5);
    }

    #[test]
    fn unbracketed_returns_error() {
        let r = bisection(|x: f64| x * x + 1.0, 0.0, 2.0, 1e-12, 10);
        assert!(r.is_err());
    }
}
