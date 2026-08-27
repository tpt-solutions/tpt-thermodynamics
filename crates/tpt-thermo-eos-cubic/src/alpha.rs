//! Temperature-dependent attractive-parameter alpha functions
//! `α(T_r, ω)`, implementing the [`AlphaFunction`] trait. These scale the pure
//! critical `a_i` into `a_i(T) = a_i(T_c) · α(T_r, ω)`.

use alloc::boxed::Box;

/// A temperature-dependent alpha function `α(T_r, ω)` for a cubic EoS.
///
/// `tr` is the reduced temperature `T / T_c` and `ω` is the Pitzer acentric
/// factor. Implementations return a dimensionless `α ≥ 0`; the default
/// [`AlphaFunction::d_alpha_d_tr`] is a central-difference fallback but every
/// built-in variant below provides an analytic derivative.
pub trait AlphaFunction: Send + Sync {
    /// `α(T_r, ω)`.
    fn alpha(&self, tr: f64, omega: f64) -> f64;

    /// `∂α/∂T_r`, used for residual enthalpy/entropy departures.
    fn d_alpha_d_tr(&self, tr: f64, omega: f64) -> f64 {
        let h = tr.abs().max(1e-4) * 1e-4;
        (self.alpha(tr + h, omega) - self.alpha(tr - h, omega)) / (2.0 * h)
    }
}

/// Soave (1972) alpha: `α = [1 + κ (1 − √T_r)]²` with
/// `κ = 0.48 + 1.574 ω − 0.176 ω²`. This is the classic PR/SRK default.
#[derive(Debug, Clone, Copy, Default)]
pub struct SoaveAlpha;

impl SoaveAlpha {
    /// The Soave `κ` coefficient for a given acentric factor.
    pub fn kappa(omega: f64) -> f64 {
        0.48 + 1.574 * omega - 0.176 * omega * omega
    }
}

impl AlphaFunction for SoaveAlpha {
    fn alpha(&self, tr: f64, omega: f64) -> f64 {
        let kappa = Self::kappa(omega);
        let m = 1.0 + kappa * (1.0 - tr.sqrt());
        m * m
    }

    fn d_alpha_d_tr(&self, tr: f64, omega: f64) -> f64 {
        let kappa = Self::kappa(omega);
        let m = 1.0 + kappa * (1.0 - tr.sqrt());
        -kappa * m / tr.sqrt()
    }
}

/// Twu (1980) alpha function: `α = T_r^{N(M−1)} · exp(L (1 − T_r^{N M}))`.
///
/// Defaults to the widely-used single-parameter set `L = 0.65392`,
/// `M = 1.22600`, `N = 0.53087`; component-specific coefficients can be
/// supplied via [`TwuAlpha::new`].
#[derive(Debug, Clone, Copy)]
pub struct TwuAlpha {
    /// `L` coefficient.
    pub l: f64,
    /// `M` coefficient.
    pub m: f64,
    /// `N` coefficient.
    pub n: f64,
}

impl Default for TwuAlpha {
    fn default() -> Self {
        Self {
            l: 0.65392,
            m: 1.22600,
            n: 0.53087,
        }
    }
}

impl TwuAlpha {
    /// Build with explicit `L`, `M`, `N` coefficients.
    pub fn new(l: f64, m: f64, n: f64) -> Self {
        Self { l, m, n }
    }
}

impl AlphaFunction for TwuAlpha {
    fn alpha(&self, tr: f64, _omega: f64) -> f64 {
        let tr = tr.max(1e-6);
        let a = self.n * (self.m - 1.0);
        let b = self.n * self.m;
        tr.powf(a) * (self.l * (1.0 - tr.powf(b))).exp()
    }

    fn d_alpha_d_tr(&self, tr: f64, omega: f64) -> f64 {
        let tr = tr.max(1e-6);
        let a = self.n * (self.m - 1.0);
        let b = self.n * self.m;
        let alpha = self.alpha(tr, omega);
        alpha * (a / tr - self.l * b * tr.powf(b - 1.0))
    }
}

/// Mathias-Copeman (1983) alpha function:
/// `α = [1 + c1 (1−√T_r) + c2 (1−√T_r)² + c3 (1−√T_r)³]²`.
///
/// `c1` defaults to the Soave `κ`; `c2 = c3 = 0` reproduces Soave exactly.
/// Component-specific `c2`/`c3` (typically fitted) improve the vapor-pressure
/// representation.
#[derive(Debug, Clone, Copy)]
pub struct MathiasCopemanAlpha {
    /// Linear coefficient (defaults to Soave `κ`).
    pub c1: f64,
    /// Quadratic coefficient.
    pub c2: f64,
    /// Cubic coefficient.
    pub c3: f64,
}

impl MathiasCopemanAlpha {
    /// Build with explicit `c1`, `c2`, `c3`.
    pub fn new(c1: f64, c2: f64, c3: f64) -> Self {
        Self { c1, c2, c3 }
    }

    /// Build with `c1` set from the acentric factor (Soave `κ`) and `c2 = c3 = 0`.
    pub fn soave_like(omega: f64) -> Self {
        Self {
            c1: SoaveAlpha::kappa(omega),
            c2: 0.0,
            c3: 0.0,
        }
    }
}

impl Default for MathiasCopemanAlpha {
    fn default() -> Self {
        Self {
            c1: 0.48,
            c2: 0.0,
            c3: 0.0,
        }
    }
}

impl AlphaFunction for MathiasCopemanAlpha {
    fn alpha(&self, tr: f64, _omega: f64) -> f64 {
        let s = 1.0 - tr.sqrt();
        let inner = 1.0 + self.c1 * s + self.c2 * s * s + self.c3 * s * s * s;
        inner * inner
    }

    fn d_alpha_d_tr(&self, tr: f64, _omega: f64) -> f64 {
        let s = 1.0 - tr.sqrt();
        let inner = 1.0 + self.c1 * s + self.c2 * s * s + self.c3 * s * s * s;
        let dinner = (self.c1 + 2.0 * self.c2 * s + 3.0 * self.c3 * s * s) * (-0.5 / tr.sqrt());
        inner * dinner * 2.0
    }
}

/// Convenience constructor returning a boxed [`SoaveAlpha`] (the PR/SRK default).
pub fn soave() -> Box<dyn AlphaFunction> {
    Box::new(SoaveAlpha)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn soave_at_critical_is_one() {
        // At T_r = 1, α = 1 by construction for any ω.
        assert!((SoaveAlpha.alpha(1.0, 0.0) - 1.0).abs() < 1e-12);
        assert!((SoaveAlpha.alpha(1.0, 0.34) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn soave_monotonic_below_tc() {
        // α decreases as temperature falls below Tc for typical ω.
        let a1 = SoaveAlpha.alpha(0.6, 0.2);
        let a2 = SoaveAlpha.alpha(0.8, 0.2);
        assert!(a1 > a2);
    }

    #[test]
    fn soave_derivative_matches_numeric() {
        let tr = 0.7;
        let d = SoaveAlpha.d_alpha_d_tr(tr, 0.25);
        let h = 1e-5;
        let num = (SoaveAlpha.alpha(tr + h, 0.25) - SoaveAlpha.alpha(tr - h, 0.25)) / (2.0 * h);
        assert!((d - num).abs() / num.abs() < 1e-4);
    }

    #[test]
    fn twu_at_critical_is_one() {
        assert!((TwuAlpha::default().alpha(1.0, 0.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn mc_soave_like_equals_soave() {
        let tr = 0.7;
        let om = 0.2;
        let mc = MathiasCopemanAlpha::soave_like(om);
        assert!((mc.alpha(tr, om) - SoaveAlpha.alpha(tr, om)).abs() < 1e-12);
    }
}
