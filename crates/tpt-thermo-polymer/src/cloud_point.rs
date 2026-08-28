//! Cloud-point (UCST/LCST) calculation for binary polymer solutions.
//!
//! Uses the Flory-Huggins free energy with a temperature-dependent interaction
//! parameter. The binodal is solved from the equal-chemical-potential conditions,
//! and the critical point (which coincides with the cloud-point maximum/minimum)
//! is located analytically from the spinodal.

use tpt_thermo_core::ThermoError;

/// Temperature dependence of the binary `χ` parameter, `χ(T) = a + b/T`.
#[derive(Debug, Clone, Copy)]
pub enum ChiTemperature {
    /// Constant `χ`.
    Constant(f64),
    /// `χ(T) = a + b/T`. `b > 0` ⇒ `χ` decreases with `T` (UCST behaviour);
    /// `b < 0` ⇒ `χ` increases with `T` (LCST behaviour).
    LinearDecreasing { a: f64, b: f64 },
}

impl ChiTemperature {
    /// Value of `χ` at temperature `t` (K).
    pub fn at(&self, t: f64) -> f64 {
        match self {
            ChiTemperature::Constant(c) => *c,
            ChiTemperature::LinearDecreasing { a, b } => a + b / t,
        }
    }

    /// `dχ/dT` (per K).
    pub fn dchi_dt(&self, t: f64) -> f64 {
        match self {
            ChiTemperature::Constant(_) => 0.0,
            ChiTemperature::LinearDecreasing { b, .. } => -b / (t * t),
        }
    }
}

/// Whether the cloud-point is an upper or lower critical solution temperature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CloudPointKind {
    /// Miscible above `T_c`, demixes on cooling.
    Ucst,
    /// Miscible below `T_c`, demixes on heating.
    Lcst,
}

/// Critical (cloud-point) result for a binary mixture.
#[derive(Debug, Clone)]
pub struct CloudPointResult {
    /// Critical temperature (K). `None` for a temperature-independent model.
    pub temperature: Option<f64>,
    /// Critical polymer (component 2) volume fraction.
    pub critical_volume_fraction: f64,
    /// Critical `χ`.
    pub critical_chi: f64,
    /// UCST or LCST.
    pub kind: Option<CloudPointKind>,
}

/// Critical point of a binary Flory-Huggins mixture.
///
/// The critical volume fraction is `φ_c = 1/(1 + √(r₂/r₁))` and the critical
/// `χ` is the spinodal value there; the critical temperature follows from the
/// `χ(T)` model. Returns `Err` if the model is not temperature-dependent.
pub fn critical_point(
    r1: f64,
    r2: f64,
    chi: &ChiTemperature,
) -> Result<CloudPointResult, ThermoError> {
    if r1 <= 0.0 || r2 <= 0.0 {
        return Err(ThermoError::InvalidInput("segment counts must be positive"));
    }
    let phi_c = 1.0 / (1.0 + (r2 / r1).sqrt());
    let chi_c = 0.5 * (1.0 / (r1 * phi_c) + 1.0 / (r2 * (1.0 - phi_c)));
    let (temperature, kind) = match chi {
        ChiTemperature::Constant(_) => (None, None),
        ChiTemperature::LinearDecreasing { a, b } => {
            let t_c = (chi_c - a).recip() * *b;
            if !t_c.is_finite() || t_c <= 0.0 {
                return Err(ThermoError::Unsupported(
                    "χ(T) does not cross χ_c at a positive temperature",
                ));
            }
            let kind = if *b > 0.0 {
                CloudPointKind::Ucst
            } else {
                CloudPointKind::Lcst
            };
            (Some(t_c), Some(kind))
        }
    };
    Ok(CloudPointResult {
        temperature,
        critical_volume_fraction: phi_c,
        critical_chi: chi_c,
        kind,
    })
}

/// Residual (per-segment) chemical potential of component 1 (solvent, `r1`) and
/// component 2 (polymer, `r2`) at polymer volume fraction `phi`.
#[inline]
fn mu1(phi: f64, r1: f64, chi: f64) -> f64 {
    (1.0 - phi).ln() + (1.0 - 1.0 / r1) * phi + chi * phi * phi
}
#[inline]
fn mu2(phi: f64, r2: f64, chi: f64) -> f64 {
    phi.ln() + (1.0 - 1.0 / r2) * (1.0 - phi) + chi * (1.0 - phi) * (1.0 - phi)
}

/// Solve the binodal (coexistence) volume fractions `(φ_dilute, φ_concentrated)`
/// for a binary mixture at the given `χ`. Returns `None` when no two-phase region
/// exists (single phase).
pub fn binodal(_t: f64, r1: f64, r2: f64, chi: f64) -> Option<(f64, f64)> {
    if chi <= 0.0 {
        return None;
    }
    let phi_c = 1.0 / (1.0 + (r2 / r1).sqrt());
    let chi_c = 0.5 * (1.0 / (r1 * phi_c) + 1.0 / (r2 * (1.0 - phi_c)));
    if chi <= chi_c {
        return None;
    }
    // The equal-chemical-potential system `μ₁(pa)=μ₁(pb)`, `μ₂(pa)=μ₂(pb)` has the
    // trivial root `pa = pb` for *every* `χ`, which is a strong attractor for Newton
    // and would collapse the continuation at the (zero-width) critical point. We exclude
    // it by construction: the dilute phase always sits below `φ_c` and the concentrated
    // phase above `φ_c`, so `pa < φ_c < pb` is enforced by clamping at every step.
    let lo_lo = 1e-9;
    let lo_hi = phi_c - 1e-9;
    let hi_lo = phi_c + 1e-9;
    let hi_hi = 1.0 - 1e-9;
    // March `χ` *down* from a widened value to the target `χ`: at the top the binodal
    // is well separated and the continuation stays on the non-trivial branch.
    let chi_max = chi + chi.max(2.0 * (chi - chi_c));
    let steps = 120;
    let mut pa = phi_c * 0.5;
    let mut pb = (phi_c + 1.0) * 0.5;
    for s in 1..=steps {
        let chi_s = chi_max + (chi - chi_max) * (s as f64 / steps as f64);
        for _ in 0..80 {
            let f1 = mu1(pa, r1, chi_s) - mu1(pb, r1, chi_s);
            let f2 = mu2(pa, r2, chi_s) - mu2(pb, r2, chi_s);
            let dmu1_dpa = -1.0 / (1.0 - pa) + (1.0 - 1.0 / r1) + 2.0 * chi_s * pa;
            let dmu1_dpb = 1.0 / (1.0 - pb) - (1.0 - 1.0 / r1) - 2.0 * chi_s * pb;
            let dmu2_dpa = 1.0 / pa - (1.0 - 1.0 / r2) - 2.0 * chi_s * (1.0 - pa);
            let dmu2_dpb = -1.0 / pb + (1.0 - 1.0 / r2) + 2.0 * chi_s * (1.0 - pb);
            let det = dmu1_dpa * dmu2_dpb - dmu1_dpb * dmu2_dpa;
            if det.abs() < 1e-14 {
                break;
            }
            let dpa = (f1 * dmu2_dpb - f2 * dmu1_dpb) / det;
            let dpb = (dmu1_dpa * f2 - dmu2_dpa * f1) / det;
            // Damp the step to stay on the binodal branch.
            let dpa = dpa.clamp(-0.1, 0.1);
            let dpb = dpb.clamp(-0.1, 0.1);
            pa = (pa - dpa).clamp(lo_lo, lo_hi);
            pb = (pb - dpb).clamp(hi_lo, hi_hi);
            if (dpa.abs() + dpb.abs()) < 1e-12 {
                break;
            }
        }
    }
    let f1 = mu1(pa, r1, chi) - mu1(pb, r1, chi);
    let f2 = mu2(pa, r2, chi) - mu2(pb, r2, chi);
    if (f1 * f1 + f2 * f2).sqrt() > 1e-6 {
        return None;
    }
    let (lo, hi) = if pa < pb { (pa, pb) } else { (pb, pa) };
    Some((lo, hi))
}

/// Sample the cloud-point curve `(T, φ_dilute)` for temperatures in
/// `[t_min, t_max]`. For each `T` the binodal is solved; points where no
/// two-phase region exists (single phase) are skipped.
pub fn cloud_point_curve(
    r1: f64,
    r2: f64,
    chi: &ChiTemperature,
    t_min: f64,
    t_max: f64,
    n: usize,
) -> alloc::vec::Vec<(f64, f64)> {
    let mut out = alloc::vec::Vec::with_capacity(n);
    if n == 0 {
        return out;
    }
    for k in 0..=n {
        let t = t_min + (t_max - t_min) * (k as f64) / (n as f64);
        let c = chi.at(t);
        if let Some((lo, _hi)) = binodal(t, r1, r2, c) {
            out.push((t, lo));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn critical_point_matches_analytic() {
        // r1 = 1, r2 = 1000 → φ_c = 1/(1+√1000) ≈ 0.03065.
        let phi_c = 1.0 / (1.0 + 1000.0_f64.sqrt());
        let chi_c = 0.5 * (1.0 / (1.0 * phi_c) + 1.0 / (1000.0 * (1.0 - phi_c)));
        let cp = critical_point(
            1.0,
            1000.0,
            &ChiTemperature::LinearDecreasing {
                a: 0.0,
                b: chi_c * 400.0,
            },
        )
        .unwrap();
        assert!((cp.critical_volume_fraction - phi_c).abs() < 1e-6);
        assert!((cp.critical_chi - chi_c).abs() < 1e-6);
        assert!((cp.temperature.unwrap() - 400.0).abs() < 1e-3);
        assert_eq!(cp.kind, Some(CloudPointKind::Ucst));
    }

    #[test]
    fn binodal_exists_below_critical_and_respects_spinodal() {
        let r1 = 1.0_f64;
        let r2 = 1000.0_f64;
        let phi_c = 1.0_f64 / (1.0_f64 + f64::sqrt(r2 / r1));
        let chi_c = 0.5 * (1.0 / (r1 * phi_c) + 1.0 / (r2 * (1.0 - phi_c)));
        let b = chi_c * 400.0; // T_c = 400 K
        let chi = ChiTemperature::LinearDecreasing { a: 0.0, b };
        // At T = 300 < T_c, χ > χ_c → demixed.
        let c = chi.at(300.0);
        let (lo, hi) = binodal(300.0, r1, r2, c).unwrap();
        assert!(lo < phi_c && hi > phi_c, "binodal must straddle φ_c");
        // At T = 500 > T_c, χ < χ_c → single phase.
        assert!(binodal(500.0, r1, r2, chi.at(500.0)).is_none());
    }

    #[test]
    fn lcast_classified_for_negative_b() {
        let cp = critical_point(
            1.0,
            1000.0,
            &ChiTemperature::LinearDecreasing {
                a: 35.0,
                b: -16.0 * 400.0,
            },
        )
        .unwrap();
        assert_eq!(cp.kind, Some(CloudPointKind::Lcst));
    }
}
