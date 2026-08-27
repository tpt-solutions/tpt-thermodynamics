//! Two-phase K-values from the fugacity-coefficient equality `K_i = φ_i^L/φ_i^V`
//! and the liquid/vapor molar-volume gap used to locate the two-phase boundary.

use crate::kprovider::{KProvider, Phase};
use crate::nonphysical;
use alloc::vec;
use alloc::vec::Vec;
use tpt_thermo_core::convergence::{ConvergenceStatus, NumericalIssueReason};
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};

/// Whether the requested point is the bubble (liquid feed) or dew (vapor feed)
/// condition; controls which phase composition is iterated to consistency.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// Liquid feed; the vapor composition is iterated.
    Bubble,
    /// Vapor feed; the liquid composition is iterated.
    Dew,
}

/// Result of iterating a phase to fugacity consistency at a two-phase point.
pub struct Equilibrium {
    /// Equilibrium K-values `K_i = y_i / x_i` (bubble) or `x_i / y_i` (dew).
    pub k: Vec<f64>,
    /// The iterated other-phase composition (`y` for bubble, `x` for dew).
    pub other: Vec<f64>,
}

/// Liquid/vapor molar-volume gap at `(T, P, z)`.
///
/// Returns the positive gap `v_vapor − v_liquid` when both phases exist (two
/// phases), or `-1.0` when the mixture is single-phase. The gap collapses to
/// zero at the two-phase boundary, giving a sign-changing (negative↔positive)
/// scalar that Brent's method roots to locate bubble/dew points.
pub fn phase_gap(
    eos: &dyn KProvider,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> f64 {
    let vl: MolarVolume = match eos.phase_volume(t, p, z, Phase::Liquid) {
        Ok(v) => v,
        Err(_) => return -1.0,
    };
    let vv: MolarVolume = match eos.phase_volume(t, p, z, Phase::Vapor) {
        Ok(v) => v,
        Err(_) => return -1.0,
    };
    let gap = vv.value - vl.value;
    if gap > 1e-9 {
        gap
    } else {
        -1.0
    }
}

/// Iterate the equilibrium (other-phase) composition to fugacity consistency at
/// a two-phase `(T, P, feed)` and return the converged K-values and composition.
///
/// For [`Kind::Bubble`] the feed is the liquid composition `x`; `y` is iterated
/// via `y_i = x_i K_i / Σ x_j K_j`. For [`Kind::Dew`] the feed is the vapor `y`;
/// `x` is iterated via `x_i = y_i / K_i / Σ y_j / K_j`.
pub fn equilibrate(
    eos: &dyn KProvider,
    t: Temperature,
    p: Pressure,
    feed: &[f64],
    kind: Kind,
) -> Result<Equilibrium, ThermoError> {
    let n = feed.len();
    if n == 0 {
        return Err(ThermoError::InvalidInput("empty composition"));
    }
    let mut liquid = feed.to_vec();
    let mut vapor = feed.to_vec();
    let mut k = vec![1.0; n];
    let tol = 1e-10;

    for _ in 0..100 {
        let v_l = eos.phase_volume(t, p, &liquid, Phase::Liquid)?;
        let v_v = eos.phase_volume(t, p, &vapor, Phase::Vapor)?;
        for i in 0..n {
            let phi_l = (eos.ln_fugacity_coefficient(t, v_l, &liquid, i)?).exp();
            let phi_v = (eos.ln_fugacity_coefficient(t, v_v, &vapor, i)?).exp();
            if phi_v <= 0.0 {
                return Err(ThermoError::Numerical(ConvergenceStatus::NumericalIssue(
                    NumericalIssueReason::NonPhysical,
                )));
            }
            k[i] = phi_l / phi_v;
        }

        match kind {
            Kind::Bubble => {
                let mut nv = vec![0.0; n];
                let mut s = 0.0_f64;
                for i in 0..n {
                    nv[i] = feed[i] * k[i];
                    s += nv[i];
                }
                if s <= 0.0 {
                    return Err(nonphysical());
                }
                for item in nv.iter_mut() {
                    *item /= s;
                }
                let mut d = 0.0_f64;
                for i in 0..n {
                    d = d.max((nv[i] - vapor[i]).abs());
                }
                vapor = nv;
                if d < tol {
                    return Ok(Equilibrium { k, other: vapor });
                }
            }
            Kind::Dew => {
                let mut nl = vec![0.0; n];
                let mut s = 0.0_f64;
                for i in 0..n {
                    nl[i] = feed[i] / k[i];
                    s += nl[i];
                }
                if s <= 0.0 {
                    return Err(nonphysical());
                }
                for item in nl.iter_mut() {
                    *item /= s;
                }
                let mut d = 0.0_f64;
                for i in 0..n {
                    d = d.max((nl[i] - liquid[i]).abs());
                }
                liquid = nl;
                if d < tol {
                    return Ok(Equilibrium { k, other: liquid });
                }
            }
        }
    }
    Err(ThermoError::Numerical(ConvergenceStatus::NotConverged))
}
