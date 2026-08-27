//! Phase-selected molar-volume recovery for an arbitrary [`EquationOfState`].
//!
//! Given `(T, P, z)` we root-solve `P(v) = target` for the positive molar
//! volumes. In the two-phase region a cubic (or SAFT) EoS yields three real
//! roots: the smallest is the liquid volume and the largest is the vapor volume.
//! A single positive root means the mixture is single-phase, in which case both
//! phases collapse to that volume (so `K_i = 1`).

use tpt_thermo_core::bisection;
use tpt_thermo_core::convergence::ConvergenceStatus;
use tpt_thermo_core::eos::EquationOfState;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use uom::si::molar_volume::cubic_meter_per_mole;

/// Which phase a molar volume is requested for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    /// Liquid (smallest positive root of `P(v)=P`).
    Liquid,
    /// Vapor (largest positive root of `P(v)=P`).
    Vapor,
}

/// Residual `P(v) − target` for the root find.
fn f<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    p: Pressure,
    z: &[f64],
    v: f64,
) -> f64 {
    match eos.pressure(t, MolarVolume::new::<cubic_meter_per_mole>(v), z) {
        Ok(pv) => pv.value - p.value,
        Err(_) => f64::NAN,
    }
}

/// Recover the molar volume of `phase` at `(T, P, z)` from the EoS pressure.
///
/// The search brackets `[1e-9, 1.0]` m³·mol⁻¹ in log space; every sign change of
/// `P(v)−P` is refined with Brent's method. Returns the smallest (`Liquid`) or
/// largest (`Vapor`) positive root found, or [`ConvergenceStatus::NotConverged`]
/// if no physical root exists (e.g. far outside the model's domain).
pub fn phase_volume<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    p: Pressure,
    z: &[f64],
    phase: Phase,
) -> Result<MolarVolume, ThermoError> {
    let vmin: f64 = 1e-9;
    let vmax: f64 = 1.0;
    let log_lo = vmin.ln();
    let log_hi = vmax.ln();
    let nstep = 600;
    let dlog = (log_hi - log_lo) / nstep as f64;

    let mut roots: alloc::vec::Vec<f64> = alloc::vec::Vec::new();
    let mut v_prev = vmin;
    let mut f_prev = f(eos, t, p, z, v_prev);
    for s in 1..=nstep {
        let v = (log_lo + dlog * s as f64).exp();
        let fv = f(eos, t, p, z, v);
        if f_prev.is_finite() && fv.is_finite() && f_prev.signum() != fv.signum() {
            // Bisection on the bracket: it converges on the bracket *width*, which is
            // robust regardless of the (large) scale of `P(v) − P`, unlike a tight
            // absolute residual tolerance.
            if let Ok(r) = bisection(
                |vv| f(eos, t, p, z, vv),
                v_prev,
                v,
                1e-12,
                200,
            ) {
                if r > 0.0 {
                    roots.push(r);
                }
            }
        }
        v_prev = v;
        f_prev = fv;
    }

    if roots.is_empty() {
        return Err(ThermoError::Numerical(ConvergenceStatus::NotConverged));
    }
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let chosen = match phase {
        Phase::Liquid => roots[0],
        Phase::Vapor => *roots.last().unwrap(),
    };
    Ok(MolarVolume::new::<cubic_meter_per_mole>(chosen))
}
