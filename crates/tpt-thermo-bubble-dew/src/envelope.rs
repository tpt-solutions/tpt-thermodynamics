//! Phase envelope: trace the bubble and dew curves of a fixed overall
//! composition `z` across a pressure sweep (a simple, robust continuation that
//! steps in pressure and re-solves each bubble/dew point).
//!
//! > The approved plan names Phase 8's arc-length continuation machinery for
//! > this trace. That crate lands later; the pressure-swept form here is the
//! > same physical curve and is structured so it can be replaced by an
//! > arc-length driver without changing callers.

use crate::equilibrium::Kind;
use crate::{BubbleDewSolver, BubblePoint, DewPoint};
use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use uom::si::thermodynamic_temperature::kelvin;

/// A single point on a phase envelope curve.
#[derive(Debug, Clone)]
pub struct EnvelopePoint {
    /// Pressure.
    pub pressure: Pressure,
    /// Temperature.
    pub temperature: Temperature,
    /// Liquid composition at this point.
    pub liquid: Vec<f64>,
    /// Vapor composition at this point.
    pub vapor: Vec<f64>,
}

/// A traced phase envelope: the bubble curve and the dew curve.
#[derive(Debug, Clone)]
pub struct Envelope {
    /// Bubble-curve points (liquid feed `z`), low-to-high pressure.
    pub bubble: Vec<EnvelopePoint>,
    /// Dew-curve points (vapor feed `z`), low-to-high pressure.
    pub dew: Vec<EnvelopePoint>,
}

fn point_from_bubble(bp: &BubblePoint) -> EnvelopePoint {
    EnvelopePoint {
        pressure: bp.pressure,
        temperature: bp.temperature,
        liquid: bp.liquid.clone(),
        vapor: bp.vapor.clone(),
    }
}

fn point_from_dew(dp: &DewPoint) -> EnvelopePoint {
    EnvelopePoint {
        pressure: dp.pressure,
        temperature: dp.temperature,
        liquid: dp.liquid.clone(),
        vapor: dp.vapor.clone(),
    }
}

/// Trace the bubble and dew curves of overall composition `z` over the given
/// pressure samples. Points that fail to converge (e.g. past the critical point)
/// are skipped, so callers should sort `pressures` and expect a possibly
/// truncated curve near the critical region.
///
/// The trace uses continuation: each pressure step seeds its root solver with
/// the previous converged temperature, which keeps the curve smooth and avoids
/// the spurious low-/high-temperature roots that a fresh full-range scan would
/// otherwise pick up at high pressure.
pub fn bubble_dew_envelope(
    solver: &BubbleDewSolver<'_>,
    z: &[f64],
    pressures: &[Pressure],
) -> Result<Envelope, ThermoError> {
    let mut bubble = Vec::new();
    let mut dew = Vec::new();
    let mut tb_guess: Option<f64> = None;
    let mut td_guess: Option<f64> = None;
    for &p in pressures {
        let bp = match tb_guess {
            Some(tg) => solver
                .boundary_temperature_guess(p, z, tg)
                .or_else(|_| solver.boundary_temperature(p, z))
                .ok()
                .map(|tb| {
                    let t = Temperature::new::<kelvin>(tb);
                    let eq = solver.equilibrium_at(t, p, z, Kind::Bubble)?;
                    Ok::<_, ThermoError>(BubblePoint {
                        temperature: t,
                        pressure: p,
                        liquid: z.to_vec(),
                        vapor: eq.other,
                        k_values: eq.k,
                        converged: true,
                    })
                }),
            None => solver.bubble_point_temperature(p, z).map(Ok).ok(),
        };
        if let Some(Ok(bp)) = bp {
            tb_guess = Some(bp.temperature.value);
            bubble.push(point_from_bubble(&bp));
        }
        let dp = match td_guess {
            Some(tg) => solver
                .boundary_temperature_dew_guess(p, z, tg)
                .or_else(|_| solver.boundary_temperature_dew(p, z))
                .ok()
                .map(|td| {
                    let t = Temperature::new::<kelvin>(td);
                    let eq = solver.equilibrium_at(t, p, z, Kind::Dew)?;
                    Ok::<_, ThermoError>(DewPoint {
                        temperature: t,
                        pressure: p,
                        liquid: eq.other,
                        vapor: z.to_vec(),
                        k_values: eq.k,
                        converged: true,
                    })
                }),
            None => solver.dew_point_temperature(p, z).map(Ok).ok(),
        };
        if let Some(Ok(dp)) = dp {
            td_guess = Some(dp.temperature.value);
            dew.push(point_from_dew(&dp));
        }
    }
    Ok(Envelope { bubble, dew })
}
