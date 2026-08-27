//! Phase envelope: trace the bubble and dew curves of a fixed overall
//! composition `z` across a pressure sweep (a simple, robust continuation that
//! steps in pressure and re-solves each bubble/dew point).
//!
//! > The approved plan names Phase 8's arc-length continuation machinery for
//! > this trace. That crate lands later; the pressure-swept form here is the
//! > same physical curve and is structured so it can be replaced by an
//! > arc-length driver without changing callers.

use crate::{BubbleDewSolver, BubblePoint, DewPoint};
use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{Pressure, Temperature};

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
pub fn bubble_dew_envelope(
    solver: &BubbleDewSolver<'_>,
    z: &[f64],
    pressures: &[Pressure],
) -> Result<Envelope, ThermoError> {
    let mut bubble = Vec::new();
    let mut dew = Vec::new();
    for &p in pressures {
        if let Ok(bp) = solver.bubble_point_temperature(p, z) {
            bubble.push(point_from_bubble(&bp));
        }
        if let Ok(dp) = solver.dew_point_temperature(p, z) {
            dew.push(point_from_dew(&dp));
        }
    }
    Ok(Envelope { bubble, dew })
}
