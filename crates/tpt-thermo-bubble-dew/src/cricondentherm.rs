//! Cricondenbar and cricondentherm of a phase envelope.
//!
//! * **Cricondenbar** — the pressure maximum of the two-phase envelope (and the
//!   temperature at which it occurs).
//! * **Cricondentherm** — the temperature maximum of the two-phase envelope (and
//!   the pressure at which it occurs).
//!
//! Both are read directly off the traced [`Envelope`] (bubble + dew points).

use crate::envelope::Envelope;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};

/// The cricondenbar and cricondentherm of an envelope.
#[derive(Debug, Clone)]
pub struct Criconden {
    /// `(P, T)` at the cricondenbar (maximum pressure).
    pub cricondenbar: (Pressure, Temperature),
    /// `(P, T)` at the cricondentherm (maximum temperature).
    pub cricondentherm: (Pressure, Temperature),
}

/// Extract the cricondenbar and cricondentherm from a traced envelope.
///
/// Returns `None` if the envelope carries no points.
pub fn cricondenbar_cricondentherm(envelope: &Envelope) -> Option<Criconden> {
    let mut pts: Vec<(f64, f64)> = envelope
        .bubble
        .iter()
        .map(|p| (p.pressure.value, p.temperature.value))
        .collect();
    pts.extend(
        envelope
            .dew
            .iter()
            .map(|p| (p.pressure.value, p.temperature.value)),
    );
    if pts.is_empty() {
        return None;
    }
    let max_p = pts
        .iter()
        .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
        .unwrap();
    let max_t = pts
        .iter()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();
    Some(Criconden {
        cricondenbar: (
            Pressure::new::<pascal>(max_p.0),
            Temperature::new::<kelvin>(max_p.1),
        ),
        cricondentherm: (
            Pressure::new::<pascal>(max_t.0),
            Temperature::new::<kelvin>(max_t.1),
        ),
    })
}
