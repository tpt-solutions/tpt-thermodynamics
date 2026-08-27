//! Multiphase classification (V / L / V-L / L-L / V-L-L) via two-direction
//! tangent-plane-distance testing.

use crate::phase_volume::PhaseVolume;
use crate::tpd::TangentPlaneDistance;
use crate::tpd::TPD_TOL;
use alloc::vec::Vec;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use tpt_thermo_core::{ComponentDatabase, EquationOfState};
use tpt_thermo_eos_cubic::cubic_solver::Phase;

/// Outcome of a multiphase classification.
#[derive(Debug, Clone)]
pub struct MultiphaseResult {
    /// Inferred number of co-existing phases: 1 (stable), 2, or 3 (V-L-L).
    pub num_phases: usize,
    /// True when no incipient phase was found in either direction.
    pub stable: bool,
    /// Minimised vapor-like trial composition (when an incipient vapor exists).
    pub vapor_trial: Option<Vec<f64>>,
    /// Minimised liquid-like trial composition (when an incipient liquid exists).
    pub liquid_trial: Option<Vec<f64>>,
}

/// Classify the phase count of `z` at `(T, P)` by testing both an incipient
/// vapor (feed taken as liquid) and an incipient liquid (feed taken as vapor).
pub fn detect_phases<E: EquationOfState + ?Sized>(
    eos: &E,
    volume: &dyn PhaseVolume,
    db: &dyn ComponentDatabase,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> MultiphaseResult {
    let calc = TangentPlaneDistance::new(eos, volume, db, t, p, z.to_vec());
    let mut vapor_unstable = false;
    let mut liquid_unstable = false;
    let mut vapor_trial = None;
    let mut liquid_trial = None;

    // Feed as liquid → look for an incipient vapor phase.
    if let Some(sol) = calc.minimize(Phase::Liquid, Phase::Vapor) {
        if sol.tpd < -TPD_TOL {
            vapor_unstable = true;
        }
        vapor_trial = Some(sol.composition);
    }
    // Feed as vapor → look for an incipient liquid phase.
    if let Some(sol) = calc.minimize(Phase::Vapor, Phase::Liquid) {
        if sol.tpd < -TPD_TOL {
            liquid_unstable = true;
        }
        liquid_trial = Some(sol.composition);
    }

    let num_phases = match (vapor_unstable, liquid_unstable) {
        (false, false) => 1,
        (true, false) | (false, true) => 2,
        (true, true) => 3,
    };

    MultiphaseResult {
        num_phases,
        stable: !(vapor_unstable || liquid_unstable),
        vapor_trial,
        liquid_trial,
    }
}
