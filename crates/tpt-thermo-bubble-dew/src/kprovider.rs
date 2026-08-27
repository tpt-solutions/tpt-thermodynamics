//! Phase-selected molar-volume access for equations of state, the minimal
//! surface [`BubbleDewSolver`] needs beyond [`EquationOfState`].

use tpt_thermo_core::eos::EquationOfState;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_eos_cubic::cubic_solver::Phase as CubicPhase;
use tpt_thermo_eos_cubic::{PengRobinson, SoaveRedlichKwong, VolumeTranslated};

/// Which phase a molar volume is requested for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    /// Liquid (smallest molar volume root).
    Liquid,
    /// Vapor (largest molar volume root).
    Vapor,
}

impl Phase {
    fn to_cubic(self) -> CubicPhase {
        match self {
            Phase::Liquid => CubicPhase::Liquid,
            Phase::Vapor => CubicPhase::Vapor,
        }
    }
}

/// An equation of state that can return phase-selected molar volumes, in
/// addition to the [`EquationOfState`] surface (pressure and fugacity
/// coefficients). This is what lets the bubble/dew solver build two-phase
/// K-values and locate the two-phase boundary.
pub trait KProvider: EquationOfState + Send + Sync {
    /// Molar volume of `phase` at `(T, P, z)`.
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Result<MolarVolume, ThermoError>;
}

impl KProvider for PengRobinson {
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Result<MolarVolume, ThermoError> {
        self.solve_phase(t, p, z, phase.to_cubic())
    }
}

impl KProvider for SoaveRedlichKwong {
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Result<MolarVolume, ThermoError> {
        self.solve_phase(t, p, z, phase.to_cubic())
    }
}

impl KProvider for VolumeTranslated {
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Result<MolarVolume, ThermoError> {
        self.solve_phase(t, p, z, phase.to_cubic())
    }
}
