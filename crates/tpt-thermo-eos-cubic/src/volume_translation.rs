//! Volume-translated (Peneloux) correction applied to a cubic EoS.
//!
//! The Peneloux translation shifts the predicted molar volume by a per-component
//! `c_i` while leaving the pressure and fugacity-coefficient relations
//! unchanged, markedly improving liquid-density predictions. It is layered on top
//! of an existing cubic EoS (Peng-Robinson by default).

use crate::alpha::soave;
use crate::cubic_solver::CubicModel;
use crate::engine::CubicEos;
use crate::mixing::VdwMixing;
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{MolarEnergy, MolarEntropy, MolarHeatCapacity, MolarVolume, Pressure, Temperature, Velocity};
use tpt_thermo_core::{EquationOfState, ThermoError, R};
use alloc::vec::Vec;

/// Peneloux volume-translation parameter `c_i` for component `(Tc, Pc, ω)`.
///
/// `c_i = −0.40768 · (R Tc / Pc) · (Z_RA − 0.29441)` with
/// `Z_RA = 0.29056 − 0.08775 ω` (so `c_i > 0`, the EoS volume is reduced).
pub fn peneloux_c(tc: f64, pc: f64, omega: f64) -> f64 {
    let z_ra = 0.29056 - 0.08775 * omega;
    -0.40768 * (R * tc / pc) * (z_ra - 0.29441)
}

/// A volume-translated cubic EoS (Peng-Robinson by default).
pub struct VolumeTranslated {
    inner: CubicEos,
}

impl VolumeTranslated {
    /// Build a volume-translated Peng-Robinson model from a database.
    pub fn peng_robinson(db: &dyn ComponentDatabase) -> Result<Self, ThermoError> {
        Self::build(CubicModel::PengRobinson, db)
    }

    /// Build a volume-translated Soave-Redlich-Kwong model from a database.
    pub fn soave_redlich_kwong(db: &dyn ComponentDatabase) -> Result<Self, ThermoError> {
        Self::build(CubicModel::SoaveRedlichKwong, db)
    }

    fn build(model: CubicModel, db: &dyn ComponentDatabase) -> Result<Self, ThermoError> {
        let n = db.num_components();
        let mut c = Vec::with_capacity(n);
        for i in 0..n {
            let tc = db.critical_temperature(i)?.value;
            let pc = db.critical_pressure(i)?.value;
            let omega = db.acentric_factor(i)?;
            c.push(peneloux_c(tc, pc, omega));
        }
        let inner = CubicEos::from_database(model, db, soave(), Box::new(VdwMixing::new(n)))?
            .with_volume_translation(c);
        Ok(Self { inner })
    }

    /// Access the underlying engine.
    pub fn engine(&self) -> &CubicEos {
        &self.inner
    }

    /// Solve for the molar volume of `phase` at `(T, P, z)`.
    pub fn solve_phase(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: crate::cubic_solver::Phase,
    ) -> Result<MolarVolume, ThermoError> {
        self.inner.solve_phase(t, p, z, phase)
    }

    /// Pure-component critical point (of the underlying cubic; the translation
    /// shifts the reported volume).
    pub fn critical_point_pure(
        &self,
        i: usize,
    ) -> Result<(Temperature, Pressure, MolarVolume), ThermoError> {
        self.inner.critical_point_pure(i)
    }
}

impl EquationOfState for VolumeTranslated {
    fn num_components(&self) -> usize {
        self.inner.num_components()
    }
    fn pressure(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Pressure, ThermoError> {
        self.inner.pressure(t, v, z)
    }
    fn ln_fugacity_coefficient(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        self.inner.ln_fugacity_coefficient(t, v, z, i)
    }
    fn molar_enthalpy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEnergy, ThermoError> {
        self.inner.molar_enthalpy(t, v, z)
    }
    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEntropy, ThermoError> {
        self.inner.molar_entropy(t, v, z)
    }
    fn molar_isobaric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        self.inner.molar_isobaric_heat_capacity(t, v, z)
    }
    fn molar_isochoric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        self.inner.molar_isochoric_heat_capacity(t, v, z)
    }
    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Velocity, ThermoError> {
        self.inner.speed_of_sound(t, v, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peneloux_c_positive_for_nonpolar() {
        // Methane-like: positive translation.
        let c = peneloux_c(190.6, 4.599e6, 0.011);
        assert!(c > 0.0);
    }

    #[test]
    fn peneloux_reduces_liquid_volume() {
        // c_i > 0 ⇒ physical volume = v_eos − c < v_eos (denser liquid).
        let c = peneloux_c(647.0, 22.06e6, 0.344);
        assert!(c > 0.0);
    }
}
