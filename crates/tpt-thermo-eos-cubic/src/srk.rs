//! Soave-Redlich-Kwong (1972) equation of state.

use crate::alpha::soave;
use crate::cubic_solver::CubicModel;
use crate::engine::CubicEos;
use crate::mixing::{CubicMixing, VdwMixing};
use alloc::boxed::Box;
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{
    MolarEnergy, MolarEntropy, MolarHeatCapacity, MolarVolume, Pressure, Temperature, Velocity,
};
use tpt_thermo_core::{EquationOfState, ThermoError};

/// Soave-Redlich-Kwong cubic equation of state:
/// `P = RT/(v−b) − aα/(v² + bv)`.
///
/// Defaults to the Soave alpha function and van der Waals one-fluid mixing.
pub struct SoaveRedlichKwong {
    inner: CubicEos,
}

impl SoaveRedlichKwong {
    /// Build from a [`ComponentDatabase`] with the default alpha and vdW mixing
    /// (binary interactions default to `0.0`).
    pub fn from_database(db: &dyn ComponentDatabase) -> Result<Self, ThermoError> {
        let inner = CubicEos::from_database(
            CubicModel::SoaveRedlichKwong,
            db,
            soave(),
            Box::new(VdwMixing::new(db.num_components())),
        )?;
        Ok(Self { inner })
    }

    /// Build from a [`ComponentDatabase`] with the default alpha and van der Waals
    /// one-fluid mixing that consumes the database's seeded/fitted `k_ij`
    /// parameters (opt-in; use [`SoaveRedlichKwong::from_database`] for the
    /// zero-BIP default).
    pub fn from_database_with_kij(db: &dyn ComponentDatabase) -> Result<Self, ThermoError> {
        let inner = CubicEos::from_database(
            CubicModel::SoaveRedlichKwong,
            db,
            soave(),
            Box::new(VdwMixing::from_database(db)),
        )?;
        Ok(Self { inner })
    }

    /// Build with an explicit mixing rule.
    pub fn with_mixing(
        db: &dyn ComponentDatabase,
        mixing: Box<dyn CubicMixing>,
    ) -> Result<Self, ThermoError> {
        let inner = CubicEos::from_database(CubicModel::SoaveRedlichKwong, db, soave(), mixing)?;
        Ok(Self { inner })
    }

    /// Build with an explicit mixing rule and alpha function.
    pub fn with_alpha_and_mixing(
        db: &dyn ComponentDatabase,
        alpha: Box<dyn crate::alpha::AlphaFunction>,
        mixing: Box<dyn CubicMixing>,
    ) -> Result<Self, ThermoError> {
        let inner = CubicEos::from_database(CubicModel::SoaveRedlichKwong, db, alpha, mixing)?;
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

    /// Pure-component critical point.
    pub fn critical_point_pure(
        &self,
        i: usize,
    ) -> Result<(Temperature, Pressure, MolarVolume), ThermoError> {
        self.inner.critical_point_pure(i)
    }
}

impl EquationOfState for SoaveRedlichKwong {
    fn num_components(&self) -> usize {
        self.inner.num_components()
    }
    fn pressure(&self, t: Temperature, v: MolarVolume, z: &[f64]) -> Result<Pressure, ThermoError> {
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
