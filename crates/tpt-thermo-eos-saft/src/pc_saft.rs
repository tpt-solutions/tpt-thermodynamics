//! PC-SAFT (perturbed-chain SAFT, Gross & Sadowski 2001).
//!
//! [`PcSaft`] is the primary SAFT model of this crate: hard-chain + dispersion
//! + association implemented over the shared [`SaftEngine`]. Build it from the
//! seed database or directly from [`SaftParameters`].

use crate::engine::{SaftEngine, SaftFlavor};
use crate::parameters::SaftParameters;
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_core::EquationOfState;
use uom::si::molar_volume::cubic_meter_per_mole;
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;

/// A perturbed-chain SAFT equation of state.
#[derive(Debug, Clone)]
pub struct PcSaft(pub(crate) SaftEngine);

impl PcSaft {
    /// Build directly from a parameter set and per-component molar masses
    /// (kg·mol⁻¹).
    pub fn new(params: SaftParameters, molar_masses: Vec<f64>) -> Self {
        Self(SaftEngine::new(params, molar_masses))
    }

    /// Build from the seed database (SAFT parameters + molar masses by name).
    pub fn from_seed_database(db: &dyn ComponentDatabase) -> Result<Self, ThermoError> {
        Ok(Self(SaftEngine::from_seed_database(db, SaftFlavor::PcSaft)?))
    }

    /// Attach a binary interaction matrix `k_ij`.
    pub fn with_kij(self, kij: Vec<Vec<f64>>) -> Self {
        Self(self.0.with_kij(kij))
    }

    /// Underlying SAFT parameters.
    pub fn parameters(&self) -> &SaftParameters {
        self.0.parameters()
    }

    /// Number of components.
    pub fn num_components(&self) -> usize {
        self.0.num_components()
    }

    /// Saturated vapor pressure and (vapor, liquid) molar volumes at `T`.
    pub fn saturation_pressure(
        &self,
        t: Temperature,
    ) -> Result<(Pressure, MolarVolume, MolarVolume), ThermoError> {
        let (p, vv, vl) = self.0.saturation_pressure(t)?;
        Ok((p, vv, vl))
    }
}

impl EquationOfState for PcSaft {
    fn num_components(&self) -> usize {
        self.0.num_components()
    }
    fn pressure(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Pressure, ThermoError> {
        self.0.pressure(t, v, z)
    }
    fn ln_fugacity_coefficient(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        self.0.ln_fugacity_coefficient(t, v, z, i)
    }
    fn molar_enthalpy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::MolarEnergy, ThermoError> {
        self.0.molar_enthalpy(t, v, z)
    }
    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::MolarEntropy, ThermoError> {
        self.0.molar_entropy(t, v, z)
    }
    fn molar_isobaric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::MolarHeatCapacity, ThermoError> {
        self.0.molar_isobaric_heat_capacity(t, v, z)
    }
    fn molar_isochoric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::MolarHeatCapacity, ThermoError> {
        self.0.molar_isochoric_heat_capacity(t, v, z)
    }
    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::Velocity, ThermoError> {
        self.0.speed_of_sound(t, v, z)
    }
}
