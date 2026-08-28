//! PC-SAFT for polymers — a thin specialization of Phase 6's PC-SAFT.
//!
//! Polymer chains are simply PC-SAFT molecules with a large segment count `m`, so
//! this module re-uses [`tpt_thermo_eos_saft::PcSaft`] unchanged and adds a
//! convenience constructor that maps a monomer (segment) parameter set and a chain
//! length to a polymer [`SaftComponent`](tpt_thermo_eos_saft::SaftComponent).
//! [`PolymerPcSaft`] therefore reduces exactly to [`PcSaft`] in the `m →` limit
//! (regression-tested below).

use alloc::vec::Vec;
use tpt_thermo_core::quantities::{
    MolarEnergy, MolarEntropy, MolarHeatCapacity, MolarVolume, Pressure, Temperature, Velocity,
};
use tpt_thermo_core::{EquationOfState, ThermoError};
use tpt_thermo_eos_saft::{PcSaft, SaftComponent, SaftParameters};

/// A polymer PC-SAFT model: a thin wrapper over [`PcSaft`].
#[derive(Debug, Clone)]
pub struct PolymerPcSaft(PcSaft);

impl PolymerPcSaft {
    /// Build from a [`SaftParameters`] table and per-component molar masses.
    pub fn new(params: SaftParameters, molar_masses: Vec<f64>) -> Self {
        Self(PcSaft::new(params, molar_masses))
    }

    /// Build a single-component polymer model from monomer (per-segment)
    /// parameters and a number of segments per chain.
    ///
    /// The polymer inherits the monomer's `σ` and `ε/k`; only the segment count
    /// scales by `segments`.
    pub fn from_monomer(
        name: &'static str,
        monomer_m: f64,
        monomer_sigma: f64,
        monomer_epsilon_k: f64,
        segments: f64,
        molar_mass: f64,
    ) -> Self {
        let comp =
            SaftComponent::pc_saft(name, monomer_m * segments, monomer_sigma, monomer_epsilon_k);
        Self(PcSaft::new(
            SaftParameters::new(alloc::vec![comp]),
            vec![molar_mass],
        ))
    }

    /// Access the underlying engine.
    pub fn inner(&self) -> &PcSaft {
        &self.0
    }
}

impl EquationOfState for PolymerPcSaft {
    fn num_components(&self) -> usize {
        self.0.num_components()
    }
    fn pressure(&self, t: Temperature, v: MolarVolume, z: &[f64]) -> Result<Pressure, ThermoError> {
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
    ) -> Result<MolarEnergy, ThermoError> {
        self.0.molar_enthalpy(t, v, z)
    }
    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEntropy, ThermoError> {
        self.0.molar_entropy(t, v, z)
    }
    fn molar_isobaric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        self.0.molar_isobaric_heat_capacity(t, v, z)
    }
    fn molar_isochoric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        self.0.molar_isochoric_heat_capacity(t, v, z)
    }
    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Velocity, ThermoError> {
        self.0.speed_of_sound(t, v, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::molar_volume::cubic_meter_per_mole;
    use uom::si::thermodynamic_temperature::kelvin;

    #[test]
    fn reduces_to_base_pc_saft() {
        // A polymer with m = 5 and a plain PC-SAFT component with m = 5 must agree.
        let polymer = PolymerPcSaft::from_monomer("poly", 1.0, 3.0, 250.0, 5.0, 0.05);
        let base = PcSaft::new(
            SaftParameters::new(vec![SaftComponent::pc_saft("poly", 5.0, 3.0, 250.0)]),
            vec![0.05],
        );
        let t = Temperature::new::<kelvin>(350.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(0.01);
        let pp = polymer.pressure(t, v, &[1.0]).unwrap();
        let pb = base.pressure(t, v, &[1.0]).unwrap();
        assert!((pp.value - pb.value).abs() / pb.value < 1e-9);
        let lnp = polymer.ln_fugacity_coefficient(t, v, &[1.0], 0).unwrap();
        let lnb = base.ln_fugacity_coefficient(t, v, &[1.0], 0).unwrap();
        assert!((lnp - lnb).abs() < 1e-9);
    }
}
