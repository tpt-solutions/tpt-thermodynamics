//! The [`EquationOfState`] trait: the central abstraction every model in the
//! workspace implements, plus a [`State`] value and a fully-working
//! [`IdealGas`] reference implementation used as living documentation/tests.

use crate::composition::Composition;
use crate::convergence::ConvergenceStatus;
use crate::error::ThermoError;
use crate::numerics::{brent, ROOT_MAX_ITER, ROOT_TOL};
use crate::quantities::{
    gas_constant, molar_entropy, MolarEnergy, MolarEntropy, MolarHeatCapacity, MolarMass,
    MolarVolume, Pressure, Temperature, Velocity,
};
use alloc::vec::Vec;
use uom::si::molar_energy::joule_per_mole;
use uom::si::molar_volume::cubic_meter_per_mole;
use uom::si::pressure::pascal as pascal_unit;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::velocity::meter_per_second;

/// A thermodynamic state of a mixture, stored canonically as `(T, P, v, z)`.
#[derive(Debug, Clone, PartialEq)]
pub struct State {
    /// Temperature.
    pub temperature: Temperature,
    /// Pressure.
    pub pressure: Pressure,
    /// Molar volume.
    pub molar_volume: MolarVolume,
    /// Mole fractions (normalised, length = number of components).
    pub composition: Vec<f64>,
}

impl State {
    /// Construct directly from primitives. `composition` must be non-empty and
    /// sum to 1.
    pub fn new(
        temperature: Temperature,
        pressure: Pressure,
        molar_volume: MolarVolume,
        composition: Vec<f64>,
    ) -> Result<Self, ThermoError> {
        if composition.is_empty() {
            return Err(ThermoError::InvalidInput("empty composition"));
        }
        let sum: f64 = composition.iter().sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(ThermoError::InvalidInput("composition does not sum to 1"));
        }
        Ok(Self {
            temperature,
            pressure,
            molar_volume,
            composition,
        })
    }

    /// Number of components.
    pub fn num_components(&self) -> usize {
        self.composition.len()
    }
}

/// Builder for [`State`] with light validation; composition is normalised on
/// build.
#[derive(Debug, Clone, Default)]
pub struct StateBuilder {
    temperature: Option<Temperature>,
    pressure: Option<Pressure>,
    molar_volume: Option<MolarVolume>,
    composition: Option<Vec<f64>>,
}

impl StateBuilder {
    /// Begin a build for `temperature`.
    pub fn new(temperature: Temperature) -> Self {
        Self {
            temperature: Some(temperature),
            ..Default::default()
        }
    }

    /// Set the pressure.
    pub fn pressure(mut self, p: Pressure) -> Self {
        self.pressure = Some(p);
        self
    }

    /// Set the molar volume.
    pub fn molar_volume(mut self, v: MolarVolume) -> Self {
        self.molar_volume = Some(v);
        self
    }

    /// Set the (un-normalised) composition; it will be normalised on build.
    pub fn composition(mut self, z: Vec<f64>) -> Self {
        self.composition = Some(z);
        self
    }

    /// Finalise.
    pub fn build(self) -> Result<State, ThermoError> {
        let t = self
            .temperature
            .ok_or(ThermoError::InvalidInput("missing temperature"))?;
        let p = self
            .pressure
            .ok_or(ThermoError::InvalidInput("missing pressure"))?;
        let v = self
            .molar_volume
            .ok_or(ThermoError::InvalidInput("missing molar volume"))?;
        let z = self
            .composition
            .ok_or(ThermoError::InvalidInput("missing composition"))?;
        let comp = Composition::from_mole_fractions(z)
            .map_err(|_| ThermoError::InvalidInput("invalid composition"))?;
        Ok(State {
            temperature: t,
            pressure: p,
            molar_volume: v,
            composition: comp.mole_fractions().to_vec(),
        })
    }
}

/// An equation of state: maps `(T, v, z)` to pressure and the derived
/// thermodynamic properties needed across the workspace.
///
/// `pressure` and `ln_fugacity_coefficient` are required (each model supplies
/// its own closure relation); the remaining methods have numerical-default
/// implementations where the spec permits, computed from `pressure`.
pub trait EquationOfState: Send + Sync {
    /// Number of components the model describes.
    fn num_components(&self) -> usize;

    /// Pressure at `(T, v, z)`.
    fn pressure(&self, t: Temperature, v: MolarVolume, z: &[f64]) -> Result<Pressure, ThermoError>;

    /// Natural log of the fugacity coefficient of component `i` at `(T, v, z)`.
    fn ln_fugacity_coefficient(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError>;

    /// Molar enthalpy at `(T, v, z)`.
    fn molar_enthalpy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEnergy, ThermoError>;

    /// Molar entropy at `(T, v, z)`.
    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEntropy, ThermoError>;

    /// Molar isobaric heat capacity at `(T, v, z)`.
    fn molar_isobaric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError>;

    /// Molar isochoric heat capacity at `(T, v, z)`.
    fn molar_isochoric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError>;

    /// Speed of sound at `(T, v, z)`.
    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Velocity, ThermoError>;

    /// Compressibility factor `Z = P v / (R T)`.
    fn compressibility_factor(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<f64, ThermoError> {
        let p = self.pressure(t, v, z)?;
        let pv = p * v;
        let rt = gas_constant() * t;
        Ok((pv / rt).value)
    }

    /// Fugacity of component `i`: `f_i = φ_i · z_i · P`.
    fn fugacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
        i: usize,
    ) -> Result<Pressure, ThermoError> {
        let phi = (self.ln_fugacity_coefficient(t, v, z, i)?).exp();
        let p = self.pressure(t, v, z)?;
        let zi = *z.get(i).ok_or(ThermoError::IndexOutOfRange(i))?;
        Ok(Pressure::new::<pascal_unit>(phi * zi * p.value))
    }

    /// Isothermal compressibility `κ_T = -(1/V)(∂V/∂P)_T` (Pa⁻¹), via central
    /// differences on [`pressure`](EquationOfState::pressure).
    fn isothermal_compressibility(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<f64, ThermoError> {
        let dpdv = self.dp_dv(t, v, z)?;
        if dpdv.abs() < 1e-30 {
            return Err(ThermoError::Numerical(ConvergenceStatus::NumericalIssue(
                crate::convergence::NumericalIssueReason::SingularJacobian,
            )));
        }
        Ok(-1.0 / (v.value * dpdv))
    }

    /// Thermal expansion coefficient `α = (1/V)(∂V/∂T)_P` (K⁻¹), derived from
    /// `(∂P/∂T)_V` and `(∂P/∂V)_T`.
    fn thermal_expansion_coefficient(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<f64, ThermoError> {
        let dpdt = self.dp_dt(t, v, z)?;
        let dpdv = self.dp_dv(t, v, z)?;
        if dpdv.abs() < 1e-30 {
            return Err(ThermoError::Numerical(ConvergenceStatus::NumericalIssue(
                crate::convergence::NumericalIssueReason::SingularJacobian,
            )));
        }
        Ok(-dpdt / (v.value * dpdv))
    }

    /// Invert the EoS to find the molar volume giving `target_pressure` at
    /// `(T, z)`, bracketed by `[v_min, v_max]`. Uses Brent's method.
    fn solve_molar_volume(
        &self,
        t: Temperature,
        target_pressure: Pressure,
        z: &[f64],
        v_min: MolarVolume,
        v_max: MolarVolume,
    ) -> Result<MolarVolume, ThermoError> {
        let p_target = target_pressure.value;
        let f = |v: f64| -> f64 {
            self.pressure(t, MolarVolume::new::<cubic_meter_per_mole>(v), z)
                .map(|p| p.value - p_target)
                .unwrap_or(f64::INFINITY)
        };
        let root = brent(f, v_min.value, v_max.value, ROOT_TOL, ROOT_MAX_ITER)
            .map_err(ThermoError::Numerical)?;
        Ok(MolarVolume::new::<cubic_meter_per_mole>(root))
    }

    /// `(∂P/∂V)_T` by central difference.
    fn dp_dv(&self, t: Temperature, v: MolarVolume, z: &[f64]) -> Result<f64, ThermoError> {
        let h = v.value.abs().max(1e-8) * 1e-6;
        let pm = self
            .pressure(t, MolarVolume::new::<cubic_meter_per_mole>(v.value - h), z)?
            .value;
        let pp = self
            .pressure(t, MolarVolume::new::<cubic_meter_per_mole>(v.value + h), z)?
            .value;
        Ok((pp - pm) / (2.0 * h))
    }

    /// `(∂P/∂T)_V` by central difference.
    fn dp_dt(&self, t: Temperature, v: MolarVolume, z: &[f64]) -> Result<f64, ThermoError> {
        let h = t.value.abs().max(1.0) * 1e-6;
        let pm = self
            .pressure(Temperature::new::<kelvin>(t.value - h), v, z)?
            .value;
        let pp = self
            .pressure(Temperature::new::<kelvin>(t.value + h), v, z)?
            .value;
        Ok((pp - pm) / (2.0 * h))
    }
}

/// Ideal-gas equation of state: `Z = 1`, `φ_i = 1`, with constant-`Cp`
/// residual-free enthalpy/entropy and a mixture mean molar mass for the speed
/// of sound. This is the reference implementation every other model is checked
/// against.
#[derive(Debug, Clone, PartialEq)]
pub struct IdealGas {
    /// Number of components.
    pub num_components: usize,
    /// Constant molar isobaric heat capacity (J·mol⁻¹·K⁻¹).
    pub cp: MolarHeatCapacity,
    /// Mean molar mass (kg·mol⁻¹), used for the speed of sound.
    pub molar_mass: MolarMass,
    /// Reference temperature for the enthalpy/entropy zero (K).
    pub t_ref: Temperature,
    /// Reference pressure for the entropy zero (Pa).
    pub p_ref: Pressure,
}

impl IdealGas {
    /// Construct for `num_components` with the given constants.
    pub fn new(
        num_components: usize,
        cp: MolarHeatCapacity,
        molar_mass: MolarMass,
        t_ref: Temperature,
        p_ref: Pressure,
    ) -> Self {
        Self {
            num_components,
            cp,
            molar_mass,
            t_ref,
            p_ref,
        }
    }
}

impl EquationOfState for IdealGas {
    fn num_components(&self) -> usize {
        self.num_components
    }

    fn pressure(&self, t: Temperature, v: MolarVolume, z: &[f64]) -> Result<Pressure, ThermoError> {
        if z.len() != self.num_components {
            return Err(ThermoError::InvalidInput("composition length mismatch"));
        }
        // P = R T / v.
        Ok(gas_constant() * t / v)
    }

    fn ln_fugacity_coefficient(
        &self,
        _t: Temperature,
        _v: MolarVolume,
        z: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        if i >= z.len() {
            return Err(ThermoError::IndexOutOfRange(i));
        }
        Ok(0.0)
    }

    fn molar_enthalpy(
        &self,
        t: Temperature,
        _v: MolarVolume,
        _z: &[f64],
    ) -> Result<MolarEnergy, ThermoError> {
        // Enthalpy reference: constant-Cp, h = Cp·(T − T_ref). Subtract raw
        // values because `Temperature − Temperature` is a `TemperatureInterval`
        // (a different uom `Kind`) that does not multiply `MolarHeatCapacity`.
        let dt = t.value - self.t_ref.value;
        Ok(MolarEnergy::new::<joule_per_mole>(self.cp.value * dt))
    }

    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEntropy, ThermoError> {
        let p = self.pressure(t, v, z)?;
        let term1 = self.cp.value * (t.value / self.t_ref.value).ln();
        let term2 = crate::R * (p.value / self.p_ref.value).ln();
        let mixing: f64 = -crate::R
            * z.iter()
                .filter(|&&xi| xi > 0.0)
                .map(|&xi| xi * xi.ln())
                .sum::<f64>();
        Ok(molar_entropy(term1 - term2 + mixing))
    }

    fn molar_isobaric_heat_capacity(
        &self,
        _t: Temperature,
        _v: MolarVolume,
        _z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        Ok(self.cp)
    }

    fn molar_isochoric_heat_capacity(
        &self,
        _t: Temperature,
        _v: MolarVolume,
        _z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        Ok(self.cp - gas_constant())
    }

    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Velocity, ThermoError> {
        let p = self.pressure(t, v, z)?;
        let rho = self.molar_mass / v; // mass density
        let gamma = self.cp.value / (self.cp.value - crate::R);
        let a2 = (p / rho).value; // m²·s⁻²
        Ok(Velocity::new::<meter_per_second>((gamma * a2).sqrt()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::molar_heat_capacity::joule_per_kelvin_mole;
    use uom::si::molar_mass::kilogram_per_mole;
    use uom::si::pressure::pascal;
    use uom::si::thermodynamic_temperature::kelvin;

    fn ideal() -> IdealGas {
        IdealGas::new(
            1,
            MolarHeatCapacity::new::<joule_per_kelvin_mole>(33.6),
            MolarMass::new::<kilogram_per_mole>(0.018),
            Temperature::new::<kelvin>(298.15),
            Pressure::new::<pascal>(1.0e5),
        )
    }

    #[test]
    fn ideal_gas_law_pv_rt() {
        let m = ideal();
        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(crate::R * 300.0 / 1.0e5);
        let p = m.pressure(t, v, &[1.0]).unwrap();
        assert!((p.value - 1.0e5).abs() / 1.0e5 < 1e-9);
        let z = m.compressibility_factor(t, v, &[1.0]).unwrap();
        assert!((z - 1.0).abs() < 1e-9);
    }

    #[test]
    fn ideal_gas_fugacity_is_pressure() {
        let m = ideal();
        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(crate::R * 300.0 / 2.0e5);
        let f = m.fugacity(t, v, &[1.0], 0).unwrap();
        let p = m.pressure(t, v, &[1.0]).unwrap();
        assert!((f.value - p.value).abs() / p.value < 1e-9);
    }

    #[test]
    fn solve_molar_volume_recovers_rt_over_p() {
        let m = ideal();
        let t = Temperature::new::<kelvin>(350.0);
        let target = Pressure::new::<pascal>(3.0e5);
        let v = m
            .solve_molar_volume(
                t,
                target,
                &[1.0],
                MolarVolume::new::<cubic_meter_per_mole>(1e-6),
                MolarVolume::new::<cubic_meter_per_mole>(1e-1),
            )
            .unwrap();
        let expected = crate::R * 350.0 / 3.0e5;
        assert!((v.value - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn compressibility_and_expansion_sanity() {
        let m = ideal();
        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(crate::R * 300.0 / 1.0e5);
        let kt = m.isothermal_compressibility(t, v, &[1.0]).unwrap();
        assert!((kt - 1.0 / 1.0e5).abs() / (1.0 / 1.0e5) < 1e-6);
        let alpha = m.thermal_expansion_coefficient(t, v, &[1.0]).unwrap();
        assert!((alpha - 1.0 / 300.0).abs() / (1.0 / 300.0) < 1e-6);
    }

    #[test]
    fn state_builder_round_trip() {
        let s = StateBuilder::new(Temperature::new::<kelvin>(300.0))
            .pressure(Pressure::new::<pascal>(1.0e5))
            .molar_volume(MolarVolume::new::<cubic_meter_per_mole>(0.024))
            .composition(alloc::vec![0.5, 0.5])
            .build()
            .unwrap();
        assert_eq!(s.num_components(), 2);
        assert!((s.composition[0] - 0.5).abs() < 1e-9);
    }
}
