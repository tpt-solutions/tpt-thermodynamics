//! Dew-point calculations: find `T` at fixed `P`, and find `P` at fixed `T`.

use crate::equilibrium::Kind;
use crate::{BubbleDewSolver, one_atm};
use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};

/// A dew-point solution.
#[derive(Debug, Clone)]
pub struct DewPoint {
    /// Dew-point temperature.
    pub temperature: Temperature,
    /// Dew-point pressure.
    pub pressure: Pressure,
    /// Incipient liquid composition `x` at the dew point.
    pub liquid: Vec<f64>,
    /// Vapor (feed) composition `y`.
    pub vapor: Vec<f64>,
    /// Equilibrium K-values `K_i = y_i / x_i`.
    pub k_values: Vec<f64>,
    /// Whether the inner fugacity iteration converged.
    pub converged: bool,
}

impl<'a> BubbleDewSolver<'a> {
    /// Dew-point temperature of a vapor composition `y` at pressure `p`:
    /// the highest temperature at which a liquid phase appears.
    pub fn dew_point_temperature(
        &self,
        p: Pressure,
        y: &[f64],
    ) -> Result<DewPoint, ThermoError> {
        let td = self.boundary_temperature(p, y, false)?;
        let t = Temperature::new::<kelvin>(td);
        let eq = self.equilibrium_at(t, p, y, Kind::Dew)?;
        Ok(DewPoint {
            temperature: t,
            pressure: p,
            liquid: eq.other,
            vapor: y.to_vec(),
            k_values: eq.k,
            converged: true,
        })
    }

    /// Dew-point pressure of a vapor composition `y` at temperature `t`:
    /// the lowest pressure at which a liquid phase appears.
    pub fn dew_point_pressure(
        &self,
        t: Temperature,
        y: &[f64],
    ) -> Result<DewPoint, ThermoError> {
        let pd = self.boundary_pressure(t, y, true)?;
        let p = Pressure::new::<pascal>(pd);
        let eq = self.equilibrium_at(t, p, y, Kind::Dew)?;
        Ok(DewPoint {
            temperature: t,
            pressure: p,
            liquid: eq.other,
            vapor: y.to_vec(),
            k_values: eq.k,
            converged: true,
        })
    }

    /// Convenience: dew-point temperature at 1 atm.
    pub fn dew_point_1atm(&self, y: &[f64]) -> Result<DewPoint, ThermoError> {
        self.dew_point_temperature(one_atm(), y)
    }
}
