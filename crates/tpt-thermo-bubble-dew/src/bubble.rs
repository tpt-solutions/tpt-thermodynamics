//! Bubble-point calculations: find `T` at fixed `P`, and find `P` at fixed `T`.

use crate::equilibrium::Kind;
use crate::{BubbleDewSolver, one_atm};
use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};

/// A bubble-point solution.
#[derive(Debug, Clone)]
pub struct BubblePoint {
    /// Bubble-point temperature.
    pub temperature: Temperature,
    /// Bubble-point pressure.
    pub pressure: Pressure,
    /// Liquid (feed) composition `x`.
    pub liquid: Vec<f64>,
    /// Incipient vapor composition `y` at the bubble point.
    pub vapor: Vec<f64>,
    /// Equilibrium K-values `K_i = y_i / x_i`.
    pub k_values: Vec<f64>,
    /// Whether the inner fugacity iteration converged.
    pub converged: bool,
}

impl<'a> BubbleDewSolver<'a> {
    /// Bubble-point temperature of a liquid composition `x` at pressure `p`:
    /// the lowest temperature at which a vapor phase appears.
    pub fn bubble_point_temperature(
        &self,
        p: Pressure,
        x: &[f64],
    ) -> Result<BubblePoint, ThermoError> {
        let tb = self.boundary_temperature(p, x, true)?;
        let t = Temperature::new::<kelvin>(tb);
        let eq = self.equilibrium_at(t, p, x, Kind::Bubble)?;
        Ok(BubblePoint {
            temperature: t,
            pressure: p,
            liquid: x.to_vec(),
            vapor: eq.other,
            k_values: eq.k,
            converged: true,
        })
    }

    /// Bubble-point pressure of a liquid composition `x` at temperature `t`:
    /// the highest pressure at which a vapor phase appears.
    pub fn bubble_point_pressure(
        &self,
        t: Temperature,
        x: &[f64],
    ) -> Result<BubblePoint, ThermoError> {
        let pb = self.boundary_pressure(t, x, false)?;
        let p = Pressure::new::<pascal>(pb);
        let eq = self.equilibrium_at(t, p, x, Kind::Bubble)?;
        Ok(BubblePoint {
            temperature: t,
            pressure: p,
            liquid: x.to_vec(),
            vapor: eq.other,
            k_values: eq.k,
            converged: true,
        })
    }

    /// Convenience: bubble-point temperature at 1 atm.
    pub fn bubble_point_1atm(&self, x: &[f64]) -> Result<BubblePoint, ThermoError> {
        self.bubble_point_temperature(one_atm(), x)
    }
}
