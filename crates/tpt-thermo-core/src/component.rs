//! Component property access, unit-safe.
//!
//! Implemented by `tpt-thermo-data` (Phase 3) over the curated seed dataset;
//! declared here so EoS crates can be generic over any source of critical
//! constants, acentric factors, and molar masses.

use crate::error::ThermoError;
use crate::quantities::{MolarMass, Pressure, Temperature};

/// Read-only access to per-component physical constants and parameters.
pub trait ComponentDatabase: Send + Sync {
    /// Number of components in the database.
    fn num_components(&self) -> usize;

    /// Canonical name of component `i`.
    fn name(&self, i: usize) -> Result<&str, ThermoError>;

    /// Critical temperature.
    fn critical_temperature(&self, i: usize) -> Result<Temperature, ThermoError>;

    /// Critical pressure.
    fn critical_pressure(&self, i: usize) -> Result<Pressure, ThermoError>;

    /// Acentric factor `ω`.
    fn acentric_factor(&self, i: usize) -> Result<f64, ThermoError>;

    /// Molar mass (kg·mol⁻¹).
    fn molar_mass(&self, i: usize) -> Result<MolarMass, ThermoError>;

    /// Optional binary interaction parameter `k_ij` for the `(i, j)` pair.
    fn binary_interaction(&self, i: usize, j: usize) -> Result<f64, ThermoError> {
        if i >= self.num_components() || j >= self.num_components() {
            return Err(ThermoError::IndexOutOfRange(if i > j { i } else { j }));
        }
        // Default: no correction.
        Ok(0.0)
    }
}
