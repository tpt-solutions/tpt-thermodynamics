//! A single component's physical constants and provenance.
//!
//! Stored in *base units* (kelvin, pascal, kg·mol⁻¹) as plain `f64` values so the
//! record is trivially (de)serialisable with `serde` (the `uom` quantities do
//! not implement `serde`). Conversion helpers provide the typed
//! [`tpt_thermo_core::quantities`] views used by EoS crates.

use serde::{Deserialize, Serialize};

/// A curated/derived component record.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ComponentRecord {
    /// Schema version for dataset evolution (not a full audit log).
    pub schema_version: u32,
    /// Canonical component name (unique within a database).
    pub name: String,
    /// Molecular formula, if known.
    #[serde(default)]
    pub formula: Option<String>,
    /// CAS registry number, if known.
    #[serde(default)]
    pub cas: Option<String>,
    /// Critical temperature, K.
    pub critical_temperature_k: f64,
    /// Critical pressure, Pa.
    pub critical_pressure_pa: f64,
    /// Pitzer acentric factor, dimensionless.
    pub acentric_factor: f64,
    /// Molar mass, kg·mol⁻¹.
    pub molar_mass_kg_per_mol: f64,
    /// Normal boiling point, K (optional).
    #[serde(default)]
    pub normal_boiling_point_k: Option<f64>,
    /// Provenance description for the values.
    #[serde(default)]
    pub source: Option<String>,
}

impl ComponentRecord {
    /// Validate the physical-constraint sanity of the record. Returns the list
    /// of human-readable issues found (empty when valid).
    pub fn validate(&self) -> Vec<String> {
        let mut issues = Vec::new();
        if self.schema_version == 0 {
            issues.push("schema_version must be >= 1".into());
        }
        if self.name.trim().is_empty() {
            issues.push("name must not be empty".into());
        }
        if !(0.0..=2000.0).contains(&self.critical_temperature_k) {
            issues.push(format!(
                "critical_temperature_k {} out of plausible range",
                self.critical_temperature_k
            ));
        }
        if !(0.0..=1.0e9).contains(&self.critical_pressure_pa) {
            issues.push(format!(
                "critical_pressure_pa {} out of plausible range",
                self.critical_pressure_pa
            ));
        }
        if !(1e-4..=1.0).contains(&self.molar_mass_kg_per_mol) {
            issues.push(format!(
                "molar_mass_kg_per_mol {} out of plausible range",
                self.molar_mass_kg_per_mol
            ));
        }
        // Acentric factor: physically in roughly [-0.5, 1.5].
        if !(-0.6..=1.6).contains(&self.acentric_factor) {
            issues.push(format!(
                "acentric_factor {} out of plausible range",
                self.acentric_factor
            ));
        }
        if let Some(tb) = self.normal_boiling_point_k {
            if tb <= 0.0 || tb > self.critical_temperature_k {
                issues.push(format!(
                    "normal_boiling_point_k {} inconsistent with critical temperature",
                    tb
                ));
            }
        }
        issues
    }
}
