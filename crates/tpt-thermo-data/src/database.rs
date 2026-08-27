//! A [`ComponentDatabase`](tpt_thermo_core::component::ComponentDatabase)
//! implementation backed by an in-memory, serde-loaded component set.

use crate::{bip::BipTable, record::ComponentRecord};
use serde::{Deserialize, Serialize};
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{MolarMass, Pressure, Temperature};
use uom::si::{molar_mass::kilogram_per_mole, pressure::pascal, thermodynamic_temperature::kelvin};

/// Helper for parsing/serializing TOML documents that use `[[components]]` at
/// the root.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct RawDb {
    components: Vec<ComponentRecord>,
}

/// An in-memory component database, optionally carrying a [`BipTable`].
#[derive(Debug, Clone)]
pub struct SeedComponentDatabase {
    records: Vec<ComponentRecord>,
    bip: BipTable,
}

impl SeedComponentDatabase {
    /// Load the embedded curated seed dataset and validate it.
    ///
    /// # Panics
    ///
    /// Panics if the embedded seed dataset fails validation (it is part of the
    /// crate and must always be valid).
    pub fn from_seed() -> Self {
        Self::from_toml_str(crate::seed::SEED_TOML).expect("embedded seed dataset failed to parse")
    }

    /// Parse a database from a TOML document of `[[components]]` records.
    pub fn from_toml_str(s: &str) -> Result<Self, String> {
        let raw: RawDb = toml::from_str(s).map_err(|e| e.to_string())?;
        let records = raw.components;
        let db = Self {
            records,
            bip: BipTable::default(),
        };
        let issues = db.validate();
        if !issues.is_empty() {
            return Err(issues.join("; "));
        }
        Ok(db)
    }

    /// Serialize the database back to canonical TOML (`[[components]]`).
    pub fn to_toml(&self) -> Result<String, String> {
        let raw = RawDb {
            components: self.records.clone(),
        };
        toml::to_string(&raw).map_err(|e| e.to_string())
    }

    /// Parse a database from a JSON document.
    pub fn from_json_str(s: &str) -> Result<Self, String> {
        let records: Vec<ComponentRecord> = serde_json::from_str(s).map_err(|e| e.to_string())?;
        let db = Self {
            records,
            bip: BipTable::default(),
        };
        let issues = db.validate();
        if !issues.is_empty() {
            return Err(issues.join("; "));
        }
        Ok(db)
    }

    /// The raw records.
    pub fn records(&self) -> &[ComponentRecord] {
        &self.records
    }

    /// Index of a component by name, if present.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.records.iter().position(|r| r.name == name)
    }

    /// Replace the binary interaction table.
    pub fn with_bip(mut self, bip: BipTable) -> Self {
        self.bip = bip;
        self
    }

    /// Validate every record and the dataset as a whole (e.g. unique names).
    /// Returns the list of issues (empty when valid).
    pub fn validate(&self) -> Vec<String> {
        let mut issues = Vec::new();
        for r in &self.records {
            issues.extend(r.validate());
        }
        let mut seen = std::collections::HashSet::new();
        for r in &self.records {
            if !seen.insert(r.name.clone()) {
                issues.push(format!("duplicate component name '{}'", r.name));
            }
        }
        issues
    }
}

impl ComponentDatabase for SeedComponentDatabase {
    fn num_components(&self) -> usize {
        self.records.len()
    }

    fn name(&self, i: usize) -> Result<&str, ThermoError> {
        self.records
            .get(i)
            .map(|r| r.name.as_str())
            .ok_or(ThermoError::IndexOutOfRange(i))
    }

    fn critical_temperature(&self, i: usize) -> Result<Temperature, ThermoError> {
        let r = self.records.get(i).ok_or(ThermoError::IndexOutOfRange(i))?;
        Ok(Temperature::new::<kelvin>(r.critical_temperature_k))
    }

    fn critical_pressure(&self, i: usize) -> Result<Pressure, ThermoError> {
        let r = self.records.get(i).ok_or(ThermoError::IndexOutOfRange(i))?;
        Ok(Pressure::new::<pascal>(r.critical_pressure_pa))
    }

    fn acentric_factor(&self, i: usize) -> Result<f64, ThermoError> {
        let r = self.records.get(i).ok_or(ThermoError::IndexOutOfRange(i))?;
        Ok(r.acentric_factor)
    }

    fn molar_mass(&self, i: usize) -> Result<MolarMass, ThermoError> {
        let r = self.records.get(i).ok_or(ThermoError::IndexOutOfRange(i))?;
        Ok(MolarMass::new::<kilogram_per_mole>(r.molar_mass_kg_per_mol))
    }

    fn binary_interaction(&self, i: usize, j: usize) -> Result<f64, ThermoError> {
        if i >= self.num_components() || j >= self.num_components() {
            return Err(ThermoError::IndexOutOfRange(if i > j { i } else { j }));
        }
        Ok(self.bip.get(i, j))
    }
}
