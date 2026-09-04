//! Binary interaction parameter (k_ij) tables.
//!
//! Phase 3 ships the structure and loader. The curated seed (`data/seed.toml`)
//! now seeds a name-keyed `[[binary_interactions]]` section of fitted PR/SRK
//! k_ij values for common pairs; [`BipTable::from_name_records`] resolves those
//! names to database indices. Every pair not explicitly listed defaults to `0.0`.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A single binary interaction parameter, referenced by the **names** of the
/// two components (so the table is robust to reordering the component list).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BinaryInteractionRecord {
    /// First component name (must match a `ComponentRecord::name`).
    pub a: String,
    /// Second component name (must match a `ComponentRecord::name`).
    pub b: String,
    /// Dimensionless binary interaction parameter `k_ij` (symmetric).
    pub k_ij: f64,
    /// Optional temperature-dependent coefficients `(a, b, c)` so that
    /// `k_ij(T) = a + b/T + c·ln(T)`. The constant `k_ij` above is used when
    /// these are absent — the cubic crate currently consumes only the constant
    /// value through the core `binary_interaction` trait method.
    #[serde(default)]
    pub td_a: Option<f64>,
    #[serde(default)]
    pub td_b: Option<f64>,
    #[serde(default)]
    pub td_c: Option<f64>,
    /// Provenance description for the value.
    #[serde(default)]
    pub source: Option<String>,
}

/// A table of binary interaction parameters, keyed by an `"i_j"` string with
/// `i <= j`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BipTable {
    /// Map from `"i_j"` (i <= j, both decimal) to the dimensionless k_ij.
    #[serde(default)]
    pub entries: HashMap<String, f64>,
}

impl BipTable {
    /// Build a table from name-keyed records, resolving each `(a, b)` pair to the
    /// components' indices via `index_of` (e.g.
    /// `SeedComponentDatabase::index_of`). Returns an error if either name is not
    /// present in the database.
    pub fn from_name_records(
        records: &[BinaryInteractionRecord],
        index_of: impl Fn(&str) -> Option<usize>,
    ) -> Result<Self, String> {
        let mut table = BipTable::default();
        for r in records {
            let (Some(i), Some(j)) = (index_of(&r.a), index_of(&r.b)) else {
                return Err(format!(
                    "binary interaction references unknown component(s): '{}' / '{}'",
                    r.a, r.b
                ));
            };
            table.set(i, j, r.k_ij);
        }
        Ok(table)
    }

    /// Look up `k_ij` for the `(i, j)` pair (symmetric). Returns `0.0` for a
    /// diagonal entry or any unset pair.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i == j {
            return 0.0;
        }
        let key = format!("{}_{}", i.min(j), i.max(j));
        *self.entries.get(&key).unwrap_or(&0.0)
    }

    /// Set `k_ij` for the `(i, j)` pair.
    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        if i == j {
            return;
        }
        let key = format!("{}_{}", i.min(j), i.max(j));
        self.entries.insert(key, value);
    }
}
