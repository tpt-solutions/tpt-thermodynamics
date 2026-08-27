//! Binary interaction parameter (k_ij) tables.
//!
//! Phase 3 ships the structure and loader; the actual fitted values are seeded
//! alongside the cubic (Phase 4) and activity (Phase 5) crates once their
//! consuming API shape is known. Until then every pair defaults to `0.0`.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A table of binary interaction parameters, keyed by an `"i_j"` string with
/// `i <= j`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BipTable {
    /// Map from `"i_j"` (i <= j, both decimal) to the dimensionless k_ij.
    #[serde(default)]
    pub entries: HashMap<String, f64>,
}

impl BipTable {
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
