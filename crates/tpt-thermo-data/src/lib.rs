//! `tpt-thermo-data` — curated component dataset, database, and BIP tables.
//!
//! This crate provides the data substrate the EoS/flash/phase crates build on:
//!
//! * [`record::ComponentRecord`] — a serde (TOML/JSON) schema in base units
//!   with physical-constraint validation.
//! * [`database::SeedComponentDatabase`] — a `ComponentDatabase`
//!   ([`tpt_thermo_core::component::ComponentDatabase`]) implementation loaded
//!   from the embedded curated seed set (`data/seed.toml`, ~24 well-known
//!   compounds) or any user-supplied TOML/JSON.
//! * [`bip::BipTable`] — binary-interaction-parameter storage (defaults to
//!   `0.0`; fitted values are seeded alongside Phases 4/5).
//!
//! Parameter estimation utilities (spec 3d) are deferred to Phase 4+.

pub mod bip;
pub mod database;
pub mod record;
pub mod seed;

pub use bip::BipTable;
pub use database::SeedComponentDatabase;
pub use record::ComponentRecord;

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_thermo_core::component::ComponentDatabase;

    #[test]
    fn seed_loads_and_is_valid() {
        let db = SeedComponentDatabase::from_seed();
        assert!(db.num_components() >= 20, "expected a sizeable seed set");
        assert!(db.validate().is_empty(), "seed should be self-consistent");
    }

    #[test]
    fn seed_values_match_literature() {
        let db = SeedComponentDatabase::from_seed();
        let water = db.index_of("water").expect("water in seed");
        let tc = db.critical_temperature(water).unwrap().value;
        assert!((tc - 647.096).abs() < 1e-3, "water Tc");
        let co2 = db.index_of("carbon dioxide").expect("co2 in seed");
        let pc = db.critical_pressure(co2).unwrap().value;
        assert!((pc - 7.377e6).abs() / 7.377e6 < 1e-3, "co2 Pc");
        let m = db.molar_mass(water).unwrap().value;
        assert!((m - 0.01801528).abs() < 1e-6, "water M");
    }

    #[test]
    fn name_lookup_and_bip() {
        let db = SeedComponentDatabase::from_seed();
        let i = db.index_of("methane").unwrap();
        assert_eq!(db.name(i).unwrap(), "methane");
        // Default k_ij is 0; diagonal is always 0.
        assert_eq!(db.binary_interaction(i, i).unwrap(), 0.0);
        assert_eq!(db.binary_interaction(0, 1).unwrap(), 0.0);
    }

    #[test]
    fn toml_round_trip() {
        let db = SeedComponentDatabase::from_seed();
        let toml_str = db.to_toml().unwrap();
        let reparsed = SeedComponentDatabase::from_toml_str(&toml_str).unwrap();
        assert_eq!(reparsed.num_components(), db.num_components());
        assert_eq!(reparsed.records()[0], db.records()[0]);
    }

    #[test]
    fn json_round_trip() {
        let db = SeedComponentDatabase::from_seed();
        let json_str = serde_json::to_string(db.records()).unwrap();
        let reparsed = SeedComponentDatabase::from_json_str(&json_str).unwrap();
        assert_eq!(reparsed.num_components(), db.num_components());
    }

    #[test]
    fn invalid_record_rejected() {
        let bad = r#"
[[components]]
schema_version = 1
name = "bad"
critical_temperature_k = -10.0
critical_pressure_pa = 1.0e6
acentric_factor = 0.0
molar_mass_kg_per_mol = 0.016
"#;
        let err = SeedComponentDatabase::from_toml_str(bad);
        assert!(err.is_err(), "negative Tc should fail validation");
    }

    #[test]
    fn duplicate_names_rejected() {
        let dup = r#"
[[components]]
schema_version = 1
name = "x"
critical_temperature_k = 300.0
critical_pressure_pa = 1.0e6
acentric_factor = 0.0
molar_mass_kg_per_mol = 0.016

[[components]]
schema_version = 1
name = "x"
critical_temperature_k = 300.0
critical_pressure_pa = 1.0e6
acentric_factor = 0.0
molar_mass_kg_per_mol = 0.016
"#;
        assert!(SeedComponentDatabase::from_toml_str(dup).is_err());
    }
}
