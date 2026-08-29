//! `tpt-thermo-data` — curated component dataset, database, and BIP tables.
//!
//! This crate provides the data substrate the EoS/flash/phase crates build on:
//!
//! * [`record::ComponentRecord`] — a serde (TOML/JSON) schema in base units
//!   with physical-constraint validation.
//! * [`database::SeedComponentDatabase`] — a `ComponentDatabase`
//! ([`tpt_thermo_core::component::ComponentDatabase`]) implementation loaded
//! from the embedded curated seed set (`data/seed.toml`, ~2300 compounds) or
//! any user-supplied TOML/JSON.
//! * [`bip::BipTable`] — binary-interaction-parameter storage. The curated seed
//!   now ships fitted PR/SRK `k_ij` values for common pairs
//!   (`[[binary_interactions]]` in `data/seed.toml`); every other pair defaults
//!   to `0.0`. These are consumed opt-in by the cubic crate.
//!
//! Parameter estimation utilities (spec 3d) are deferred to Phase 4+.

pub mod bip;
pub mod database;
pub mod record;
pub mod seed;

pub use bip::{BinaryInteractionRecord, BipTable};
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
        // Diagonal k_ij is always 0.
        assert_eq!(db.binary_interaction(i, i).unwrap(), 0.0);
        // An unseeded pair defaults to 0.
        let neopentane = db.index_of("neopentane").unwrap();
        assert_eq!(db.binary_interaction(i, neopentane).unwrap(), 0.0);
    }

    #[test]
    fn expanded_seed_covers_common_chemicals() {
        // The deferred-scope expansion pushed the curated set well beyond the
        // original ~58 compounds; assert a representative slice resolves.
        let db = SeedComponentDatabase::from_seed();
        assert!(
            db.num_components() >= 2000,
            "expanded seed should exceed 2000 compounds"
        );
        for name in [
            "butene",
            "pentanol",
            "ethyl acetate",
            "pyridine",
            "xylene",
            "methylcyclohexane",
            "dichloromethane",
            "chlorodifluoromethane",
            "thiophene",
            "n-eicosane",
            "pyrene",
            "benzoic acid",
            "bromobenzene",
            "r245fa",
        ] {
            assert!(
                db.index_of(name).is_some(),
                "missing expanded compound {name}"
            );
        }
    }

    #[test]
    fn seeded_binary_interactions_resolve() {
        let db = SeedComponentDatabase::from_seed();
        let co2 = db.index_of("carbon dioxide").unwrap();
        let methane = db.index_of("methane").unwrap();
        // Seeded PR k_ij for CO2-methane per Poling et al. 2001.
        assert!((db.binary_interaction(co2, methane).unwrap() - 0.10).abs() < 1e-9);
        // Symmetric.
        assert!((db.binary_interaction(methane, co2).unwrap() - 0.10).abs() < 1e-9);
        // Water-methane is a strongly non-ideal pair.
        let water = db.index_of("water").unwrap();
        assert!((db.binary_interaction(water, methane).unwrap() - 0.49).abs() < 1e-9);
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
