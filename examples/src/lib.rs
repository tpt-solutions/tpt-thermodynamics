//! Workspace example crate for `tpt-thermodynamics`.
//!
//! This crate is a scratch pad for runnable, documented example code that
//! exercises the `tpt-thermo-*` crates as each Phase 2-13 lands. It is **not
//! published** and is excluded from the packaged workspace surface.
//!
//! Examples are added per phase so that every public API has at least one
//! executable, tested showcase.

/// Phase 2 examples: composition conversion and the ideal-gas EoS toy.
pub mod phase2 {
    use tpt_thermo_core::{
        composition::Composition,
        eos::{EquationOfState, IdealGas},
        quantities::{MolarHeatCapacity, MolarMass, MolarVolume, Pressure, Temperature},
    };
    use uom::si::{
        molar_heat_capacity::joule_per_kelvin_mole, molar_mass::kilogram_per_mole,
        molar_volume::cubic_meter_per_mole, pressure::pascal as pascal_unit,
        thermodynamic_temperature::kelvin,
    };

    /// Convert a mole-fraction composition to mass fractions and back, using
    /// the `tpt-thermo-core` [`Composition`] helper.
    pub fn mole_to_mass_round_trip() -> (Vec<f64>, Vec<f64>) {
        // 50/50 mole methane (16 g/mol) / water-like (18 g/mol).
        let c = Composition::from_mole_fractions(vec![0.5, 0.5]).unwrap();
        let mm = [
            MolarMass::new::<kilogram_per_mole>(0.016),
            MolarMass::new::<kilogram_per_mole>(0.018),
        ];
        let mass = c.mass_fractions(&mm).unwrap();
        let back = Composition::from_mass_fractions(mass.clone(), &mm).unwrap();
        (mass, back.mole_fractions().to_vec())
    }

    /// Ideal-gas P-V-T: recover pressure from `P = R T / v` for a 300 K sample
    /// at `R·300/1e5` m³·mol⁻¹ (so the expected pressure is 1 bar).
    pub fn ideal_gas_pressure() -> f64 {
        let m = IdealGas::new(
            1,
            MolarHeatCapacity::new::<joule_per_kelvin_mole>(33.6),
            MolarMass::new::<kilogram_per_mole>(0.028),
            Temperature::new::<kelvin>(298.15),
            Pressure::new::<pascal_unit>(1.0e5),
        );
        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(tpt_thermo_core::R * 300.0 / 1.0e5);
        m.pressure(t, v, &[1.0]).unwrap().value
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn composition_round_trips() {
            let (mass, back) = mole_to_mass_round_trip();
            assert!((mass[0] - 0.016 / 0.034).abs() < 1e-12);
            assert!((back[0] - 0.5).abs() < 1e-12);
        }

        #[test]
        fn ideal_gas_matches_rt_over_v() {
            let p = ideal_gas_pressure();
            assert!((p - 1.0e5).abs() / 1.0e5 < 1e-9);
        }
    }
}

/// Phase 4 examples: cubic-equation-of-state P-V-T behaviour with Peng-Robinson
/// over the curated seed set.
pub mod phase4 {
    use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
    use tpt_thermo_data::SeedComponentDatabase;
    use tpt_thermo_eos_cubic::cubic_solver::Phase;
    use tpt_thermo_eos_cubic::PengRobinson;
    use uom::si::{
        molar_volume::cubic_meter_per_mole, pressure::pascal, thermodynamic_temperature::kelvin,
    };

    fn db() -> SeedComponentDatabase {
        SeedComponentDatabase::from_seed()
    }

    fn unit(i: usize) -> Vec<f64> {
        let n = db().num_components();
        let mut z = vec![0.0; n];
        z[i] = 1.0;
        z
    }

    /// A representative two-phase (subcritical) pressure for `component` at `t_k`:
    /// the midpoint of the pressure interval over which the cubic shows three
    /// real compressibility roots (liquid + unstable + vapor).
    fn two_phase_pressure(component: &str, t_k: f64) -> f64 {
        let eos = PengRobinson::from_database(&db()).unwrap();
        let idx = db().index_of(component).unwrap();
        let t = Temperature::new::<kelvin>(t_k);
        let z = unit(idx);
        let mut start = None;
        let mut end = None;
        for k in 1..2000 {
            let p = 1.0e3 * (k as f64);
            let three = eos
                .engine()
                .z_roots(t, Pressure::new::<pascal>(p), &z)
                .len()
                == 3;
            if three && start.is_none() {
                start = Some(p);
            }
            if start.is_some() && !three {
                end = Some(p);
                break;
            }
        }
        let start = start.expect("two-phase region exists below Pc");
        let end = end.unwrap_or(0.95 * db().critical_pressure(idx).unwrap().value);
        0.5 * (start + end)
    }

    /// Saturated liquid and vapor molar volumes (m³/mol) of `component` at `t_k`,
    /// obtained from the cubic's physically-meaningful root selection.
    pub fn saturated_volumes(component: &str, t_k: f64) -> (f64, f64) {
        let eos = PengRobinson::from_database(&db()).unwrap();
        let idx = db().index_of(component).unwrap();
        let t = Temperature::new::<kelvin>(t_k);
        let z = unit(idx);
        let p = Pressure::new::<pascal>(two_phase_pressure(component, t_k));
        let v_l = eos.solve_phase(t, p, &z, Phase::Liquid).unwrap().value;
        let v_v = eos.solve_phase(t, p, &z, Phase::Vapor).unwrap().value;
        (v_l, v_v)
    }

    /// Isothermal P-V-T round trip: `P(v_liquid)` recovers the pressure used to
    /// obtain `v_liquid`.
    pub fn pressure_round_trip(component: &str, t_k: f64) -> (f64, f64) {
        let eos = PengRobinson::from_database(&db()).unwrap();
        let idx = db().index_of(component).unwrap();
        let t = Temperature::new::<kelvin>(t_k);
        let z = unit(idx);
        let p = Pressure::new::<pascal>(two_phase_pressure(component, t_k));
        let v_l = eos.solve_phase(t, p, &z, Phase::Liquid).unwrap();
        let p_back = eos.pressure(t, v_l, &z).unwrap().value;
        (p.value, p_back)
    }

    /// Pure-component critical point `(Tc, Pc, vc)` via the engine's 2D Newton
    /// solver on the criticality conditions.
    pub fn critical(component: &str) -> (f64, f64, f64) {
        let eos = PengRobinson::from_database(&db()).unwrap();
        let idx = db().index_of(component).unwrap();
        let (tc, pc, vc) = eos.critical_point_pure(idx).unwrap();
        (tc.value, pc.value, vc.value)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn liquid_denser_than_vapor() {
            let (vl, vv) = saturated_volumes("methane", 150.0);
            assert!(vl < vv, "liquid {vl} must be denser than vapor {vv}");
        }

        #[test]
        fn pv_round_trip_close() {
            let (p, p_back) = pressure_round_trip("methane", 150.0);
            assert!((p - p_back).abs() / p < 1e-6, "P {p} vs recovered {p_back}");
        }

        #[test]
        fn critical_is_physical() {
            let (tc, pc, vc) = critical("methane");
            assert!(tc > 0.0 && pc > 0.0 && vc > 0.0);
        }
    }
}

/// Phase 9 examples: bubble/dew points, the phase envelope, azeotrope and
/// criconden detection, driven by the Peng-Robinson cubic EoS over the curated
/// seed set.
pub mod phase9 {
    use tpt_thermo_bubble_dew::{bubble_dew_envelope, BubbleDewSolver, KProvider};
    use tpt_thermo_core::quantities::Pressure;
    use tpt_thermo_data::SeedComponentDatabase;
    use tpt_thermo_eos_cubic::PengRobinson;
    use uom::si::pressure::atmosphere;

    fn benz_tol_db() -> SeedComponentDatabase {
        // Minimal 2-component seed for benzene/toluene (curated constants).
        let toml = r#"
[[components]]
schema_version = 1
name = "benzene"
critical_temperature_k = 562.05
critical_pressure_pa = 4894000.0
acentric_factor = 0.210
molar_mass_kg_per_mol = 0.07811184

[[components]]
schema_version = 1
name = "toluene"
critical_temperature_k = 591.79
critical_pressure_pa = 4108000.0
acentric_factor = 0.257
molar_mass_kg_per_mol = 0.09213842
"#;
        SeedComponentDatabase::from_toml_str(toml).unwrap()
    }

    /// Trace the benzene/toluene bubble and dew curves at 0.5–2 atm and return
    /// the number of points on each.
    pub fn benzene_toluene_envelope() -> (usize, usize) {
        let db = benz_tol_db();
        let eos = PengRobinson::from_database(&db).unwrap();
        let solver = BubbleDewSolver::new(&eos as &dyn KProvider, &db);
        let z = vec![0.5, 0.5];
        let pressures: Vec<Pressure> = (1..=15)
            .map(|i| Pressure::new::<atmosphere>(0.5 + 0.1 * i as f64))
            .collect();
        let env = bubble_dew_envelope(&solver, &z, &pressures).unwrap();
        (env.bubble.len(), env.dew.len())
    }

    /// Bubble and dew temperatures of an equimolar benzene/toluene liquid at 1 atm.
    pub fn benzene_toluene_bubble_dew_1atm() -> (f64, f64) {
        let db = benz_tol_db();
        let eos = PengRobinson::from_database(&db).unwrap();
        let solver = BubbleDewSolver::new(&eos as &dyn KProvider, &db);
        let x = vec![0.5, 0.5];
        let p = Pressure::new::<atmosphere>(1.0);
        let tb = solver
            .bubble_point_temperature(p, &x)
            .unwrap()
            .temperature
            .value;
        let td = solver
            .dew_point_temperature(p, &x)
            .unwrap()
            .temperature
            .value;
        (tb, td)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn envelope_has_points() {
            let (b, d) = benzene_toluene_envelope();
            assert!(b > 0 && d > 0);
        }

        #[test]
        fn bubble_below_dew() {
            let (tb, td) = benzene_toluene_bubble_dew_1atm();
            assert!(tb < td, "bubble {tb} must be below dew {td}");
            // Between the pure-component normal boiling points.
            assert!(tb > 353.0 && tb < 384.0);
        }
    }
}
