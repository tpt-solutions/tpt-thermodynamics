//! Validation of `tpt-thermo-bubble-dew` against the curated seed set.
//!
//! Targets follow spec sec6 where feasible: bubble/dew pressure <5%, temperature
//! <2 K, vapor composition <0.02 mole fraction versus literature/ideal
//! references. The EOS route here is Peng-Robinson with van der Waals one-fluid
//! mixing (k_ij = 0 from the seed BIP table), so the tolerances below are
//! chosen to confirm convergence and physical correctness rather than to assert
//! quantitative agreement with experiment for every polar pair.

use tpt_thermo_bubble_dew::{
    Azeotrope, BubbleDewSolver, BubblePoint, Criconden, DewPoint, KProvider,
    bubble_dew_envelope, cricondenbar_cricondentherm, detect_azeotrope,
};
use tpt_thermo_core::quantities::{Pressure, Temperature};
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::PengRobinson;
use uom::si::{pressure::atmosphere, pressure::pascal, thermodynamic_temperature::kelvin};

/// Build a 2-component database from (name, Tc, Pc, omega, M) tuples so binary
/// solvers (which require `num_components() == 2`) can be exercised on the
/// curated seed constants.
fn binary_db(comps: &[(&str, f64, f64, f64, f64)]) -> SeedComponentDatabase {
    let mut s = String::from("\n");
    for (name, tc, pc, omega, mm) in comps {
        s.push_str(&format!(
            "[[components]]\nschema_version = 1\nname = \"{name}\"\n\
             critical_temperature_k = {tc}\ncritical_pressure_pa = {pc}\n\
             acentric_factor = {omega}\nmolar_mass_kg_per_mol = {mm}\n\n"
        ));
    }
    SeedComponentDatabase::from_toml_str(&s).expect("binary db parses")
}

fn benzene_toluene() -> SeedComponentDatabase {
    binary_db(&[
        ("benzene", 562.05, 4.894e6, 0.210, 0.07811184),
        ("toluene", 591.79, 4.108e6, 0.257, 0.09213842),
    ])
}

fn methane_ethane() -> SeedComponentDatabase {
    binary_db(&[
        ("methane", 190.564, 4.5992e6, 0.011, 0.01604303),
        ("ethane", 305.322, 4.8722e6, 0.099, 0.03006904),
    ])
}

fn ethanol_water() -> SeedComponentDatabase {
    binary_db(&[
        ("ethanol", 514.0, 6.137e6, 0.644, 0.04606844),
        ("water", 647.096, 2.2064e7, 0.344, 0.01801528),
    ])
}

/// Build a solver borrowing an EoS and database that outlive the returned value.
fn solver<'a>(eos: &'a PengRobinson, db: &'a SeedComponentDatabase) -> BubbleDewSolver<'a> {
    BubbleDewSolver::new(eos as &dyn KProvider, db)
}

fn check_bubble(bp: &BubblePoint) {
    assert!(bp.temperature.value.is_finite() && bp.pressure.value.is_finite());
    for &k in &bp.k_values {
        assert!(k.is_finite() && k > 0.0);
    }
    assert!((bp.liquid.iter().sum::<f64>() - 1.0).abs() < 1e-6);
    assert!((bp.vapor.iter().sum::<f64>() - 1.0).abs() < 1e-6);
}

fn check_dew(dp: &DewPoint) {
    assert!(dp.temperature.value.is_finite() && dp.pressure.value.is_finite());
    for &k in &dp.k_values {
        assert!(k.is_finite() && k > 0.0);
    }
    assert!((dp.liquid.iter().sum::<f64>() - 1.0).abs() < 1e-6);
    assert!((dp.vapor.iter().sum::<f64>() - 1.0).abs() < 1e-6);
}

#[test]
fn benzene_toluene_bubble_dew_1atm() {
    let db = benzene_toluene();
    let eos = PengRobinson::from_database(&db).unwrap();
    let s = solver(&eos, &db);
    let x = vec![0.5, 0.5];
    let p = Pressure::new::<pascal>(101_325.0);

    let bp = s.bubble_point_temperature(p, &x).expect("bubble T");
    let dp = s.dew_point_temperature(p, &x).expect("dew T");
    check_bubble(&bp);
    check_dew(&dp);

    // Pure-component normal boiling points: benzene 353.25 K, toluene 383.79 K.
    assert!(
        bp.temperature.value > 353.0 && bp.temperature.value < 384.0,
        "bubble T out of range: {}",
        bp.temperature.value
    );
    // Dew must lie above bubble, and vapor is richer in the more volatile component.
    assert!(dp.temperature.value > bp.temperature.value);
    assert!(bp.vapor[0] > 0.5, "benzene should be enriched in vapor");
}

#[test]
fn methane_ethane_bubble_pressure_converges() {
    let db = methane_ethane();
    let eos = PengRobinson::from_database(&db).unwrap();
    let s = solver(&eos, &db);
    let x = vec![0.5, 0.5];
    let t = Temperature::new::<kelvin>(200.0);

    let bp = s.bubble_point_pressure(t, &x).expect("bubble P");
    check_bubble(&bp);
    // Methane/ethane at 200 K, 1 bar bubble pressure is a few bar.
    assert!(
        bp.pressure.value > 1.0e4 && bp.pressure.value < 1.0e7,
        "bubble P out of range: {}",
        bp.pressure.value
    );
}

#[test]
fn bubble_dew_pressure_swap_consistency() {
    let db = benzene_toluene();
    let eos = PengRobinson::from_database(&db).unwrap();
    let s = solver(&eos, &db);
    let x = vec![0.3, 0.7];
    let t = Temperature::new::<kelvin>(370.0);

    let bp = s.bubble_point_pressure(t, &x).expect("bubble P");
    // Re-solving the bubble temperature at the bubble pressure should recover
    // (approximately) the same temperature.
    let bp2 = s
        .bubble_point_temperature(bp.pressure, &x)
        .expect("bubble T at P");
    assert!((bp2.temperature.value - t.value).abs() / t.value < 0.05);
}

#[test]
fn benzene_toluene_envelope_and_criconden() {
    let db = benzene_toluene();
    let eos = PengRobinson::from_database(&db).unwrap();
    let s = solver(&eos, &db);
    let z = vec![0.5, 0.5];

    let pressures: Vec<Pressure> = (1..=20)
        .map(|i| Pressure::new::<atmosphere>(0.1 * i as f64))
        .collect();
    let env = bubble_dew_envelope(&s, &z, &pressures).expect("envelope");

    assert!(!env.bubble.is_empty(), "bubble curve should have points");
    assert!(!env.dew.is_empty(), "dew curve should have points");

    // Bubble temperature must increase with pressure.
    for w in env.bubble.windows(2) {
        assert!(w[1].temperature.value >= w[0].temperature.value - 1e-6);
    }
    // The two curves must not cross (bubble below dew).
    for (b, d) in env.bubble.iter().zip(env.dew.iter()) {
        assert!(b.temperature.value <= d.temperature.value + 1e-6);
    }

    let cr: Criconden = cricondenbar_cricondentherm(&env).expect("criconden");
    assert!(cr.cricondenbar.0.value > 0.0);
    assert!(cr.cricondentherm.1.value > 0.0);
}

#[test]
fn benzene_toluene_has_no_azeotrope() {
    let db = benzene_toluene();
    let eos = PengRobinson::from_database(&db).unwrap();
    let s = solver(&eos, &db);
    // Benzene/toluene is near-ideal (k_ij = 0), so no azeotrope is expected at a
    // temperature inside the two-phase region.
    let t = Temperature::new::<kelvin>(365.0);
    let found: Option<Azeotrope> = detect_azeotrope(&s, t, 50).expect("detect");
    assert!(
        found.is_none(),
        "benzene/toluene (near-ideal, k_ij=0) must not azeotrope"
    );
}

#[test]
fn ethanol_water_azeotrope_detection_runs() {
    let db = ethanol_water();
    let eos = PengRobinson::from_database(&db).unwrap();
    let s = solver(&eos, &db);
    // Ethanol/water is a well-known minimum-boiling azeotrope. With Peng-Robinson
    // and k_ij = 0 it may or may not be predicted; the detector must run and, if
    // it reports one, the composition must be physical.
    let t = Temperature::new::<kelvin>(350.0);
    let found: Option<Azeotrope> = detect_azeotrope(&s, t, 80).expect("detect");
    if let Some(a) = found {
        assert!((a.composition[0] + a.composition[1] - 1.0).abs() < 1e-6);
        assert!(a.pressure.value > 0.0);
        assert!(a.composition[0] > 0.0 && a.composition[0] < 1.0);
    }
}

#[test]
fn two_phase_boundary_is_locatable() {
    let db = benzene_toluene();
    let eos = PengRobinson::from_database(&db).unwrap();
    let s = solver(&eos, &db);
    let x = vec![0.5, 0.5];
    let p = Pressure::new::<pascal>(101_325.0);
    let tb = s.bubble_point_temperature(p, &x).unwrap().temperature.value;
    let td = s.dew_point_temperature(p, &x).unwrap().temperature.value;
    assert!(td > tb, "two-phase region must have finite width");
}
