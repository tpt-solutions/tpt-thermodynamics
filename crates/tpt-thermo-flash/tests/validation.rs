//! Validation tests for the flash solvers against the cubic EoS + seed dataset.

use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use tpt_thermo_core::EquationOfState;
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::PengRobinson;
use tpt_thermo_flash::{flash_pt, FlashCalculator};
use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};

fn methane_ethane() -> (SeedComponentDatabase, PengRobinson, Vec<f64>) {
    let db = SeedComponentDatabase::from_seed();
    let eos = PengRobinson::from_database(&db).unwrap();
    let methane = db.index_of("methane").unwrap();
    let ethane = db.index_of("ethane").unwrap();
    let mut z = vec![0.0_f64; db.num_components()];
    z[methane] = 0.7;
    z[ethane] = 0.3;
    (db, eos, z)
}

#[test]
fn pt_flash_two_phase_split() {
    let (db, eos, z) = methane_ethane();
    let calc = FlashCalculator::with_db(&eos, &db);
    // 200 K, 3 MPa: both components are below/near their criticals for this
    // mixture, giving a genuine vapor–liquid split.
    let t = Temperature::new::<kelvin>(200.0);
    let p = Pressure::new::<pascal>(3.0e6);
    let res = calc.flash_pt(t, p, &z).unwrap();
    assert!(res.converged, "PT flash should converge");
    assert_eq!(res.phase_flag, tpt_thermo_flash::pt::PhaseFlag::TwoPhase);
    assert!((0.0..=1.0).contains(&res.vapor_fraction));
    // Vapor is richer in the lighter component (methane).
    let methane = db.index_of("methane").unwrap();
    assert!(
        res.vapor_composition[methane] > res.liquid_composition[methane],
        "methane should enrich the vapor"
    );
}

#[test]
fn pt_flash_single_phase_at_high_t() {
    let (db, eos, z) = methane_ethane();
    let calc = FlashCalculator::with_db(&eos, &db);
    // Well above the criconden: single vapor phase.
    let t = Temperature::new::<kelvin>(450.0);
    let p = Pressure::new::<pascal>(1.0e6);
    let res = calc.flash_pt(t, p, &z).unwrap();
    assert!(res.converged);
    assert_eq!(res.phase_flag, tpt_thermo_flash::pt::PhaseFlag::SinglePhase);
}

#[test]
fn free_function_matches_calculator() {
    let (db, eos, z) = methane_ethane();
    let t = Temperature::new::<kelvin>(300.0);
    let p = Pressure::new::<pascal>(3.0e6);
    let r1 = flash_pt(&eos, Some(&db), t, p, &z).unwrap();
    let calc = FlashCalculator::with_db(&eos, &db);
    let r2 = calc.flash_pt(t, p, &z).unwrap();
    assert!((r1.vapor_fraction - r2.vapor_fraction).abs() < 1e-9);
}

#[test]
fn ph_flash_matches_target_enthalpy() {
    let (db, eos, z) = methane_ethane();
    let calc = FlashCalculator::with_db(&eos, &db);
    // First do a PT flash to get a target enthalpy.
    let t0 = Temperature::new::<kelvin>(300.0);
    let p = Pressure::new::<pascal>(2.0e6);
    let pt = calc.flash_pt(t0, p, &z).unwrap();
    let h: f64 = (1.0 - pt.vapor_fraction)
        * eos
            .molar_enthalpy(t0, pt.liquid_volume, &pt.liquid_composition)
            .unwrap()
            .value
        + pt.vapor_fraction
            * eos
                .molar_enthalpy(t0, pt.vapor_volume, &pt.vapor_composition)
                .unwrap()
                .value;
    let h_target =
        tpt_thermo_core::quantities::MolarEnergy::new::<uom::si::molar_energy::joule_per_mole>(h);
    let ph = calc.flash_ph(h_target, p, &z).unwrap();
    // Recompute enthalpy of the PH result and confirm it matches the target.
    let h2: f64 = (1.0 - ph.vapor_fraction)
        * eos
            .molar_enthalpy(t0, ph.liquid_volume, &ph.liquid_composition)
            .unwrap()
            .value
        + ph.vapor_fraction
            * eos
                .molar_enthalpy(t0, ph.vapor_volume, &ph.vapor_composition)
                .unwrap()
                .value;
    assert!(
        (h2 - h).abs() / h.abs() < 1e-3,
        "PH flash enthalpy should match target"
    );
}

#[test]
fn lle_isoactivity_splits_nonideal_binary() {
    use tpt_thermo_eos_activity::{parameters::TdParam, Nrtl};
    use tpt_thermo_flash::lle::lle_isoactivity;
    use uom::si::pressure::pascal;

    // Strongly non-ideal (partially-miscible) binary: large asymmetric tau terms
    // force a liquid–liquid split.
    let nrtl = Nrtl::binary(
        TdParam::new(6.0, 0.0, 0.0),
        TdParam::new(-3.0, 0.0, 0.0),
        0.3,
    )
    .unwrap();
    let t = Temperature::new::<kelvin>(298.15);
    let p = Pressure::new::<pascal>(1.0e5);
    let z = vec![0.5_f64, 0.5];
    let res = lle_isoactivity(&nrtl, t, p, &z).unwrap();
    assert!(res.converged, "LLE should converge");
    // The two co-existing phases must differ (otherwise no split occurred).
    let mut max_diff = 0.0_f64;
    for (a, b) in res.phase_i.iter().zip(res.phase_ii.iter()) {
        max_diff = max_diff.max((a - b).abs());
    }
    assert!(max_diff > 1e-3, "LLE should produce two distinct phases");
}

#[test]
fn parallel_batch_matches_serial() {
    use tpt_thermo_flash::flash_pt_batch;
    let (db, eos, _z) = methane_ethane();
    let methane = db.index_of("methane").unwrap();
    let ethane = db.index_of("ethane").unwrap();
    // Build a feed table sweeping the methane mole fraction.
    let mut feeds = Vec::new();
    for k in 0..12 {
        let mut zk = vec![0.0_f64; db.num_components()];
        let x = 0.1 + 0.07 * k as f64;
        zk[methane] = x;
        zk[ethane] = 1.0 - x;
        feeds.push(zk);
    }
    let t = Temperature::new::<kelvin>(220.0);
    let p = Pressure::new::<pascal>(3.0e6);
    let serial = flash_pt_batch(&eos, Some(&db), t, p, &feeds).unwrap();
    let parallel =
        tpt_thermo_flash::flash_pt_batch_parallel(&eos, Some(&db), t, p, &feeds).unwrap();
    assert_eq!(serial.len(), parallel.len());
    for (s, par) in serial.iter().zip(parallel.iter()) {
        assert!(
            (s.vapor_fraction - par.vapor_fraction).abs() < 1e-9,
            "batch mismatch"
        );
    }
}

/// Spec sec6 breadth expansion (seed of the 20+ multicomponent flash target).
///
/// A natural-gas-like five-component mixture (methane/ethane/propane/n-butane/
/// n-pentane) at a temperature/pressure inside its two-phase envelope must flash
/// to a converged VLE split whose phase compositions close the overall material
/// balance `z = (1−β)·x + β·y`, and whose vapor is enriched in the light
/// component. This is one representative multicomponent system; the full breadth
/// set is a mechanical extension of the same harness.
#[test]
fn multicomponent_flash_material_balance() {
    let full = SeedComponentDatabase::from_seed();
    let comps = ["methane", "ethane", "propane", "n-butane", "n-pentane"];
    let fracs = [0.6_f64, 0.2, 0.1, 0.07, 0.03];
    let mut z = vec![0.0_f64; full.num_components()];
    for (c, f) in comps.iter().zip(fracs.iter()) {
        z[full.index_of(c).unwrap()] = *f;
    }
    let eos = PengRobinson::from_database(&full).unwrap();
    let calc = FlashCalculator::with_db(&eos, &full);

    // 250 K, 5 MPa: inside the VL envelope for this mixture.
    let t = Temperature::new::<kelvin>(250.0);
    let p = Pressure::new::<pascal>(5.0e6);
    let res = calc.flash_pt(t, p, &z).unwrap();
    assert!(res.converged, "multicomponent flash should converge");
    assert_eq!(res.phase_flag, tpt_thermo_flash::pt::PhaseFlag::TwoPhase);
    assert!(
        (0.0..1.0).contains(&res.vapor_fraction),
        "vapor fraction must be in (0, 1)"
    );

    // Component-wise material balance closure.
    let beta = res.vapor_fraction;
    let max_err = z
        .iter()
        .enumerate()
        .filter(|(_, &zi)| zi > 0.0)
        .map(|(i, &zi)| {
            let recon = (1.0 - beta) * res.liquid_composition[i] + beta * res.vapor_composition[i];
            (recon - zi).abs()
        })
        .fold(0.0_f64, f64::max);
    assert!(
        max_err < 1.0e-6,
        "material balance closure error {max_err:.2e} exceeds tolerance"
    );

    // Vapor must be enriched in the lightest component (methane).
    let methane = full.index_of("methane").unwrap();
    assert!(
        res.vapor_composition[methane] > res.liquid_composition[methane],
        "methane should enrich the vapor"
    );
}
