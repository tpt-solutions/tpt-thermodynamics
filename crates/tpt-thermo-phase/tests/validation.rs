//! Phase-8 validation against the curated seed set (spec sec6, seed subset).
//!
//! These are integration tests exercising the stability test, multiphase
//! classification, mixture critical point, and continuation end-to-end on the
//! Peng-Robinson EoS over the seed dataset.

use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use tpt_thermo_core::StabilityTest;
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::PengRobinson;
use tpt_thermo_phase::phase_volume::PhaseVolume;
use tpt_thermo_phase::{
    critical_locus_binary, detect_phases, mixture_critical_point, CriticalGuess, StabilityAnalyzer,
};
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;

/// Build a `z` vector of `db.num_components()` with the given `(index, fraction)`
/// pairs.
fn composition(db: &SeedComponentDatabase, fracs: &[(usize, f64)]) -> Vec<f64> {
    let mut z = vec![0.0_f64; db.num_components()];
    for (i, f) in fracs {
        z[*i] = *f;
    }
    z
}

#[test]
fn pure_water_superheated_is_stable() {
    let db = SeedComponentDatabase::from_seed();
    let eos = PengRobinson::from_database(&db).unwrap();
    let water = db.index_of("water").unwrap();
    let z = composition(&db, &[(water, 1.0)]);
    let vol = &eos as &dyn PhaseVolume;
    let ana = StabilityAnalyzer::new(&eos, vol, &db);
    let t = Temperature::new::<kelvin>(600.0);
    let p = Pressure::new::<pascal>(1.0e5);
    let res = ana.test(t, p, &z).unwrap();
    assert!(res.stable, "superheated pure water must be stable");
}

#[test]
fn methane_butane_binary_is_two_phase() {
    let db = SeedComponentDatabase::from_seed();
    let eos = PengRobinson::from_database(&db).unwrap();
    let methane = db.index_of("methane").unwrap();
    let nbutane = db.index_of("n-butane").unwrap();
    let z = composition(&db, &[(methane, 0.5), (nbutane, 0.5)]);
    let vol = &eos as &dyn PhaseVolume;
    let ana = StabilityAnalyzer::new(&eos, vol, &db);
    // Inside the two-phase envelope: a light/heavy binary at moderate T, P.
    let t = Temperature::new::<kelvin>(300.0);
    let p = Pressure::new::<pascal>(2.0e6);
    let res = ana.test(t, p, &z).unwrap();
    assert!(!res.stable, "methane/n-butane at 300K, 2MPa should be unstable (two-phase)");
}

#[test]
fn multiphase_classification_runs() {
    let db = SeedComponentDatabase::from_seed();
    let eos = PengRobinson::from_database(&db).unwrap();
    let methane = db.index_of("methane").unwrap();
    let nbutane = db.index_of("n-butane").unwrap();
    let z = composition(&db, &[(methane, 0.5), (nbutane, 0.5)]);
    let vol = &eos as &dyn PhaseVolume;
    let r = detect_phases(&eos, vol, &db, Temperature::new::<kelvin>(300.0), Pressure::new::<pascal>(2.0e6), &z);
    assert!(r.num_phases >= 1);
    assert!(!r.stable);
}

#[test]
fn mixture_critical_point_within_pure_bounds() {
    let db = SeedComponentDatabase::from_seed();
    let eos = PengRobinson::from_database(&db).unwrap();
    let methane = db.index_of("methane").unwrap();
    let ethane = db.index_of("ethane").unwrap();
    let z = composition(&db, &[(methane, 0.5), (ethane, 0.5)]);
    let guess = CriticalGuess::from_database(&db, &z);
    let (tc, pc, _v) = mixture_critical_point(&eos, &z, guess).unwrap();
    let tcm = db.critical_temperature(methane).unwrap().value;
    let tce = db.critical_temperature(ethane).unwrap().value;
    let pcm = db.critical_pressure(methane).unwrap().value;
    let pce = db.critical_pressure(ethane).unwrap().value;
    assert!(tc.value > tcm.min(tce) - 2.0 && tc.value < tce.max(tcm) + 2.0);
    assert!(pc.value > pcm.min(pce) * 0.5 && pc.value < pce.max(pcm) * 1.5);
}

#[test]
fn critical_locus_traces_interior() {
    let db = SeedComponentDatabase::from_seed();
    let eos = PengRobinson::from_database(&db).unwrap();
    let locus = critical_locus_binary(&eos, &db, 8);
    // Interior (non-pure) points converge; pure-component endpoints are a
    // Jacobian singularity of the raw (v,T) system and are intentionally taken
    // from the database / engine instead.
    assert!(locus.len() >= 3, "expect interior locus points to converge");
    let tcs: Vec<f64> = locus.iter().map(|(_, tc, _)| tc.value).collect();
    let tcm = db.critical_temperature(db.index_of("methane").unwrap()).unwrap().value;
    let tce = db.critical_temperature(db.index_of("ethane").unwrap()).unwrap().value;
    // Critical temperature must stay bracketed by the pure-component values.
    assert!(
        *tcs.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap() > tcm - 2.0
            && *tcs.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap() < tce + 2.0
    );
}

#[test]
fn excess_gibbs_is_finite_for_feed() {
    let db = SeedComponentDatabase::from_seed();
    let eos = PengRobinson::from_database(&db).unwrap();
    let methane = db.index_of("methane").unwrap();
    let z = composition(&db, &[(methane, 1.0)]);
    let vol = &eos as &dyn PhaseVolume;
    let ana = StabilityAnalyzer::new(&eos, vol, &db);
    let g = ana
        .excess_gibbs(
            Temperature::new::<kelvin>(300.0),
            Pressure::new::<pascal>(1.0e6),
            &z,
        )
        .unwrap();
    assert!(g.value.is_finite());
}
