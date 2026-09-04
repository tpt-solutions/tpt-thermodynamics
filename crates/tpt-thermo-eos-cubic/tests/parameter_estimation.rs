//! Parameter-estimation (spec 3d) regression tests.

use tpt_thermo_core::quantities::Temperature;
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::mixing::VdwMixing;
use tpt_thermo_eos_cubic::{bubble_pressure, fit_binary_kij, PengRobinson};
use uom::si::thermodynamic_temperature::kelvin;

#[test]
fn bubble_pressure_finite_for_propane_nbutane() {
    let full = SeedComponentDatabase::from_seed();
    let propane = full.index_of("propane").unwrap();
    let nbutane = full.index_of("n-butane").unwrap();
    let db = full.subset(&[propane, nbutane]).unwrap();
    let eos = PengRobinson::from_database(&db).unwrap();
    let t = Temperature::new::<kelvin>(350.0);
    let x = vec![0.5_f64, 0.5];
    let pb = bubble_pressure(&eos, t, &x).unwrap();
    // Expect a plausible subcritical bubble pressure (monotone in T for this pair).
    assert!(
        pb.value > 1.0e4 && pb.value < 1.0e7,
        "bubble P = {}",
        pb.value
    );
}

#[test]
fn bubble_pressure_monotone_in_temperature() {
    let full = SeedComponentDatabase::from_seed();
    let propane = full.index_of("propane").unwrap();
    let nbutane = full.index_of("n-butane").unwrap();
    let db = full.subset(&[propane, nbutane]).unwrap();
    let eos = PengRobinson::from_database(&db).unwrap();
    let mut prev = 0.0_f64;
    for &t_k in &[300.0_f64, 320.0, 340.0, 350.0] {
        let t = Temperature::new::<kelvin>(t_k);
        let x = vec![0.5_f64, 0.5];
        let pb = bubble_pressure(&eos, t, &x).unwrap();
        assert!(pb.value > prev, "bubble P should increase with T");
        prev = pb.value;
    }
}

/// Validate the full fit pipeline: generate synthetic isothermal bubble-pressure
/// VLE data with a *known* `k_ij` via an explicit mixing rule, then recover a
/// `k_ij` by least squares and confirm the fitted model reproduces the data.
///
/// Propane/n-butane is a well-conditioned non-associating pair for which the
/// bubble-pressure solver converges; for near-ideal binaries the objective is
/// mildly flat, so we assert data reproducibility (within 10 kPa) rather than
/// exact recovery of the input `k_ij`.
#[test]
fn fit_binary_kij_reproduces_bubble_pressures() {
    let full = SeedComponentDatabase::from_seed();
    let propane = full.index_of("propane").unwrap();
    let nbutane = full.index_of("n-butane").unwrap();
    let db = full.subset(&[propane, nbutane]).unwrap();

    let true_k = 0.05_f64;
    let gen_eos = PengRobinson::with_mixing(
        &db,
        Box::new(VdwMixing::from_matrix(vec![
            vec![0.0, true_k],
            vec![true_k, 0.0],
        ])),
    )
    .unwrap();
    let mut points = Vec::new();
    for &t_k in &[320.0_f64, 340.0, 350.0] {
        let t = Temperature::new::<kelvin>(t_k);
        let x = vec![0.5_f64, 0.5];
        let pb = bubble_pressure(&gen_eos, t, &x).unwrap();
        points.push((t, pb, 0.5));
    }

    let (k_fit, rms) = fit_binary_kij(&db, &points).unwrap();
    assert!(rms < 1.0e4, "rms residual {rms} Pa should be small");

    let fit_eos = PengRobinson::with_mixing(
        &db,
        Box::new(VdwMixing::from_matrix(vec![
            vec![0.0, k_fit],
            vec![k_fit, 0.0],
        ])),
    )
    .unwrap();
    for (t, pb, x1) in &points {
        let x = vec![*x1, 1.0 - *x1];
        let pb_fit = bubble_pressure(&fit_eos, *t, &x).unwrap();
        assert!(
            (pb_fit.value - pb.value).abs() < 1.0e4,
            "fitted k={k_fit} reproduces Pb within 10 kPa (got {}, expected {})",
            pb_fit.value,
            pb.value
        );
    }
}

/// The flash-based bubble solver must converge (not return `NotConverged`) for
/// associating / strongly non-ideal binaries that the previous plain-SS residual
/// could not bracket. We assert only that a finite, physically-plausible bubble
/// pressure is returned — PR/SRK with zero `k_ij` is not quantitatively accurate
/// for water, so this is a robustness test, not an accuracy test.
///
/// Note: binaries in which a component is *supercritical* at the test temperature
/// (e.g. water/methane, CO₂/methane) still require a stability-tested flash and
/// are a tracked follow-up; these subcritical–subcritical associating pairs are
/// the cases the flash-based incipient solver now handles robustly.
#[test]
fn bubble_pressure_converges_for_associating_binaries() {
    let full = SeedComponentDatabase::from_seed();

    let cases: &[(&str, &str, f64, f64)] = &[
        // (a, b, x_a, T [K]) — both components subcritical at T
        ("water", "ethanol", 0.5, 350.0),
        ("water", "methanol", 0.5, 350.0),
        ("ethanol", "water", 0.5, 373.15),
        ("methanol", "ethanol", 0.5, 350.0),
        ("ethanol", "benzene", 0.5, 350.0),
    ];
    for &(a, b, xa, t_k) in cases {
        let ia = full.index_of(a).expect(a);
        let ib = full.index_of(b).expect(b);
        let db = full.subset(&[ia, ib]).unwrap();
        let eos = PengRobinson::from_database(&db).unwrap();
        let t = Temperature::new::<kelvin>(t_k);
        let x = vec![xa, 1.0 - xa];
        let pb = bubble_pressure(&eos, t, &x)
            .unwrap_or_else(|_| panic!("bubble_pressure diverged for {a}/{b} @ {t_k} K"));
        assert!(
            pb.value.is_finite() && pb.value > 1.0e2 && pb.value < 1.0e8,
            "bubble P for {a}/{b} = {} Pa out of plausible band",
            pb.value
        );
    }
}

/// Near-critical subcritical binaries (both components below their critical
/// temperature, but at a temperature close to a component's critical point)
/// exercise the flash-based solver in the regime where the old residual was
/// non-smooth. These converge with the GDEM-accelerated incipient-phase solve.
#[test]
fn bubble_pressure_converges_for_near_critical_binaries() {
    let full = SeedComponentDatabase::from_seed();
    let pairs: &[(&str, &str, f64)] = &[
        // (a, b, T [K]) — T within ~10% of a component's critical temperature
        ("ethane", "propane", 300.0),
        ("carbon dioxide", "ethane", 280.0),
        ("benzene", "toluene", 400.0),
    ];
    for &(a, b, t_k) in pairs {
        let ia = full.index_of(a).expect(a);
        let ib = full.index_of(b).expect(b);
        let db = full.subset(&[ia, ib]).unwrap();
        let eos = PengRobinson::from_database(&db).unwrap();
        let t = Temperature::new::<kelvin>(t_k);
        let x = vec![0.5_f64, 0.5];
        let pb = bubble_pressure(&eos, t, &x)
            .unwrap_or_else(|_| panic!("bubble_pressure diverged for {a}/{b} @ {t_k} K"));
        assert!(
            pb.value.is_finite() && pb.value > 1.0e2 && pb.value < 1.0e8,
            "bubble P for {a}/{b} @ {t_k} K = {} Pa out of plausible band",
            pb.value
        );
    }
}

/// End-to-end self-consistency for an associating binary: generate synthetic
/// isothermal bubble-pressure VLE with a *known* `k_ij`, recover it by least
/// squares, and confirm the fitted model reproduces the data. This exercises the
/// full flash-based `bubble_pressure` + `fit_binary_kij` pipeline on a
/// strongly non-ideal pair (ethanol/water), which previously did not converge.
#[test]
fn fit_binary_kij_roundtrip_water_ethanol() {
    let full = SeedComponentDatabase::from_seed();
    let iw = full.index_of("water").unwrap();
    let ie = full.index_of("ethanol").unwrap();
    let db = full.subset(&[iw, ie]).unwrap();

    let true_k = -0.05_f64;
    let gen_eos = PengRobinson::with_mixing(
        &db,
        Box::new(VdwMixing::from_matrix(vec![
            vec![0.0, true_k],
            vec![true_k, 0.0],
        ])),
    )
    .unwrap();
    let mut points = Vec::new();
    for &t_k in &[340.0_f64, 350.0, 360.0] {
        let t = Temperature::new::<kelvin>(t_k);
        let x = vec![0.5_f64, 0.5];
        let pb = bubble_pressure(&gen_eos, t, &x).unwrap();
        points.push((t, pb, 0.5));
    }

    let (k_fit, rms) = fit_binary_kij(&db, &points).unwrap();
    assert!(rms < 1.0e4, "rms residual {rms} Pa should be small");

    let fit_eos = PengRobinson::with_mixing(
        &db,
        Box::new(VdwMixing::from_matrix(vec![
            vec![0.0, k_fit],
            vec![k_fit, 0.0],
        ])),
    )
    .unwrap();
    for (t, pb, x1) in &points {
        let x = vec![*x1, 1.0 - *x1];
        let pb_fit = bubble_pressure(&fit_eos, *t, &x).unwrap();
        assert!(
            (pb_fit.value - pb.value).abs() < 1.0e4,
            "fitted k={k_fit} reproduces Pb within 10 kPa (got {}, expected {})",
            pb_fit.value,
            pb.value
        );
    }
}
