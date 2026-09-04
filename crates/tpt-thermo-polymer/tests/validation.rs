//! Integration validation for `tpt-thermo-polymer` against analytical references.

use tpt_thermo_core::quantities::{MolarEnergy, MolarVolume, Pressure, Temperature};
use tpt_thermo_core::EquationOfState;
use tpt_thermo_polymer::{
    chi_from_osmotic_pressure, cloud_point, critical_point, flory_melting_depression,
    most_probable, schulz_zimm, ChiTemperature, FloryHuggins, SanchezLacombe,
};
use uom::si::molar_energy::joule_per_mole;
use uom::si::molar_volume::cubic_meter_per_mole;
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;

#[test]
fn flory_huggins_cloud_point_matches_analytic() {
    // Symmetric-ish polymer/solvent: r1 = 1, r2 = 1000.
    let r1 = 1.0_f64;
    let r2 = 1000.0_f64;
    let phi_c = 1.0_f64 / (1.0 + f64::sqrt(r2 / r1));
    let chi_c = 0.5 * (1.0 / (r1 * phi_c) + 1.0 / (r2 * (1.0 - phi_c)));
    let cp = critical_point(
        r1,
        r2,
        &ChiTemperature::LinearDecreasing {
            a: 0.0,
            b: chi_c * 400.0,
        },
    )
    .unwrap();
    assert!((cp.critical_volume_fraction - phi_c).abs() < 1e-6);
    assert!((cp.temperature.unwrap() - 400.0).abs() < 1e-3);
    // Binodal at T = 300 K straddles φ_c.
    let c = ChiTemperature::LinearDecreasing {
        a: 0.0,
        b: chi_c * 400.0,
    }
    .at(300.0);
    let (lo, hi) = cloud_point::binodal(300.0, r1, r2, c).unwrap();
    assert!(lo < phi_c && hi > phi_c);
}

#[test]
fn flory_huggins_activity_coefficient_known_limit() {
    // r1 = r2 = 1 reduces to regular solution: ln γ₁ = χ x₂².
    let m = FloryHuggins::new_scalar(vec![1.0, 1.0], 0.8);
    let lng = m.ln_gamma(&[0.25, 0.75])[0];
    assert!((lng - 0.8 * 0.75 * 0.75).abs() < 1e-9);
}

#[test]
fn sanchez_lacombe_ideal_limit() {
    let m = SanchezLacombe::new(vec![1.0], vec![1e-5], vec![0.0]);
    let t = Temperature::new::<kelvin>(300.0);
    // Low density ⇒ P ≈ RT/v (ideal gas).
    let v = MolarVolume::new::<cubic_meter_per_mole>(0.5);
    let p = m.pressure(t, v, &[1.0]).unwrap();
    let expected = tpt_thermo_core::R * 300.0 / 0.5;
    assert!((p.value - expected).abs() / expected < 1e-3);
    // Fugacity coefficient ≈ 1 at low density.
    let lng = m.ln_fugacity_coefficient(t, v, &[1.0], 0).unwrap();
    assert!(lng.abs() < 1e-2);
}

#[test]
fn molecular_weight_distributions() {
    let mp = most_probable(20.0, 400);
    assert!((mp.number_average() - 20.0).abs() < 1e-4);
    assert!((mp.dispersity() - (2.0 - 1.0 / 20.0)).abs() < 1e-3);
    let sz = schulz_zimm(100.0, 1.4, 800);
    assert!((sz.number_average() - 100.0).abs() < 1e-3);
    assert!((sz.dispersity() - 1.4).abs() < 1e-2);
}

#[test]
fn crystallization_depression_direction() {
    let tm0 = Temperature::new::<kelvin>(450.0);
    let dh = MolarEnergy::new::<joule_per_mole>(1.2e7);
    let tm_dilute = flory_melting_depression(tm0, dh, 0.5, 1.0e6);
    let tm_pure = flory_melting_depression(tm0, dh, 1.0 - 1e-9, 1.0e6);
    assert!(tm_dilute.value < tm_pure.value);
    assert!((tm_pure.value - 450.0).abs() < 1e-3);
}

#[test]
fn osmotic_pressure_recovers_chi() {
    let phi = 0.05_f64;
    let r2 = 800.0_f64;
    let chi = 0.3_f64;
    let t = Temperature::new::<kelvin>(300.0);
    let vstar = MolarVolume::new::<cubic_meter_per_mole>(1.0e-5);
    let lhs = -(1.0 / phi) * ((1.0 - phi).ln() + (1.0 - 1.0 / r2) * phi + chi * phi * phi);
    let pi = Pressure::new::<pascal>(lhs * tpt_thermo_core::R * t.value / vstar.value);
    let recovered = chi_from_osmotic_pressure(phi, pi, t, vstar, r2).unwrap();
    assert!((recovered - chi).abs() < 1e-9);
}
