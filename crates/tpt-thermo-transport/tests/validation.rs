//! Validation tests for the transport correlations against the seed dataset.

use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use tpt_thermo_data::SeedComponentDatabase;
use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};

#[test]
fn chung_gas_viscosity_nitrogen() {
    use tpt_thermo_core::quantities::DynamicViscosity;
    use uom::si::dynamic_viscosity::pascal_second;
    let db = SeedComponentDatabase::from_seed();
    let n2 = db.index_of("nitrogen").unwrap();
    let mut z = vec![0.0_f64; db.num_components()];
    z[n2] = 1.0;
    let t = Temperature::new::<kelvin>(300.0);
    let rho = tpt_thermo_transport::conductivity::ideal_molar_density(
        t,
        Pressure::new::<pascal>(1.0e5),
    );
    let eta = tpt_thermo_transport::viscosity::chung_gas_viscosity(&db, t, rho, &z).unwrap();
    let v = eta.value; // Pa·s
    // N2 at 300 K, 1 atm: ≈ 1.78e-5 Pa·s (Chung is within ~10%).
    assert!(
        (1.4e-5..2.2e-5).contains(&v),
        "N2 viscosity out of range: {v:e}"
    );
    let _ = DynamicViscosity::new::<pascal_second>(v);
}

#[test]
fn fuller_diffusivity_n2_co2() {
    use uom::si::diffusion_coefficient::square_meter_per_second;
    let db = SeedComponentDatabase::from_seed();
    let n2 = db.index_of("nitrogen").unwrap();
    let co2 = db.index_of("carbon dioxide").unwrap();
    let t = Temperature::new::<kelvin>(298.15);
    let p = Pressure::new::<pascal>(1.0e5);
    let d = tpt_thermo_transport::diffusivity::fuller_schettler_giddings(&db, t, p, n2, co2)
        .unwrap();
    let v = d.value; // m²/s
    // N2–CO2 at 298 K, 1 atm: ≈ 1.6e-5 m²/s.
    assert!(
        (1.0e-5..2.3e-5).contains(&v),
        "N2-CO2 diffusivity out of range: {v:e}"
    );
    let _ = tpt_thermo_core::quantities::DiffusionCoefficient::new::<square_meter_per_second>(v);
}

#[test]
fn lucas_liquid_viscosity_positive_and_plausible() {
    let db = SeedComponentDatabase::from_seed();
    let water = db.index_of("water").unwrap();
    let mut z = vec![0.0_f64; db.num_components()];
    z[water] = 1.0;
    let t = Temperature::new::<kelvin>(300.0);
    let p = Pressure::new::<pascal>(1.0e5);
    let eta = tpt_thermo_transport::viscosity::lucas_liquid_viscosity(&db, t, p, &z).unwrap();
    // Smoke band: liquid viscosity must be positive and physically plausible
    // (documented as an approximate corresponding-states estimate).
    assert!((1e-4..10.0).contains(&eta.value), "eta = {}", eta.value);
}

#[test]
fn gas_thermal_conductivity_order_of_magnitude() {
    let db = SeedComponentDatabase::from_seed();
    let n2 = db.index_of("nitrogen").unwrap();
    let mut z = vec![0.0_f64; db.num_components()];
    z[n2] = 1.0;
    let t = Temperature::new::<kelvin>(300.0);
    let rho = tpt_thermo_transport::conductivity::ideal_molar_density(
        t,
        Pressure::new::<pascal>(1.0e5),
    );
    let lambda =
        tpt_thermo_transport::conductivity::chung_gas_thermal_conductivity(&db, t, rho, &z)
            .unwrap();
    // N2 gasconductivity at 300 K ≈ 0.026 W·m⁻¹·K⁻¹; Eucken closure is within ~30%.
    assert!(
        (0.008..0.05).contains(&lambda.value),
        "lambda = {}",
        lambda.value
    );
}
