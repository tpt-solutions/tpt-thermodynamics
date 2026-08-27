//! Cross-crate integration test: `tpt-thermo-eos-cubic`'s Huron-Vidal mixing
//! rule consuming an activity model from this crate via the core
//! [`ExcessGibbsModel`] trait object.

use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_core::EquationOfState;
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_activity::parameters::TdParam;
use tpt_thermo_eos_activity::Nrtl;
use tpt_thermo_eos_cubic::cubic_solver::Phase;
use tpt_thermo_eos_cubic::{HuronVidal, HvVariant, PengRobinson};
use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};

#[test]
fn huron_vidal_consumes_nrtl_excess_model() {
    let db = SeedComponentDatabase::from_seed();
    // A binary NRTL parameterisation (indices are arbitrary; the database only
    // supplies critical data for the cubic EoS here).
    let nrtl = Nrtl::binary(
        TdParam::new(0.5, 100.0, 0.0),
        TdParam::new(-0.2, 50.0, 0.0),
        0.3,
    )
    .unwrap();
    // Box it as the core trait object the cubic mixing rule is generic over.
    let hv = HuronVidal::new(HvVariant::Mhv1, Box::new(nrtl));
    let eos = PengRobinson::with_mixing(&db, Box::new(hv)).unwrap();

    let t = Temperature::new::<kelvin>(350.0);
    let p = Pressure::new::<pascal>(1.0e6);
    let v: MolarVolume = eos.solve_phase(t, p, &[0.5, 0.5], Phase::Vapor).unwrap();
    assert!(v.value > 0.0);

    // Fugacity must be finite and consistent with the coupled excess model. The
    // coupling path (CubicMixing -> ExcessGibbsModel::reduced_excess_gibbs /
    // ln_gamma) is exercised inside `ln_fugacity_coefficient`.
    let ln_phi = eos.ln_fugacity_coefficient(t, v, &[0.5, 0.5], 0).unwrap();
    assert!(ln_phi.is_finite());
}
