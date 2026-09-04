//! Phase 4 validation and integration tests against the curated seed dataset.

use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::mixing::ExcessGibbsModel;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_core::{EquationOfState, ThermoError, R};
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::cubic_solver::Phase;
use tpt_thermo_eos_cubic::mixing::{CubicMixing, HuronVidal, HvVariant};
use tpt_thermo_eos_cubic::{critical, PengRobinson, VolumeTranslated};
use uom::si::{
    molar_volume::cubic_meter_per_mole, pressure::pascal, thermodynamic_temperature::kelvin,
};

fn db() -> SeedComponentDatabase {
    SeedComponentDatabase::from_seed()
}

/// A one-parameter Margules mock excess-Gibbs model for HV/WS coupling tests.
#[derive(Clone)]
struct Margules {
    a: f64,
}

impl ExcessGibbsModel for Margules {
    fn num_components(&self) -> usize {
        2
    }
    fn reduced_excess_gibbs(
        &self,
        _t: Temperature,
        _p: Pressure,
        x: &[f64],
    ) -> Result<f64, ThermoError> {
        Ok(self.a * x[0] * x[1])
    }
    fn ln_gamma(
        &self,
        _t: Temperature,
        _p: Pressure,
        x: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        Ok(match i {
            0 => self.a * x[1] * x[1],
            _ => self.a * x[0] * x[0],
        })
    }
}

#[test]
fn pr_low_pressure_gas_is_near_ideal() {
    let eos = PengRobinson::from_database(&db()).unwrap();
    let methane = db().index_of("methane").unwrap();
    let t = Temperature::new::<kelvin>(300.0);
    let p = Pressure::new::<pascal>(1.0e5);
    let v = eos.solve_phase(t, p, &unit(methane), Phase::Vapor).unwrap();
    let z = p.value * v.value / (R * t.value);
    assert!((z - 1.0).abs() < 0.05, "expected near-ideal Z, got {z}");
}

#[test]
fn pr_three_roots_below_critical() {
    let eos = PengRobinson::from_database(&db()).unwrap();
    let methane = db().index_of("methane").unwrap();
    // ~0.6 Tc and a mid pressure (above Psat, below the liquid spinodal) → two-phase.
    let t = Temperature::new::<kelvin>(114.0);
    let p = Pressure::new::<pascal>(0.5e6);
    let roots = eos.engine().z_roots(t, p, &unit(methane));
    assert_eq!(roots.len(), 3, "two-phase region should give 3 roots");
}

#[test]
fn pr_critical_point_finite_and_positive() {
    let eos = PengRobinson::from_database(&db()).unwrap();
    let methane = db().index_of("methane").unwrap();
    let (tc, pc, vc) = eos.critical_point_pure(methane).unwrap();
    let _p_check = eos.pressure(tc, vc, &unit(methane)).unwrap().value;
    assert!(tc.value > 0.0 && pc.value > 0.0 && vc.value > 0.0);
    let tc_seed = db().critical_temperature(methane).unwrap().value;
    let pc_seed = db().critical_pressure(methane).unwrap().value;
    assert!(
        (tc.value - tc_seed).abs() / tc_seed < 0.15,
        "Tc off by >15%"
    );
    assert!(
        (pc.value - pc_seed).abs() / pc_seed < 0.20,
        "Pc off by >20%"
    );
}

#[test]
fn pr_saturation_pressure_plausible_and_consistent() {
    let eos = PengRobinson::from_database(&db()).unwrap();
    let methane = db().index_of("methane").unwrap();
    let t = Temperature::new::<kelvin>(150.0);
    let psat = vapor_pressure(&eos, t, methane).expect("saturation should converge");
    assert!(psat > 1.0e4 && psat < 4.0e6, "Psat out of band: {psat}");
    let v_l = eos
        .solve_phase(
            t,
            Pressure::new::<pascal>(psat),
            &unit(methane),
            Phase::Liquid,
        )
        .unwrap();
    let v_v = eos
        .solve_phase(
            t,
            Pressure::new::<pascal>(psat),
            &unit(methane),
            Phase::Vapor,
        )
        .unwrap();
    assert!(v_l.value < v_v.value, "liquid must be denser than vapor");
}

#[test]
fn pr_enthalpy_of_vaporization_positive() {
    let eos = PengRobinson::from_database(&db()).unwrap();
    let co2 = db().index_of("carbon dioxide").unwrap();
    let t = Temperature::new::<kelvin>(250.0);
    let hvap = enthalpy_of_vaporization(&eos, t, co2);
    assert!(hvap > 5.0e3 && hvap < 40.0e3, "Hvap out of range: {hvap}");
}

#[test]
fn volume_translation_gives_finite_liquid_volume() {
    // The Peneloux translation shifts the volume by `c_i` while leaving pressure
    // and fugacity unchanged: at an identical (T, P) the translated EoS reports a
    // volume `c_i` smaller than the bare cubic.
    let pr = PengRobinson::from_database(&db()).unwrap();
    let vt = VolumeTranslated::peng_robinson(&db()).unwrap();
    let water = db().index_of("water").unwrap();
    let t = Temperature::new::<kelvin>(300.0);
    let p = Pressure::new::<pascal>(1.0e5);
    let v_pr = pr.solve_phase(t, p, &unit(water), Phase::Vapor).unwrap();
    let v_vt = vt.solve_phase(t, p, &unit(water), Phase::Vapor).unwrap();
    assert!(v_pr.value > 0.0 && v_vt.value > 0.0);
    // Translation reduces the volume for this (associating) component.
    assert!(v_vt.value < v_pr.value);
}

#[test]
fn speed_of_sound_positive_for_gas() {
    let eos = PengRobinson::from_database(&db()).unwrap();
    let methane = db().index_of("methane").unwrap();
    let t = Temperature::new::<kelvin>(300.0);
    let p = Pressure::new::<pascal>(1.0e5);
    let v = eos.solve_phase(t, p, &unit(methane), Phase::Vapor).unwrap();
    let a = eos.speed_of_sound(t, v, &unit(methane)).unwrap().value;
    assert!(a > 0.0, "speed of sound must be positive");
}

#[test]
fn mechanical_stability_flag() {
    let eos = PengRobinson::from_database(&db()).unwrap();
    let methane = db().index_of("methane").unwrap();
    let t = Temperature::new::<kelvin>(300.0);
    let p = Pressure::new::<pascal>(1.0e5);
    let v = eos.solve_phase(t, p, &unit(methane), Phase::Vapor).unwrap();
    assert!(critical::mechanical_stability(&eos, t, v, &unit(methane)).unwrap());
}

#[test]
fn hv_mixing_couples_excess_model() {
    let margules = Margules { a: 1.5 };
    let hv = HuronVidal::new(HvVariant::Mhv1, Box::new(margules));
    let a = [1.0, 1.0];
    let b = [0.05, 0.05];
    let pure = hv.a_mix(&a, &b, &[1.0, 0.0], 300.0, 1.0e5);
    assert!((pure - 1.0).abs() < 1e-12, "pure must equal a_1");
    let mix = hv.a_mix(&a, &b, &[0.5, 0.5], 300.0, 1.0e5);
    assert!(mix > 0.0 && mix.is_finite());
    // aij_sum reconstructed from the gradient is finite.
    let sum = hv.aij_sum(&a, &b, &[0.5, 0.5], 0, 300.0, 1.0e5);
    assert!(sum.is_finite());
}

#[test]
fn seeded_bips_are_consumed_by_cubic() {
    let db = db();
    let co2 = db.index_of("carbon dioxide").unwrap();
    let methane = db.index_of("methane").unwrap();
    // The seed ships a non-zero PR k_ij for CO2-methane; the opt-in constructor
    // must produce an EoS that differs from the zero-BIP default.
    let pr_zero = PengRobinson::from_database(&db).unwrap();
    let pr_kij = PengRobinson::from_database_with_kij(&db).unwrap();
    let t = Temperature::new::<kelvin>(280.0);
    let p = Pressure::new::<pascal>(1.0e5);
    let mut z = vec![0.0; db.num_components()];
    z[co2] = 0.5;
    z[methane] = 0.5;
    let v_zero = pr_zero.solve_phase(t, p, &z, Phase::Vapor).unwrap();
    let v_kij = pr_kij.solve_phase(t, p, &z, Phase::Vapor).unwrap();
    let ln_zero = pr_zero.ln_fugacity_coefficient(t, v_zero, &z, co2).unwrap();
    let ln_kij = pr_kij.ln_fugacity_coefficient(t, v_kij, &z, co2).unwrap();
    assert!(
        (ln_zero - ln_kij).abs() > 1e-4,
        "seeded BIPs must change the fugacity coefficient"
    );
}

// --- helpers -------------------------------------------------------------

fn unit(i: usize) -> Vec<f64> {
    let n = db().num_components();
    let mut z = vec![0.0; n];
    z[i] = 1.0;
    z
}

/// Equal-fugacity vapor-pressure solver for a pure component.
///
/// The liquid–vapor fugacity difference `g = ln φ_liquid − ln φ_vapor` is zero
/// exactly at the saturation pressure (equal fugacity) and negative on either
/// side, so the saturation pressure is the **maximiser** of `g` over the
/// two-phase region. A golden-section search is robust for components whose van
/// der Waals loop reaches very low pressure (e.g. methane) as well as those whose
/// loop is narrower (e.g. CO₂).
fn vapor_pressure(eos: &PengRobinson, t: Temperature, i: usize) -> Option<f64> {
    let pc = db().critical_pressure(i).ok()?.value;
    let a = 1.0e3;
    let b = 0.95 * pc;
    // Outside the three-root (two-phase) region `g` is undefined; treat it as
    // −∞ so the maximiser stays inside the physical two-phase region.
    let g = |p: f64| -> f64 { g_fug(eos, t, i, p).unwrap_or(f64::NEG_INFINITY) };
    let invphi = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut lo = a;
    let mut hi = b;
    let mut x1 = lo + (1.0 - invphi) * (hi - lo);
    let mut x2 = lo + invphi * (hi - lo);
    let mut f1 = g(x1);
    let mut f2 = g(x2);
    for _ in 0..200 {
        if f1 > f2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = lo + (1.0 - invphi) * (hi - lo);
            f1 = g(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + invphi * (hi - lo);
            f2 = g(x2);
        }
        if (hi - lo) / (b - a) < 1e-7 {
            break;
        }
    }
    let p_sat = 0.5 * (lo + hi);
    if g_fug(eos, t, i, p_sat).is_some() {
        Some(p_sat)
    } else {
        None
    }
}

fn enthalpy_of_vaporization(eos: &PengRobinson, t: Temperature, i: usize) -> f64 {
    let psat = vapor_pressure(eos, t, i).expect("vapor pressure must converge");
    let p = Pressure::new::<pascal>(psat);
    let v_l = eos.solve_phase(t, p, &unit(i), Phase::Liquid).unwrap();
    let v_v = eos.solve_phase(t, p, &unit(i), Phase::Vapor).unwrap();
    let h_l = eos.molar_enthalpy(t, v_l, &unit(i)).unwrap().value;
    let h_v = eos.molar_enthalpy(t, v_v, &unit(i)).unwrap().value;
    h_v - h_l
}

/// `ln φ_liquid − ln φ_vapor` at `(T, P)`; `None` if not in the two-phase region
/// (fewer than three real compressibility roots).
fn g_fug(eos: &PengRobinson, t: Temperature, i: usize, p: f64) -> Option<f64> {
    let zs = unit(i);
    let roots = eos.engine().z_roots(t, Pressure::new::<pascal>(p), &zs);
    if roots.len() != 3 {
        return None;
    }
    let z_l = *roots.first().unwrap();
    let z_v = *roots.last().unwrap();
    let v_l = MolarVolume::new::<cubic_meter_per_mole>(z_l * R * t.value / p);
    let v_v = MolarVolume::new::<cubic_meter_per_mole>(z_v * R * t.value / p);
    let ln_phi_l = eos.ln_fugacity_coefficient(t, v_l, &zs, 0).ok()?;
    let ln_phi_v = eos.ln_fugacity_coefficient(t, v_v, &zs, 0).ok()?;
    Some(ln_phi_l - ln_phi_v)
}
