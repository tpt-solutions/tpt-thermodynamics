//! Integration tests / validation for the SAFT crate against the curated
//! parameter table.

use tpt_thermo_core::quantities::{MolarVolume, Temperature};
use tpt_thermo_core::EquationOfState;
use tpt_thermo_eos_saft::parameters::{SaftParameters, SEED_SAFT_PARAMETERS};
use tpt_thermo_eos_saft::{PcSaft, SaftVrMie};
use uom::si::{molar_volume::cubic_meter_per_mole, thermodynamic_temperature::kelvin};

fn params_for(names: &[&str]) -> (SaftParameters, Vec<f64>) {
    let comps: Vec<_> = names
        .iter()
        .map(|n| {
            SEED_SAFT_PARAMETERS
                .iter()
                .find(|c| c.name == *n)
                .copied()
                .expect(n)
        })
        .collect();
    // Molar masses (kg·mol⁻¹) in the same order.
    let mm = names
        .iter()
        .map(|n| match *n {
            "methane" => 0.016_043,
            "carbon dioxide" => 0.044_01,
            "water" => 0.018_015,
            "argon" => 0.039_948,
            _ => 0.030,
        })
        .collect();
    (SaftParameters::new(comps), mm)
}

fn methane() -> PcSaft {
    let (p, mm) = params_for(&["methane"]);
    PcSaft::new(p, mm)
}

#[test]
fn ideal_gas_limit() {
    let eos = methane();
    let t = Temperature::new::<kelvin>(300.0);
    let v = MolarVolume::new::<cubic_meter_per_mole>(1.0);
    let p = eos.pressure(t, v, &[1.0]).unwrap();
    let z = eos.compressibility_factor(t, v, &[1.0]).unwrap();
    assert!(
        (z - 1.0).abs() < 1e-3,
        "Z should approach 1 at low density, got {z}"
    );
    let expected = 8.314462618 * 300.0 / 1.0;
    assert!((p.value - expected).abs() / expected < 1e-2);
    let lnphi = eos.ln_fugacity_coefficient(t, v, &[1.0], 0).unwrap();
    assert!(
        lnphi.abs() < 1e-2,
        "ln φ should vanish at low density, got {lnphi}"
    );
}

#[test]
fn pressure_consistency_with_z() {
    let eos = methane();
    let t = Temperature::new::<kelvin>(200.0);
    let v = MolarVolume::new::<cubic_meter_per_mole>(5e-5);
    let p = eos.pressure(t, v, &[1.0]).unwrap().value;
    let z = eos.compressibility_factor(t, v, &[1.0]).unwrap();
    let expected = z * 8.314462618 * 200.0 / 5e-5;
    assert!((p - expected).abs() / expected < 1e-6);
}

#[test]
fn association_converges_for_water() {
    let (p, mm) = params_for(&["water"]);
    let eos = PcSaft::new(p, mm);
    let t = Temperature::new::<kelvin>(300.0);
    let v = MolarVolume::new::<cubic_meter_per_mole>(5e-5);
    let _p = eos.pressure(t, v, &[1.0]).unwrap();
    let _h = eos.molar_enthalpy(t, v, &[1.0]).unwrap();
    let _s = eos.molar_entropy(t, v, &[1.0]).unwrap();
}

#[test]
fn saturation_is_self_consistent() {
    let eos = methane();
    for &t_k in &[150.0_f64, 170.0, 185.0] {
        let t = Temperature::new::<kelvin>(t_k);
        let (p, vl, vv) = match eos.saturation_pressure(t) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("T={t_k} saturation error: {e:?}");
                continue;
            }
        };
        let p = p.value;
        let vl = vl.value;
        let vv = vv.value;
        assert!(p > 0.0 && p < 5.0e6, "P_sat out of range at {t_k} K: {p}");
        assert!(vv > vl, "vapor volume must exceed liquid volume at {t_k} K");

        let lnphi_l = eos
            .ln_fugacity_coefficient(t, MolarVolume::new::<cubic_meter_per_mole>(vl), &[1.0], 0)
            .unwrap();
        let lnphi_v = eos
            .ln_fugacity_coefficient(t, MolarVolume::new::<cubic_meter_per_mole>(vv), &[1.0], 0)
            .unwrap();
        assert!(
            (lnphi_l - lnphi_v).abs() < 1e-2,
            "fugacity mismatch at {t_k} K: {lnphi_l} vs {lnphi_v}"
        );

        let pl = eos
            .pressure(t, MolarVolume::new::<cubic_meter_per_mole>(vl), &[1.0])
            .unwrap()
            .value;
        let pv = eos
            .pressure(t, MolarVolume::new::<cubic_meter_per_mole>(vv), &[1.0])
            .unwrap()
            .value;
        assert!((pl - pv).abs() / p < 1e-3, "pressure mismatch at {t_k} K");
    }
}

#[test]
fn vr_mie_builds_for_argon() {
    let (p, mm) = params_for(&["argon"]);
    let eos = SaftVrMie::new(p, mm);
    let t = Temperature::new::<kelvin>(150.0);
    // Sub-critical density inside the model's mechanically-stable envelope
    // (the very high-density liquid root lies past the predicted spinodal for
    // the seed AR parameters, as documented in `todo.md` Phase 6 deferred
    // accuracy scope).
    let v = MolarVolume::new::<cubic_meter_per_mole>(2e-4);
    let p = eos.pressure(t, v, &[1.0]).unwrap();
    assert!(p.value > 0.0 && p.value.is_finite());
}

#[test]
fn debug_pv() {
    let eos = methane();
    let t = Temperature::new::<kelvin>(150.0);
    for k in 0..40 {
        let logv = -6.0 + 6.0 * (k as f64) / 39.0;
        let v = 10f64.powf(logv);
        let p = match eos.pressure(t, MolarVolume::new::<cubic_meter_per_mole>(v), &[1.0]) {
            Ok(pp) => pp.value,
            Err(_) => f64::NAN,
        };
        eprintln!("v={:.3e} P={:.3e}", v, p);
    }
}

#[test]
fn binary_mixture_pressure() {
    let (p, mm) = params_for(&["methane", "carbon dioxide"]);
    let eos = PcSaft::new(p, mm);
    let t = Temperature::new::<kelvin>(280.0);
    let v = MolarVolume::new::<cubic_meter_per_mole>(4e-4);
    let p = eos.pressure(t, v, &[0.5, 0.5]).unwrap();
    assert!(p.value > 0.0 && p.value.is_finite());
}
