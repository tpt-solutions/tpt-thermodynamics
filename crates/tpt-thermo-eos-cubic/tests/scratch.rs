use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_core::EquationOfState;
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::PengRobinson;
use uom::si::{molar_volume::cubic_meter_per_mole, pressure::pascal, thermodynamic_temperature::kelvin};

fn db() -> SeedComponentDatabase { SeedComponentDatabase::from_seed() }
fn unit(i: usize) -> Vec<f64> {
    let n = db().num_components();
    let mut z = vec![0.0; n];
    z[i] = 1.0; z
}

#[test]
fn inspect_methane() {
    let eos = PengRobinson::from_database(&db()).unwrap();
    let m = db().index_of("methane").unwrap();
    let t = Temperature::new::<kelvin>(150.0);
    let vc = eos.critical_point_pure(m).unwrap().2.value;
    eprintln!("vc={}", vc);
    for p in [1.0e3, 5.0e5, 1.0e6, 1.2e6, 1.4e6, 1.6e6, 2.0e6] {
        let r = eos.engine().z_roots(t, Pressure::new::<pascal>(p), &unit(m));
        let info = if r.len() == 3 {
            let zl = *r.first().unwrap();
            let zv = *r.last().unwrap();
            let vl = zl * 8.314 * t.value / p;
            let vv = zv * 8.314 * t.value / p;
            let lnl = eos.ln_fugacity_coefficient(t, MolarVolume::new::<cubic_meter_per_mole>(vl), &unit(m), 0).unwrap();
            let lnv = eos.ln_fugacity_coefficient(t, MolarVolume::new::<cubic_meter_per_mole>(vv), &unit(m), 0).unwrap();
            format!("vl={:.3e} vv={:.3e} g={:.4}", vl, vv, lnl - lnv)
        } else {
            format!("1root z={:.4}", r[0])
        };
        eprintln!("p={:.3e} n={} {}", p, r.len(), info);
    }
}
