use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_core::EquationOfState;
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::PengRobinson;
use tpt_thermo_eos_cubic::cubic_solver::Phase;
use uom::si::{molar_volume::cubic_meter_per_mole, pressure::pascal, thermodynamic_temperature::kelvin};

fn db() -> SeedComponentDatabase { SeedComponentDatabase::from_seed() }
fn unit(i: usize) -> Vec<f64> {
    let n = db().num_components();
    let mut z = vec![0.0; n];
    z[i] = 1.0; z
}

#[test]
fn inspect_vp() {
    let eos = PengRobinson::from_database(&db()).unwrap();
    let co2 = db().index_of("carbon dioxide").unwrap();
    let t = Temperature::new::<kelvin>(250.0);
    let pc = db().critical_pressure(co2).unwrap().value;
    for k in 1..=20 {
        let p = pc * (k as f64) / 20.0;
        let roots = eos.engine().z_roots(t, Pressure::new::<pascal>(p), &unit(co2));
        let g = if roots.len() == 3 {
            let zl = *roots.first().unwrap();
            let zv = *roots.last().unwrap();
            let vl = zl * 8.314 * t.value / p;
            let vv = zv * 8.314 * t.value / p;
            let lnl = eos.ln_fugacity_coefficient(t, MolarVolume::new::<cubic_meter_per_mole>(vl), &unit(co2), 0).unwrap();
            let lnv = eos.ln_fugacity_coefficient(t, MolarVolume::new::<cubic_meter_per_mole>(vv), &unit(co2), 0).unwrap();
            Some(lnl - lnv)
        } else { None };
        eprintln!("p={:.3e} nroots={} g={:?}", p, roots.len(), g);
    }
}
