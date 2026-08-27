use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::cubic_solver::{compressibility_roots, CubicModel};
use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};

fn db() -> SeedComponentDatabase { SeedComponentDatabase::from_seed() }
fn unit(i: usize) -> Vec<f64> {
    let n = db().num_components();
    let mut z = vec![0.0; n];
    z[i] = 1.0; z
}

fn disc(a: f64, b: f64, c: f64, d: f64) -> f64 {
    let b = b / a; let c = c / a; let d = d / a;
    let p = c - b * b / 3.0;
    let q = 2.0 * b * b * b / 27.0 - b * c / 3.0 + d;
    q * q / 4.0 + p * p * p / 27.0
}

#[test]
fn inspect_disc() {
    let eos = tpt_thermo_eos_cubic::PengRobinson::from_database(&db()).unwrap();
    let m = db().index_of("methane").unwrap;
    let t = Temperature::new::<kelvin>(150.0);
    let (amix, bmix, _a, _b) = eos.engine().mix_params(t.value, &unit(m));
    for p in [1.0e3, 1.0e6] {
        let a = amix * p / (8.314 * 8.314 * t.value * t.value);
        let b = bmix * p / (8.314 * t.value);
        let (c2, c1, c0) = (b - 1.0, a - 2.0 * b - 3.0 * b * b, b * b + b * b * b - a * b);
        let d = disc(1.0, c2, c1, c0);
        let r = compressibility_roots(CubicModel::PengRobinson, a, b);
        eprintln!("p={:.3e} A={:.3e} B={:.3e} disc={:.3e} roots={:?}", p, a, b, d, r);
    }
}
