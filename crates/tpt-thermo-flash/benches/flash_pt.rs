//! Criterion benchmark for the PT flash (target: < 1 ms for a 10-component feed).

use criterion::{criterion_group, criterion_main, Criterion};
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::PengRobinson;
use tpt_thermo_flash::FlashCalculator;
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;

fn bench_pt(c: &mut Criterion) {
    let db = SeedComponentDatabase::from_seed();
    let eos = PengRobinson::from_database(&db).unwrap();
    let calc = FlashCalculator::with_db(&eos, &db);
    let n = db.num_components().min(10);
    let mut z = vec![0.0_f64; db.num_components()];
    for (i, zi) in z.iter_mut().take(n).enumerate() {
        *zi = 1.0 / n as f64;
        let _ = i;
    }
    let t = Temperature::new::<kelvin>(300.0);
    let p = Pressure::new::<pascal>(2.0e6);
    c.bench_function("flash_pt_10_component", |b| {
        b.iter(|| calc.flash_pt(t, p, &z).unwrap())
    });
}

criterion_group!(benches, bench_pt);
criterion_main!(benches);
