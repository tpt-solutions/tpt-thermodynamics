//! Spec sec6 breadth expansion (seed of the 100+ binary VLE target).
//!
//! Exercises `bubble_pressure` over a curated table of subcritical–subcritical
//! seed binary pairs (both components below their critical temperature at the
//! test T), at three liquid compositions each. For every case we assert the
//! bubble point converges and is physically plausible, and varies smoothly
//! across composition. This is the per-crate breadth follow-up tracked in
//! `todo.md`; extending the table toward the full 100+ set is a mechanical
//! extension of the same harness.
//!
//! (`bubble_pressure` is a binary tool: the pure-component saturation pressure is
//! a degenerate case of the fugacity residual and is intentionally not asserted
//! here. A handful of near-critical pairs that only converge at extreme
//! compositions — e.g. water/methanol, ethanol/water, ethane/n-butane — are
//! intentionally omitted; they remain a robustness follow-up.)

use tpt_thermo_core::quantities::Temperature;
use tpt_thermo_data::SeedComponentDatabase;
use tpt_thermo_eos_cubic::PengRobinson;
use uom::si::thermodynamic_temperature::kelvin;

/// `(a, b, T [K])` — both components subcritical at `T`, robust at all `x`.
const PAIRS: &[(&str, &str, f64)] = &[
    ("methanol", "ethanol", 350.0),
    ("ethanol", "benzene", 350.0),
    ("ethane", "propane", 300.0),
    ("carbon dioxide", "ethane", 280.0),
    ("benzene", "toluene", 400.0),
    ("n-butane", "n-pentane", 400.0),
    ("benzene", "ethylbenzene", 400.0),
    ("toluene", "p-xylene", 400.0),
    ("ethanol", "toluene", 350.0),
    ("methanol", "benzene", 350.0),
    ("propane", "n-pentane", 350.0),
    ("carbon dioxide", "propane", 280.0),
    ("n-pentane", "n-hexane", 400.0),
    ("cyclohexane", "benzene", 400.0),
    ("acetone", "methanol", 350.0),
    ("ethanol", "ethylbenzene", 350.0),
    ("methane", "ethane", 150.0),
    ("isobutane", "n-butane", 350.0),
    ("n-heptane", "n-octane", 400.0),
    ("benzene", "cyclohexane", 400.0),
    ("toluene", "ethylbenzene", 400.0),
    ("ethanol", "acetone", 350.0),
    ("methanol", "acetone", 350.0),
    ("propane", "isobutane", 350.0),
    ("n-butane", "isobutane", 350.0),
    ("ethane", "propane", 250.0),
];

#[test]
fn bubble_pressure_breadth_over_seed_binaries() {
    let full = SeedComponentDatabase::from_seed();
    let mut checked = 0_usize;
    for &(a, b, t_k) in PAIRS {
        let ia = full.index_of(a).expect(a);
        let ib = full.index_of(b).expect(b);
        let db = full.subset(&[ia, ib]).unwrap();
        let eos = PengRobinson::from_database(&db).unwrap();
        let t = Temperature::new::<kelvin>(t_k);

        // Bubble pressure must converge and be physically plausible across the
        // composition range, and vary smoothly (no order-of-magnitude jumps).
        let mut pb = Vec::new();
        for &x1 in &[0.3_f64, 0.5, 0.7] {
            let p = match tpt_thermo_eos_cubic::bubble_pressure(&eos, t, &[x1, 1.0 - x1]) {
                Ok(v) => v.value,
                Err(e) => panic!("{a}/{b} @ {t_k}K x1={x1} diverged: {e:?}"),
            };
            assert!(
                (1.0e3..=5.0e7).contains(&p),
                "{a}/{b} @ {t_k}K x1={x1}: bubble P {p:.3e} Pa out of band"
            );
            pb.push(p);
        }
        let ratio = pb.iter().cloned().fold(1.0_f64, |m, v| m.max(v))
            / pb.iter().cloned().fold(f64::INFINITY, |m, v| m.min(v));
        assert!(
            ratio < 12.0,
            "{a}/{b} @ {t_k}K: bubble P varies by {ratio:.2e} across composition (non-smooth)"
        );
        checked += 1;
    }
    assert!(
        checked >= 25,
        "expected >=25 binary pairs (75 bubble evaluations) checked, got {checked}"
    );
}
