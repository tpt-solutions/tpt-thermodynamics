//! Trial-composition initialisation strategies for the stability test: Wilson
//! K-values (analytic, from critical constants) and structured grids.

use alloc::vec;
use alloc::vec::Vec;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use tpt_thermo_core::ComponentDatabase;

/// Wilson (1969) approximate K-values `K_i ≈ (P_c/P)·exp[5.37(1+ω_i)(1−T_c/T)]`.
///
/// These seed the tangent-plane-distance successive-substitution iteration.
pub fn wilson_k_values(db: &dyn ComponentDatabase, t: Temperature, p: Pressure) -> Vec<f64> {
    let nc = db.num_components();
    let mut k = Vec::with_capacity(nc);
    for i in 0..nc {
        let tc = db.critical_temperature(i).map(|x| x.value).unwrap_or(300.0);
        let pc = db.critical_pressure(i).map(|x| x.value).unwrap_or(1.0e6);
        let w = db.acentric_factor(i).unwrap_or(0.0);
        let ki = (pc / p.value.max(1.0)) * (5.37 * (1.0 + w) * (1.0 - tc / t.value.max(1.0))).exp();
        k.push(ki);
    }
    k
}

/// Unit-vector (pure-component) trial compositions — one per component.
pub fn pure_component_trials(nc: usize) -> Vec<Vec<f64>> {
    (0..nc)
        .map(|i| {
            let mut v = vec![0.0_f64; nc];
            v[i] = 1.0;
            v
        })
        .collect()
}

/// A structured grid of binary splits between every pair of components, plus the
/// pure-component vertices. `divisions` controls the number of interior points
/// per pair; kept small (the TPD iteration does the real work from here).
pub fn regular_grid_trials(nc: usize, divisions: usize) -> Vec<Vec<f64>> {
    let mut out = pure_component_trials(nc);
    if nc >= 2 {
        for i in 0..nc {
            for j in (i + 1)..nc {
                for d in 1..divisions {
                    let f = d as f64 / divisions as f64;
                    let mut v = vec![0.0_f64; nc];
                    v[i] = f;
                    v[j] = 1.0 - f;
                    out.push(v);
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_thermo_data::SeedComponentDatabase;
    use uom::si::pressure::pascal;
    use uom::si::thermodynamic_temperature::kelvin;

    #[test]
    fn wilson_more_volatile_at_higher_t() {
        let db = SeedComponentDatabase::from_seed();
        let t = Temperature::new::<kelvin>(300.0);
        let p = Pressure::new::<pascal>(1.0e5);
        let k = wilson_k_values(&db, t, p);
        let methane = db.index_of("methane").unwrap();
        let water = db.index_of("water").unwrap();
        // Methane (light) should have a larger K than water (heavy) at 300 K, 1 bar.
        assert!(
            k[methane] > k[water],
            "light component should be more volatile"
        );
    }

    #[test]
    fn grid_trials_normalised() {
        let trials = regular_grid_trials(3, 4);
        for tr in &trials {
            let s: f64 = tr.iter().sum();
            assert!((s - 1.0).abs() < 1e-12);
        }
        assert!(trials.len() >= 3);
    }
}
