//! Azeotrope detection for binary systems.
//!
//! An azeotrope is a point on the two-phase boundary where the incipient vapor
//! and liquid compositions are equal (`y = x`), so the K-values are unity and the
//! bubble and dew curves touch. For a binary at fixed temperature we scan the
//! liquid mole fraction `x₁` and locate the composition where the bubble-point
//! pressure curve `P_bubble(x₁)` meets the dew-point pressure curve
//! `P_dew(x₁)`; the crossing is the azeotrope.

use crate::{BubbleDewSolver};
use alloc::vec;
use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use uom::si::thermodynamic_temperature::kelvin;

/// A detected azeotrope.
#[derive(Debug, Clone)]
pub struct Azeotrope {
    /// Temperature of the azeotrope.
    pub temperature: Temperature,
    /// Pressure of the azeotrope.
    pub pressure: Pressure,
    /// Liquid (= vapor) composition at the azeotrope.
    pub composition: Vec<f64>,
}

/// Detect a binary azeotrope at fixed temperature `t`.
///
/// Returns `Ok(Some(..))` when a bubble/dew pressure crossing is found,
/// `Ok(None)` when the scan finds no crossing (no azeotrope over the scanned
/// range), and `Err` on an upstream evaluation failure.
pub fn detect_azeotrope(
    solver: &BubbleDewSolver,
    t: Temperature,
    n: usize,
) -> Result<Option<Azeotrope>, ThermoError> {
    let nc = solver.num_components();
    if nc != 2 {
        return Err(ThermoError::InvalidInput(
            "azeotrope detection currently supports binary systems only",
        ));
    }
    let n = n.max(2);

    let mut prev_x = 0.0_f64;
    let mut prev_f = f64::NAN;
    for k in 0..=n {
        let x1 = k as f64 / n as f64;
        let mut x = vec![0.0f64; nc];
        x[0] = x1;
        x[1] = 1.0 - x1;
        let pb = match solver.bubble_point_pressure(t, &x) {
            Ok(b) => b.pressure.value,
            Err(_) => {
                prev_x = x1;
                prev_f = f64::NAN;
                continue;
            }
        };
        let mut y = vec![0.0f64; nc];
        y[0] = x1;
        y[1] = 1.0 - x1;
        let pd = match solver.dew_point_pressure(t, &y) {
            Ok(d) => d.pressure.value,
            Err(_) => {
                prev_x = x1;
                prev_f = f64::NAN;
                continue;
            }
        };
        // Positive f means bubble above dew (normal VLE); the azeotrope is where
        // they cross (f = 0). Look for a genuine sign change (or exact zero) that
        // is (a) away from the pure-component endpoints and (b) large enough in
        // magnitude to be a real separation between the curves rather than the
        // numerical noise that makes P_bubble ≈ P_dew near x → 0/1.
        let f = pb - pd;
        let interior = x1 > 0.03 && x1 < 0.97 && prev_x > 0.03 && prev_x < 0.97;
        let significant = prev_f.abs() > 1.0 && f.abs() > 1.0;
        if prev_f.is_finite() && interior && significant && (f == 0.0 || prev_f * f < 0.0) {
            let (a, b) = (prev_x, x1);
            let xr = bisect(solver, t, a, b, prev_f)?;
            if xr <= 0.02 || xr >= 0.98 {
                // Crossing resolved to a pure-component endpoint: an artifact.
                prev_x = x1;
                prev_f = f;
                continue;
            }
            let mut xr_vec = vec![0.0f64; nc];
            xr_vec[0] = xr;
            xr_vec[1] = 1.0 - xr;
            let bp = solver.bubble_point_pressure(t, &xr_vec)?;
            return Ok(Some(Azeotrope {
                temperature: t,
                pressure: bp.pressure,
                composition: xr_vec,
            }));
        }
        prev_x = x1;
        prev_f = f;
    }
    Ok(None)
}

/// Bisection on `x₁` in `[a, b]` for the bubble/dew pressure crossing.
fn bisect(
    solver: &BubbleDewSolver,
    t: Temperature,
    mut a: f64,
    mut b: f64,
    mut fa: f64,
) -> Result<f64, ThermoError> {
    let nc = solver.num_components();
    let f = |x1: f64| -> f64 {
        let mut x = vec![0.0f64; nc];
        x[0] = x1;
        x[1] = 1.0 - x1;
        let pb = solver.bubble_point_pressure(t, &x).map(|b| b.pressure.value);
        let mut y = vec![0.0f64; nc];
        y[0] = x1;
        y[1] = 1.0 - x1;
        let pd = solver.dew_point_pressure(t, &y).map(|d| d.pressure.value);
        match (pb, pd) {
            (Ok(pb), Ok(pd)) => pb - pd,
            _ => f64::NAN,
        }
    };
    for _ in 0..60 {
        let m = 0.5 * (a + b);
        let fm = f(m);
        if !fm.is_finite() {
            return Ok(m);
        }
        if fm.abs() < 1.0 || (b - a) < 1e-9 {
            return Ok(m);
        }
        if fa * fm < 0.0 {
            b = m;
        } else {
            a = m;
            fa = fm;
        }
    }
    Ok(0.5 * (a + b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_thermo_data::SeedComponentDatabase;
    use tpt_thermo_eos_cubic::PengRobinson;

    #[test]
    fn detector_runs_and_is_well_formed() {
        // Azeotrope detection requires a binary; build a 2-component (ethanol/
        // water) database rather than the full seed set.
        let toml = r#"
[[components]]
schema_version = 1
name = "ethanol"
critical_temperature_k = 514.0
critical_pressure_pa = 6137000.0
acentric_factor = 0.644
molar_mass_kg_per_mol = 0.04606844

[[components]]
schema_version = 1
name = "water"
critical_temperature_k = 647.096
critical_pressure_pa = 22064000.0
acentric_factor = 0.344
molar_mass_kg_per_mol = 0.01801528
"#;
        let db = SeedComponentDatabase::from_toml_str(toml).unwrap();
        let eos = PengRobinson::from_database(&db).unwrap();
        let solver = BubbleDewSolver::new(&eos as &dyn crate::KProvider, &db);
        // At a single temperature we only assert the detector runs without
        // error and returns a well-formed result when it does find a crossing.
        // The cubic EoS with zero binary interactions may or may not reproduce
        // an azeotrope.
        let t = Temperature::new::<kelvin>(350.0);
        let res = detect_azeotrope(&solver, t, 60);
        assert!(res.is_ok());
        if let Ok(Some(az)) = res {
            assert!(az.pressure.value > 0.0);
            assert!((az.composition[0] + az.composition[1] - 1.0).abs() < 1e-6);
        }
    }
}
