//! Mixture critical-point calculation (Heidemann–Rahal condition: a horizontal
//! inflection on the `P–v` isotherm at fixed composition) and critical-locus
//! tracing for binaries.

use alloc::vec;
use alloc::vec::Vec;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_core::R;
use tpt_thermo_core::{ComponentDatabase, EquationOfState, ThermoError};
use uom::si::molar_volume::cubic_meter_per_mole;
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;

/// Initial guess for [`mixture_critical_point`], typically from pure-component
/// critical constants (mole-fraction averaged).
pub struct CriticalGuess {
    /// Guess temperature.
    pub t: Temperature,
    /// Guess pressure.
    pub p: Pressure,
}

impl CriticalGuess {
    /// Mole-fraction-averaged critical `T`/`P` from the database as a guess.
    pub fn from_database(db: &dyn ComponentDatabase, z: &[f64]) -> Self {
        let mut tm = 0.0_f64;
        let mut pm = 0.0_f64;
        for (i, zi) in z.iter().enumerate() {
            tm += zi * db.critical_temperature(i).map(|x| x.value).unwrap_or(0.0);
            pm += zi * db.critical_pressure(i).map(|x| x.value).unwrap_or(0.0);
        }
        CriticalGuess {
            t: Temperature::new::<kelvin>(tm.max(1.0)),
            p: Pressure::new::<pascal>(pm.max(1.0)),
        }
    }
}

/// Solve `J d = -F` for the 2×2 `J` with Levenberg–Marquardt damping, so the
/// (near-)singular Jacobian at a critical point still yields a usable step.
fn damped_solve(j: &[[f64; 2]; 2], f: &[f64; 2]) -> Option<(f64, f64)> {
    for lambda in [0.0_f64, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1] {
        let a = [[j[0][0] + lambda, j[0][1]], [j[1][0], j[1][1] + lambda]];
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        if det.abs() > 1e-30 {
            let dv = (-f[0] * a[1][1] + f[1] * a[0][1]) / det;
            let dt = (f[0] * a[1][0] - f[1] * a[0][0]) / det;
            return Some((dv, dt));
        }
    }
    None
}

/// Second derivative `(∂²P/∂v²)_T` via central differences on [`EquationOfState::dp_dv`].
fn d2p_dv2<E: EquationOfState + ?Sized>(eos: &E, t: Temperature, v: MolarVolume, z: &[f64]) -> f64 {
    let h = v.value.abs().max(1e-8) * 1e-4;
    let a = eos
        .dp_dv(t, MolarVolume::new::<cubic_meter_per_mole>(v.value - h), z)
        .unwrap_or(f64::NAN);
    let b = eos
        .dp_dv(t, MolarVolume::new::<cubic_meter_per_mole>(v.value + h), z)
        .unwrap_or(f64::NAN);
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        (b - a) / (2.0 * h)
    }
}

/// Mixture critical point at fixed composition `z`.
///
/// Solves the Heidemann–Rahal conditions `(∂P/∂v)_T = 0` and `(∂²P/∂v²)_T = 0`
/// for `(v, T)` by 2-D Newton iteration; the critical pressure then follows from
/// `P = pressure(v_c, T_c, z)`.
pub fn mixture_critical_point<E: EquationOfState + ?Sized>(
    eos: &E,
    z: &[f64],
    guess: CriticalGuess,
) -> Result<(Temperature, Pressure, MolarVolume), ThermoError> {
    let mut t = guess.t;
    // Start from a dense volume (typical critical compressibility ~0.25–0.3).
    let v0 = 0.25 * R * t.value / guess.p.value.max(1.0);
    let mut v = MolarVolume::new::<cubic_meter_per_mole>(v0.max(1e-6));
    for _ in 0..400 {
        let pv = eos.dp_dv(t, v, z)?;
        let pvv = d2p_dv2(eos, t, v, z);
        if !pv.is_finite() || !pvv.is_finite() {
            return Err(ThermoError::Numerical(
                tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                    tpt_thermo_core::NumericalIssueReason::NonPhysical,
                ),
            ));
        }
        // Dimensionless residuals (the absolute ∂²P/∂v² is enormous in Pa·m⁻⁶,
        // so scale by the ideal-gas pressure R T / v).
        let pmag = R * t.value / v.value.max(1e-12);
        let r1 = pv * v.value / pmag;
        let r2 = pvv * v.value * v.value / pmag;
        if r1.abs() < 1e-6 && r2.abs() < 1e-6 {
            let p = eos.pressure(t, v, z)?;
            return Ok((t, p, v));
        }
        let hv = v.value.abs().max(1e-8) * 1e-4;
        let ht = t.value.abs().max(1.0) * 1e-4;
        let pv_vp = eos
            .dp_dv(t, MolarVolume::new::<cubic_meter_per_mole>(v.value + hv), z)
            .unwrap_or(f64::NAN);
        let pv_vm = eos
            .dp_dv(t, MolarVolume::new::<cubic_meter_per_mole>(v.value - hv), z)
            .unwrap_or(f64::NAN);
        let pv_tp = eos
            .dp_dv(Temperature::new::<kelvin>(t.value + ht), v, z)
            .unwrap_or(f64::NAN);
        let pv_tm = eos
            .dp_dv(Temperature::new::<kelvin>(t.value - ht), v, z)
            .unwrap_or(f64::NAN);
        let pvv_vp = d2p_dv2(
            eos,
            t,
            MolarVolume::new::<cubic_meter_per_mole>(v.value + hv),
            z,
        );
        let pvv_vm = d2p_dv2(
            eos,
            t,
            MolarVolume::new::<cubic_meter_per_mole>(v.value - hv),
            z,
        );
        let pvv_tp = d2p_dv2(eos, Temperature::new::<kelvin>(t.value + ht), v, z);
        let pvv_tm = d2p_dv2(eos, Temperature::new::<kelvin>(t.value - ht), v, z);
        if [pv_vp, pv_vm, pv_tp, pv_tm, pvv_vp, pvv_vm, pvv_tp, pvv_tm]
            .iter()
            .any(|x| x.is_nan())
        {
            return Err(ThermoError::Numerical(
                tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                    tpt_thermo_core::NumericalIssueReason::NonPhysical,
                ),
            ));
        }
        let j = [
            [(pv_vp - pv_vm) / (2.0 * hv), (pv_tp - pv_tm) / (2.0 * ht)],
            [
                (pvv_vp - pvv_vm) / (2.0 * hv),
                (pvv_tp - pvv_tm) / (2.0 * ht),
            ],
        ];
        // The critical point is a horizontal inflection, so the (v,T) Jacobian
        // is singular there; regularise with Levenberg–Marquardt damping.
        let f = [pv, pvv];
        let (dv, dt) = damped_solve(&j, &f).ok_or({
            ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::NumericalIssueReason::SingularJacobian,
            ))
        })?;
        let mut step = 1.0_f64;
        loop {
            let nv = v.value + step * dv;
            let nt = t.value + step * dt;
            if nv > 1e-7 && nt > 1.0 {
                v = MolarVolume::new::<cubic_meter_per_mole>(nv);
                t = Temperature::new::<kelvin>(nt);
                break;
            }
            step *= 0.5;
            if step < 1e-6 {
                return Err(ThermoError::Numerical(
                    tpt_thermo_core::ConvergenceStatus::NotConverged,
                ));
            }
        }
    }
    Err(ThermoError::Numerical(
        tpt_thermo_core::ConvergenceStatus::NotConverged,
    ))
}

/// Trace the critical locus of a binary (`components i0`/`i1`) as the mole
/// fraction `z_{i0}` sweeps `0 → 1` in `n` steps, returning
/// `(z_{i0}, T_c, P_c)` points that converged. The composition is built at the
/// full database length (other components set to zero).
pub fn critical_locus_binary<E: EquationOfState + ?Sized>(
    eos: &E,
    db: &dyn ComponentDatabase,
    i0: usize,
    i1: usize,
    n: usize,
) -> Vec<(f64, Temperature, Pressure)> {
    let mut out = Vec::new();
    if n == 0 {
        return out;
    }
    for k in 0..=n {
        let z1 = k as f64 / n as f64;
        let mut z = vec![0.0_f64; db.num_components()];
        z[i0] = z1;
        z[i1] = 1.0 - z1;
        let guess = CriticalGuess::from_database(db, &z);
        if let Ok((t, p, _v)) = mixture_critical_point(eos, &z, guess) {
            out.push((z1, t, p));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_thermo_core::component::ComponentDatabase;
    use tpt_thermo_data::SeedComponentDatabase;
    use tpt_thermo_eos_cubic::PengRobinson;

    #[test]
    fn binary_critical_between_pure_endpoints() {
        let db = SeedComponentDatabase::from_seed();
        let eos = PengRobinson::from_database(&db).unwrap();
        let methane = db.index_of("methane").unwrap();
        let ethane = db.index_of("ethane").unwrap();
        let mut z = vec![0.0; db.num_components()];
        z[methane] = 0.5;
        z[ethane] = 0.5;
        let guess = CriticalGuess::from_database(&db, &z);
        let (tc, pc, _v) = mixture_critical_point(&eos, &z, guess).unwrap();
        let tcm = db.critical_temperature(methane).unwrap().value;
        let tce = db.critical_temperature(ethane).unwrap().value;
        let pcm = db.critical_pressure(methane).unwrap().value;
        let pce = db.critical_pressure(ethane).unwrap().value;
        assert!(tc.value > tcm.min(tce) - 1.0 && tc.value < tce.max(tcm) + 1.0);
        assert!(pc.value > pcm.min(pce) * 0.5 && pc.value < pce.max(pcm) * 1.5);
    }
}
