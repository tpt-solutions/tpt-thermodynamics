//! The four remaining flash specifications (PH, TV, TS, PU, PV).
//!
//! Each is an outer root-find on the "free" variable (`T` for PH/PU/PV, `P` for
//! TS/TV) around an inner [`flash_pt`](crate::pt::flash_pt), matching the target
//! molar enthalpy / volume / entropy / internal energy. The mixture property is
//! reconstructed from the inner PT split's phase volumes and compositions.

use tpt_thermo_core::convergence::ConvergenceStatus;
use tpt_thermo_core::eos::EquationOfState;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{MolarEnergy, MolarEntropy, MolarVolume, Pressure, Temperature};
use uom::si::molar_energy::joule_per_mole;
use uom::si::molar_volume::cubic_meter_per_mole;
use uom::si::thermodynamic_temperature::kelvin;

use crate::pt::{flash_pt_impl, FlashCalculator, FlashResult};
use crate::FlashError;

/// Default iteration budget for the outer loops.
const OUTER_MAX_ITER: usize = 80;
/// Default target tolerance for the outer residual (J·mol⁻¹ or m³·mol⁻¹).
const OUTER_TOL: f64 = 1e-5;

/// Mixture molar enthalpy (J·mol⁻¹) from a PT split. The EoS returns the
/// *mixture* molar property directly, so the two-phase blend is the phase-fraction
/// weighted sum of the liquid- and vapor-phase mixture enthalpies.
fn mix_enthalpy<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    r: &FlashResult,
) -> Result<f64, ThermoError> {
    let hl = eos
        .molar_enthalpy(t, r.liquid_volume, &r.liquid_composition)?
        .value;
    let hv = eos
        .molar_enthalpy(t, r.vapor_volume, &r.vapor_composition)?
        .value;
    Ok((1.0 - r.vapor_fraction) * hl + r.vapor_fraction * hv)
}

/// Mixture molar entropy (J·mol⁻¹·K⁻¹) from a PT split.
fn mix_entropy<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    r: &FlashResult,
) -> Result<f64, ThermoError> {
    let sl = eos
        .molar_entropy(t, r.liquid_volume, &r.liquid_composition)?
        .value;
    let sv = eos
        .molar_entropy(t, r.vapor_volume, &r.vapor_composition)?
        .value;
    Ok((1.0 - r.vapor_fraction) * sl + r.vapor_fraction * sv)
}

/// Mixture molar volume (m³·mol⁻¹) from a PT split.
fn mix_volume(r: &FlashResult) -> f64 {
    r.mixture_molar_volume().value
}

/// Mixture molar internal energy (J·mol⁻¹): `U = H − P·V`.
fn mix_internal_energy<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    p: Pressure,
    r: &FlashResult,
) -> Result<f64, ThermoError> {
    let h = mix_enthalpy(eos, t, r)?;
    Ok(h - p.value * mix_volume(r))
}

/// Bisection root-find on a scalar (`T` or `P`) whose `eval` returns the residual
/// and the inner PT split. Brackets are expanded if the initial bracket does not
/// contain a sign change.
fn solve_outer<F>(mut var_lo: f64, mut var_hi: f64, mut eval: F) -> Result<FlashResult, FlashError>
where
    F: FnMut(f64) -> Result<(f64, FlashResult), FlashError>,
{
    let (mut ra, mut best_a) = eval(var_lo)?;
    let (mut rb, mut best_b) = eval(var_hi)?;
    // Expand the bracket until a sign change is found.
    for _ in 0..20 {
        if ra * rb <= 0.0 {
            break;
        }
        let new_hi = var_hi * 1.5;
        let (rb_new, bk) = eval(new_hi)?;
        var_hi = new_hi;
        rb = rb_new;
        best_b = bk;
    }
    for _ in 0..OUTER_MAX_ITER {
        if ra * rb <= 0.0 {
            let mid = 0.5 * (var_lo + var_hi);
            let (rm, fr) = eval(mid)?;
            if rm.abs() < OUTER_TOL {
                return Ok(fr);
            }
            if (var_hi - var_lo).abs() < 1e-9 {
                return Ok(fr);
            }
            if ra * rm <= 0.0 {
                var_hi = mid;
                rb = rm;
                best_b = fr;
            } else {
                var_lo = mid;
                ra = rm;
                best_a = fr;
            }
        } else {
            // No sign change found in the searchable bracket: return the closer split.
            return if ra.abs() < rb.abs() {
                Ok(best_a)
            } else {
                Ok(best_b)
            };
        }
    }
    if ra.abs() < rb.abs() {
        Ok(best_a)
    } else {
        Ok(best_b)
    }
}

// ---------------------------------------------------------------------------
// PH flash
// ---------------------------------------------------------------------------

/// PH flash implementation: specified molar enthalpy `h` at pressure `p`.
pub(crate) fn flash_ph_impl<E: EquationOfState + ?Sized>(
    calc: &FlashCalculator<'_, E>,
    h: MolarEnergy,
    p: Pressure,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    let eos = calc.eos_ref();
    let db = calc.db_opt();
    let nc = calc.comps();
    let target = h.value;
    let (t_lo, t_hi) = temperature_bracket(db, nc, z);
    let p_clone = p;
    solve_outer(t_lo, t_hi, |tk| {
        let t = Temperature::new::<kelvin>(tk);
        let r = flash_pt_impl(
            eos,
            db,
            nc,
            t,
            p_clone,
            z,
            crate::pt::PT_MAX_ITER,
            crate::pt::PT_TOL,
        )?;
        let res = mix_enthalpy(eos, t, &r).map_err(FlashError::Thermo)?;
        Ok((res - target, r))
    })
}

// ---------------------------------------------------------------------------
// PU flash
// ---------------------------------------------------------------------------

/// PU flash implementation: specified molar internal energy `u` at pressure `p`.
pub(crate) fn flash_pu_impl<E: EquationOfState + ?Sized>(
    calc: &FlashCalculator<'_, E>,
    u: MolarEnergy,
    p: Pressure,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    let eos = calc.eos_ref();
    let db = calc.db_opt();
    let nc = calc.comps();
    let target = u.value;
    let (t_lo, t_hi) = temperature_bracket(db, nc, z);
    let p_clone = p;
    solve_outer(t_lo, t_hi, |tk| {
        let t = Temperature::new::<kelvin>(tk);
        let r = flash_pt_impl(
            eos,
            db,
            nc,
            t,
            p_clone,
            z,
            crate::pt::PT_MAX_ITER,
            crate::pt::PT_TOL,
        )?;
        let res = mix_internal_energy(eos, t, p_clone, &r).map_err(FlashError::Thermo)?;
        Ok((res - target, r))
    })
}

// ---------------------------------------------------------------------------
// PV flash
// ---------------------------------------------------------------------------

/// PV flash implementation: specified molar volume `v` at pressure `p` (free
/// variable is temperature).
pub(crate) fn flash_pv_impl<E: EquationOfState + ?Sized>(
    calc: &FlashCalculator<'_, E>,
    p: Pressure,
    v: MolarVolume,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    let eos = calc.eos_ref();
    let db = calc.db_opt();
    let nc = calc.comps();
    let target = v.value;
    let (t_lo, t_hi) = temperature_bracket(db, nc, z);
    let p_clone = p;
    solve_outer(t_lo, t_hi, |tk| {
        let t = Temperature::new::<kelvin>(tk);
        let r = flash_pt_impl(
            eos,
            db,
            nc,
            t,
            p_clone,
            z,
            crate::pt::PT_MAX_ITER,
            crate::pt::PT_TOL,
        )?;
        Ok((mix_volume(&r) - target, r))
    })
}

// ---------------------------------------------------------------------------
// TS flash (vary P, fixed T)
// ---------------------------------------------------------------------------

/// TS flash implementation: specified molar entropy `s` at temperature `t`.
pub(crate) fn flash_ts_impl<E: EquationOfState + ?Sized>(
    calc: &FlashCalculator<'_, E>,
    t: Temperature,
    s: MolarEntropy,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    let eos = calc.eos_ref();
    let db = calc.db_opt();
    let nc = calc.comps();
    let target = s.value;
    let (p_lo, p_hi) = pressure_bracket(db, nc, z);
    let t_clone = t;
    solve_outer(p_lo, p_hi, |pk| {
        let p = Pressure::new::<uom::si::pressure::pascal>(pk);
        let r = flash_pt_impl(
            eos,
            db,
            nc,
            t_clone,
            p,
            z,
            crate::pt::PT_MAX_ITER,
            crate::pt::PT_TOL,
        )?;
        let res = mix_entropy(eos, t_clone, &r).map_err(FlashError::Thermo)?;
        Ok((res - target, r))
    })
}

// ---------------------------------------------------------------------------
// TV flash (vary P, fixed T)
// ---------------------------------------------------------------------------

/// TV flash implementation: specified molar volume `v` at temperature `t`.
pub(crate) fn flash_tv_impl<E: EquationOfState + ?Sized>(
    calc: &FlashCalculator<'_, E>,
    t: Temperature,
    v: MolarVolume,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    let eos = calc.eos_ref();
    let db = calc.db_opt();
    let nc = calc.comps();
    let target = v.value;
    let (p_lo, p_hi) = pressure_bracket(db, nc, z);
    let t_clone = t;
    solve_outer(p_lo, p_hi, |pk| {
        let p = Pressure::new::<uom::si::pressure::pascal>(pk);
        let r = flash_pt_impl(
            eos,
            db,
            nc,
            t_clone,
            p,
            z,
            crate::pt::PT_MAX_ITER,
            crate::pt::PT_TOL,
        )?;
        Ok((mix_volume(&r) - target, r))
    })
}

// ---------------------------------------------------------------------------
// Brackets
// ---------------------------------------------------------------------------

fn temperature_bracket(
    db: Option<&dyn tpt_thermo_core::component::ComponentDatabase>,
    nc: usize,
    z: &[f64],
) -> (f64, f64) {
    let mut tmin = f64::INFINITY;
    let mut tmax = 0.0_f64;
    if let Some(d) = db {
        for i in 0..nc {
            if i < z.len() && z[i] > 0.0 {
                if let Ok(tc) = d.critical_temperature(i) {
                    tmin = tmin.min(tc.value);
                    tmax = tmax.max(tc.value);
                }
            }
        }
    }
    if !tmin.is_finite() {
        tmin = 200.0;
        tmax = 1000.0;
    }
    (0.4 * tmin, 1.2 * tmax)
}

fn pressure_bracket(
    db: Option<&dyn tpt_thermo_core::component::ComponentDatabase>,
    nc: usize,
    z: &[f64],
) -> (f64, f64) {
    let mut pmin = f64::INFINITY;
    if let Some(d) = db {
        for i in 0..nc {
            if i < z.len() && z[i] > 0.0 {
                if let Ok(pc) = d.critical_pressure(i) {
                    pmin = pmin.min(pc.value);
                }
            }
        }
    }
    if !pmin.is_finite() {
        pmin = 5.0e6;
    }
    (1.0, (2.0 * pmin).max(5.0e6))
}

// ---------------------------------------------------------------------------
// Convenience free functions (build a calculator internally).
// ---------------------------------------------------------------------------

/// PH flash (free function). See [`FlashCalculator::flash_ph`].
pub fn flash_ph<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn tpt_thermo_core::component::ComponentDatabase>,
    h: MolarEnergy,
    p: Pressure,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    match db {
        Some(d) => FlashCalculator::with_db(eos, d).flash_ph(h, p, z),
        None => FlashCalculator::new(eos).flash_ph(h, p, z),
    }
}

/// TV flash (free function). See [`FlashCalculator::flash_tv`].
pub fn flash_tv<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn tpt_thermo_core::component::ComponentDatabase>,
    t: Temperature,
    v: MolarVolume,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    match db {
        Some(d) => FlashCalculator::with_db(eos, d).flash_tv(t, v, z),
        None => FlashCalculator::new(eos).flash_tv(t, v, z),
    }
}

/// TS flash (free function). See [`FlashCalculator::flash_ts`].
pub fn flash_ts<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn tpt_thermo_core::component::ComponentDatabase>,
    t: Temperature,
    s: MolarEntropy,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    match db {
        Some(d) => FlashCalculator::with_db(eos, d).flash_ts(t, s, z),
        None => FlashCalculator::new(eos).flash_ts(t, s, z),
    }
}

/// PU flash (free function). See [`FlashCalculator::flash_pu`].
pub fn flash_pu<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn tpt_thermo_core::component::ComponentDatabase>,
    u: MolarEnergy,
    p: Pressure,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    match db {
        Some(d) => FlashCalculator::with_db(eos, d).flash_pu(u, p, z),
        None => FlashCalculator::new(eos).flash_pu(u, p, z),
    }
}

/// PV flash (free function). See [`FlashCalculator::flash_pv`].
pub fn flash_pv<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn tpt_thermo_core::component::ComponentDatabase>,
    p: Pressure,
    v: MolarVolume,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    match db {
        Some(d) => FlashCalculator::with_db(eos, d).flash_pv(p, v, z),
        None => FlashCalculator::new(eos).flash_pv(p, v, z),
    }
}

/// Used to keep [`ConvergenceStatus`] referenced for downstream reporting.
#[allow(dead_code)]
fn _assert_status_in_scope(_s: ConvergenceStatus) {}

#[allow(dead_code)]
fn _assert_joule(_: MolarEnergy) -> f64 {
    MolarEnergy::new::<joule_per_mole>(1.0).value
}

#[allow(dead_code)]
fn _assert_vol(_: MolarVolume) -> f64 {
    MolarVolume::new::<cubic_meter_per_mole>(1.0).value
}
