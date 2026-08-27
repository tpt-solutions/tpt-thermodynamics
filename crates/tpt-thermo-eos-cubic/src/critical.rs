//! Critical-point, spinodal, and mechanical-stability helpers.

use crate::engine::CubicEos;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_core::{EquationOfState, ThermoError};
use uom::si::molar_volume::cubic_meter_per_mole;

/// Pure-component critical point `(T_c, P_c, v_c)` of a [`CubicEos`].
pub fn critical_point(
    eos: &CubicEos,
    i: usize,
) -> Result<(Temperature, Pressure, MolarVolume), ThermoError> {
    eos.critical_point_pure(i)
}

/// Mechanical stability test at `(T, v, z)`: a single phase is mechanically
/// stable iff `(∂P/∂v)_T < 0`.
pub fn mechanical_stability<E: EquationOfState>(
    eos: &E,
    t: Temperature,
    v: MolarVolume,
    z: &[f64],
) -> Result<bool, ThermoError> {
    let dpdv = eos.dp_dv(t, v, z)?;
    Ok(dpdv < 0.0)
}

/// Locate the liquid and vapor spinodal volumes (where `(∂P/∂v)_T = 0`) at
/// subcritical `T` for a pure component, by scanning a wide molar-volume range
/// for sign changes of `(∂P/∂v)_T`. Returns `None` if only a single
/// (supercritical) spinodal exists.
pub fn spinodal_roots<E: EquationOfState>(
    eos: &E,
    t: Temperature,
    z: &[f64],
) -> Result<Option<(MolarVolume, MolarVolume)>, ThermoError> {
    let steps = 4000usize;
    let v_lo = 1e-6;
    let v_hi = 1.0;
    let mut crossings = alloc::vec::Vec::new();
    let mut prev_v = v_lo;
    let mut prev_dp = eos.dp_dv(t, MolarVolume::new::<cubic_meter_per_mole>(v_lo), z)?;
    for k in 1..=steps {
        let v = v_lo + (v_hi - v_lo) * (k as f64) / (steps as f64);
        let dp = eos.dp_dv(t, MolarVolume::new::<cubic_meter_per_mole>(v), z)?;
        if prev_dp <= 0.0 && dp > 0.0 {
            // dp/dv crosses from negative to positive: liquid spinodal (inner).
            crossings.push(prev_v);
        } else if prev_dp >= 0.0 && dp < 0.0 {
            // Crosses back negative: vapor spinodal (outer).
            crossings.push(v);
        }
        prev_v = v;
        prev_dp = dp;
    }
    if crossings.len() >= 2 {
        Ok(Some((
            MolarVolume::new::<cubic_meter_per_mole>(crossings[0]),
            MolarVolume::new::<cubic_meter_per_mole>(crossings[1]),
        )))
    } else {
        Ok(None)
    }
}
