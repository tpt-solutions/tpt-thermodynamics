//! K-value initialisation (Wilson correlation) for the flash iterations.

use alloc::vec::Vec;
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::Pressure;

/// Wilson correlation for initial K-values:
/// `K_i = (Pc_i / P) · exp(5.3727·(1 + ω_i)·(1 − Tc_i/T))`.
///
/// Uses the critical constants and acentric factors from `db`. Returns one K-value
/// per component. `z` is only used to size the output and need not be normalised.
pub fn wilson_k_values(
    t_kelvin: f64,
    p: Pressure,
    db: &dyn ComponentDatabase,
    z: &[f64],
) -> Result<Vec<f64>, ThermoError> {
    let nc = db.num_components();
    if z.len() != nc {
        return Err(ThermoError::InvalidInput("feed length mismatch"));
    }
    let tk = if t_kelvin > 0.0 { t_kelvin } else { 1.0 };
    let p_val = p.value.max(1.0);
    let mut k = alloc::vec![1.0_f64; nc];
    for (i, ki) in k.iter_mut().enumerate() {
        let tc = db.critical_temperature(i)?.value;
        let pc = db.critical_pressure(i)?.value;
        let omega = db.acentric_factor(i)?;
        let exponent = 5.3727 * (1.0 + omega) * (1.0 - tc / tk);
        *ki = (pc / p_val) * exponent.exp();
    }
    Ok(k)
}
