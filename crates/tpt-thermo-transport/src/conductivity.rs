//! Thermal-conductivity correlations: Chung et al. (1988) gas and a
//! corresponding-states liquid estimate, with Filippov mixing.

use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{Pressure, Temperature, ThermalConductivity};
use uom::si::thermal_conductivity::watt_per_meter_kelvin;

use crate::viscosity::chung_gas_viscosity;

/// Ideal-gas molar isochoric heat capacity estimate (J·mol⁻¹·K⁻¹): a polyatomic
/// default of `2.5·R`, used by the Eucken-style conductivity closure.
fn cv_molar_estimate() -> f64 {
    2.5 * tpt_thermo_core::R
}

/// Dilute-gas thermal conductivity via Chung et al. (1988) / Eucken closure.
///
/// Computes the Chung gas viscosity, then applies the modified-Eucken relation
/// `λ = (η / M)·(Cv + 5/4·R)`. Returns W·m⁻¹·K⁻¹.
pub fn chung_gas_thermal_conductivity(
    db: &dyn ComponentDatabase,
    t: Temperature,
    molar_density: f64,
    z: &[f64],
) -> Result<ThermalConductivity, ThermoError> {
    let eta = chung_gas_viscosity(db, t, molar_density, z)?;
    let m_mix: f64 = z
        .iter()
        .zip(0..db.num_components())
        .map(|(xi, i)| xi * db.molar_mass(i).map(|m| m.value).unwrap_or(0.0))
        .sum();
    if m_mix <= 0.0 {
        return Err(ThermoError::InvalidInput("zero mixture molar mass"));
    }
    let lambda = (eta.value / m_mix) * (cv_molar_estimate() + 1.25 * tpt_thermo_core::R);
    Ok(ThermalConductivity::new::<watt_per_meter_kelvin>(lambda))
}

/// Liquid thermal conductivity (Poling et al. corresponding-states form).
///
/// `λ_liq = (Tc^(5/6) / (M^0.5 · Pc^(2/3))) · 0.0227 · (1 − Tr)^0.38` (SI, with
/// `Tc` K, `M` g·mol⁻¹, `Pc` bar), mole-fraction averaged over components for a
/// mixture. Order-of-magnitude estimate for non-associating liquids.
pub fn liquid_thermal_conductivity(
    db: &dyn ComponentDatabase,
    t: Temperature,
    z: &[f64],
) -> Result<ThermalConductivity, ThermoError> {
    let tk = t.value;
    let n = db.num_components();
    if z.len() != n {
        return Err(ThermoError::InvalidInput("feed length mismatch"));
    }
    let mut sum = 0.0_f64;
    let mut total = 0.0_f64;
    for i in 0..n {
        if z[i] <= 0.0 {
            continue;
        }
        let tc = db.critical_temperature(i)?.value;
        let pc_bar = db.critical_pressure(i)?.value / 1.0e5;
        let m_g = db.molar_mass(i)?.value * 1000.0;
        let tr = tk / tc;
        if tr >= 1.0 || tc <= 0.0 || pc_bar <= 0.0 {
            return Err(ThermoError::Numerical(
                tpt_thermo_core::convergence::ConvergenceStatus::NumericalIssue(
                    tpt_thermo_core::convergence::NumericalIssueReason::OutOfDomain,
                ),
            ));
        }
        let lambda_i =
            (tc.powf(5.0 / 6.0) / (m_g.sqrt() * pc_bar.powf(2.0 / 3.0))) * 0.227
                * (1.0 - tr).powf(0.38);
        sum += z[i] * lambda_i;
        total += z[i];
    }
    if total <= 0.0 {
        return Err(ThermoError::InvalidInput("no positive mole fractions"));
    }
    Ok(ThermalConductivity::new::<watt_per_meter_kelvin>(sum / 1.0e3))
}

/// Filippov (binary) mixture thermal conductivity (W·m⁻¹·K⁻¹).
pub fn filippov_thermal_conductivity(x1: f64, lambda1: f64, lambda2: f64) -> f64 {
    x1 * lambda1 + (1.0 - x1) * lambda2 - 0.72 * x1 * (1.0 - x1) * (lambda1 - lambda2)
}

/// Convenience: ideal-gas molar density (mol·m⁻³) at `(T, P)`.
pub fn ideal_molar_density(t: Temperature, p: Pressure) -> f64 {
    p.value / (tpt_thermo_core::R * t.value)
}
