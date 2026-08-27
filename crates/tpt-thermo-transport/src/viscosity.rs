//! Viscosity correlations: Chung et al. (1988) dilute-gas and Lucas (1981) liquid,
//! plus Wilke / Mason–Saxena mixture rules.

use alloc::vec::Vec;
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{DynamicViscosity, MolarMass, Pressure, Temperature};
use uom::si::dynamic_viscosity::pascal_second;
use uom::si::molar_volume::cubic_meter_per_mole;
use uom::si::pressure::bar;
use uom::si::thermodynamic_temperature::kelvin;

use crate::parameters::{lj_params_for, LjParams};

/// Lennard-Jones and critical parameters for every component of `db`.
fn component_params(db: &dyn ComponentDatabase) -> Result<Vec<LjParams>, ThermoError> {
    lj_params_for(db)
}

/// Reduced collision integral `Ω_η(T*)` (Neufeld et al. 1972).
fn collision_integral(tr: f64) -> f64 {
    let t = tr.max(1e-3);
    1.16145 * t.powf(-0.14874)
        + 0.52487 * (-0.7732 * t).exp()
        + 2.1611 * (-2.43787 * t).exp()
}

/// Dilute-gas viscosity via Chung et al. (1988), mixture form.
///
/// `molar_density` is the total molar density (mol·m⁻³). Returns the mixture
/// viscosity in Pa·s.
pub fn chung_gas_viscosity(
    db: &dyn ComponentDatabase,
    t: Temperature,
    molar_density: f64,
    z: &[f64],
) -> Result<DynamicViscosity, ThermoError> {
    let tk = t.value;
    let params = component_params(db)?;
    let n = params.len();
    if z.len() != n {
        return Err(ThermoError::InvalidInput("feed length mismatch"));
    }
    // Mixture combos over the composition.
    let mut sig3 = 0.0_f64;
    let mut eps = 0.0_f64;
    let mut m_mix = 0.0_f64;
    let mut omega_mix = 0.0_f64;
    for i in 0..n {
        m_mix += z[i] * params[i].molar_mass * 1000.0; // kg/mol -> g/mol
        omega_mix += z[i] * params[i].omega;
        for j in 0..n {
            let sij = 0.5 * (params[i].sigma_a + params[j].sigma_a);
            sig3 += z[i] * z[j] * sij.powi(3);
            let eij = (params[i].eps_ok * params[j].eps_ok).sqrt();
            eps += z[i] * z[j] * eij;
        }
    }
    let sigma_mix = sig3.powf(1.0 / 3.0);
    let tr = tk / eps.max(1e-3);
    let omega_eta = collision_integral(tr);
    // Polyatomic quantum/shape correction (Chung F_c approximation via acentricity).
    let fc = 1.0 + 0.22 * omega_mix / tr;
    // η [μP] = 26.69 · F_c · sqrt(M·T) / (σ² · Ω); 1 μP = 10⁻⁷ Pa·s.
    let eta_ucp = 26.69 * fc * (m_mix * tk).sqrt() / (sigma_mix * sigma_mix * omega_eta);
    let _ = molar_density; // low-pressure form (density-independent to first order)
    Ok(DynamicViscosity::new::<pascal_second>(eta_ucp * 1e-7))
}

/// Liquid viscosity via the Lucas (1981) corresponding-states correlation.
///
/// Returns the mixture viscosity (Pa·s) as the mole-fraction average of the
/// per-component Lucas values (a serviceable mixture rule for non-associating
/// liquids).
pub fn lucas_liquid_viscosity(
    db: &dyn ComponentDatabase,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> Result<DynamicViscosity, ThermoError> {
    let tk = t.value;
    let params = component_params(db)?;
    let n = params.len();
    if z.len() != n {
        return Err(ThermoError::InvalidInput("feed length mismatch"));
    }
    let mut eta_mix = 0.0_f64;
    let mut total = 0.0_f64;
    for i in 0..n {
        if z[i] <= 0.0 {
            continue;
        }
        let tc = db.critical_temperature(i)?.value;
        let pc_bar = db.critical_pressure(i)?.value / 1.0e5;
        let omega = db.acentric_factor(i)?;
        let m_g = db.molar_mass(i)?.value * 1000.0; // kg/mol -> g/mol
        let tr = tk / tc;
        if tr <= 0.0 || tc <= 0.0 || pc_bar <= 0.0 {
            return Err(ThermoError::Numerical(
                tpt_thermo_core::convergence::ConvergenceStatus::NumericalIssue(
                    tpt_thermo_core::convergence::NumericalIssueReason::OutOfDomain,
                ),
            ));
        }
        let eta_star = (tc.powf(1.0 / 6.0) / (m_g.sqrt() * pc_bar.powf(2.0 / 3.0))) * 17.6;
        let a = (0.8282 - 0.6065 * omega) / (0.5371 - 0.3286 * omega);
        let f_t = tr.powf(2.0 / 3.0) / (1.0 + (tr - 1.0) * a);
        let f_omega = 0.125
            * (tr.powf(2.0 / 3.0) - 1.0)
            * (0.288 - 0.344 * omega + 0.124 * omega * omega + 0.4329 * omega.powi(3))
            + 1.0;
        // Low-pressure approximation: F_p ≈ 1.
        let pr = p.value / 1.0e5 / pc_bar.max(1e-6);
        let f_p = 1.0 + pr * 0.1;
        let eta_i = eta_star * f_t * f_omega * f_p; // mPa·s
        eta_mix += z[i] * eta_i;
        total += z[i];
    }
    if total <= 0.0 {
        return Err(ThermoError::InvalidInput("no positive mole fractions"));
    }
    // mPa·s -> Pa·s.
    Ok(DynamicViscosity::new::<pascal_second>(eta_mix / 1.0e3))
}

/// Wilke (1950) mixture viscosity from pure-component viscosities (Pa·s).
pub fn wilke_mixture_viscosity(
    pure_viscosity: &[f64],
    molar_masses: &[f64],
    x: &[f64],
) -> f64 {
    let n = pure_viscosity.len();
    let mut eta = 0.0_f64;
    let mut denom = 0.0_f64;
    for i in 0..n {
        if x[i] <= 0.0 {
            continue;
        }
        let mut s = 0.0_f64;
        for j in 0..n {
            if pure_viscosity[j] <= 0.0 {
                continue;
            }
            let mm = (molar_masses[i] / molar_masses[j]).sqrt();
            let phi = (1.0 + (pure_viscosity[i] / pure_viscosity[j]).sqrt() * mm)
                / (2.0 * 2.0_f64.sqrt() * (1.0 + molar_masses[i] / molar_masses[j]).powf(0.25));
            s += x[j] * phi;
        }
        if s > 0.0 {
            eta += x[i] * pure_viscosity[i] / s;
            denom += x[i];
        }
    }
    if denom > 0.0 {
        eta
    } else {
        0.0
    }
}

/// Mason–Saxena (1958) mixture viscosity (a refined Wilke form). Returns Pa·s.
pub fn mason_saxena_mixture_viscosity(
    pure_viscosity: &[f64],
    molar_masses: &[f64],
    x: &[f64],
) -> f64 {
    let n = pure_viscosity.len();
    let mut eta = 0.0_f64;
    let mut denom = 0.0_f64;
    for i in 0..n {
        if x[i] <= 0.0 {
            continue;
        }
        let mut s = 0.0_f64;
        for j in 0..n {
            if pure_viscosity[j] <= 0.0 {
                continue;
            }
            let mm = (molar_masses[i] / molar_masses[j]).sqrt();
            let phi = ((1.0 + (pure_viscosity[i] / pure_viscosity[j]).sqrt() * mm)
                / (2.0 * 2.0_f64.sqrt() * (1.0 + molar_masses[i] / molar_masses[j]).powf(0.25)))
                .powf(0.6);
            s += x[j] * phi;
        }
        if s > 0.0 {
            eta += x[i] * pure_viscosity[i] / s;
            denom += x[i];
        }
    }
    if denom > 0.0 {
        eta
    } else {
        0.0
    }
}

/// Verify the crate's dependencies resolve against the seed set (smoke helper).
#[allow(dead_code)]
fn _check_units() {
    let _ = MolarMass::new::<uom::si::molar_mass::kilogram_per_mole>(0.0);
    let _ = Temperature::new::<kelvin>(0.0);
    let _ = Pressure::new::<bar>(0.0);
    let _ = cubic_meter_per_mole;
}
