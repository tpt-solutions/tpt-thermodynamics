//! Diffusivity correlations: Fuller–Schettler–Giddings (gas) and Vignes / Darken
//! (liquid).

use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{DiffusionCoefficient, Pressure, Temperature};
use uom::si::diffusion_coefficient::square_meter_per_second;

/// Fuller–Schettler–Giddings diffusion volume `Σv` (cm³·mol⁻¹)^(1/3)-scale) for a
/// named component. Covers the seed set; unknown names fall back to a rough
/// constant.
pub fn diffusion_volume(name: &str) -> f64 {
    match name.trim().to_lowercase().as_str() {
        "hydrogen" | "h2" => 7.07,
        "helium" | "he" => 2.88,
        "nitrogen" | "n2" => 17.9,
        "oxygen" | "o2" => 16.6,
        "argon" | "ar" => 16.2,
        "carbon monoxide" | "co" => 18.9,
        "carbon dioxide" | "co2" => 26.9,
        "water" | "h2o" => 12.7,
        "ammonia" | "nh3" => 11.47,
        "hydrogen sulfide" | "h2s" => 27.52,
        "hydrogen chloride" | "hcl" => 21.81,
        "methane" => 25.14,
        "ethane" => 45.66,
        "propane" => 66.18,
        "n-butane" => 86.7,
        "n-pentane" => 107.2,
        "n-hexane" => 127.7,
        "n-heptane" => 148.2,
        "n-octane" => 168.8,
        "ethylene" | "ethene" => 41.04,
        "propylene" | "propene" => 61.56,
        "benzene" => 89.06,
        "toluene" => 129.8,
        "methanol" => 33.56,
        "ethanol" => 51.77,
        _ => 20.0,
    }
}

/// Fuller–Schettler–Giddings binary gas diffusivity `D_AB` (m²·s⁻¹).
pub fn fuller_schettler_giddings(
    db: &dyn ComponentDatabase,
    t: Temperature,
    p: Pressure,
    i: usize,
    j: usize,
) -> Result<DiffusionCoefficient, ThermoError> {
    let tk = t.value;
    let p_atm = p.value / 1.01325e5;
    let mi = db.molar_mass(i)?.value * 1000.0; // kg/mol -> g/mol
    let mj = db.molar_mass(j)?.value * 1000.0;
    let ni = db.name(i)?;
    let nj = db.name(j)?;
    let vi = diffusion_volume(ni);
    let vj = diffusion_volume(nj);
    let s = (vi.powf(1.0 / 3.0) + vj.powf(1.0 / 3.0)).powi(2);
    let d_cm2_s =
        (1.0133e-3 * tk.powf(1.75) * (1.0 / mi + 1.0 / mj).sqrt()) / (p_atm * s);
    Ok(DiffusionCoefficient::new::<square_meter_per_second>(
        d_cm2_s * 1.0e-4,
    ))
}

/// Vignes (1966) liquid interdiffusion: `D = D12^{x2} · D21^{x1}`.
pub fn vignes_liquid_binary(d12: f64, d21: f64, x1: f64) -> f64 {
    if d12 <= 0.0 || d21 <= 0.0 {
        return 0.0;
    }
    d12.powf(1.0 - x1) * d21.powf(x1)
}

/// Darken (1948) interdiffusion: `D = x2·D1^* + x1·D2^*`.
pub fn darken_liquid_binary(d1_star: f64, d2_star: f64, x1: f64) -> f64 {
    (1.0 - x1) * d1_star + x1 * d2_star
}

/// Convenience wrapper returning the binary gas diffusivity in m²·s⁻¹ as a plain
/// `f64` (used by callers that do not need the quantity type).
pub fn fuller_schettler_giddings_value(
    db: &dyn ComponentDatabase,
    t: Temperature,
    p: Pressure,
    i: usize,
    j: usize,
) -> Result<f64, ThermoError> {
    Ok(fuller_schettler_giddings(db, t, p, i, j)?.value)
}

#[allow(dead_code)]
fn _unused(_: &[f64]) {}
