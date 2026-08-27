//! Solid–liquid equilibrium (SLE): ideal melting-point-depression solubility.

use tpt_thermo_core::quantities::Temperature;
use tpt_thermo_core::R;

/// Ideal (melting-point-depression) solubility of a solute at temperature `t`,
/// given its melting temperature `t_melt` and fusion enthalpy `dh_fus`.
///
/// `x = exp[ −ΔH_fus/R·(1/T − 1/T_m) ]`, clamped to `[0, 1]`. Returns the solid
/// mole fraction `x_solid = 1 − x`? No — this returns the *liquid-phase* solute
/// mole fraction `x` (the solubility). At `T = T_m` it is 1 (fully miscible);
/// below `T_m` it falls below 1.
pub fn solid_liquid_solubility(t_melt: Temperature, dh_fus: tpt_thermo_core::quantities::MolarEnergy, t: Temperature) -> f64 {
    let tm = t_melt.value;
    let tt = t.value;
    if tt <= 0.0 || tm <= 0.0 {
        return 0.0;
    }
    // x = exp[ −ΔH_fus/R · (1/T − 1/T_m) ]; below T_m this is < 1.
    let arg = (-dh_fus.value / (R * tt)) * (1.0 / tt - 1.0 / tm);
    arg.exp().clamp(0.0, 1.0)
}

/// Temperature-dependent solubility extension including a heat-capacity change
/// `dc_p` (J·mol⁻¹·K⁻¹) between solid and liquid, via the integrated van't Hoff
/// relation with a `ΔC_p` correction term.
pub fn solid_liquid_solubility_tdependent(
    t_melt: Temperature,
    dh_fus: tpt_thermo_core::quantities::MolarEnergy,
    dc_p: f64,
    t: Temperature,
) -> f64 {
    let tm = t_melt.value;
    let tt = t.value;
    if tt <= 0.0 || tm <= 0.0 {
        return 0.0;
    }
    let dhm = dh_fus.value + dc_p * (tt - tm);
    let term = (-dhm / (R * tt)) * (1.0 / tt - 1.0 / tm)
        + (dc_p / R) * ((tt / tm).ln() - (1.0 - tm / tt));
    term.exp().clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::molar_energy::joule_per_mole;
    use uom::si::thermodynamic_temperature::kelvin;

    #[test]
    fn solubility_unity_at_melt_then_decreases() {
        let tm = Temperature::new::<kelvin>(350.0);
        let dh = tpt_thermo_core::quantities::MolarEnergy::new::<joule_per_mole>(10.0e3);
        let at_melt = solid_liquid_solubility(tm, dh, tm);
        let below = solid_liquid_solubility(tm, dh, Temperature::new::<kelvin>(300.0));
        let above = solid_liquid_solubility(tm, dh, Temperature::new::<kelvin>(400.0));
        assert!((at_melt - 1.0).abs() < 1e-12, "miscible at melting point");
        assert!(below < 1.0 && below > 0.0, "subsaturated below melt");
        assert!((above - 1.0).abs() < 1e-12, "miscible above melt (clamped)");
    }
}
