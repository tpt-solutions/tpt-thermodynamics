//! Crystallization: Flory melting-point depression.
//!
//! For a semicrystalline polymer diluted by a solvent (or diluent), the melting
//! point `T_m` is depressed relative to the pure-polymer melting point `T_m^0`
//! according to the Flory equation:
//!
//! ```text
//! 1/T_m − 1/T_m^0 = −(R/Δh_f) [ ln φ₂ + (1 − 1/r)·(1 − φ₂) ]
//! ```
//!
//! where `φ₂` is the polymer volume fraction, `r` the degree of polymerization,
//! `Δh_f` the (molar) heat of fusion per repeat unit, and `R` the gas constant.

use tpt_thermo_core::quantities::{MolarEnergy, Temperature};
use tpt_thermo_core::R;

/// Flory melting-point depression. Returns the depressed melting temperature.
pub fn flory_melting_depression(
    pure_melting_point: Temperature,
    heat_of_fusion_per_repeat_unit: MolarEnergy,
    polymer_volume_fraction: f64,
    degree_of_polymerization: f64,
) -> Temperature {
    let phi = polymer_volume_fraction.clamp(1e-9, 1.0 - 1e-9);
    let r = degree_of_polymerization.max(1.0);
    let bracket = phi.ln() + (1.0 - 1.0 / r) * (1.0 - phi);
    let inv_tm =
        1.0 / pure_melting_point.value - (R / heat_of_fusion_per_repeat_unit.value) * bracket;
    Temperature::new::<uom::si::thermodynamic_temperature::kelvin>(1.0 / inv_tm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::molar_energy::joule_per_mole;
    use uom::si::thermodynamic_temperature::kelvin;

    #[test]
    fn depresses_with_dilution() {
        let tm0 = Temperature::new::<kelvin>(400.0);
        let dh = MolarEnergy::new::<joule_per_mole>(1.0e7);
        let tm = flory_melting_depression(tm0, dh, 0.9, 1000.0);
        assert!(tm.value < 400.0, "melting point should be depressed");
        // As r → ∞, the (1 − 1/r) term → 1; verify monotonic with φ.
        let tm_low = flory_melting_depression(tm0, dh, 0.5, 1.0e6);
        let tm_high = flory_melting_depression(tm0, dh, 0.95, 1.0e6);
        assert!(
            tm_low.value < tm_high.value,
            "more polymer ⇒ less depression"
        );
    }

    #[test]
    fn pure_solvent_limit() {
        let tm0 = Temperature::new::<kelvin>(400.0);
        let dh = MolarEnergy::new::<joule_per_mole>(1.0e7);
        // φ → 1 (no diluent) ⇒ no depression.
        let tm = flory_melting_depression(tm0, dh, 1.0 - 1e-9, 1000.0);
        assert!((tm.value - 400.0).abs() < 1e-3);
    }
}
