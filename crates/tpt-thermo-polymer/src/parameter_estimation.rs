//! Fitting a Flory-Huggins `χ` parameter from experimental data.
//!
//! Currently supports recovery of `χ` from a single osmotic-pressure
//! measurement (the most common polymer-characterization route). Extension to
//! tie-line / activity data is tracked as Deferred Scope.

use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_core::ThermoError;
use tpt_thermo_core::R;

/// Recover the binary `χ` from an osmotic-pressure measurement of a polymer
/// (component 2, segment count `r2`) solution.
///
/// Uses the Flory-Huggins osmotic-pressure relation (solvent `r1 = 1`):
///
/// ```text
/// π·v* /(RT) = −(1/φ) [ ln(1−φ) + (1 − 1/r₂)·φ + χ·φ² ]
/// ```
///
/// where `φ` is the polymer volume fraction, `v*` the solvent characteristic
/// volume, and `π` the osmotic pressure.
pub fn chi_from_osmotic_pressure(
    polymer_volume_fraction: f64,
    osmotic_pressure: Pressure,
    t: Temperature,
    solvent_characteristic_volume: MolarVolume,
    r2: f64,
) -> Result<f64, ThermoError> {
    let phi = polymer_volume_fraction;
    if !(0.0..1.0).contains(&phi) {
        return Err(ThermoError::InvalidInput(
            "polymer volume fraction must be in (0,1)",
        ));
    }
    if r2 <= 0.0 {
        return Err(ThermoError::InvalidInput("segment count must be positive"));
    }
    let lhs = (osmotic_pressure.value * solvent_characteristic_volume.value) / (R * t.value);
    let num = -lhs * phi - (1.0 - phi).ln() - (1.0 - 1.0 / r2) * phi;
    Ok(num / (phi * phi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::molar_volume::cubic_meter_per_mole;
    use uom::si::pressure::pascal;
    use uom::si::thermodynamic_temperature::kelvin;

    #[test]
    fn round_trip_chi() {
        // Pick χ = 0.4, reconstruct π, then recover χ.
        let phi = 0.1_f64;
        let r2 = 1000.0;
        let chi = 0.4;
        let t = Temperature::new::<kelvin>(300.0);
        let vstar = MolarVolume::new::<cubic_meter_per_mole>(1.0e-5);
        let lhs = -(1.0 / phi) * ((1.0 - phi).ln() + (1.0 - 1.0 / r2) * phi + chi * phi * phi);
        let pi = Pressure::new::<pascal>(lhs * R * t.value / vstar.value);
        let recovered = chi_from_osmotic_pressure(phi, pi, t, vstar, r2).unwrap();
        assert!((recovered - chi).abs() < 1e-9);
    }
}
