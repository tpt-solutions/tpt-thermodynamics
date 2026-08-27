//! Workspace example crate for `tpt-thermodynamics`.
//!
//! This crate is a scratch pad for runnable, documented example code that
//! exercises the `tpt-thermo-*` crates as each Phase 2-13 lands. It is **not
//! published** and is excluded from the packaged workspace surface.
//!
//! Examples are added per phase so that every public API has at least one
//! executable, tested showcase.

/// Phase 2 examples: composition conversion and the ideal-gas EoS toy.
pub mod phase2 {
    use tpt_thermo_core::{
        composition::Composition,
        eos::{EquationOfState, IdealGas},
        quantities::{MolarHeatCapacity, MolarMass, MolarVolume, Pressure, Temperature},
    };
    use uom::si::{
        molar_heat_capacity::joule_per_kelvin_mole, molar_mass::kilogram_per_mole,
        molar_volume::cubic_meter_per_mole, pressure::pascal as pascal_unit,
        thermodynamic_temperature::kelvin,
    };

    /// Convert a mole-fraction composition to mass fractions and back, using
    /// the `tpt-thermo-core` [`Composition`] helper.
    pub fn mole_to_mass_round_trip() -> (Vec<f64>, Vec<f64>) {
        // 50/50 mole methane (16 g/mol) / water-like (18 g/mol).
        let c = Composition::from_mole_fractions(vec![0.5, 0.5]).unwrap();
        let mm = [
            MolarMass::new::<kilogram_per_mole>(0.016),
            MolarMass::new::<kilogram_per_mole>(0.018),
        ];
        let mass = c.mass_fractions(&mm).unwrap();
        let back = Composition::from_mass_fractions(mass.clone(), &mm).unwrap();
        (mass, back.mole_fractions().to_vec())
    }

    /// Ideal-gas P-V-T: recover pressure from `P = R T / v` for a 300 K sample
    /// at `R·300/1e5` m³·mol⁻¹ (so the expected pressure is 1 bar).
    pub fn ideal_gas_pressure() -> f64 {
        let m = IdealGas::new(
            1,
            MolarHeatCapacity::new::<joule_per_kelvin_mole>(33.6),
            MolarMass::new::<kilogram_per_mole>(0.028),
            Temperature::new::<kelvin>(298.15),
            Pressure::new::<pascal_unit>(1.0e5),
        );
        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(tpt_thermo_core::R * 300.0 / 1.0e5);
        m.pressure(t, v, &[1.0]).unwrap().value
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn composition_round_trips() {
            let (mass, back) = mole_to_mass_round_trip();
            assert!((mass[0] - 0.016 / 0.034).abs() < 1e-12);
            assert!((back[0] - 0.5).abs() < 1e-12);
        }

        #[test]
        fn ideal_gas_matches_rt_over_v() {
            let p = ideal_gas_pressure();
            assert!((p - 1.0e5).abs() / 1.0e5 < 1e-9);
        }
    }
}
