//! PC-SAFT for polymers: a thin specialisation of [`tpt_thermo_eos_saft::PcSaft`].
//!
//! Polymer segments are described by the same SAFT parameters (`m`, `σ`, `ε/k`) used
//! by [`tpt_thermo_eos_saft`]; this module just provides a convenience constructor
//! that builds a `PcSaft` directly from a polymer-segment specification and proves
//! (via the [`tests`]) that it reduces exactly to plain PC-SAFT for an equivalent
//! pure component.

use tpt_thermo_core::EquationOfState;
use tpt_thermo_eos_saft::parameters::{SaftComponent, SaftParameters};
use tpt_thermo_eos_saft::PcSaft;

/// A single polymer (or segment) SAFT specification.
#[derive(Debug, Clone, Copy)]
pub struct PolymerSaftSpec {
    /// Component name (matched against the seed database where relevant).
    pub name: &'static str,
    /// Number of segments per chain `m`.
    pub m: f64,
    /// Segment diameter `σ` (Å).
    pub sigma: f64,
    /// Segment energy `ε/k` (K).
    pub epsilon_k: f64,
}

/// Build a [`PcSaft`] from polymer-segment specs and the corresponding molar masses.
pub fn build_pc_saft_polymer(specs: &[PolymerSaftSpec], molar_masses: &[f64]) -> PcSaft {
    assert_eq!(specs.len(), molar_masses.len(), "spec/mass length mismatch");
    let comps: Vec<SaftComponent> = specs
        .iter()
        .map(|s| SaftComponent::pc_saft(s.name, s.m, s.sigma, s.epsilon_k))
        .collect();
    PcSaft::new(SaftParameters::new(comps), molar_masses.to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
    use tpt_thermo_eos_saft::parameters::SEED_SAFT_PARAMETERS;
    use uom::si::{molar_volume::cubic_meter_per_mole, pressure::pascal, thermodynamic_temperature::kelvin};

    #[test]
    fn reduces_to_plain_pc_saft() {
        // A methane-like polymer segment should reproduce plain PC-SAFT built from the
        // seed table for the same pure component.
        let seed_methane = SEED_SAFT_PARAMETERS
            .iter()
            .find(|c| c.name == "methane")
            .copied()
            .unwrap();
        let polymer = build_pc_saft_polymer(
            &[PolymerSaftSpec {
                name: "methane",
                m: seed_methane.m,
                sigma: seed_methane.sigma,
                epsilon_k: seed_methane.epsilon_k,
            }],
            &[0.016_043],
        );
        let plain = PcSaft::new(
            SaftParameters::new(vec![seed_methane]),
            vec![0.016_043],
        );

        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(0.025);
        let p_poly = polymer.pressure(t, v, &[1.0]).unwrap();
        let p_plain = plain.pressure(t, v, &[1.0]).unwrap();
        assert!((p_poly.value - p_plain.value).abs() / p_plain.value < 1e-12);

        let phi_poly = polymer.ln_fugacity_coefficient(t, v, &[1.0], 0).unwrap();
        let phi_plain = plain.ln_fugacity_coefficient(t, v, &[1.0], 0).unwrap();
        assert!((phi_poly - phi_plain).abs() < 1e-12);
    }

    #[test]
    fn polymer_eos_is_usable() {
        // A long-chain polymer (large m) should be buildable and evaluate.
        let eos = build_pc_saft_polymer(
            &[PolymerSaftSpec {
                name: "polymer-a",
                m: 50.0,
                sigma: 4.0,
                epsilon_k: 300.0,
            }],
            &[0.5],
        );
        let t = Temperature::new::<kelvin>(400.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(0.05);
        let p = eos.pressure(t, v, &[1.0]).unwrap();
        assert!(p.value.is_finite());
        let _ = Pressure::new::<pascal>(0.0);
    }
}
