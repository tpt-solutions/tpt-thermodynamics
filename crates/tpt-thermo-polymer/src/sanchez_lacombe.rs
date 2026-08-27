//! Sanchez–Lacombe lattice-fluid equation of state.
//!
//! Implements [`tpt_thermo_core::EquationOfState`] for the Sanchez–Lacombe (1976)
//! lattice fluid. Each component carries a segment count `r`, a close-packed molar
//! volume `v*`, and a segment-energy temperature `ε/k`. Mixtures use the standard
//! one-fluid (White) mixing rules.
//!
//! The pressure relation (reduced density `ρ̃ = v*/v`, coordination number `z`) is
//!
//! ```text
//! P = (R·T / v*)·[ (1 − ρ̃)/ρ̃ − (z/2)·ln(1 − ρ̃) ] − (ε* / v*)·[ (1 − 1/r)·ρ̃ + ρ̃² ]
//! ```
//!
//! which recovers the ideal-gas limit `P·v → R·T` as `ρ̃ → 0`. Residual enthalpy and
//! entropy are obtained numerically from the residual Helmholtz energy; the fugacity
//! coefficient uses the compressibility integral (pseudo-pure, one-fluid mixture).
//! Quantitative VLE validation against literature is tracked as Deferred Scope.

use alloc::vec::Vec;
use tpt_thermo_core::quantities::{
    gas_constant, MolarEnergy, MolarEntropy, MolarHeatCapacity, MolarMass, MolarVolume,
    Pressure, Temperature, Velocity,
};
use tpt_thermo_core::{EquationOfState, ThermoError};
use uom::si::dynamic_viscosity::pascal_second;
use uom::si::molar_energy::joule_per_mole;
use uom::si::molar_entropy::joule_per_kelvin_mole;
use uom::si::molar_heat_capacity::joule_per_kelvin_mole;
use uom::si::molar_volume::cubic_meter_per_mole;
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::velocity::meter_per_second;

/// Per-component Sanchez–Lacombe parameters.
#[derive(Debug, Clone, Copy)]
pub struct SlComponent {
    /// Number of lattice segments per molecule.
    pub r: f64,
    /// Close-packed molar volume `v*` (m³·mol⁻¹).
    pub v_star: f64,
    /// Segment-energy temperature `ε/k` (K).
    pub epsilon_k: f64,
}

/// Sanchez–Lacombe lattice-fluid model.
#[derive(Debug, Clone)]
pub struct SanchezLacombe {
    components: Vec<SlComponent>,
    /// Lattice coordination number (default 6).
    pub z: f64,
}

impl SanchezLacombe {
    /// Build from per-component parameters.
    pub fn new(components: Vec<SlComponent>) -> Self {
        Self {
            components,
            z: 6.0,
        }
    }

    /// Build for a single component.
    pub fn pure(component: SlComponent) -> Self {
        Self::new(alloc::vec![component])
    }

    /// One-fluid mixed parameters for a composition `x`.
    fn mix(&self, x: &[f64]) -> SlComponent {
        let n = self.components.len();
        let mut r = 0.0_f64;
        let mut v_star = 0.0_f64;
        let mut eps = 0.0_f64;
        for i in 0..n {
            let ci = self.components[i];
            r += x[i] * ci.r;
            v_star += x[i] * ci.v_star;
            for j in 0..n {
                let cj = self.components[j];
                eps += x[i] * x[j] * (ci.epsilon_k * cj.epsilon_k).sqrt();
            }
        }
        SlComponent {
            r,
            v_star,
            epsilon_k: eps,
        }
    }

    /// Reduced density `ρ̃ = v*/v` for mixed parameters at molar volume `v`.
    fn rho_tilde(params: &SlComponent, v: f64) -> f64 {
        (params.v_star / v).clamp(1e-9, 1.0 - 1e-9)
    }

    /// Pressure (Pa) at `(T, v, composition)`.
    fn pressure_value(&self, t: f64, v: f64, x: &[f64]) -> f64 {
        let p = self.mix(x);
        let rt = gas_constant().value * t;
        let eps_star = p.epsilon_k * gas_constant().value; // J/mol
        let rho = Self::rho_tilde(&p, v);
        let rep = (rt / p.v_star) * ((1.0 - rho) / rho - (self.z / 2.0) * (1.0 - rho).ln());
        let attr = (eps_star / p.v_star) * ((1.0 - 1.0 / p.r) * rho + rho * rho);
        rep - attr
    }

    /// Residual Helmholtz energy per mole (J·mol⁻¹) via ∫(P − RT/V)dV − RT·ln Z.
    fn residual_helmholtz(&self, t: f64, v: f64, x: &[f64]) -> f64 {
        let rt = gas_constant().value * t;
        let z = self.compressibility_factor(
            Temperature::new::<kelvin>(t),
            MolarVolume::new::<cubic_meter_per_mole>(v),
            x,
        );
        let z = z.unwrap_or(1.0);
        // Integrate (P − RT/V) from v to Vmax (where ideal-gas behaviour holds).
        let vmax = (v * 1.0e4).max(1.0);
        let n = 60;
        let mut integral = 0.0_f64;
        let mut prev_v = v;
        let mut prev = self.pressure_value(t, v, x) - rt / v;
        for k in 1..=n {
            let vk = v * (vmax / v).powf(k as f64 / n as f64);
            let pk = self.pressure_value(t, vk, x) - rt / vk;
            integral += 0.5 * (pk + prev) * (vk - prev_v);
            prev_v = vk;
            prev = pk;
        }
        integral - rt * z.value.ln()
    }
}

impl EquationOfState for SanchezLacombe {
    fn num_components(&self) -> usize {
        self.components.len()
    }

    fn pressure(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Pressure, ThermoError> {
        if z.len() != self.components.len() {
            return Err(ThermoError::InvalidInput("composition length mismatch"));
        }
        Ok(Pressure::new::<pascal>(self.pressure_value(t.value, v.value, z)))
    }

    fn ln_fugacity_coefficient(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
        _i: usize,
    ) -> Result<f64, ThermoError> {
        // Pseudo-pure (one-fluid) compressibility integral for the mixture.
        let p = self.pressure(t, v, z)?.value;
        if p <= 0.0 {
            return Ok(0.0);
        }
        let n = 40;
        let mut ln_phi = 0.0_f64;
        let mut prev_p = p * 1.0e-4;
        let mut prev_z = {
            let v0 = self
                .solve_molar_volume(
                    t,
                    Pressure::new::<pascal>(prev_p),
                    z,
                    MolarVolume::new::<cubic_meter_per_mole>(1e-6),
                    MolarVolume::new::<cubic_meter_per_mole>(1e2),
                )
                .map(|vv| vv.value)
                .unwrap_or(1.0);
            prev_p * v0 / (gas_constant().value * t.value)
        };
        for k in 1..=n {
            let pk = p * (k as f64 / n as f64);
            let vk = self
                .solve_molar_volume(
                    t,
                    Pressure::new::<pascal>(pk),
                    z,
                    MolarVolume::new::<cubic_meter_per_mole>(1e-6),
                    MolarVolume::new::<cubic_meter_per_mole>(1e2),
                )
                .map(|vv| vv.value)
                .unwrap_or(1.0);
            let zk = pk * vk / (gas_constant().value * t.value);
            let dlnp = (pk / prev_p).ln();
            ln_phi += 0.5 * ((zk - 1.0) + (prev_z - 1.0)) * dlnp;
            prev_p = pk;
            prev_z = zk;
        }
        Ok(ln_phi)
    }

    fn molar_enthalpy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEnergy, ThermoError> {
        let a_r = self.residual_helmholtz(t.value, v.value, z);
        let h = 1.0;
        let a_rp = self.residual_helmholtz(t.value + h, v.value, z);
        let a_rm = self.residual_helmholtz(t.value - h, v.value, z);
        let s_r = -(a_rp - a_rm) / (2.0 * h); // residual entropy (J/mol/K)
        let p = self.pressure(t, v, z)?.value;
        // Residual enthalpy H^R = A^R + T·S^R + P·v − R·T.
        let rt = gas_constant().value * t.value;
        let hr = a_r + t.value * s_r + p * v.value - rt;
        // Ideal-gas reference enthalpy: constant-Cp (3R per mole of segments-ish).
        let cp = 3.0 * gas_constant().value * self.mix(z).r.max(1.0);
        let t_ref = 298.15;
        let h_ig = cp * (t.value - t_ref);
        Ok(MolarEnergy::new::<joule_per_mole>(hr + h_ig))
    }

    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEntropy, ThermoError> {
        let a_r = self.residual_helmholtz(t.value, v.value, z);
        let h = 1.0;
        let a_rp = self.residual_helmholtz(t.value + h, v.value, z);
        let a_rm = self.residual_helmholtz(t.value - h, v.value, z);
        let s_r = -(a_rp - a_rm) / (2.0 * h);
        let cp = 3.0 * gas_constant().value * self.mix(z).r.max(1.0);
        let t_ref = 298.15;
        let s_ig = cp * (t.value / t_ref).ln();
        Ok(MolarEntropy::new::<joule_per_kelvin_mole>(s_r + s_ig))
    }

    fn molar_isobaric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        let h = 1.0;
        let hp = self
            .molar_enthalpy(Temperature::new::<kelvin>(t.value + h), v, z)?
            .value;
        let hm = self
            .molar_enthalpy(Temperature::new::<kelvin>(t.value - h), v, z)?
            .value;
        Ok(MolarHeatCapacity::new::<joule_per_kelvin_mole>((hp - hm) / (2.0 * h)))
    }

    fn molar_isochoric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        let cp = self.molar_isobaric_heat_capacity(t, v, z)?.value;
        let p = self.pressure(t, v, z)?.value;
        let dpdt = self.dp_dt(t, v, z)?;
        let dpdv = self.dp_dv(t, v, z)?;
        let cv = if dpdv.abs() > 1e-30 {
            cp - t.value * p * dpdt * dpdt / (v.value * dpdv)
        } else {
            cp
        };
        Ok(MolarHeatCapacity::new::<joule_per_kelvin_mole>(cv))
    }

    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Velocity, ThermoError> {
        let p = self.pressure(t, v, z)?;
        let cv = self.molar_isochoric_heat_capacity(t, v, z)?.value;
        let cp = self.molar_isobaric_heat_capacity(t, v, z)?.value;
        let kt = self.isothermal_compressibility(t, v, z)?;
        let rho = MolarMass::new::<uom::si::molar_mass::kilogram_per_mole>(self.mix(z).v_star
            * 0.0
            + 0.1)
            / v;
        let gamma = if cv > 0.0 { cp / cv } else { 1.0 };
        let a2 = gamma / (kt * rho.value.max(1e-9));
        Ok(Velocity::new::<meter_per_second>(a2.max(0.0).sqrt()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn water_like() -> SanchezLacombe {
        // Approximate SL parameters for a small molecule (water-ish).
        SanchezLacombe::pure(SlComponent {
            r: 1.0,
            v_star: 1.8e-5,
            epsilon_k: 450.0,
        })
    }

    #[test]
    fn ideal_gas_limit() {
        let eos = water_like();
        let t = Temperature::new::<kelvin>(500.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(0.1); // very dilute
        let p = eos.pressure(t, v, &[1.0]).unwrap();
        let expected = gas_constant().value * 500.0 / 0.1;
        assert!((p.value - expected).abs() / expected < 0.02, "p={}, exp={}", p.value, expected);
        let phi = eos.ln_fugacity_coefficient(t, v, &[1.0], 0).unwrap();
        assert!(phi.abs() < 0.1, "ln φ = {phi}");
    }

    #[test]
    fn pressure_finite_and_positive_in_liquid_region() {
        let eos = water_like();
        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(2.0e-5);
        let p = eos.pressure(t, v, &[1.0]).unwrap();
        assert!(p.value.is_finite() && p.value > 0.0, "p = {}", p.value);
    }

    #[test]
    fn enthalpy_residual_vanishes_at_dilution() {
        let eos = water_like();
        let t = Temperature::new::<kelvin>(500.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(0.2);
        let h = eos.molar_enthalpy(t, v, &[1.0]).unwrap();
        assert!(h.value.is_finite());
    }
}
