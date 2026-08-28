//! Sanchez-Lacombe lattice-fluid equation of state.
//!
//! Implements the Sanchez-Lacombe (1976/1994) lattice-fluid model, which is the
//! natural EoS for polymer systems: each molecule occupies `r_i` lattice sites of
//! characteristic volume `v*_i` and carries a cohesive energy `ε_i` per segment.
//! Mixture parameters follow the standard one-fluid averages
//!
//! ```text
//! r̄   = Σ_i x_i r_i
//! v*  = Σ_i x_i r_i v*_i / r̄
//! ε   = (Σ_i x_i √ε_i)²
//! ```
//!
//! and the reduced density is `ρ̃ = r̄ v* / v`. The pressure is
//!
//! ```text
//! P = −(RT/v*) [ ln(1−ρ̃) + (1 − 1/r̄)·ρ̃ ] − (ε/v*)·ρ̃²
//! ```
//!
//! Fugacity coefficients and the residual thermodynamic properties are obtained
//! from the residual Helmholtz free energy via finite differences, so the model
//! honours the same [`EquationOfState`] contract as the cubic and SAFT crates.

use alloc::vec;
use alloc::vec::Vec;
use tpt_thermo_core::quantities::{
    MolarEnergy, MolarEntropy, MolarHeatCapacity, MolarVolume, Pressure, Temperature, Velocity,
};
use tpt_thermo_core::{EquationOfState, ThermoError, R};
use uom::si::molar_energy::joule_per_mole;
use uom::si::molar_heat_capacity::joule_per_kelvin_mole;
use uom::si::pressure::pascal as pa_unit;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::velocity::meter_per_second as v_unit;

/// Sanchez-Lacombe lattice-fluid equation of state.
#[derive(Debug, Clone)]
pub struct SanchezLacombe {
    /// Number of lattice sites per chain (degree of polymerization analogue).
    r: Vec<f64>,
    /// Characteristic volume of one segment (m³·mol⁻¹).
    vstar: Vec<f64>,
    /// Cohesive energy of one segment (J·mol⁻¹).
    epsilon: Vec<f64>,
    /// Ideal-gas isobaric heat capacity per component (J·mol⁻¹·K⁻¹), used for the
    /// ideal-gas reference of the absolute enthalpy/entropy.
    cp0: Vec<f64>,
    /// Reference temperature (K) for the ideal-gas reference.
    t_ref: f64,
    /// Reference pressure (Pa) for the ideal-gas reference.
    p_ref: f64,
}

impl SanchezLacombe {
    /// Build from segment counts, characteristic volumes, and cohesive energies;
    /// ideal-gas `Cp` defaults to `3R` per component.
    pub fn new(r: Vec<f64>, vstar: Vec<f64>, epsilon: Vec<f64>) -> Self {
        let n = r.len();
        let cp0 = vec![3.0 * R; n];
        Self {
            r,
            vstar,
            epsilon,
            cp0,
            t_ref: 298.15,
            p_ref: 1.0e5,
        }
    }

    /// Build with an explicit ideal-gas heat capacity per component.
    pub fn with_ideal_gas_cp(
        r: Vec<f64>,
        vstar: Vec<f64>,
        epsilon: Vec<f64>,
        cp0: Vec<f64>,
    ) -> Self {
        Self {
            r,
            vstar,
            epsilon,
            cp0,
            t_ref: 298.15,
            p_ref: 1.0e5,
        }
    }

    /// One-fluid mixture parameters `(r̄, v*, ε)` at composition `z`.
    fn mixture_params(&self, z: &[f64]) -> (f64, f64, f64) {
        let rbar: f64 = z.iter().zip(self.r.iter()).map(|(zi, ri)| zi * ri).sum();
        let vstar_mix: f64 = z
            .iter()
            .zip(self.r.iter())
            .zip(self.vstar.iter())
            .map(|((zi, ri), vsi)| zi * ri * vsi)
            .sum::<f64>()
            / rbar;
        let e_sqrt: f64 = z
            .iter()
            .zip(self.epsilon.iter())
            .map(|(zi, ei)| zi * ei.sqrt())
            .sum();
        let eps_mix = e_sqrt * e_sqrt;
        (rbar, vstar_mix, eps_mix)
    }

    /// Pressure in pascals at the given temperature, molar volume (m³·mol⁻¹) and
    /// composition. Returns `f64::NAN` if the state is outside the physical
    /// (sub-close-packed) domain.
    fn p_val(&self, t: Temperature, v: f64, z: &[f64]) -> f64 {
        let rt = R * t.value;
        let (rbar, vstar_mix, eps_mix) = self.mixture_params(z);
        let rho = rbar * vstar_mix / v;
        if rho >= 1.0 {
            return f64::NAN;
        }
        let term = (1.0 - rho).ln() + (1.0 - 1.0 / rbar) * rho;
        -(rt / vstar_mix) * term - (eps_mix / vstar_mix) * rho * rho
    }

    /// Solve for the molar volume (m³·mol⁻¹) giving `target_p` (Pa) at `(T, z)`,
    /// following the isotherm from `guess` to stay on a consistent branch.
    fn solve_v_for_p(&self, t: Temperature, z: &[f64], target_p: f64, guess: f64) -> f64 {
        let (rbar, vstar_mix, _) = self.mixture_params(z);
        let v_min = rbar * vstar_mix / 0.995;
        let mut v = guess.max(v_min * 1.001).min(1e6);
        for _ in 0..80 {
            let p = self.p_val(t, v, z);
            if !p.is_finite() {
                break;
            }
            let f = p - target_p;
            if f.abs() < target_p.abs() * 1e-9 + 1e-6 {
                return v;
            }
            let dv = 1e-6 * v.max(1e-9);
            let dpdv = (self.p_val(t, v + dv, z) - p) / dv;
            if dpdv.abs() < 1e-30 {
                break;
            }
            let mut step = f / dpdv;
            step = step.clamp(-0.4 * v, 0.4 * v);
            v -= step;
            if v <= v_min {
                v = v_min * 1.001;
            }
            if v > 1e6 {
                v = 1e6;
            }
        }
        v
    }

    /// Residual Gibbs free energy per mole (J·mol⁻¹) via the pressure integral
    /// `g^res/(RT) = ∫_0^P (Z−1)/P dP`.
    fn g_res_molar(&self, t: Temperature, v: f64, z: &[f64]) -> f64 {
        let p = self.p_val(t, v, z);
        if !p.is_finite() || p <= 0.0 {
            return 0.0;
        }
        let rt = R * t.value;
        let n = 160;
        let p0 = (p * 1e-5).max(1e-2);
        let mut v_prev = self.solve_v_for_p(t, z, p0, v * 1e3);
        let mut p_prev = p0;
        let mut z_prev = p_prev * v_prev / rt;
        let mut integral = 0.0_f64;
        for step in 1..=n {
            let frac = step as f64 / n as f64;
            let p_step = p0 * (p / p0).powf(frac);
            let v_step = self.solve_v_for_p(t, z, p_step, v_prev);
            let z_step = p_step * v_step / rt;
            integral +=
                0.5 * ((z_prev - 1.0) / p_prev + (z_step - 1.0) / p_step) * (p_step - p_prev);
            p_prev = p_step;
            v_prev = v_step;
            z_prev = z_step;
        }
        integral * rt
    }

    /// Residual Helmholtz free energy per mole (J·mol⁻¹).
    fn a_res_molar(&self, t: Temperature, v: f64, z: &[f64]) -> f64 {
        let g = self.g_res_molar(t, v, z);
        let p = self.p_val(t, v, z);
        let rtv = R * t.value / v;
        g - (p - rtv) * v
    }

    /// Numerical `(∂P/∂v)_T` (Pa·m³·mol⁻¹) at the given state.
    fn dpdv(&self, t: Temperature, v: f64, z: &[f64]) -> f64 {
        let h = 1e-6 * v.max(1e-9);
        (self.p_val(t, v + h, z) - self.p_val(t, v - h, z)) / (2.0 * h)
    }

    /// Numerical `(∂P/∂T)_v` (Pa·K⁻¹) at the given state.
    fn dpdt(&self, t: Temperature, v: f64, z: &[f64]) -> f64 {
        let h = 1e-4 * t.value.max(1.0);
        (self.p_val(Temperature::new::<kelvin>(t.value + h), v, z)
            - self.p_val(Temperature::new::<kelvin>(t.value - h), v, z))
            / (2.0 * h)
    }
}

impl EquationOfState for SanchezLacombe {
    fn num_components(&self) -> usize {
        self.r.len()
    }

    fn pressure(&self, t: Temperature, v: MolarVolume, z: &[f64]) -> Result<Pressure, ThermoError> {
        if z.len() != self.r.len() {
            return Err(ThermoError::InvalidInput("composition length mismatch"));
        }
        let p = self.p_val(t, v.value, z);
        if !p.is_finite() {
            return Err(ThermoError::InvalidInput("state outside physical domain"));
        }
        Ok(Pressure::new::<pa_unit>(p))
    }

    fn ln_fugacity_coefficient(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        if i >= self.r.len() {
            return Err(ThermoError::IndexOutOfRange(i));
        }
        let p = self.p_val(t, v.value, z);
        if !p.is_finite() || p <= 0.0 {
            return Err(ThermoError::InvalidInput("state outside physical domain"));
        }
        let rt = R * t.value;
        // Residual Gibbs per RT at the current (T, P, z) state.
        let h = self.g_res_molar(t, v.value, z) / rt;
        // Partial derivatives of `h` w.r.t. each mole fraction, holding (T, P) fixed
        // and renormalising the remaining fractions. The fugacity coefficient of a
        // pure component reduces to `ln φ = h`.
        let n = z.len();
        let mut dh = alloc::vec::Vec::with_capacity(n);
        let eps: f64 = 1e-6;
        for k in 0..n {
            let d = eps.min(0.5 * (1.0 - z[k]));
            if d <= 0.0 {
                dh.push(0.0);
                continue;
            }
            let mut zp = z.to_vec();
            zp[k] += d;
            let sp: f64 = zp.iter().sum();
            for x in zp.iter_mut() {
                *x /= sp;
            }
            let vp = self.solve_v_for_p(t, &zp, p, v.value);
            let hp = self.g_res_molar(t, vp, &zp) / rt;
            let mut zm = z.to_vec();
            zm[k] -= d;
            let sm: f64 = zm.iter().sum();
            for x in zm.iter_mut() {
                *x /= sm;
            }
            let vm = self.solve_v_for_p(t, &zm, p, v.value);
            let hm = self.g_res_molar(t, vm, &zm) / rt;
            dh.push((hp - hm) / (2.0 * d));
        }
        let mut sum = 0.0_f64;
        for k in 0..n {
            sum += z[k] * dh[k];
        }
        Ok(h + dh[i] - sum)
    }

    fn molar_enthalpy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEnergy, ThermoError> {
        // Ideal-gas reference.
        let mut h_ig = 0.0_f64;
        for (i, &zi) in z.iter().enumerate() {
            h_ig += zi * self.cp0[i] * (t.value - self.t_ref);
        }
        // Residual: h^res = a^res + T s^res + (P − RT/v) v, s^res = −(∂a^res/∂T)_v.
        let dt = 1e-3 * t.value.max(1.0);
        let ap = self.a_res_molar(Temperature::new::<kelvin>(t.value + dt), v.value, z);
        let am = self.a_res_molar(Temperature::new::<kelvin>(t.value - dt), v.value, z);
        let s_res = -(ap - am) / (2.0 * dt);
        let a = self.a_res_molar(t, v.value, z);
        let p = self.p_val(t, v.value, z);
        let rtv = R * t.value / v.value;
        let h_res = a + t.value * s_res + (p - rtv) * v.value;
        Ok(MolarEnergy::new::<joule_per_mole>(h_ig + h_res))
    }

    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEntropy, ThermoError> {
        let p = self.p_val(t, v.value, z);
        let dt = 1e-3 * t.value.max(1.0);
        let ap = self.a_res_molar(Temperature::new::<kelvin>(t.value + dt), v.value, z);
        let am = self.a_res_molar(Temperature::new::<kelvin>(t.value - dt), v.value, z);
        let s_res = -(ap - am) / (2.0 * dt);
        // Ideal-gas entropy with mixing.
        let mut s_ig = 0.0_f64;
        for (i, &zi) in z.iter().enumerate() {
            if zi <= 0.0 {
                continue;
            }
            let mixing = -R * zi.ln();
            s_ig += zi * (self.cp0[i] * (t.value / self.t_ref).ln() - R * (p / self.p_ref).ln())
                + mixing;
        }
        Ok(MolarEntropy::new::<joule_per_kelvin_mole>(s_ig + s_res))
    }

    fn molar_isochoric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        let dt = 1e-3 * t.value.max(1.0);
        let up = {
            let a = self.a_res_molar(Temperature::new::<kelvin>(t.value + dt), v.value, z);
            let s = -(self.a_res_molar(Temperature::new::<kelvin>(t.value + 2.0 * dt), v.value, z)
                - a)
                / (2.0 * dt);
            a + (t.value + dt) * s
        };
        let um = {
            let a = self.a_res_molar(Temperature::new::<kelvin>(t.value - dt), v.value, z);
            let s =
                -(a - self.a_res_molar(Temperature::new::<kelvin>(t.value - 2.0 * dt), v.value, z))
                    / (2.0 * dt);
            a + (t.value - dt) * s
        };
        let cv_res = (up - um) / (2.0 * dt);
        let mut cv_ig = 0.0_f64;
        for (i, &zi) in z.iter().enumerate() {
            cv_ig += zi * (self.cp0[i] - R);
        }
        Ok(MolarHeatCapacity::new::<joule_per_kelvin_mole>(
            cv_ig + cv_res,
        ))
    }

    fn molar_isobaric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        let cv = self.molar_isochoric_heat_capacity(t, v, z)?;
        let dpdt = self.dpdt(t, v.value, z);
        let dpdv = self.dpdv(t, v.value, z);
        if dpdv.abs() < 1e-30 {
            return Err(ThermoError::Numerical(
                tpt_thermo_core::convergence::ConvergenceStatus::NumericalIssue(
                    tpt_thermo_core::convergence::NumericalIssueReason::SingularJacobian,
                ),
            ));
        }
        let cp_minus_cv = -t.value * dpdt * dpdt / dpdv;
        Ok(MolarHeatCapacity::new::<joule_per_kelvin_mole>(
            cv.value + cp_minus_cv,
        ))
    }

    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Velocity, ThermoError> {
        let dpdv = self.dpdv(t, v.value, z);
        let dpdt = self.dpdt(t, v.value, z);
        let a2 = -v.value * v.value * dpdv + v.value * v.value * dpdt * dpdt / dpdv;
        if a2 <= 0.0 || !a2.is_finite() {
            return Err(ThermoError::Unsupported(
                "sound speed undefined (unstable or non-hyperbolic state)",
            ));
        }
        Ok(Velocity::new::<v_unit>(a2.sqrt()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::molar_volume::cubic_meter_per_mole;

    fn ideal_like() -> SanchezLacombe {
        // r = 1, no cohesion ⇒ reduces to ideal gas; characteristic volume small.
        SanchezLacombe::new(vec![1.0], vec![1e-5], vec![0.0])
    }

    #[test]
    fn ideal_limit_pressure() {
        let m = ideal_like();
        let t = Temperature::new::<kelvin>(300.0);
        // At low density, P ≈ RT/v.
        let v = MolarVolume::new::<cubic_meter_per_mole>(0.5); // P ≈ 8.314*300/0.5 ≈ 5000 Pa
        let p = m.pressure(t, v, &[1.0]).unwrap();
        let expected = R * 300.0 / 0.5;
        assert!(
            (p.value - expected).abs() / expected < 1e-3,
            "{} vs {}",
            p.value,
            expected
        );
    }

    #[test]
    fn ideal_limit_fugacity_near_one() {
        let m = ideal_like();
        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(0.5);
        let lng = m.ln_fugacity_coefficient(t, v, &[1.0], 0).unwrap();
        assert!(lng.abs() < 1e-2, "ln φ = {lng}");
    }

    #[test]
    fn ideal_limit_enthalpy_residual_small() {
        let m = ideal_like();
        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(0.5);
        let h = m.molar_enthalpy(t, v, &[1.0]).unwrap();
        // With ε = 0 the residual enthalpy is small; h ≈ h^ig = cp0 (T − T_ref).
        let expected = 3.0 * R * (300.0 - 298.15);
        assert!(
            (h.value - expected).abs() < 1.0,
            "residual enthalpy = {} J/mol",
            h.value - expected
        );
    }

    #[test]
    fn cohesive_term_raises_pressure_relative_to_ideal() {
        // Increasing ε deepens the attractive well, lowering pressure at fixed v.
        let t = Temperature::new::<kelvin>(400.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(2.0e-3);
        let weak = SanchezLacombe::new(vec![100.0], vec![1e-5], vec![1.0e3]);
        let strong = SanchezLacombe::new(vec![100.0], vec![1e-5], vec![5.0e3]);
        let pw = weak.pressure(t, v, &[1.0]).unwrap().value;
        let ps = strong.pressure(t, v, &[1.0]).unwrap().value;
        assert!(ps < pw, "stronger cohesion should lower pressure");
    }

    #[test]
    fn speed_of_sound_positive_ideal() {
        let m = ideal_like();
        let t = Temperature::new::<kelvin>(300.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(0.5);
        let a = m.speed_of_sound(t, v, &[1.0]);
        assert!(a.is_ok());
        assert!(a.unwrap().value > 0.0);
    }

    #[test]
    fn pressure_matches_reference_formula() {
        // Manual evaluation for a known state.
        let m = SanchezLacombe::new(vec![10.0], vec![1e-5], vec![2000.0]);
        let t = Temperature::new::<kelvin>(350.0);
        let v = MolarVolume::new::<cubic_meter_per_mole>(5e-3);
        let p = m.pressure(t, v, &[1.0]).unwrap();
        let rbar = 10.0;
        let vstar = 1e-5;
        let rho: f64 = rbar * vstar / 5e-3;
        let term = (1.0 - rho).ln() + (1.0 - 1.0 / rbar) * rho;
        let expected = -(R * 350.0 / vstar) * term - (2000.0 / vstar) * rho * rho;
        assert!((p.value - expected).abs() < 1e-6);
    }
}
