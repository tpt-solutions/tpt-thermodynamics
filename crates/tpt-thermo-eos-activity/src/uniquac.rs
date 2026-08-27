//! UNIQUAC — the Universal Quasi-Chemical model (Abrams & Prausnitz, 1975).
//!
//! Requires per-component structural parameters `r_i` (volume) and `q_i`
//! (surface area), plus an asymmetric energy-parameter matrix `a_ij` (K) in the
//! reduced form `τ_ij = exp(-(a_ij − a_ii)/T)`:
//!
//! ```text
//! φ_i = x_i r_i / Σ_j x_j r_j        (segment/volume fraction)
//! θ_i = x_i q_i / Σ_j x_j q_j        (surface fraction)
//! l_i = z/2·(r_i − q_i) − (r_i − 1)  (z = 10)
//! g^E/(RT) = Σ_i x_i[ ln(φ_i/x_i) + z/2·q_i·ln(θ_i/φ_i) + l_i − φ_i/x_i·Σ_j x_j l_j ]
//!            + Σ_i x_i q_i[ 1 − ln(Σ_k θ_k τ_ki) − Σ_j θ_j τ_ij / (Σ_k θ_k τ_kj) ]
//! ln γ_i = ln(φ_i/x_i) + z/2·q_i·ln(θ_i/φ_i) + l_i − φ_i/x_i·Σ_j x_j l_j
//!          + q_i[ 1 − ln(Σ_k θ_k τ_ki) − Σ_j θ_j τ_ij / (Σ_k θ_k τ_kj) ]
//! ```

use crate::parameters::{self, InteractionMatrix};
use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::mixing::ExcessGibbsModel;
use tpt_thermo_core::quantities::{Pressure, Temperature};

/// Lattice coordination number used by UNIQUAC.
const Z: f64 = 10.0;

/// Per-component structural parameters `(r_i, q_i)`.
#[derive(Debug, Clone)]
pub struct StructuralParams(pub Vec<(f64, f64)>);

impl StructuralParams {
    /// Build from a list of `(r_i, q_i)` pairs.
    pub fn new(params: Vec<(f64, f64)>) -> Self {
        Self(params)
    }

    /// `r_i` (volume parameter) for component `i`.
    pub fn r(&self, i: usize) -> f64 {
        self.0[i].0
    }

    /// `q_i` (surface-area parameter) for component `i`.
    pub fn q(&self, i: usize) -> f64 {
        self.0[i].1
    }

    /// Number of components.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// UNIQUAC activity model.
#[derive(Debug, Clone)]
pub struct Uniquac {
    n: usize,
    struc: StructuralParams,
    a: InteractionMatrix,
}

impl Uniquac {
    /// Build from structural parameters and an energy-parameter matrix (K).
    pub fn new(struc: StructuralParams, a: InteractionMatrix) -> Result<Self, ThermoError> {
        let n = struc.len();
        if a.len() != n {
            return Err(ThermoError::InvalidInput("UNIQUAC matrix size mismatch"));
        }
        Ok(Self { n, struc, a })
    }

    /// Convenience constructor for a binary system from `(r1,q1)`, `(r2,q2)` and
    /// the two energy differences `Δa_12 = a_12 − a_11`, `Δa_21 = a_21 − a_22`.
    pub fn binary(
        r1: f64,
        q1: f64,
        r2: f64,
        q2: f64,
        da12: f64,
        da21: f64,
    ) -> Result<Self, ThermoError> {
        let mut m = InteractionMatrix::zeros(2);
        m.set(0, 1, da12);
        m.set(1, 0, da21);
        Self::new(StructuralParams::new(vec![(r1, q1), (r2, q2)]), m)
    }

    /// Number of components.
    pub fn num_components(&self) -> usize {
        self.n
    }

    fn tau(&self, i: usize, j: usize, tk: f64) -> f64 {
        if i == j {
            1.0
        } else {
            let da = self.a.get(i, j) - self.a.get(i, i);
            libm::exp(-da / tk)
        }
    }

    fn fractions(&self, x: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let mut phi = vec![0.0f64; self.n];
        let mut theta = vec![0.0f64; self.n];
        let mut rs = 0.0;
        let mut qs = 0.0;
        for i in 0..self.n {
            rs += x[i] * self.struc.r(i);
            qs += x[i] * self.struc.q(i);
        }
        for i in 0..self.n {
            phi[i] = x[i] * self.struc.r(i) / rs;
            theta[i] = x[i] * self.struc.q(i) / qs;
        }
        (phi, theta)
    }

    /// Reduced excess Gibbs energy `g^E/(R T)`.
    pub fn reduced_excess_gibbs_at(&self, t: Temperature, x: &[f64]) -> Result<f64, ThermoError> {
        parameters::check_composition(x)?;
        if x.len() != self.n {
            return Err(ThermoError::InvalidInput("composition length mismatch"));
        }
        let tk = parameters::tk(t);
        let (phi, theta) = self.fractions(x);
        let mut l = vec![0.0f64; self.n];
        for i in 0..self.n {
            l[i] = 0.5 * Z * (self.struc.r(i) - self.struc.q(i)) - (self.struc.r(i) - 1.0);
        }
        let mut sum_lx = 0.0;
        for j in 0..self.n {
            sum_lx += x[j] * l[j];
        }
        let mut total = 0.0;
        for i in 0..self.n {
            let mut sum_tau = 0.0;
            for k in 0..self.n {
                sum_tau += theta[k] * self.tau(k, i, tk);
            }
            let mut resid_inner = 0.0;
            for j in 0..self.n {
                let mut denom = 0.0;
                for k in 0..self.n {
                    denom += theta[k] * self.tau(k, j, tk);
                }
                resid_inner += theta[j] * self.tau(i, j, tk) / denom;
            }
            let comb = libm::log(phi[i] / x[i])
                + 0.5 * Z * self.struc.q(i) * libm::log(theta[i] / phi[i])
                + l[i]
                - (phi[i] / x[i]) * sum_lx;
            total += x[i] * (comb + self.struc.q(i) * (1.0 - libm::log(sum_tau) - resid_inner));
        }
        Ok(total)
    }

    /// Natural log of the activity coefficient of component `i`.
    pub fn ln_gamma_at(&self, t: Temperature, x: &[f64], i: usize) -> Result<f64, ThermoError> {
        parameters::check_composition(x)?;
        if x.len() != self.n || i >= self.n {
            return Err(ThermoError::InvalidInput("composition/ index mismatch"));
        }
        let tk = parameters::tk(t);
        let (phi, theta) = self.fractions(x);
        let mut l = vec![0.0f64; self.n];
        for k in 0..self.n {
            l[k] = 0.5 * Z * (self.struc.r(k) - self.struc.q(k)) - (self.struc.r(k) - 1.0);
        }
        let mut sum_lx = 0.0;
        for j in 0..self.n {
            sum_lx += x[j] * l[j];
        }
        let mut sum_tau_i = 0.0;
        for k in 0..self.n {
            sum_tau_i += theta[k] * self.tau(k, i, tk);
        }
        let mut resid = 0.0;
        for j in 0..self.n {
            let mut denom = 0.0;
            for k in 0..self.n {
                denom += theta[k] * self.tau(k, j, tk);
            }
            resid += theta[j] * self.tau(i, j, tk) / denom;
        }
        let comb = libm::log(phi[i] / x[i])
            + 0.5 * Z * self.struc.q(i) * libm::log(theta[i] / phi[i])
            + l[i]
            - (phi[i] / x[i]) * sum_lx;
        Ok(comb + self.struc.q(i) * (1.0 - libm::log(sum_tau_i) - resid))
    }

    /// Activity coefficient (not log) of component `i`.
    pub fn gamma_at(&self, t: Temperature, x: &[f64], i: usize) -> Result<f64, ThermoError> {
        Ok(libm::exp(self.ln_gamma_at(t, x, i)?))
    }
}

impl ExcessGibbsModel for Uniquac {
    fn num_components(&self) -> usize {
        self.n
    }

    fn reduced_excess_gibbs(
        &self,
        t: Temperature,
        _p: Pressure,
        x: &[f64],
    ) -> Result<f64, ThermoError> {
        self.reduced_excess_gibbs_at(t, x)
    }

    fn ln_gamma(
        &self,
        t: Temperature,
        _p: Pressure,
        x: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        self.ln_gamma_at(t, x, i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::thermodynamic_temperature::kelvin;

    // Water (1) / ethanol (2) — literature-ish r,q and energy differences (K).
    fn water_ethanol() -> Uniquac {
        Uniquac::binary(0.92, 1.40, 2.11, 1.97, -229.1, 119.7).unwrap()
    }

    #[test]
    fn pure_component_gives_zero_ln_gamma() {
        let m = water_ethanol();
        let t = Temperature::new::<kelvin>(298.15);
        assert!(m.ln_gamma_at(t, &[1.0, 0.0], 0).unwrap().abs() < 1e-12);
        assert!(m.ln_gamma_at(t, &[0.0, 1.0], 1).unwrap().abs() < 1e-12);
    }

    #[test]
    fn gibbs_duhem_consistency_binary() {
        let m = water_ethanol();
        let t = Temperature::new::<kelvin>(298.15);
        for x1 in [0.1f64, 0.3, 0.5, 0.7, 0.9] {
            let x = [x1, 1.0 - x1];
            let ge = m.reduced_excess_gibbs_at(t, &x).unwrap();
            let gd =
                x[0] * m.ln_gamma_at(t, &x, 0).unwrap() + x[1] * m.ln_gamma_at(t, &x, 1).unwrap();
            assert!(
                (ge - gd).abs() < 1e-9,
                "g^E/RT != Σ x lnγ at x1={x1}: ge={ge}, gd={gd}"
            );
        }
    }

    #[test]
    fn activity_coefficient_is_positive() {
        let m = water_ethanol();
        let t = Temperature::new::<kelvin>(298.15);
        let g = m.gamma_at(t, &[0.4, 0.6], 0).unwrap();
        assert!(g.is_finite() && g > 0.0);
    }
}
