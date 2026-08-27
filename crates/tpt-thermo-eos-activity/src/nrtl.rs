//! NRTL — the Non-Random Two-Liquid model (Renon & Prausnitz, 1968).
//!
//! Implements [`tpt_thermo_core::mixing::ExcessGibbsModel`] for the symmetric
//! (or asymmetric) *τ* / *α* parameterisation:
//!
//! ```text
//! τ_ij  = a_ij + b_ij/T + c_ij·ln(T)         (see [`TdParam`](crate::parameters::TdParam))
//! G_ij  = exp(-α_ij · τ_ij)
//! g^E/(RT) = Σ_i x_i · ( Σ_j x_j G_ji τ_ji / Σ_k x_k G_ki )
//! ln γ_i = Σ_j [ x_j G_ji τ_ji / Σ_k x_k G_ki ]
//!          + Σ_j [ x_j G_ij / Σ_k x_k G_kj · ( τ_ij − Σ_m x_m G_mj τ_mj / Σ_k x_k G_kj ) ]
//! ```
//!
//! *α_ij* is a non-temperature-dependent randomness parameter (typically
//! 0.2–0.47); for `i == j` it is zero.

use crate::parameters::{self, TdMatrix, TdParam};
use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::mixing::ExcessGibbsModel;
use tpt_thermo_core::quantities::{Pressure, Temperature};

/// NRTL activity model.
#[derive(Debug, Clone)]
pub struct Nrtl {
    n: usize,
    /// τ_ij parameters (temperature dependent), full `n×n`.
    tau: TdMatrix,
    /// α_ij randomness parameters, full `n×n` (α_ii = 0).
    alpha: Vec<Vec<f64>>,
}

impl Nrtl {
    /// Build from a full `n×n` τ-matrix and α-matrix.
    ///
    /// `tau[i][j]` is the τ_ij parameter and `alpha[i][j]` the α_ij randomness
    /// parameter. Diagonal entries are unused and forced to zero on evaluation.
    pub fn new(tau: TdMatrix, alpha: Vec<Vec<f64>>) -> Result<Self, ThermoError> {
        let n = tau.len();
        if alpha.len() != n || alpha.iter().any(|row| row.len() != n) {
            return Err(ThermoError::InvalidInput("NRTL α matrix must match τ matrix size"));
        }
        Ok(Self { n, tau, alpha })
    }

    /// Convenience constructor for a binary system from the two τ pairs and α.
    pub fn binary(tau12: TdParam, tau21: TdParam, alpha: f64) -> Result<Self, ThermoError> {
        let mut tau_m = TdMatrix::zeros(2);
        tau_m.set(0, 1, tau12);
        tau_m.set(1, 0, tau21);
        let alpha_m = vec![vec![0.0, alpha], vec![alpha, 0.0]];
        Self::new(tau_m, alpha_m)
    }

    /// Number of components.
    pub fn num_components(&self) -> usize {
        self.n
    }

    /// The reduced excess Gibbs energy `g^E/(R T)` at `(T, x)`.
    pub fn reduced_excess_gibbs_at(&self, t: Temperature, x: &[f64]) -> Result<f64, ThermoError> {
        parameters::check_composition(x)?;
        if x.len() != self.n {
            return Err(ThermoError::InvalidInput("composition length mismatch"));
        }
        let mut g = vec![0.0f64; self.n * self.n]; // G_ij = exp(-α_ij τ_ij)
        for i in 0..self.n {
            for j in 0..self.n {
                let a = if i == j { 0.0 } else { self.alpha[i][j] };
                let tau = if i == j { 0.0 } else { self.tau.value_at(i, j, t) };
                g[i * self.n + j] = libm::exp(-a * tau);
            }
        }
        let mut total = 0.0;
        for i in 0..self.n {
            let mut denom = 0.0;
            for k in 0..self.n {
                denom += x[k] * g[k * self.n + i];
            }
            let mut num = 0.0;
            for j in 0..self.n {
                let tau = if j == i { 0.0 } else { self.tau.value_at(j, i, t) };
                num += x[j] * g[j * self.n + i] * tau;
            }
            total += x[i] * num / denom;
        }
        Ok(total)
    }

    /// Natural log of the activity coefficient of component `i` at `(T, x)`.
    pub fn ln_gamma_at(&self, t: Temperature, x: &[f64], i: usize) -> Result<f64, ThermoError> {
        parameters::check_composition(x)?;
        if x.len() != self.n || i >= self.n {
            return Err(ThermoError::InvalidInput("composition/ index mismatch"));
        }
        let mut g = vec![0.0f64; self.n * self.n];
        for a in 0..self.n {
            for b in 0..self.n {
                let alpha = if a == b { 0.0 } else { self.alpha[a][b] };
                let tau = if a == b { 0.0 } else { self.tau.value_at(a, b, t) };
                g[a * self.n + b] = libm::exp(-alpha * tau);
            }
        }
        // First term: Σ_j x_j G_ji τ_ji / Σ_k x_k G_ki
        let mut denom_i = 0.0;
        for k in 0..self.n {
            denom_i += x[k] * g[k * self.n + i];
        }
        let mut first = 0.0;
        for j in 0..self.n {
            let tau = if j == i { 0.0 } else { self.tau.value_at(j, i, t) };
            first += x[j] * g[j * self.n + i] * tau;
        }
        first /= denom_i;

        // Second term: Σ_j [ x_j G_ij / Σ_k x_k G_kj · (τ_ij − Σ_m x_m G_mj τ_mj / Σ_k x_k G_kj) ]
        let mut second = 0.0;
        for j in 0..self.n {
            let mut denom_j = 0.0;
            for k in 0..self.n {
                denom_j += x[k] * g[k * self.n + j];
            }
            let mut num = 0.0;
            for m in 0..self.n {
                let tau = if m == j { 0.0 } else { self.tau.value_at(m, j, t) };
                num += x[m] * g[m * self.n + j] * tau;
            }
            let tau_ij = if i == j { 0.0 } else { self.tau.value_at(i, j, t) };
            second += (x[j] * g[i * self.n + j] / denom_j) * (tau_ij - num / denom_j);
        }
        Ok(first + second)
    }

    /// Activity coefficient (not log) of component `i`.
    pub fn gamma_at(&self, t: Temperature, x: &[f64], i: usize) -> Result<f64, ThermoError> {
        Ok(libm::exp(self.ln_gamma_at(t, x, i)?))
    }
}

impl ExcessGibbsModel for Nrtl {
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

    fn binary() -> Nrtl {
        // Symmetric dummy parameters; α = 0.3.
        let tau12 = TdParam::new(0.5, 100.0, 0.0);
        let tau21 = TdParam::new(-0.2, 50.0, 0.0);
        Nrtl::binary(tau12, tau21, 0.3).unwrap()
    }

    #[test]
    fn pure_component_gives_zero_ln_gamma() {
        let m = binary();
        let t = Temperature::new::<kelvin>(333.15);
        assert!(m.ln_gamma_at(t, &[1.0, 0.0], 0).unwrap().abs() < 1e-12);
        assert!(m.ln_gamma_at(t, &[0.0, 1.0], 1).unwrap().abs() < 1e-12);
    }

    #[test]
    fn gibbs_duhem_consistency_binary() {
        let m = binary();
        let t = Temperature::new::<kelvin>(333.15);
        // Identity check: g^E/(R T) must equal Σ_i x_i ln γ_i at every
        // composition for a self-consistent model.
        for x1 in [0.1f64, 0.3, 0.5, 0.7, 0.9] {
            let x = [x1, 1.0 - x1];
            let ge = m.reduced_excess_gibbs_at(t, &x).unwrap();
            let gd = x[0] * m.ln_gamma_at(t, &x, 0).unwrap()
                + x[1] * m.ln_gamma_at(t, &x, 1).unwrap();
            assert!((ge - gd).abs() < 1e-9, "g^E/RT != Σ x lnγ at x1={x1}: ge={ge}, gd={gd}");
        }
    }

    #[test]
    fn infinite_dilution_limit_matches_one_term() {
        // At x_i → 1 the first term of ln γ_j should be exp(τ_ij)·G_ij form.
        let m = binary();
        let t = Temperature::new::<kelvin>(333.15);
        let tau12 = m.tau.value_at(0, 1, t);
        let a = m.alpha[0][1];
        let tau21 = m.tau.value_at(1, 0, t);
        let g21 = libm::exp(-a * tau21);
        // Correct infinite-dilution limit at x1 = 1: ln γ_2^∞ = τ_12 + G_21 τ_21.
        let expected = tau12 + g21 * tau21;
        let got = m.ln_gamma_at(t, &[1.0 - 1e-9, 1e-9], 1).unwrap();
        assert!((got - expected).abs() < 1e-6, "got {got}, expected {expected}");
    }
}
