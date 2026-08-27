//! Wilson's equation (1964) — a local-composition model built from component
//! molar volumes and binary energy parameters.
//!
//! ```text
//! Λ_ij = (V_j / V_i) · exp(-(a_ij − a_ii)/T)        (a in K; Λ_ii = 1)
//! g^E/(RT) = − Σ_i x_i · ln( Σ_j x_j Λ_ij )
//! ln γ_i = 1 − ln( Σ_j x_j Λ_ij ) − Σ_k x_k Λ_ki / ( Σ_j x_j Λ_kj )
//! ```
//!
//! Note Wilson's `g^E` has no liquid-liquid miscibility gap term, so it cannot
//! describe partially-miscible (LLE) systems — appropriate for the VLE systems
//! this crate targets.

use crate::parameters::{self, InteractionMatrix, Volumes};
use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::mixing::ExcessGibbsModel;
use tpt_thermo_core::quantities::{Pressure, Temperature};

/// Wilson activity model.
#[derive(Debug, Clone)]
pub struct Wilson {
    n: usize,
    /// Energy parameters `a_ij` in kelvin (asymmetric; `a_ii` arbitrary, used in
    /// the `a_ij − a_ii` difference).
    a: InteractionMatrix,
    /// Component molar volumes (m³·mol⁻¹).
    volumes: Volumes,
}

impl Wilson {
    /// Build from an energy-parameter matrix (K) and per-component molar volumes.
    pub fn new(a: InteractionMatrix, volumes: Volumes) -> Result<Self, ThermoError> {
        let n = a.len();
        if volumes.len() != n {
            return Err(ThermoError::InvalidInput("Wilson volumes length mismatch"));
        }
        Ok(Self { n, a, volumes })
    }

    /// Convenience constructor for a binary system from the two energy
    /// differences `Δa_12 = a_12 − a_11`, `Δa_21 = a_21 − a_22` and the two
    /// molar volumes `V_1, V_2`.
    pub fn binary(da12: f64, da21: f64, v1: f64, v2: f64) -> Result<Self, ThermoError> {
        let mut m = InteractionMatrix::zeros(2);
        m.set(0, 1, da12);
        m.set(1, 0, da21);
        // Diagonal contributes only via the difference, so leave at 0.
        Self::new(m, Volumes::new(vec![v1, v2]))
    }

    /// Number of components.
    pub fn num_components(&self) -> usize {
        self.n
    }

    /// Build the `Λ_ij` matrix at temperature `t`.
    fn lambda_matrix(&self, t: Temperature) -> Vec<Vec<f64>> {
        let tk = parameters::tk(t);
        let mut lam = vec![vec![0.0f64; self.n]; self.n];
        for i in 0..self.n {
            for j in 0..self.n {
                if i == j {
                    lam[i][j] = 1.0;
                } else {
                    let ratio = self.volumes.get(j) / self.volumes.get(i);
                    let da = self.a.get(i, j) - self.a.get(i, i);
                    lam[i][j] = ratio * libm::exp(-da / tk);
                }
            }
        }
        lam
    }

    /// Reduced excess Gibbs energy `g^E/(R T)`.
    pub fn reduced_excess_gibbs_at(&self, t: Temperature, x: &[f64]) -> Result<f64, ThermoError> {
        parameters::check_composition(x)?;
        if x.len() != self.n {
            return Err(ThermoError::InvalidInput("composition length mismatch"));
        }
        let lam = self.lambda_matrix(t);
        let mut total = 0.0;
        for i in 0..self.n {
            let mut s = 0.0;
            for j in 0..self.n {
                s += x[j] * lam[i][j];
            }
            total -= x[i] * libm::log(s);
        }
        Ok(total)
    }

    /// Natural log of the activity coefficient of component `i`.
    pub fn ln_gamma_at(&self, t: Temperature, x: &[f64], i: usize) -> Result<f64, ThermoError> {
        parameters::check_composition(x)?;
        if x.len() != self.n || i >= self.n {
            return Err(ThermoError::InvalidInput("composition/ index mismatch"));
        }
        let lam = self.lambda_matrix(t);
        let mut s_i = 0.0;
        for j in 0..self.n {
            s_i += x[j] * lam[i][j];
        }
        let mut second = 0.0;
        for k in 0..self.n {
            let mut sk = 0.0;
            for j in 0..self.n {
                sk += x[j] * lam[k][j];
            }
            second += x[k] * lam[k][i] / sk;
        }
        Ok(1.0 - libm::log(s_i) - second)
    }

    /// Activity coefficient of component `i` (not log).
    pub fn gamma_at(&self, t: Temperature, x: &[f64], i: usize) -> Result<f64, ThermoError> {
        Ok(libm::exp(self.ln_gamma_at(t, x, i)?))
    }
}

impl ExcessGibbsModel for Wilson {
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

    fn binary() -> Wilson {
        Wilson::binary(-200.0, 400.0, 1.0e-4, 1.5e-4).unwrap()
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
    fn lambda_diagonal_is_one() {
        let m = binary();
        let t = Temperature::new::<kelvin>(333.15);
        let lam = m.lambda_matrix(t);
        assert!((lam[0][0] - 1.0).abs() < 1e-12);
        assert!((lam[1][1] - 1.0).abs() < 1e-12);
    }
}
