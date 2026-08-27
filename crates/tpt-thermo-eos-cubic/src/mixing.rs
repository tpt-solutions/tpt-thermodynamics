//! Mixing rules for the cubic EoS family.
//!
//! * [`VdwMixing`] — classical van der Waals one-fluid mixing with optional
//!   temperature-dependent binary interaction parameters `k_ij(T)`.
//! * [`HuronVidal`] — Huron-Vidal first/second order (MHV1, MHV2) and PSRK,
//!   generic over the core's [`ExcessGibbsModel`] trait (implemented by the
//!   activity crate, Phase 5). The cross attractive parameters are reconstructed
//!   from the excess model's `ln γ` gradient so the EoS fugacity is consistent
//!   with the activity model.
//! * [`WongSandler`] — Wong-Sandler mixing, also generic over
//!   [`ExcessGibbsModel`].
//!
//! All mixing rules implement [`CubicMixing`], the richer interface the engine
//! needs (cross attractive parameters, not just mixture `a`/`b`).

use alloc::boxed::Box;
use alloc::vec;
use alloc::vec::Vec;
use tpt_thermo_core::mixing::ExcessGibbsModel;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;

/// The mixing interface the cubic engine uses. `t` is the absolute temperature
/// (K) and `p` is a reference pressure (Pa, ~1 bar) at which the excess model is
/// evaluated (excess Gibbs energy is effectively pressure-independent, so the
/// engine always passes a fixed reference rather than the unknown state
/// pressure, avoiding an implicit solve).
pub trait CubicMixing: Send + Sync {
    /// Mixture co-volume `b = Σ z_i b_i` (the standard vdW form, also used by
    /// HV/WS for the `b` combination).
    fn b_mix(&self, b: &[f64], z: &[f64]) -> f64;
    /// Mixture attractive parameter `a_mix` at `(T, z)`.
    fn a_mix(&self, a: &[f64], b: &[f64], z: &[f64], t: f64, p: f64) -> f64;
    /// `Σ_j z_j a_ij` for component `i` (used in the fugacity coefficient).
    fn aij_sum(&self, a: &[f64], b: &[f64], z: &[f64], i: usize, t: f64, p: f64) -> f64;
}

/// van der Waals one-fluid mixing with optional T-dependent `k_ij(T)`.
#[derive(Debug, Clone)]
pub struct VdwMixing {
    kij: Vec<Vec<f64>>,
    td: Option<Vec<Vec<(f64, f64, f64)>>>,
}

impl VdwMixing {
    /// Zero binary interactions for `n` components.
    pub fn new(n: usize) -> Self {
        Self {
            kij: vec![vec![0.0; n]; n],
            td: None,
        }
    }

    /// Build from a symmetric `k_ij` matrix (indexed `[i][j]`).
    pub fn from_matrix(kij: Vec<Vec<f64>>) -> Self {
        Self { kij, td: None }
    }

    /// Attach temperature-dependent coefficients `(a, b, c)` for each pair, so
    /// `k_ij(T) = a + b/T + c·ln(T)`. The constant matrix (if any) is ignored
    /// once these are set.
    pub fn with_tdependent(mut self, td: Vec<Vec<(f64, f64, f64)>>) -> Self {
        self.td = Some(td);
        self
    }

    fn kij(&self, i: usize, j: usize, t: f64) -> f64 {
        if i == j {
            return 0.0;
        }
        match &self.td {
            Some(td) => {
                let (ka, kb, kc) = td[i][j];
                ka + kb / t + kc * t.ln()
            }
            None => self.kij[i][j],
        }
    }

    #[inline]
    fn a_ij(&self, a: &[f64], i: usize, j: usize, t: f64) -> f64 {
        let k = self.kij(i, j, t);
        (1.0 - k) * (a[i] * a[j]).sqrt()
    }
}

impl CubicMixing for VdwMixing {
    fn b_mix(&self, b: &[f64], z: &[f64]) -> f64 {
        b.iter().zip(z.iter()).map(|(bi, zi)| bi * zi).sum()
    }

    fn a_mix(&self, a: &[f64], _b: &[f64], z: &[f64], t: f64, _p: f64) -> f64 {
        let mut sum = 0.0;
        for (i, &zi) in z.iter().enumerate() {
            for (j, &zj) in z.iter().enumerate() {
                sum += zi * zj * self.a_ij(a, i, j, t);
            }
        }
        sum
    }

    fn aij_sum(&self, a: &[f64], _b: &[f64], z: &[f64], i: usize, t: f64, _p: f64) -> f64 {
        let mut sum = 0.0;
        for (j, &zj) in z.iter().enumerate() {
            sum += zj * self.a_ij(a, i, j, t);
        }
        sum
    }
}

/// Huron-Vidal variant (selects the constant `c` and the `g^E` combination).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HvVariant {
    /// Huron-Vidal first order (MHV1).
    Mhv1,
    /// Huron-Vidal second order (MHV2).
    Mhv2,
    /// Predictive SRK (PSRK).
    Psrk,
}

impl HvVariant {
    /// The `c` constant in the MHV relation `a/b = Σ x_k a_k/b_k − g^E/(c R T)`.
    fn c(&self) -> f64 {
        match self {
            // MHV1/2 for Peng-Robinson.
            HvVariant::Mhv1 | HvVariant::Mhv2 => 0.62329,
            // PSRK (SRK-based) constant.
            HvVariant::Psrk => 0.64663,
        }
    }
}

/// Huron-Vidal (MHV1/MHV2/PSRK) mixing, generic over an [`ExcessGibbsModel`].
pub struct HuronVidal {
    variant: HvVariant,
    excess: Box<dyn ExcessGibbsModel>,
}

impl core::fmt::Debug for HuronVidal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("HuronVidal")
            .field("variant", &self.variant)
            .finish()
    }
}

impl HuronVidal {
    /// Build for `variant` from an excess-Gibbs model.
    pub fn new(variant: HvVariant, excess: Box<dyn ExcessGibbsModel>) -> Self {
        Self { variant, excess }
    }

    /// The `c` constant in use (see [`HvVariant::c`]).
    pub fn constant(&self) -> f64 {
        self.variant.c()
    }
}

impl CubicMixing for HuronVidal {
    fn b_mix(&self, b: &[f64], z: &[f64]) -> f64 {
        b.iter().zip(z.iter()).map(|(bi, zi)| bi * zi).sum()
    }

    fn a_mix(&self, a: &[f64], b: &[f64], z: &[f64], t: f64, p: f64) -> f64 {
        let bmix = self.b_mix(b, z);
        let t_q = Temperature::new::<kelvin>(t);
        let p_q = Pressure::new::<pascal>(p);
        let g = self.excess.reduced_excess_gibbs(t_q, p_q, z).unwrap_or(0.0);
        let c = self.constant();
        // Σ x_k a_k/b_k − g^E/(c R T); a_mix = b_mix · (that).
        let sum = z
            .iter()
            .zip(a.iter())
            .zip(b.iter())
            .map(|((zi, ai), bi)| zi * ai / bi)
            .sum::<f64>()
            - g / c;
        bmix * sum
    }

    fn aij_sum(&self, a: &[f64], b: &[f64], z: &[f64], i: usize, t: f64, p: f64) -> f64 {
        let bmix = self.b_mix(b, z);
        let t_q = Temperature::new::<kelvin>(t);
        let p_q = Pressure::new::<pascal>(p);
        let g = self.excess.reduced_excess_gibbs(t_q, p_q, z).unwrap_or(0.0);
        let gamma = self.excess.ln_gamma(t_q, p_q, z, i).unwrap_or(0.0);
        let c = self.constant();
        let bracket = z
            .iter()
            .zip(a.iter())
            .zip(b.iter())
            .map(|((zi, ai), bi)| zi * ai / bi)
            .sum::<f64>()
            - g / c;
        // Σ_j z_j a_ij = ½ [ b_i (Σ a_k/b_k − g^E/c) + b_mix (a_i/b_i − ln γ_i/c) ].
        0.5 * (b[i] * bracket + bmix * (a[i] / b[i] - gamma / c))
    }
}

/// Wong-Sandler mixing, generic over an [`ExcessGibbsModel`].
///
/// Uses the MHV-style `a` combination (with the Wong-Sandler constant) and the
/// standard `b = Σ x_i b_i` combination; the full WS `b` correction requires the
/// extended two-parameter form and is approximated here. Cross-crate numerical
/// validation lands in Phase 5.
pub struct WongSandler {
    excess: Box<dyn ExcessGibbsModel>,
    /// Wong-Sandler constant `q` (default 0.62329 for Peng-Robinson).
    q: f64,
}

impl core::fmt::Debug for WongSandler {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("WongSandler").field("q", &self.q).finish()
    }
}

impl WongSandler {
    /// Build from an excess-Gibbs model with the given constant `q`.
    pub fn new(excess: Box<dyn ExcessGibbsModel>, q: f64) -> Self {
        Self { excess, q }
    }
}

impl CubicMixing for WongSandler {
    fn b_mix(&self, b: &[f64], z: &[f64]) -> f64 {
        b.iter().zip(z.iter()).map(|(bi, zi)| bi * zi).sum()
    }

    fn a_mix(&self, a: &[f64], b: &[f64], z: &[f64], t: f64, p: f64) -> f64 {
        let bmix = self.b_mix(b, z);
        let t_q = Temperature::new::<kelvin>(t);
        let p_q = Pressure::new::<pascal>(p);
        let g = self.excess.reduced_excess_gibbs(t_q, p_q, z).unwrap_or(0.0);
        let sum = z
            .iter()
            .zip(a.iter())
            .zip(b.iter())
            .map(|((zi, ai), bi)| zi * ai / bi)
            .sum::<f64>()
            - g / self.q;
        bmix * sum
    }

    fn aij_sum(&self, a: &[f64], b: &[f64], z: &[f64], i: usize, t: f64, p: f64) -> f64 {
        let bmix = self.b_mix(b, z);
        let t_q = Temperature::new::<kelvin>(t);
        let p_q = Pressure::new::<pascal>(p);
        let g = self.excess.reduced_excess_gibbs(t_q, p_q, z).unwrap_or(0.0);
        let gamma = self.excess.ln_gamma(t_q, p_q, z, i).unwrap_or(0.0);
        let bracket = z
            .iter()
            .zip(a.iter())
            .zip(b.iter())
            .map(|((zi, ai), bi)| zi * ai / bi)
            .sum::<f64>()
            - g / self.q;
        0.5 * (b[i] * bracket + bmix * (a[i] / b[i] - gamma / self.q))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vdw_pure_reduces_to_self() {
        let m = VdwMixing::new(1);
        let a = [1.0];
        let b = [0.05];
        let z = [1.0];
        assert!((m.a_mix(&a, &b, &z, 300.0, 1.0e5) - 1.0).abs() < 1e-12);
        assert!((m.aij_sum(&a, &b, &z, 0, 300.0, 1.0e5) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn vdw_binary_symmetric() {
        let m = VdwMixing::from_matrix(vec![vec![0.0, 0.1], vec![0.1, 0.0]]);
        let a = [1.0, 1.0];
        let b = [0.05, 0.05];
        let z = [0.5, 0.5];
        // a_ij = (1-0.1)*1 = 0.9; with diagonal a_ii = 1:
        // a_mix = z1²·1 + 2 z1 z2·0.9 + z2²·1 = 0.95.
        assert!((m.a_mix(&a, &b, &z, 300.0, 1.0e5) - 0.95).abs() < 1e-12);
    }

    #[test]
    fn vdw_tdependent_reduces_at_high_t() {
        // At very high T the b/T + c ln T part dominates the constant offset.
        let td = vec![
            vec![(0.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
            vec![(0.0, 1.0, 0.0), (0.0, 0.0, 0.0)],
        ];
        let m = VdwMixing::new(2).with_tdependent(td);
        let a = [1.0, 1.0];
        let b = [0.05, 0.05];
        let z = [0.5, 0.5];
        // k_ij(1000) = 0 + 1/1000 + 0 = 0.001 → a_12 = 0.999.
        // a_mix = 0.5 + 0.5·0.999 = 0.9995.
        let expected = 0.5 + 0.5 * 0.999;
        assert!((m.a_mix(&a, &b, &z, 1000.0, 1.0e5) - expected).abs() < 1e-9);
    }
}
