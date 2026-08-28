//! Flory–Huggins activity model for polymer solutions.
//!
//! The mean-field Flory–Huggins (FH) theory expresses the activity coefficient of
//! component `i` from its volume (segment) fraction `φ_i` and its segment count
//! `r_i`:
//!
//! ```text
//! ln γ_i = ln φ_i + (1 − 1/r_i)·(1 − φ_i) + Σ_{j≠i} χ_ij·φ_j²
//! ```
//!
//! which reduces to the classic binary form `ln γ_1 = ln φ_1 + (1 − 1/r_1)·φ_2 +
//! χ·φ_2²`. A single scalar `χ` treats every cross-interaction equally (`χ_ij = χ`
//! for `i ≠ j`); a full `χ` matrix may be supplied for multicomponent generality.

/// Flory–Huggins model parameters.
#[derive(Debug, Clone)]
pub struct FloryHuggins {
    /// Segment counts `r_i` (degree of polymerisation; `r = 1` for the solvent).
    pub r: Vec<f64>,
    /// Cross interaction `χ_ij` matrix (`χ_ii` is ignored / forced to zero).
    pub chi: Vec<Vec<f64>>,
}

impl FloryHuggins {
    /// Build from segment counts `r` and a single scalar `χ` applied to every
    /// off-diagonal pair.
    pub fn new_scalar(r: Vec<f64>, chi: f64) -> Self {
        let n = r.len();
        let chi = (0..n)
            .map(|i| (0..n).map(|j| if i == j { 0.0 } else { chi }).collect())
            .collect();
        Self { r, chi }
    }

    /// Build from an explicit `χ` matrix.
    pub fn new_matrix(r: Vec<f64>, mut chi: Vec<Vec<f64>>) -> Self {
        let n = r.len();
        for i in 0..n {
            if chi[i].len() != n {
                chi[i] = vec![0.0; n];
            }
            chi[i][i] = 0.0;
        }
        Self { r, chi }
    }

    /// Segment (volume) fractions `φ_i = x_i·r_i / Σ x_j·r_j` for mole fractions `x`.
    pub fn volume_fractions(&self, x: &[f64]) -> Vec<f64> {
        let r_sum: f64 = x.iter().zip(&self.r).map(|(xi, ri)| xi * ri).sum();
        if r_sum <= 0.0 {
            return vec![0.0; x.len()];
        }
        x.iter()
            .zip(&self.r)
            .map(|(xi, ri)| xi * ri / r_sum)
            .collect()
    }

    /// Natural log of the activity coefficient of every component at mole fractions
    /// `x` (normalised internally).
    pub fn ln_gamma(&self, x: &[f64]) -> Vec<f64> {
        let n = self.r.len();
        assert_eq!(x.len(), n, "composition length mismatch");
        let sum: f64 = x.iter().sum();
        let x = if sum.abs() > 0.0 {
            x.iter().map(|v| v / sum).collect::<Vec<_>>()
        } else {
            x.to_vec()
        };
        let phi = self.volume_fractions(&x);
        (0..n)
            .map(|i| {
                let combinatorial = (phi[i] / x[i]).ln() + (1.0 - 1.0 / self.r[i]) * (1.0 - phi[i]);
                let interaction: f64 = (0..n)
                    .filter(|&j| j != i)
                    .map(|j| self.chi[i][j] * phi[j] * phi[j])
                    .sum();
                combinatorial + interaction
            })
            .collect()
    }

    /// Activity coefficient of component `i` at mole fractions `x`.
    pub fn gamma(&self, x: &[f64], i: usize) -> f64 {
        self.ln_gamma(x)[i].exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_reduces_to_classic_form() {
        // Polymer (r=1000) + solvent (r=1), χ=0.4. With 1% polymer by mole the
        // *volume* fraction of solvent is ~0.09 (polymer dominates volume).
        let fh = FloryHuggins::new_scalar(vec![1000.0, 1.0], 0.4);
        let x = [0.01, 0.99]; // 1% polymer by mole
        let phi = fh.volume_fractions(&x);
        let phi_solvent = (x[1] * 1.0) / (x[0] * 1000.0 + x[1] * 1.0);
        assert!(
            (phi[1] - phi_solvent).abs() < 1e-6,
            "φ_solvent = {}",
            phi[1]
        );
        let ln_g = fh.ln_gamma(&x);
        // Multicomponent-generalised FH: ln γ_solvent = ln(φ_2/x_2) + χ·φ_1²
        // (r_2 = 1 ⇒ combinatorial middle term vanishes).
        let expected = (phi[1] / x[1]).ln() + 0.4 * phi[0] * phi[0];
        assert!(
            (ln_g[1] - expected).abs() < 1e-9,
            "got {} exp {}",
            ln_g[1],
            expected
        );
        assert!(ln_g[0] > 0.0, "polymer γ must be > 1 (poor solvent) here");
    }

    #[test]
    fn ideal_limit() {
        // χ = 0 and equal segment counts → ideal mixing (γ = 1 everywhere).
        let fh = FloryHuggins::new_scalar(vec![1.0, 1.0], 0.0);
        let x = [0.3, 0.7];
        let phi = fh.volume_fractions(&x);
        let ln_g = fh.ln_gamma(&x);
        for i in 0..2 {
            // For r_i = 1, ln γ_i = ln(φ_i/x_i); with equal r, φ_i = x_i → ln γ = 0.
            assert!(phi[i].abs() > 0.0);
            assert!(ln_g[i].abs() < 1e-12, "ln γ[{}] = {}", i, ln_g[i]);
        }
    }
}
