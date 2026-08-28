//! Molecular-weight distributions and their moments.
//!
//! Implements the two work-horse discrete distributions used for polymer
//! characterization — the Flory *most-probable* (step-growth) distribution and the
//! Schulz–Zimm (living-polymer / broad) distribution — together with helpers to
//! compute number- and weight-average degree of polymerization and the dispersity
//! `Đ = M_w / M_n`.

use std::vec::Vec;

/// A molecular-weight distribution over chain lengths `r = 1..=max_r`.
#[derive(Debug, Clone)]
pub struct MolecularWeightDistribution {
    /// Number fraction `n_r` of chains of length `r` (index `r-1`).
    pub number_fraction: Vec<f64>,
    /// Weight fraction `w_r` of chains of length `r` (index `r-1`).
    pub weight_fraction: Vec<f64>,
}

impl MolecularWeightDistribution {
    /// Number-average degree of polymerization `x_n = Σ r n_r`.
    pub fn number_average(&self) -> f64 {
        self.number_fraction
            .iter()
            .enumerate()
            .map(|(i, &nr)| (i as f64 + 1.0) * nr)
            .sum()
    }

    /// Weight-average degree of polymerization `x_w = Σ r w_r`.
    pub fn weight_average(&self) -> f64 {
        self.weight_fraction
            .iter()
            .enumerate()
            .map(|(i, &wr)| (i as f64 + 1.0) * wr)
            .sum()
    }

    /// Dispersity `Đ = x_w / x_n`.
    pub fn dispersity(&self) -> f64 {
        self.weight_average() / self.number_average()
    }
}

/// Flory most-probable (step-growth) distribution for a number-average degree of
/// polymerization `x_n` (equivalently extent of reaction `p = 1 - 1/x_n`).
///
/// `n_r = (1-p)² p^{r-1}`; `x_n = 1/(1-p)`, `x_w = (1+p)/(1-p)`, `Đ = 1+p`.
pub fn most_probable(x_n: f64, max_r: usize) -> MolecularWeightDistribution {
    let p = 1.0 - 1.0 / x_n;
    let mut nf = vec![0.0; max_r];
    let mut total = 0.0_f64;
    for r in 1..=max_r {
        let v = (1.0 - p) * (1.0 - p) * p.powi(r as i32 - 1);
        nf[r - 1] = v;
        total += v;
    }
    for nr in nf.iter_mut() {
        *nr /= total;
    }
    let mut wf = vec![0.0; max_r];
    for (i, &nr) in nf.iter().enumerate() {
        let r = i as f64 + 1.0;
        wf[i] = r * nr / (x_n); // Σ r·n_r/n = x_n, so normalize by x_n
    }
    MolecularWeightDistribution {
        number_fraction: nf,
        weight_fraction: wf,
    }
}

/// Schulz–Zimm distribution with a target number-average `x_n` and dispersity
/// `d = x_w / x_n` (`d > 1`). The shape parameter is `k = 1/(d-1)`.
///
/// `n_r ∝ r^{k-1} exp(-k r / x_n)`.
pub fn schulz_zimm(x_n: f64, d: f64, max_r: usize) -> MolecularWeightDistribution {
    assert!(d > 1.0, "dispersity must exceed 1");
    let k = 1.0 / (d - 1.0);
    let mut raw = vec![0.0; max_r];
    let mut sum = 0.0_f64;
    for r in 1..=max_r {
        let v = (r as f64).powf(k - 1.0) * (-k * (r as f64) / x_n).exp();
        raw[r - 1] = v;
        sum += v;
    }
    let nf: Vec<f64> = raw.iter().map(|&v| v / sum).collect();
    let x_n_actual: f64 = nf
        .iter()
        .enumerate()
        .map(|(i, &nr)| (i as f64 + 1.0) * nr)
        .sum();
    let mut wf = vec![0.0; max_r];
    for (i, &nr) in nf.iter().enumerate() {
        let r = i as f64 + 1.0;
        wf[i] = r * nr / x_n_actual;
    }
    MolecularWeightDistribution {
        number_fraction: nf,
        weight_fraction: wf,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn most_probable_moments() {
        let d = most_probable(10.0, 200);
        assert!((d.number_average() - 10.0).abs() < 1e-6);
        // Đ = 1 + p = 1 + (1 - 1/x_n) = 2 - 1/10 = 1.9.
        assert!((d.dispersity() - 1.9).abs() < 1e-3);
    }

    #[test]
    fn schulz_zimm_recovers_target() {
        let d = schulz_zimm(50.0, 1.5, 500);
        assert!(
            (d.number_average() - 50.0).abs() < 1.0,
            "number_average = {}",
            d.number_average()
        );
        assert!(
            (d.dispersity() - 1.5).abs() < 1e-2,
            "dispersity = {}",
            d.dispersity()
        );
    }
}
