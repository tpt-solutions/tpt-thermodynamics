//! Molecular-weight distributions and their averages.
//!
//! Two canonical discrete distributions are provided:
//!
//! * [`schulz_zimm`] — the Schulz–Zimm (most common in polymer physics), parameterised
//!   by the dispersity `Đ = M_w/M_n`.
//! * [`most_probable`] — the Flory "most-probable" distribution (condensation
//!   polymerisation), `Đ = 2`.
//!
//! Both return a probability mass over chain lengths `n = 1..=n_max` (number of repeat
//! units), together with the number-/weight-average molecular weights and dispersity.

use alloc::vec::Vec;

/// Normalised Schulz–Zimm distribution over chain lengths `n = 1..=n_max`.
///
/// Returns the number-fraction `f(n)` such that `Σ f(n) = 1`, with mean chain length
/// `n_mean` and dispersity `d = M_w/M_n` (must be `> 1`).
pub fn schulz_zimm(n_mean: f64, d: f64, n_max: usize) -> Vec<f64> {
    assert!(n_mean > 0.0 && d > 1.0, "need n_mean > 0 and Đ > 1");
    let k = 1.0 / (d - 1.0);
    let a = k.powf(k + 1.0);
    let b = (k / n_mean).powf(k + 1.0);
    let mut f = vec![0.0_f64; n_max];
    let mut total = 0.0_f64;
    for (idx, fi) in f.iter_mut().enumerate() {
        let n = (idx + 1) as f64;
        let v = a * b * n.powf(k) * (-b * n).exp();
        *fi = v;
        total += v;
    }
    if total > 0.0 {
        for fi in f.iter_mut() {
            *fi /= total;
        }
    }
    f
}

/// Normalised "most-probable" (Flory) distribution over chain lengths
/// `n = 1..=n_max`, parameterised by the extent of reaction `p` (`0 < p < 1`).
///
/// `f(n) = (1 − p)·p^{n−1}`, with `n_mean = 1/(1−p)` and `Đ = 2`.
pub fn most_probable(p: f64, n_max: usize) -> Vec<f64> {
    assert!(p > 0.0 && p < 1.0, "need 0 < p < 1");
    let mut f = vec![0.0_f64; n_max];
    let mut total = 0.0_f64;
    for (idx, fi) in f.iter_mut().enumerate() {
        let n = (idx + 1) as f64;
        let v = (1.0 - p) * p.powf(n - 1.0);
        *fi = v;
        total += v;
    }
    if total > 0.0 {
        for fi in f.iter_mut() {
            *fi /= total;
        }
    }
    f
}

/// Averaged molecular-weight descriptors of a discrete distribution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MwdMoments {
    /// Number-average chain length `M_n` (in repeat-unit count).
    pub number_average: f64,
    /// Weight-average chain length `M_w`.
    pub weight_average: f64,
    /// Dispersity `Đ = M_w / M_n`.
    pub dispersity: f64,
}

/// Compute number/weight averages and dispersity from a number-fraction `f(n)`
/// (length `n = 1..=f.len()`).
pub fn moments(f: &[f64]) -> MwdMoments {
    let mut n_sum = 0.0_f64;
    let mut nw = 0.0_f64;
    let mut w_sum = 0.0_f64;
    for (idx, &fi) in f.iter().enumerate() {
        let n = (idx + 1) as f64;
        n_sum += fi * n;
        let w = fi * n;
        nw += w * n;
        w_sum += w;
    }
    let mn = if n_sum > 0.0 { n_sum } else { 1.0 };
    let mw = if w_sum > 0.0 { nw / w_sum } else { 1.0 };
    MwdMoments {
        number_average: mn,
        weight_average: mw,
        dispersity: mw / mn,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn most_probable_has_dispersity_two() {
        let f = most_probable(0.5, 200);
        let m = moments(&f);
        assert!((m.number_average - 2.0).abs() < 1e-2, "Mn = {}", m.number_average);
        assert!((m.dispersity - 2.0).abs() < 1e-2, "Đ = {}", m.dispersity);
    }

    #[test]
    fn schulz_zimm_matches_target_dispersity() {
        let f = schulz_zimm(500.0, 1.5, 4000);
        let m = moments(&f);
        assert!((m.number_average - 500.0).abs() / 500.0 < 1e-2, "Mn = {}", m.number_average);
        assert!((m.dispersity - 1.5).abs() / 1.5 < 1e-2, "Đ = {}", m.dispersity);
    }

    #[test]
    fn distributions_normalised() {
        let f = most_probable(0.9, 1000);
        let s: f64 = f.iter().sum();
        assert!((s - 1.0).abs() < 1e-9);
        let g = schulz_zimm(10.0, 1.2, 500);
        let s2: f64 = g.iter().sum();
        assert!((s2 - 1.0).abs() < 1e-9);
    }
}
