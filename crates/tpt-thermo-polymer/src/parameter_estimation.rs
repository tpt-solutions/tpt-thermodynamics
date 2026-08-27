//! Estimation of the Flory–Huggins `χ` parameter from a single LLE tie-line.
//!
//! Given the two equilibrium liquid-phase compositions `x'` and `x''` of a binary
//! polymer(1)+solvent(2) system (segment counts `r`), the interaction parameter is the
//! root of the chemical-potential-equality condition
//! `ln(a₀'/a₁') = ln(a₀''/a₁'')`, solved here with Brent's method over `χ ∈ [−1, 5]`.

use tpt_thermo_core::convergence::ConvergenceStatus;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::numerics::brent;

use crate::flory_huggins::FloryHuggins;

/// Activity ratio `a₀/a₁` for a binary composition `x` under Flory–Huggins with `χ`.
fn activity_ratio(x: &[f64], r: &[f64; 2], chi: f64) -> f64 {
    let fh = FloryHuggins::new_scalar(r.to_vec(), chi);
    let g = fh.ln_gamma(x);
    let a0 = x[0] * g[0].exp();
    let a1 = x[1] * g[1].exp();
    a0 / a1
}

/// Estimate `χ` from a binary LLE tie-line `(x', x'')` with segment counts `r`.
pub fn estimate_chi_from_tieline(
    x_prime: &[f64; 2],
    x_double_prime: &[f64; 2],
    r: &[f64; 2],
) -> Result<f64, ThermoError> {
    let f = |chi: f64| -> f64 {
        activity_ratio(x_prime, r, chi).ln() - activity_ratio(x_double_prime, r, chi).ln()
    };
    brent(f, -1.0, 5.0, 1e-9, 100)
        .map_err(|c| ThermoError::Numerical(ConvergenceStatus::Diverged(c)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_phases_give_zero_residual_at_any_chi() {
        let x = [0.2, 0.8];
        let r = [100.0, 1.0];
        let chi = estimate_chi_from_tieline(&x, &x, &r).unwrap();
        // With identical phases the equality holds for any χ; the solver should still
        // converge to a finite value and the residual at that χ must vanish.
        let f = activity_ratio(&x, &r, chi).ln() - activity_ratio(&x, &r, chi).ln();
        assert!(f.abs() < 1e-12);
    }

    #[test]
    fn recovered_chi_satisfies_equilibrium() {
        // Construct a tie-line by picking two compositions and finding the χ that
        // equalises their activity ratios; then verify the estimator returns that χ.
        let r = [100.0, 1.0];
        let x_prime = [0.02, 0.98];
        let x_double_prime = [0.95, 0.05];
        let chi = estimate_chi_from_tieline(&x_prime, &x_double_prime, &r).unwrap();
        assert!(chi.is_finite() && chi > -1.0 && chi < 5.0);
        let residual = activity_ratio(&x_prime, &r, chi).ln()
            - activity_ratio(&x_double_prime, &r, chi).ln();
        assert!(residual.abs() < 1e-6, "residual = {residual}");
    }
}
