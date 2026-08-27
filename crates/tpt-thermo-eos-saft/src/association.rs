//! The association (hydrogen-bonding) contribution to SAFT.
//!
//! Implements the Gross & Sadowski (2001) association term with the standard
//! equivalent-site simplification:
//!
//! ```text
//! a^assoc/(RT) = Σ_i x_i N^A_i [ ln X_i + (1 - X_i)/2 ]
//! X_i = 1 / (1 + ρ Σ_j x_j N^A_j X_j Δ_ij)
//! Δ_ij = g_ij(σ_ij) σ_ij³ κ^AB_ij [ exp(ε^AB_ij/(kT)) - 1 ]
//! ```
//!
//! `X_i` is the fraction of *non-bonded* sites on component `i`. It is solved
//! for each associating component by a Newton-Raphson iteration on the coupled
//! fixed-point system; the routine returns a [`ConvergenceStatus`] so callers
//! can report non-convergence instead of panicking.

use crate::parameters::SaftParameters;
use tpt_thermo_core::convergence::ConvergenceStatus;
use tpt_thermo_core::numerics::ROOT_MAX_ITER;

/// Result of the association solve.
#[derive(Debug, Clone)]
pub struct AssociationResult {
    /// Fraction of non-bonded sites per component (`0.0` for non-associating
    /// species).
    pub x_nondim: Vec<f64>,
    /// Reduced residual Helmholtz of association, `a^assoc/(RT)`.
    pub ares: f64,
}

/// Number of association sites per component (0 for non-associating).
fn n_sites(params: &SaftParameters, i: usize) -> usize {
    params
        .component(i)
        .association
        .map(|a| a.scheme.num_sites())
        .unwrap_or(0)
}

/// Cross-association strength `Δ_ij` (m³) given the mixture packing fraction
/// `zeta_3` (needs only the temperature-dependent hard-sphere diameters, which
/// already live inside `params` as the temperature-independent `σ`).
fn delta_ij(params: &SaftParameters, i: usize, j: usize, t: f64, zeta_3: f64) -> f64 {
    let ci = params.component(i);
    let cj = params.component(j);
    let ai = match ci.association {
        Some(a) => a,
        None => return 0.0,
    };
    let aj = match cj.association {
        Some(a) => a,
        None => return 0.0,
    };
    let sigma_ij = 0.5 * (ci.sigma + cj.sigma) * 1.0e-10; // Å → m
    let eps_ab = 0.5 * (ai.epsilon_ab_k + aj.epsilon_ab_k); // K, arithmetic mean
    let kappa = (ai.kappa_ab * aj.kappa_ab).sqrt(); // geometric mean
                                                    // Radial distribution at contact (mixture packing fraction).
    let g = (1.0 - 0.5 * zeta_3) / (1.0 - zeta_3).powi(3);
    let boltz = (eps_ab / t).exp() - 1.0;
    g * sigma_ij.powi(3) * kappa * boltz
}

/// Solve the association fixed-point system at `(T, ρ, x)`.
///
/// `rho` is the molecular number density (m⁻³) and `zeta_3` the mixture packing
/// fraction (computed by the hard-chain reference). Returns the non-bonded
/// fractions and the reduced association Helmholtz energy.
pub fn solve_association(
    params: &SaftParameters,
    x: &[f64],
    rho: f64,
    t: f64,
    zeta_3: f64,
) -> Result<AssociationResult, ConvergenceStatus> {
    let n = params.num_components();
    let sites: Vec<usize> = (0..n).map(|i| n_sites(params, i)).collect();
    let any_assoc = sites.iter().any(|&s| s > 0);

    if !any_assoc || rho <= 0.0 {
        return Ok(AssociationResult {
            x_nondim: vec![0.0; n],
            ares: 0.0,
        });
    }

    // Precompute Δ_ij only for associating pairs.
    let mut delta = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        if sites[i] == 0 {
            continue;
        }
        for j in 0..n {
            if sites[j] == 0 {
                continue;
            }
            delta[i][j] = delta_ij(params, i, j, t, zeta_3);
        }
    }

    // Initial guess: all sites non-bonded.
    let mut xv = vec![1.0_f64; n];
    let tol = 1.0e-12;
    let mut converged = false;
    for _ in 0..ROOT_MAX_ITER {
        // S_i = ρ Σ_j x_j N^A_j X_j Δ_ij
        let mut s = vec![0.0_f64; n];
        for i in 0..n {
            if sites[i] == 0 {
                continue;
            }
            let mut acc = 0.0;
            for j in 0..n {
                if sites[j] == 0 {
                    continue;
                }
                acc += x[j] * sites[j] as f64 * xv[j] * delta[i][j];
            }
            s[i] = rho * acc;
        }
        // Fixed-point residual: F_i = X_i - 1/(1 + S_i).
        let mut max_f: f64 = 0.0;
        let mut f = vec![0.0_f64; n];
        for i in 0..n {
            if sites[i] == 0 {
                continue;
            }
            let denom = 1.0 + s[i];
            f[i] = xv[i] - 1.0 / denom;
            max_f = max_f.max(f[i].abs());
        }
        if max_f < tol {
            converged = true;
            break;
        }

        // Newton step: dF_i/dX_k = δ_ik + (1/(1+S_i)²) ρ x_k N^A_k Δ_ik.
        let mut jac = vec![vec![0.0_f64; n]; n];
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n {
            if sites[i] == 0 {
                continue;
            }
            let denom = 1.0 + s[i];
            for k in 0..n {
                if sites[k] == 0 {
                    continue;
                }
                let mut jik = 0.0;
                if i == k {
                    jik += 1.0;
                }
                jik += (1.0 / denom / denom) * rho * x[k] * sites[k] as f64 * delta[i][k];
                jac[i][k] = jik;
            }
            rhs[i] = -f[i];
        }
        let dx = solve_linear(&jac, &rhs, n).map_err(|_| {
            ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::convergence::NumericalIssueReason::SingularJacobian,
            )
        })?;
        let mut max_step: f64 = 0.0;
        for i in 0..n {
            if sites[i] == 0 {
                continue;
            }
            xv[i] += dx[i];
            max_step = max_step.max(dx[i].abs());
            if !xv[i].is_finite() {
                return Err(ConvergenceStatus::Diverged(
                    tpt_thermo_core::convergence::DivergenceReason::NonFinite,
                ));
            }
            xv[i] = xv[i].clamp(1.0e-12, 1.0);
        }
        if max_step < tol {
            converged = true;
            break;
        }
    }

    if !converged {
        return Err(ConvergenceStatus::Diverged(
            tpt_thermo_core::convergence::DivergenceReason::MaxIterations,
        ));
    }

    // Reduced association Helmholtz energy.
    let mut ares = 0.0;
    for i in 0..n {
        if sites[i] == 0 {
            continue;
        }
        ares += x[i] * sites[i] as f64 * (xv[i].ln() + 0.5 * (1.0 - xv[i]));
    }

    Ok(AssociationResult { x_nondim: xv, ares })
}

/// Solve a small dense linear system `A x = b` by Gaussian elimination with
/// partial pivoting.
#[allow(clippy::needless_range_loop)]
fn solve_linear(a: &[Vec<f64>], b: &[f64], n: usize) -> Result<Vec<f64>, ()> {
    let mut m = vec![vec![0.0_f64; n]; n];
    let mut rhs = b.to_vec();
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a[i][j];
        }
    }
    for col in 0..n {
        // Partial pivot.
        let mut piv = col;
        let mut max = m[col][col].abs();
        for r in (col + 1)..n {
            if m[r][col].abs() > max {
                max = m[r][col].abs();
                piv = r;
            }
        }
        if max < 1.0e-30 {
            return Err(());
        }
        m.swap(col, piv);
        rhs.swap(col, piv);
        let diag = m[col][col];
        for r in (col + 1)..n {
            let factor = m[r][col] / diag;
            for c in col..n {
                m[r][c] -= factor * m[col][c];
            }
            rhs[r] -= factor * rhs[col];
        }
    }
    let mut x = vec![0.0_f64; n];
    for r in (0..n).rev() {
        let mut acc = rhs[r];
        for c in (r + 1)..n {
            acc -= m[r][c] * x[c];
        }
        x[r] = acc / m[r][r];
    }
    Ok(x)
}

/// `g_ij(σ_ij)` radial distribution at contact for a given packing fraction.
#[allow(dead_code)]
pub(crate) fn contact_rdf(zeta_3: f64) -> f64 {
    (1.0 - 0.5 * zeta_3) / (1.0 - zeta_3).powi(3)
}
