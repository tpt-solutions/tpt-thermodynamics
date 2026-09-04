//! Parameter-estimation utilities (spec 3d follow-up).
//!
//! These routines fit an equation-of-state parameter to experimental data. The
//! headline utility is [`fit_binary_kij`]: given isothermal bubble-pressure VLE
//! data for a binary, it recovers the PR/SRK binary interaction parameter `k_ij`
//! by minimising the sum of squared bubble-pressure residuals (using the core's
//! [`brent_minimize`](tpt_thermo_core::numerics::brent_minimize)).
//!
//! [`bubble_pressure`] is a self-contained isothermal bubble-point solver: it
//! locates the pressure at which an incipient vapor phase appears for a given
//! liquid composition. It is **flash-based**: at every trial pressure it runs the
//! Michelsen incipient-phase calculation (successive substitution on the
//! K-values `K_i = φ_iᴸ/φ_iⱽ` with GDEM acceleration, seeded from Wilson
//! K-values), which is exactly the inner loop of a full PT flash. The fugacity
//! residual `Σ_i K_i z_i − 1` is then a smooth function of pressure that crosses
//! zero at the bubble point — robust for associating and near-critical binaries
//! where the previous plain-successive-substitution residual did not converge
//! smoothly.
//!
//! Known limitation (tracked in `todo.md`): for binaries containing a
//! supercritical light component (e.g. water/methane, CO₂/methane) the
//! incipient-phase solve can still fail to bracket the bubble, because a
//! fully stability-tested flash (tangent-plane-distance) is required to reject
//! the spurious two-phase solutions the bare successive-substitution flash
//! converges to there. That stability test is a repo-wide follow-up (the
//! `tpt-thermo-flash` PT flash itself lacks it).

use alloc::boxed::Box;
use alloc::vec;
use alloc::vec::Vec;
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::convergence::ConvergenceStatus;
use tpt_thermo_core::eos::EquationOfState;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::numerics::{brent, brent_minimize};
use tpt_thermo_core::quantities::{Pressure, Temperature};
use uom::si::pressure::pascal;

use crate::cubic_solver::Phase;
use crate::mixing::VdwMixing;
use crate::pr::PengRobinson;

/// Scalar dominant-eigenvalue (GDEM) acceleration memory, ported from the flash
/// crate's `acceleration::gdem_step`. The parameter-estimation utilities need it
/// to converge incipient-phase K-values without depending on `tpt-thermo-flash`
/// (which depends on this crate, so a direct dependency would be cyclic).
struct GdemState {
    prev_d: Option<Vec<f64>>,
}

impl GdemState {
    fn new() -> Self {
        Self { prev_d: None }
    }

    /// Return the accelerated `Kⁿ⁺¹` given `Kⁿ` (`k_old`) and the freshly evaluated
    /// `Kⁿ⁺¹` (`k_new`). Extrapolates past the slow successive-substitution manifold
    /// when the dominant eigenvalue `λ ∈ (0, 1)`; otherwise plain substitution.
    fn step(&mut self, k_old: &[f64], k_new: &[f64]) -> Vec<f64> {
        let d_new: Vec<f64> = k_new.iter().zip(k_old.iter()).map(|(a, b)| a - b).collect();
        let out = match &self.prev_d {
            Some(d_old) if norm2(d_old) > 1e-30 => {
                let lambda = dot(&d_new, d_old) / norm2(d_old);
                if (0.0..1.0).contains(&lambda) {
                    let g = (1.0 / (1.0 - lambda)).clamp(-3.0, 3.0);
                    k_new
                        .iter()
                        .zip(d_new.iter())
                        .map(|(kn, d)| kn + g * d)
                        .collect()
                } else {
                    k_new.to_vec()
                }
            }
            _ => k_new.to_vec(),
        };
        let out = if out.iter().all(|v| v.is_finite()) {
            out
        } else {
            k_new.to_vec()
        };
        self.prev_d = Some(d_new);
        out
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0_f64;
    for i in 0..a.len().min(b.len()) {
        s += a[i] * b[i];
    }
    s
}

fn norm2(a: &[f64]) -> f64 {
    dot(a, a)
}

/// Wilson correlation for initial K-values at `(T, P)` using the pure-component
/// critical data carried by the EoS. Seeding with Wilson rather than the naive
/// `y = x` start is what lets the incipient-phase iteration converge for
/// strongly non-ideal (associating / near-critical) binaries.
fn wilson_k(eos: &PengRobinson, t: Temperature, p: f64, z: &[f64]) -> Vec<f64> {
    let n = z.len();
    let mut k = vec![1.0_f64; n];
    let tk = t.value.max(1.0);
    let pv = p.max(1.0);
    for i in 0..n {
        if z[i] <= 0.0 {
            continue;
        }
        if let Some((tc, pc, omega)) = eos.engine().component_critical(i) {
            if tc > 0.0 && pc > 0.0 {
                let exponent = 5.3727 * (1.0 + omega) * (1.0 - tc / tk);
                k[i] = (pc / pv) * exponent.exp();
            }
        }
    }
    k
}

/// Classify a single (real) compressibility root: vapor-like roots are large
/// (`Z ≳ 0.5`), liquid-like roots are small. Used only in the genuine
/// single-phase region (fewer than two real roots) to decide which side of the
/// two-phase envelope a pressure lies on.
fn is_vapor_like(roots: &[f64]) -> bool {
    let z = roots
        .iter()
        .cloned()
        .filter(|r| *r > 1e-6)
        .fold(0.0_f64, f64::max);
    z > 0.5
}

/// Flash-based incipient-phase K-values for a liquid composition `z` at pressure
/// `p` (Pa). Runs the Michelsen incipient-phase iteration:
///
/// * `K_i = φ_iᴸ(z) / φ_iⱽ(y)` with the liquid evaluated at the fixed liquid `z`
///   and the vapor at the trial composition `y`;
/// * `y_i = K_i z_i / Σ_j K_j z_j` (Rachford–Rice split at `β → 0`);
/// * successive substitution accelerated by [`GdemState`].
///
/// Returns `Some((K, residual))` with `residual = Σ_i K_i z_i − 1` (the fugacity
/// residual that vanishes at the bubble point) when a vapor root exists, or
/// `None` when no vapor root is available at this pressure (i.e. the system is
/// single-phase, above the bubble point).
fn incipient_k(eos: &PengRobinson, t: Temperature, z: &[f64], p: f64) -> Option<(Vec<f64>, f64)> {
    let n = z.len();
    let vl = eos
        .solve_phase(t, Pressure::new::<pascal>(p), z, Phase::Liquid)
        .ok()?;
    let mut k = wilson_k(eos, t, p, z);
    let mut y: Vec<f64> = {
        let s = k.iter().zip(z.iter()).map(|(ki, zi)| ki * zi).sum::<f64>();
        if s <= 0.0 {
            return None;
        }
        (0..n).map(|i| k[i] * z[i] / s).collect()
    };

    let mut mem = GdemState::new();
    for _ in 0..300 {
        let vv = eos
            .solve_phase(t, Pressure::new::<pascal>(p), &y, Phase::Vapor)
            .ok()?;
        let mut knew = vec![0.0_f64; n];
        let mut ok = true;
        for i in 0..n {
            let fl = match eos.ln_fugacity_coefficient(t, vl, z, i) {
                Ok(v) => v,
                Err(_) => {
                    ok = false;
                    break;
                }
            };
            let fv = match eos.ln_fugacity_coefficient(t, vv, &y, i) {
                Ok(v) => v,
                Err(_) => {
                    ok = false;
                    break;
                }
            };
            if !fl.is_finite() || !fv.is_finite() {
                ok = false;
                break;
            }
            knew[i] = (fl - fv).exp();
        }
        if !ok {
            return None;
        }
        let k_acc = mem.step(&k, &knew);
        let s = k_acc
            .iter()
            .zip(z.iter())
            .map(|(ki, zi)| ki * zi)
            .sum::<f64>();
        if !s.is_finite() || s <= 0.0 {
            return None;
        }
        let mut ynew = vec![0.0_f64; n];
        let mut maxd = 0.0_f64;
        for i in 0..n {
            ynew[i] = k_acc[i] * z[i] / s;
            maxd = maxd.max((ynew[i] - y[i]).abs());
        }
        k = k_acc;
        y = ynew;
        if maxd < 1e-11 {
            break;
        }
    }

    let s = k.iter().zip(z.iter()).map(|(ki, zi)| ki * zi).sum::<f64>();
    if !s.is_finite() {
        return None;
    }
    Some((k, s - 1.0))
}

/// Fugacity residual `Σ_i K_i z_i − 1` for a liquid composition `z` at pressure
/// `p` (Pa), where `K_i = φ_iᴸ/φ_iⱽ` with the liquid root at `z` and the vapor
/// root converged by the incipient-phase flash.
///
/// Outside the two-phase region (only one real root) a bracket sign is returned:
/// `+1` below the envelope (vapor-like, two-phase side) and `−1` above it
/// (liquid-like, single-phase side), so the bubble point is the unique
/// `Σ K_i z_i − 1 = 0` crossing as pressure rises.
fn bubble_residual(eos: &PengRobinson, t: Temperature, z: &[f64], p: f64) -> f64 {
    let roots = eos.engine().z_roots(t, Pressure::new::<pascal>(p), z);
    if roots.len() < 2 {
        return if is_vapor_like(&roots) { 1.0 } else { -1.0 };
    }
    match incipient_k(eos, t, z, p) {
        Some((_k, r)) => r,
        // No vapor root at this pressure: the system is single-phase (above the
        // bubble point). Return the negative bracket sign so Brent's invariant
        // `f(lo) > 0 > f(hi)` is preserved.
        None => -1.0,
    }
}

/// Isothermal bubble pressure (Pa) of a liquid mixture `z` at temperature `t`
/// using the Peng–Robinson EoS.
///
/// The bubble point is the highest pressure at the given `t` and liquid
/// composition `z` for which a vapor phase can exist. The solver brackets it
/// between the liquid-single-phase region (above the bubble point, residual
/// negative) and the two-phase region (below it, residual positive), then drives
/// the flash-based fugacity residual to zero with Brent's method.
pub fn bubble_pressure(
    eos: &PengRobinson,
    t: Temperature,
    z: &[f64],
) -> Result<Pressure, ThermoError> {
    let n = z.len();
    let pc_max = (0..n)
        .map(|i| {
            eos.critical_point_pure(i)
                .ok()
                .map(|c| c.1.value)
                .unwrap_or(1.0e9)
        })
        .fold(0.0_f64, f64::max);

    // `hi` starts safely above the bubble point (liquid single phase, above the
    // cricondenbar which cannot exceed the maximum pure critical pressure). `lo`
    // is found by stepping down from `hi` until the residual becomes positive —
    // i.e. just below the bubble point, inside the two-phase region. The bubble
    // is then the unique `Σ K_i z_i − 1 = 0` crossing in `[lo, hi]`.
    let hi = (1.5 * pc_max).clamp(1.0e4, 1.0e8);
    let mut p = hi;
    let mut lo = 1.0e3_f64;
    let mut found = false;
    for _ in 0..600 {
        let r = bubble_residual(eos, t, z, p);
        if r > 0.0 {
            lo = p;
            found = true;
            break;
        }
        p *= 0.9;
        if p < 1.0e1 {
            break;
        }
    }

    if !found || lo >= hi {
        return Err(ThermoError::Numerical(ConvergenceStatus::NotConverged));
    }

    match brent(|p: f64| bubble_residual(eos, t, z, p), lo, hi, 1e-9, 400).ok() {
        Some(p) => Ok(Pressure::new::<pascal>(p)),
        None => Err(ThermoError::Numerical(ConvergenceStatus::NotConverged)),
    }
}

/// Fit a single binary interaction parameter `k_ij` to isothermal bubble-pressure
/// VLE data `points = [(T, P_exp, x1)]` for a **two-component** database
/// (components 0 and 1). Returns `(k_ij, rms_residual_pa)`.
///
/// The Peng–Robinson EoS is rebuilt for each trial `k_ij`; the objective is the
/// sum of squared bubble-pressure residuals, minimised over `k_ij ∈ [-0.5, 0.5]`.
pub fn fit_binary_kij(
    db: &dyn ComponentDatabase,
    points: &[(Temperature, Pressure, f64)],
) -> Result<(f64, f64), ThermoError> {
    let objective = |k: f64| -> f64 {
        let mixing = VdwMixing::from_matrix(vec![vec![0.0, k], vec![k, 0.0]]);
        let eos = match PengRobinson::with_mixing(db, Box::new(mixing)) {
            Ok(e) => e,
            Err(_) => return f64::INFINITY,
        };
        let mut ssr = 0.0_f64;
        for (t, p_exp, x1) in points {
            let mut z = vec![0.0_f64; 2];
            z[0] = *x1;
            z[1] = 1.0 - *x1;
            match bubble_pressure(&eos, *t, &z) {
                Ok(p_calc) => {
                    let d = p_calc.value - p_exp.value;
                    ssr += d * d;
                }
                Err(_) => return f64::INFINITY,
            }
        }
        ssr
    };
    let (k_fit, ssr) = brent_minimize(objective, -0.5_f64, 0.5_f64, 1e-7, 200);
    Ok((k_fit, ssr.sqrt()))
}
