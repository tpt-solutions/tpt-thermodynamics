//! Tangent-plane-distance (TPD) phase-stability test (Michelsen 1982), with
//! successive-substitution + Newton-Raphson refinement, and the
//! [`StabilityAnalyzer`] that implements [`tpt_thermo_core::StabilityTest`].

use crate::linalg::solve_linear;
use crate::phase_volume::PhaseVolume;
use crate::trial_compositions;
use alloc::vec;
use alloc::vec::Vec;
use tpt_thermo_core::quantities::{MolarEnergy, Pressure, Temperature};
use tpt_thermo_core::R;
use tpt_thermo_core::{
    ComponentDatabase, EquationOfState, StabilityResult, StabilityTest, ThermoError,
};
use tpt_thermo_eos_cubic::cubic_solver::Phase;
use uom::si::molar_energy::joule_per_mole;

/// Tolerance on the minimum TPD below which a phase is declared unstable.
pub(crate) const TPD_TOL: f64 = 1e-8;
const SS_MAX_ITER: usize = 200;
const SS_TOL: f64 = 1e-10;
const NEWTON_MAX_ITER: usize = 50;

/// A tangent-plane-distance evaluator for a fixed `(T, P, z)` and EoS.
pub struct TangentPlaneDistance<'a, E: EquationOfState + ?Sized> {
    eos: &'a E,
    volume: &'a dyn PhaseVolume,
    db: &'a dyn ComponentDatabase,
    t: Temperature,
    p: Pressure,
    z: Vec<f64>,
}

impl<'a, E: EquationOfState + ?Sized> TangentPlaneDistance<'a, E> {
    /// Build an evaluator. `db` supplies critical constants for Wilson K-values.
    pub fn new(
        eos: &'a E,
        volume: &'a dyn PhaseVolume,
        db: &'a dyn ComponentDatabase,
        t: Temperature,
        p: Pressure,
        z: Vec<f64>,
    ) -> Self {
        Self {
            eos,
            volume,
            db,
            t,
            p,
            z,
        }
    }

    /// Fugacity coefficients of a trial composition at the requested phase volume.
    fn ln_phi(&self, w: &[f64], phase: Phase) -> Option<Vec<f64>> {
        let v = self.volume.phase_volume(self.t, self.p, w, phase)?;
        (0..w.len())
            .map(|i| self.eos.ln_fugacity_coefficient(self.t, v, w, i).ok())
            .collect()
    }

    /// Tangent-plane distance of trial composition `w` given reference lnφ.
    pub fn tpd(&self, w: &[f64], ln_phi_ref: &[f64], trial_phase: Phase) -> Option<f64> {
        let ln_phi = self.ln_phi(w, trial_phase)?;
        let mut d = 0.0_f64;
        for i in 0..w.len() {
            if w[i] <= 0.0 {
                continue;
            }
            d += w[i] * (w[i].ln() + ln_phi[i] - self.z[i].ln() - ln_phi_ref[i]);
        }
        Some(d)
    }

    /// Wilson K-values for this system.
    pub fn wilson_k(&self) -> Vec<f64> {
        trial_compositions::wilson_k_values(self.db, self.t, self.p)
    }

    /// Minimise the TPD for a reference phase over an opposite-phase trial, via
    /// Michelsen successive substitution followed by a Newton-Raphson polish on
    /// the stationarity conditions.
    pub fn minimize(&self, ref_phase: Phase, trial_phase: Phase) -> Option<TpdSolution> {
        let ln_phi_ref = self.ln_phi(&self.z, ref_phase)?;
        let k0 = self.wilson_k();
        let mut k = k0;
        let mut w = composition_from_k(&self.z, &k, trial_phase);
        let mut converged = false;
        for _ in 0..SS_MAX_ITER {
            let w_new = composition_from_k(&self.z, &k, trial_phase);
            let diff = max_abs_diff(&w, &w_new);
            w = w_new;
            let ln_phi_t = self.ln_phi(&w, trial_phase)?;
            for i in 0..w.len() {
                k[i] = (ln_phi_ref[i] - ln_phi_t[i]).exp();
            }
            if diff < SS_TOL {
                converged = true;
                break;
            }
        }
        // Newton polish only when every component is present (avoids -inf logs).
        let w_ref = if w.iter().all(|v| *v > 1e-12) {
            self.newton_refine(&w, &ln_phi_ref, trial_phase)
        } else {
            None
        };
        let w_final = w_ref.unwrap_or(w);
        let tpd = self.tpd(&w_final, &ln_phi_ref, trial_phase).unwrap_or(0.0);
        Some(TpdSolution {
            composition: w_final,
            tpd,
            converged,
        })
    }

    /// Newton-Raphson refinement of the stationarity conditions in reduced
    /// (logit) variables. Returns `None` if it fails to improve on the SS point.
    fn newton_refine(
        &self,
        w0: &[f64],
        ln_phi_ref: &[f64],
        trial_phase: Phase,
    ) -> Option<Vec<f64>> {
        let n = w0.len();
        if n <= 1 {
            return Some(w0.to_vec());
        }
        let mut x: Vec<f64> = (0..n - 1)
            .map(|i| (w0[i] / w0[n - 1].max(1e-12)).ln())
            .collect();
        let mut w = w_from_x(&x);
        let mut tpd_best = self.tpd(&w, ln_phi_ref, trial_phase)?;
        for _ in 0..NEWTON_MAX_ITER {
            let ln_phi = self.ln_phi(&w, trial_phase)?;
            let g: Vec<f64> = (0..n)
                .map(|i| w[i].ln() + ln_phi[i] - self.z[i].ln() - ln_phi_ref[i])
                .collect();
            let mut r = vec![0.0_f64; n - 1];
            for i in 0..n - 1 {
                r[i] = g[i] - g[n - 1];
            }
            let rnorm2: f64 = r.iter().map(|v| v * v).sum();
            if rnorm2 < 1e-20 {
                break;
            }
            let mut j = vec![vec![0.0_f64; n - 1]; n - 1];
            let h = 1e-6;
            for c in 0..n - 1 {
                let mut xp = x.clone();
                xp[c] += h;
                let wp = w_from_x(&xp);
                let ln_phi_p = self.ln_phi(&wp, trial_phase)?;
                let gp: Vec<f64> = (0..n)
                    .map(|i| wp[i].ln() + ln_phi_p[i] - self.z[i].ln() - ln_phi_ref[i])
                    .collect();
                for i in 0..n - 1 {
                    j[i][c] = (gp[i] - gp[n - 1] - r[i]) / h;
                }
            }
            let b: Vec<f64> = r.iter().map(|v| -v).collect();
            let dx = solve_linear(&j, &b)?;
            // Damped line search: accept the step only if it lowers the TPD.
            let mut step = 1.0_f64;
            let mut improved = false;
            for _ in 0..20 {
                let mut xn = x.clone();
                for i in 0..n - 1 {
                    xn[i] = x[i] + step * dx[i];
                }
                let wn = w_from_x(&xn);
                if let Some(tpdn) = self.tpd(&wn, ln_phi_ref, trial_phase) {
                    if tpdn < tpd_best - 1e-12 {
                        x = xn;
                        w = wn;
                        tpd_best = tpdn;
                        improved = true;
                        break;
                    }
                }
                step *= 0.5;
                if step < 1e-6 {
                    break;
                }
            }
            if !improved {
                break;
            }
        }
        Some(w)
    }
}

/// Result of a TPD minimisation.
#[derive(Debug, Clone)]
pub struct TpdSolution {
    /// Minimised trial (incipient-phase) composition.
    pub composition: Vec<f64>,
    /// Minimum tangent-plane distance (negative ⇒ unstable).
    pub tpd: f64,
    /// Whether successive substitution converged.
    pub converged: bool,
}

/// Build an incipient-phase trial composition from K-values.
///
/// * `Vapor` trial (feed as liquid): `W_i = z_i K_i / (1 + β(K_i − 1))`.
/// * `Liquid` trial (feed as vapor): `W_i = z_i / (1 + β(K_i − 1))`.
fn composition_from_k(z: &[f64], k: &[f64], trial_phase: Phase) -> Vec<f64> {
    let beta = 1.0_f64;
    let n = z.len();
    let mut w = vec![0.0_f64; n];
    for i in 0..n {
        let denom = (1.0 + beta * (k[i] - 1.0)).max(1e-12);
        w[i] = match trial_phase {
            Phase::Vapor => z[i] * k[i] / denom,
            Phase::Liquid => z[i] / denom,
        };
    }
    normalize(&mut w);
    w
}

/// Map reduced logit variables `x_0..x_{n-2}` to a normalised composition
/// `w_i = e^{x_i}/S` (`i < n-1`), `w_{n-1} = 1/S`, `S = 1 + Σ e^{x_j}`.
fn w_from_x(x: &[f64]) -> Vec<f64> {
    let n = x.len() + 1;
    let mut w = vec![0.0_f64; n];
    let mut s = 0.0_f64;
    for i in 0..n - 1 {
        let e = x[i].exp();
        w[i] = e;
        s += e;
    }
    s += 1.0;
    for i in 0..n - 1 {
        w[i] /= s;
    }
    w[n - 1] = 1.0 / s;
    w
}

fn normalize(w: &mut [f64]) {
    let s: f64 = w.iter().sum();
    if s > 0.0 {
        for v in w.iter_mut() {
            *v /= s;
        }
    }
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

/// Phase-stability analyzer: implements [`tpt_thermo_core::StabilityTest`] for
/// any EoS that also satisfies [`PhaseVolume`].
pub struct StabilityAnalyzer<'a, E: EquationOfState + ?Sized> {
    eos: &'a E,
    volume: &'a dyn PhaseVolume,
    db: &'a dyn ComponentDatabase,
}

impl<'a, E: EquationOfState + Send + Sync + ?Sized> StabilityAnalyzer<'a, E> {
    /// Build an analyzer over `eos` (also providing phase volumes) and `db`.
    pub fn new(eos: &'a E, volume: &'a dyn PhaseVolume, db: &'a dyn ComponentDatabase) -> Self {
        Self { eos, volume, db }
    }

    /// Full multiphase classification at `(T, P, z)` (see [`crate::multiphase`]).
    pub fn analyze(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
    ) -> crate::multiphase::MultiphaseResult {
        crate::multiphase::detect_phases(self.eos, self.volume, self.db, t, p, z)
    }
}

impl<'a, E: EquationOfState + Send + Sync + ?Sized> StabilityTest for StabilityAnalyzer<'a, E> {
    fn test(
        &self,
        t: Temperature,
        p: Pressure,
        composition: &[f64],
    ) -> Result<StabilityResult, ThermoError> {
        if composition.len() != self.eos.num_components() {
            return Err(ThermoError::InvalidInput("composition length mismatch"));
        }
        let calc =
            TangentPlaneDistance::new(self.eos, self.volume, self.db, t, p, composition.to_vec());
        let mut stable = true;
        let mut trials = Vec::new();
        for (ref_phase, trial_phase) in
            [(Phase::Vapor, Phase::Liquid), (Phase::Liquid, Phase::Vapor)]
        {
            if let Some(sol) = calc.minimize(ref_phase, trial_phase) {
                if sol.tpd < -TPD_TOL {
                    stable = false;
                }
                trials.push(sol.composition);
            }
        }
        Ok(StabilityResult {
            stable,
            trial_compositions: trials,
            found_second_phase: !stable,
        })
    }

    fn excess_gibbs(
        &self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
    ) -> Result<MolarEnergy, ThermoError> {
        // Residual (departure) Gibbs energy: g^R = R T Σ x_i ln φ_i, evaluated at
        // an available phase volume for the trial composition.
        let phase = if self.volume.phase_volume(t, p, x, Phase::Vapor).is_some() {
            Phase::Vapor
        } else {
            Phase::Liquid
        };
        let v = self.volume.phase_volume(t, p, x, phase).ok_or({
            ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::NumericalIssueReason::NonPhysical,
            ))
        })?;
        let terms: Result<Vec<f64>, ThermoError> = (0..x.len())
            .map(|i| {
                self.eos
                    .ln_fugacity_coefficient(t, v, x, i)
                    .map(|ln| x[i] * ln)
            })
            .collect();
        let g: f64 = terms?.into_iter().sum();
        Ok(MolarEnergy::new::<joule_per_mole>(R * t.value * g))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::phase_volume::BrentPhaseVolume;
    use tpt_thermo_core::component::ComponentDatabase;
    use tpt_thermo_core::quantities::MolarVolume;
    use tpt_thermo_data::SeedComponentDatabase;
    use tpt_thermo_eos_cubic::PengRobinson;
    use uom::si::molar_volume::cubic_meter_per_mole;
    use uom::si::pressure::pascal;
    use uom::si::thermodynamic_temperature::kelvin;

    #[test]
    fn pure_feed_is_stable() {
        let db = SeedComponentDatabase::from_seed();
        let eos = PengRobinson::from_database(&db).unwrap();
        let water = db.index_of("water").unwrap();
        let mut z = vec![0.0; db.num_components()];
        z[water] = 1.0;
        let vol = &eos as &dyn PhaseVolume;
        let ana = StabilityAnalyzer::new(&eos, vol, &db);
        let t = Temperature::new::<kelvin>(600.0);
        let p = Pressure::new::<pascal>(1.0e5);
        let res = ana.test(t, p, &z).unwrap();
        assert!(res.stable, "superheated pure water should be stable");
    }

    #[test]
    fn brent_volume_satisfies_pressure() {
        // Sanity: a volume from BrentPhaseVolume should reproduce P via the EoS.
        let db = SeedComponentDatabase::from_seed();
        let eos = PengRobinson::from_database(&db).unwrap();
        let methane = db.index_of("methane").unwrap();
        let mut z = vec![0.0; db.num_components()];
        z[methane] = 1.0;
        let t = Temperature::new::<kelvin>(300.0);
        let p = Pressure::new::<pascal>(5.0e5);
        let v = BrentPhaseVolume::new(&eos)
            .phase_volume(t, p, &z, Phase::Vapor)
            .unwrap();
        let pp = eos
            .pressure(t, MolarVolume::new::<cubic_meter_per_mole>(v.value), &z)
            .unwrap();
        assert!((pp.value - p.value).abs() / p.value < 1e-6);
    }
}
