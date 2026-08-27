//! Isothermal/isobaric (Rachford–Rice) PT flash, the foundation the other flash
//! specifications build on.

use alloc::vec::Vec;
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::convergence::{ConvergenceStatus, NumericalIssueReason};
use tpt_thermo_core::eos::EquationOfState;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};

use crate::acceleration::{gdem_step, AccelerationMemory};
use crate::initialization::wilson_k_values;
use crate::phase_volume::{phase_volume, Phase};
use crate::rachford_rice::{rachford_rice, RachfordRiceResult};
use crate::FlashError;

/// Default successive-substitution iteration budget.
pub const PT_MAX_ITER: usize = 200;
/// Default relative K-value convergence tolerance.
pub const PT_TOL: f64 = 1e-6;

/// Whether a flash result is two-phase or single-phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseFlag {
    /// Vapor–liquid two-phase.
    TwoPhase,
    /// Single phase (all liquid or all vapor).
    SinglePhase,
}

/// Result of a flash calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct FlashResult {
    /// Vapor mole fraction `β`.
    pub vapor_fraction: f64,
    /// Liquid-phase composition `x` (length = number of components).
    pub liquid_composition: Vec<f64>,
    /// Vapor-phase composition `y` (length = number of components).
    pub vapor_composition: Vec<f64>,
    /// Liquid-phase molar volume.
    pub liquid_volume: MolarVolume,
    /// Vapor-phase molar volume.
    pub vapor_volume: MolarVolume,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Whether the iteration converged to tolerance.
    pub converged: bool,
    /// Two-phase vs single-phase classification.
    pub phase_flag: PhaseFlag,
}

impl FlashResult {
    /// Average molar volume of the equilibrium mixture:
    /// `V = (1 − β)·vᴸ + β·vⱽ`.
    pub fn mixture_molar_volume(&self) -> MolarVolume {
        use uom::si::molar_volume::cubic_meter_per_mole;
        let v = (1.0 - self.vapor_fraction) * self.liquid_volume.value
            + self.vapor_fraction * self.vapor_volume.value;
        MolarVolume::new::<cubic_meter_per_mole>(v)
    }
}

/// Driver for flash calculations over a specific equation of state.
pub struct FlashCalculator<'a, E: EquationOfState + ?Sized> {
    eos: &'a E,
    db: Option<&'a dyn ComponentDatabase>,
    nc: usize,
}

impl<'a, E: EquationOfState + ?Sized> FlashCalculator<'a, E> {
    /// Build a calculator for `eos` without a component database (K-values start
    /// from `K_i = 1`).
    pub fn new(eos: &'a E) -> Self {
        Self {
            eos,
            db: None,
            nc: eos.num_components(),
        }
    }

    /// Build a calculator for `eos` with a component database used to seed Wilson
    /// K-values.
    pub fn with_db(eos: &'a E, db: &'a dyn ComponentDatabase) -> Self {
        Self {
            eos,
            db: Some(db),
            nc: eos.num_components(),
        }
    }

    /// Number of components the underlying EoS describes.
    pub fn num_components(&self) -> usize {
        self.nc
    }

    /// Borrow the underlying equation of state.
    #[allow(dead_code)]
    pub(crate) fn eos_ref(&self) -> &E {
        self.eos
    }

    /// Borrow the component database, if any.
    #[allow(dead_code)]
    pub(crate) fn db_opt(&self) -> Option<&dyn ComponentDatabase> {
        self.db
    }

    /// Number of components (crate-internal accessor).
    #[allow(dead_code)]
    pub(crate) fn comps(&self) -> usize {
        self.nc
    }

    /// K-values `K_i = φ_iᴸ/φ_iⱽ` evaluated at `(T, P)` for liquid composition `x`
    /// and vapor composition `y`.
    pub fn k_values(
        &self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
        y: &[f64],
    ) -> Result<Vec<f64>, ThermoError> {
        let vl = phase_volume(self.eos, t, p, x, Phase::Liquid)?;
        let vv = phase_volume(self.eos, t, p, y, Phase::Vapor)?;
        let n = x.len();
        let mut k = alloc::vec![0.0_f64; n];
        for i in 0..n {
            let phil = (self.eos.ln_fugacity_coefficient(t, vl, x, i)?).exp();
            let phiv = (self.eos.ln_fugacity_coefficient(t, vv, y, i)?).exp();
            k[i] = if phiv > 0.0 { phil / phiv } else { 1.0 };
        }
        Ok(k)
    }

    /// Isothermal/isobaric flash at `(T, P, z)`.
    pub fn flash_pt(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
    ) -> Result<FlashResult, FlashError> {
        flash_pt_impl(self.eos, self.db, self.nc, t, p, z, PT_MAX_ITER, PT_TOL)
    }

    /// PH flash: specified molar enthalpy `h`, pressure `P`.
    pub fn flash_ph(
        &self,
        h: tpt_thermo_core::quantities::MolarEnergy,
        p: Pressure,
        z: &[f64],
    ) -> Result<FlashResult, FlashError> {
        crate::variants::flash_ph_impl(self, h, p, z)
    }

    /// TV flash: specified molar volume `v`, temperature `T`.
    pub fn flash_tv(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<FlashResult, FlashError> {
        crate::variants::flash_tv_impl(self, t, v, z)
    }

    /// TS flash: specified molar entropy `s`, temperature `T`.
    pub fn flash_ts(
        &self,
        t: Temperature,
        s: tpt_thermo_core::quantities::MolarEntropy,
        z: &[f64],
    ) -> Result<FlashResult, FlashError> {
        crate::variants::flash_ts_impl(self, t, s, z)
    }

    /// PU flash: specified molar internal energy `u`, pressure `P`.
    pub fn flash_pu(
        &self,
        u: tpt_thermo_core::quantities::MolarEnergy,
        p: Pressure,
        z: &[f64],
    ) -> Result<FlashResult, FlashError> {
        crate::variants::flash_pu_impl(self, u, p, z)
    }

    /// PV flash: specified molar volume `v`, pressure `P`.
    pub fn flash_pv(
        &self,
        p: Pressure,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<FlashResult, FlashError> {
        crate::variants::flash_pv_impl(self, p, v, z)
    }
}

fn relative_change(a: &[f64], b: &[f64]) -> f64 {
    let mut m = 0.0_f64;
    for i in 0..a.len().min(b.len()) {
        let denom = a[i].abs().max(1e-12);
        m = m.max((a[i] - b[i]).abs() / denom);
    }
    m
}

/// Core PT flash implementation, shared by [`FlashCalculator::flash_pt`] and the
/// convenience free function.
pub fn flash_pt_impl<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn ComponentDatabase>,
    nc: usize,
    t: Temperature,
    p: Pressure,
    z: &[f64],
    max_iter: usize,
    tol: f64,
) -> Result<FlashResult, FlashError> {
    if z.len() != nc || z.is_empty() {
        return Err(FlashError::InvalidFeed);
    }
    let sum_z: f64 = z.iter().sum();
    if (sum_z - 1.0).abs() > 1e-6 {
        return Err(FlashError::InvalidFeed);
    }

    let mut k = match db {
        Some(d) => wilson_k_values(t.value, p, d, z)
            .map_err(|_| FlashError::NotConverged(ConvergenceStatus::NotConverged))?,
        None => alloc::vec![1.0_f64; nc],
    };

    let mut mem = AccelerationMemory::new(nc);
    let mut iterations = 0;

    for it in 0..max_iter {
        iterations = it + 1;
        let rr = rachford_rice(&k, z).map_err(FlashError::Thermo)?;
        let k_new = {
            let res = eos_k_values(eos, t, p, &rr.x, &rr.y);
            match res {
                Ok(v) => v,
                Err(_) => {
                    // Non-physical K update; stop and return what we have.
                    return build_result(eos, t, p, z, &rr, false, iterations);
                }
            }
        };
        if relative_change(&k, &k_new) < tol {
            k = k_new;
            let rr_final = rachford_rice(&k, z).map_err(FlashError::Thermo)?;
            return build_result(eos, t, p, z, &rr_final, true, iterations);
        }
        let (k_next, mem_next) = gdem_step(&k, &k_new, mem);
        k = k_next;
        mem = mem_next;
    }

    // Did not meet tolerance within the budget.
    let rr = rachford_rice(&k, z).map_err(FlashError::Thermo)?;
    build_result(eos, t, p, z, &rr, false, iterations)
}

/// Compute K-values for an arbitrary EoS (mirrors [`FlashCalculator::k_values`]).
fn eos_k_values<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    p: Pressure,
    x: &[f64],
    y: &[f64],
) -> Result<Vec<f64>, ThermoError> {
    let vl = phase_volume(eos, t, p, x, Phase::Liquid)?;
    let vv = phase_volume(eos, t, p, y, Phase::Vapor)?;
    let n = x.len();
    let mut k = alloc::vec![0.0_f64; n];
    for i in 0..n {
        let phil = (eos.ln_fugacity_coefficient(t, vl, x, i)?).exp();
        let phiv = (eos.ln_fugacity_coefficient(t, vv, y, i)?).exp();
        k[i] = if phiv > 0.0 { phil / phiv } else { 1.0 };
    }
    Ok(k)
}

/// Assemble a [`FlashResult`] from a converged Rachford–Rice split, recovering the
/// equilibrium phase volumes (single-phase collapses to one volume).
fn build_result<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    p: Pressure,
    z: &[f64],
    rr: &RachfordRiceResult,
    converged: bool,
    iterations: usize,
) -> Result<FlashResult, FlashError> {
    let beta = rr.beta.clamp(0.0, 1.0);
    let single_phase = beta <= 1e-12 || beta >= 1.0 - 1e-12;

    let (x, y) = if single_phase {
        // The whole feed is the single phase; the "other" phase is identical.
        (z.to_vec(), z.to_vec())
    } else {
        (rr.x.clone(), rr.y.clone())
    };

    let lv = match phase_volume(eos, t, p, &x, Phase::Liquid) {
        Ok(v) => v,
        Err(_) => {
            return Err(FlashError::Thermo(ThermoError::Numerical(
                ConvergenceStatus::NumericalIssue(NumericalIssueReason::NonPhysical),
            )))
        }
    };
    // For a single phase the vapor root may not exist; reuse the liquid volume.
    let vv = if single_phase {
        lv
    } else {
        phase_volume(eos, t, p, &y, Phase::Vapor).unwrap_or(lv)
    };

    Ok(FlashResult {
        vapor_fraction: beta,
        liquid_composition: x,
        vapor_composition: y,
        liquid_volume: lv,
        vapor_volume: vv,
        iterations,
        converged,
        phase_flag: if single_phase {
            PhaseFlag::SinglePhase
        } else {
            PhaseFlag::TwoPhase
        },
    })
}

/// Convenience free function: PT flash with an optional component database.
pub fn flash_pt<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn ComponentDatabase>,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    let nc = eos.num_components();
    flash_pt_impl(eos, db, nc, t, p, z, PT_MAX_ITER, PT_TOL)
}
