//! Shared cubic-EoS engine for the van der Waals / Redlich-Kwong family.
//!
//! [`CubicEos`] implements the full [`EquationOfState`] surface
//! (`pressure`, fugacity, residual enthalpy/entropy, heat capacities, speed of
//! sound) for any combination of [`CubicModel`], alpha function, and mixing
//! rule. The [`PengRobinson`](crate::PengRobinson),
//! [`SoaveRedlichKwong`](crate::SoaveRedlichKwong), and
//! [`VolumeTranslated`](crate::VolumeTranslated) types wrap it with fixed
//! model/mixing choices.

use crate::alpha::AlphaFunction;
use crate::cubic_solver::{CubicModel, Phase, compressibility_roots, select_root};
use crate::mixing::CubicMixing;
use alloc::boxed::Box;
use alloc::vec;
use alloc::vec::Vec;
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::quantities::{
    molar_entropy, MolarEnergy, MolarEntropy, MolarHeatCapacity, MolarVolume, Pressure,
    Temperature, Velocity,
};
use tpt_thermo_core::{EquationOfState, R, ThermoError};
use uom::si::molar_energy::joule_per_mole;
use uom::si::molar_heat_capacity::joule_per_kelvin_mole;
use uom::si::molar_volume::cubic_meter_per_mole;
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::velocity::meter_per_second;

const SQRT2: f64 = core::f64::consts::SQRT_2;
/// Reference pressure (Pa) at which excess-Gibbs mixing rules are evaluated; the
/// excess Gibbs energy is effectively pressure-independent.
const P_REF: f64 = 1.0e5;

#[derive(Debug, Clone)]
struct PureParams {
    tc: f64,
    #[allow(dead_code)]
    pc: f64,
    omega: f64,
    bi: f64,
    ai0: f64,
}

impl PureParams {
    fn a(&self, alpha: &dyn AlphaFunction, t: f64) -> f64 {
        self.ai0 * alpha.alpha(t / self.tc, self.omega)
    }
}

/// A cubic equation of state (Peng-Robinson / SRK family) with a chosen alpha
/// function and mixing rule.
pub struct CubicEos {
    model: CubicModel,
    pure: Vec<PureParams>,
    alpha: Box<dyn AlphaFunction>,
    mixing: Box<dyn CubicMixing>,
    volume_translation: Vec<f64>,
    molar_masses: Vec<f64>,
}

impl core::fmt::Debug for CubicEos {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("CubicEos")
            .field("model", &self.model)
            .field("n", &self.pure.len())
            .field("volume_translation", &self.volume_translation)
            .finish()
    }
}

impl CubicEos {
    /// Build from a [`ComponentDatabase`], a model, an alpha function, and a
    /// mixing rule. Pure critical parameters come from the database; the mixing
    /// rule carries any binary interactions.
    pub fn from_database(
        model: CubicModel,
        db: &dyn ComponentDatabase,
        alpha: Box<dyn AlphaFunction>,
        mixing: Box<dyn CubicMixing>,
    ) -> Result<Self, ThermoError> {
        let n = db.num_components();
        let mut pure = Vec::with_capacity(n);
        let mut molar_masses = Vec::with_capacity(n);
        for i in 0..n {
            let tc = db.critical_temperature(i)?.value;
            let pc = db.critical_pressure(i)?.value;
            let omega = db.acentric_factor(i)?;
            let bi = model.b_factor() * R * tc / pc;
            let ai0 = model.a_factor() * R * R * tc * tc / pc;
            pure.push(PureParams { tc, pc, omega, bi, ai0 });
            molar_masses.push(db.molar_mass(i)?.value);
        }
        Ok(Self {
            model,
            pure,
            alpha,
            mixing,
            volume_translation: vec![0.0; n],
            molar_masses,
        })
    }

    /// Attach per-component Peneloux volume-translation parameters (m³/mol).
    pub fn with_volume_translation(mut self, c: Vec<f64>) -> Self {
        if c.len() == self.pure.len() {
            self.volume_translation = c;
        }
        self
    }

    /// The underlying cubic model.
    pub fn model(&self) -> CubicModel {
        self.model
    }

    /// Number of components.
    pub fn num_components(&self) -> usize {
        self.pure.len()
    }

    fn c_mix(&self, z: &[f64]) -> f64 {
        self.volume_translation
            .iter()
            .zip(z.iter())
            .map(|(ci, zi)| ci * zi)
            .sum()
    }

    fn a_i_vec(&self, t: f64) -> Vec<f64> {
        self.pure.iter().map(|p| p.a(self.alpha.as_ref(), t)).collect()
    }

    fn b_i_vec(&self) -> Vec<f64> {
        self.pure.iter().map(|p| p.bi).collect()
    }

    /// Mixture `a`, `b`, and the component `a_i`, `b_i` vectors at `(T, z)`.
    pub(crate) fn mix_params(
        &self,
        t: f64,
        z: &[f64],
    ) -> (f64, f64, Vec<f64>, Vec<f64>) {
        let a = self.a_i_vec(t);
        let b = self.b_i_vec();
        let amix = self.mixing.a_mix(&a, &b, z, t, P_REF);
        let bmix = self.mixing.b_mix(&b, z);
        (amix, bmix, a, b)
    }

    fn log_term(&self, zc: f64, t: f64, bmix: f64, p: f64) -> (f64, f64) {
        // Returns `(ln-term, B)` for the fugacity / residual-departure formulas,
        // where `B = b_mix·P/(R T)` is the dimensionless co-volume. The
        // attractive coefficient is supplied separately by [`attractive_coeff`].
        let b = bmix * p / (R * t);
        let logterm = match self.model {
            CubicModel::PengRobinson => {
                ((zc + (1.0 + SQRT2) * b) / (zc + (1.0 - SQRT2) * b)).ln()
            }
            CubicModel::SoaveRedlichKwong => (1.0 + b / zc).ln(),
        };
        (logterm, b)
    }

    /// Attractive-term coefficient `1/(2√2·b_mix·R T)` (PR) or `1/(b_mix·R T)`
    /// (SRK), multiplied by the mixture `a_mix` (fugacity) or its temperature
    /// derivative (residual enthalpy / entropy).
    fn attractive_coeff(&self, bmix: f64, t: f64) -> f64 {
        match self.model {
            CubicModel::PengRobinson => 1.0 / (2.0 * SQRT2 * bmix * R * t),
            CubicModel::SoaveRedlichKwong => 1.0 / (bmix * R * t),
        }
    }

    /// Compressibility-factor roots at `(T, P, z)`.
    pub fn z_roots(&self, t: Temperature, p: Pressure, z: &[f64]) -> Vec<f64> {
        let (amix, bmix, _, _) = self.mix_params(t.value, z);
        let a = amix * p.value / (R * R * t.value * t.value);
        let b = bmix * p.value / (R * t.value);
        compressibility_roots(self.model, a, b)
    }

    /// Solve for the molar volume of `phase` at `(T, P, z)`, accounting for any
    /// volume translation (returns the physical molar volume).
    pub fn solve_phase(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Result<MolarVolume, ThermoError> {
        let roots = self.z_roots(t, p, z);
        let zroot = select_root(&roots, phase)
            .ok_or_else(|| ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NotConverged))?;
        if zroot <= 0.0 {
            return Err(ThermoError::Numerical(
                tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                    tpt_thermo_core::NumericalIssueReason::NonPhysical,
                ),
            ));
        }
        let v_eos = zroot * R * t.value / p.value;
        let v_phys = v_eos - self.c_mix(z);
        Ok(MolarVolume::new::<cubic_meter_per_mole>(v_phys))
    }

    /// `(∂P/∂v)_T` and `(∂²P/∂v²)_T` at `(T, v_eos, z)` for a pure component,
    /// used by critical/spinodal detection.
    pub(crate) fn pv_pure(&self, i: usize, t: f64, v_eos: f64) -> (f64, f64) {
        let p = self.pure[i].clone();
        let a = p.a(self.alpha.as_ref(), t);
        let b = p.bi;
        let (u, w) = (self.model.u(), self.model.w());
        let denom = v_eos * v_eos + u * b * v_eos + w * b * b;
        let pv = -R * t / (v_eos - b).powi(2) + a * (2.0 * v_eos + u * b) / denom.powi(2);
        let h = v_eos.abs().max(1e-8) * 1e-6;
        let pvp = self.pv_at(i, t, v_eos + h);
        let pvm = self.pv_at(i, t, v_eos - h);
        (pv, (pvp - pvm) / (2.0 * h))
    }

    fn pv_at(&self, i: usize, t: f64, v_eos: f64) -> f64 {
        let p = self.pure[i].clone();
        let a = p.a(self.alpha.as_ref(), t);
        let b = p.bi;
        let (u, w) = (self.model.u(), self.model.w());
        let denom = v_eos * v_eos + u * b * v_eos + w * b * b;
        -R * t / (v_eos - b).powi(2) + a * (2.0 * v_eos + u * b) / denom.powi(2)
    }

    /// Pure-component critical point `(T_c, P_c, v_c)` via 2D Newton on
    /// `(∂P/∂v)_T = 0` and `(∂²P/∂v²)_T = 0`, solved in scaled variables
    /// `x = v/b` and `y = T/T_c` so the residuals are O(1) and the search stays
    /// in the physical neighbourhood of the (analytic) critical point.
    pub fn critical_point_pure(
        &self,
        i: usize,
    ) -> Result<(Temperature, Pressure, MolarVolume), ThermoError> {
        let tc0 = self.pure[i].tc;
        let bi = self.pure[i].bi;
        let scale_pv = R * tc0 / (bi * bi);
        let scale_pvv = R * tc0 / (bi * bi * bi);
        let mut x = 4.0_f64; // v / b
        let mut y = 0.95_f64; // T / Tc
        let mut converged = false;
        for _ in 0..200 {
            let v = x * bi;
            let t = y * tc0;
            let (pv, pvv) = self.pv_pure(i, t, v);
            let f1 = pv / scale_pv;
            let f2 = pvv / scale_pvv;
            if f1.abs() < 1e-9 && f2.abs() < 1e-9 {
                converged = true;
                break;
            }
            let hx = 1e-4;
            let hy = 1e-4;
            let (pv_xp, pvv_xp) = self.pv_pure(i, y * tc0, (x + hx) * bi);
            let (pv_xm, pvv_xm) = self.pv_pure(i, y * tc0, (x - hx) * bi);
            let (pv_yp, pvv_yp) = self.pv_pure(i, (y + hy) * tc0, x * bi);
            let (pv_ym, pvv_ym) = self.pv_pure(i, (y - hy) * tc0, x * bi);
            let df1_dx = (pv_xp / scale_pv - pv_xm / scale_pv) / (2.0 * hx);
            let df2_dx = (pvv_xp / scale_pvv - pvv_xm / scale_pvv) / (2.0 * hx);
            let df1_dy = (pv_yp / scale_pv - pv_ym / scale_pv) / (2.0 * hy);
            let df2_dy = (pvv_yp / scale_pvv - pvv_ym / scale_pvv) / (2.0 * hy);
            let det = df1_dx * df2_dy - df1_dy * df2_dx;
            if det.abs() < 1e-30 {
                break;
            }
            let dx = (f1 * df2_dy - f2 * df1_dy) / det;
            let dy = (df1_dx * f2 - df2_dx * f1) / det;
            x -= dx;
            y -= dy;
            // Keep the iterate in the physical window around the critical point.
            x = x.clamp(1.05, 30.0);
            y = y.clamp(0.3, 1.05);
            if dx.abs() < 1e-8 && dy.abs() < 1e-8 {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NotConverged));
        }
        let v_eos = x * bi;
        let t = y * tc0;
        // Build a proper pure-component composition vector of length `n`.
        let mut z_pure = alloc::vec::Vec::with_capacity(self.pure.len());
        for k in 0..self.pure.len() {
            z_pure.push(if k == i { 1.0 } else { 0.0 });
        }
        let (amix, bmix, _, _) = self.mix_params(t, &z_pure);
        let denom = v_eos * v_eos + self.model.u() * bmix * v_eos + self.model.w() * bmix * bmix;
        let p = R * t / (v_eos - bmix) - amix / denom;
        let v_phys = v_eos - self.volume_translation[i];
        Ok((
            Temperature::new::<kelvin>(t),
            Pressure::new::<pascal>(p),
            MolarVolume::new::<cubic_meter_per_mole>(v_phys),
        ))
    }
}

impl EquationOfState for CubicEos {
    fn num_components(&self) -> usize {
        self.pure.len()
    }

    fn pressure(&self, t: Temperature, v: MolarVolume, z: &[f64]) -> Result<Pressure, ThermoError> {
        if z.len() != self.pure.len() {
            return Err(ThermoError::InvalidInput("composition length mismatch"));
        }
        let v_eos = v.value + self.c_mix(z);
        let (amix, bmix, _, _) = self.mix_params(t.value, z);
        let (u, w) = (self.model.u(), self.model.w());
        let denom = v_eos * v_eos + u * bmix * v_eos + w * bmix * bmix;
        if denom <= 0.0 || v_eos <= bmix {
            return Err(ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::NumericalIssueReason::NonPhysical,
            )));
        }
        let p = R * t.value / (v_eos - bmix) - amix / denom;
        Ok(Pressure::new::<pascal>(p))
    }

    fn ln_fugacity_coefficient(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        if i >= z.len() {
            return Err(ThermoError::IndexOutOfRange(i));
        }
        let v_eos = v.value + self.c_mix(z);
        let (amix, bmix, a, b) = self.mix_params(t.value, z);
        let (u, w) = (self.model.u(), self.model.w());
        let denom = v_eos * v_eos + u * bmix * v_eos + w * bmix * bmix;
        if denom <= 0.0 || v_eos <= bmix {
            return Err(ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::NumericalIssueReason::NonPhysical,
            )));
        }
        let p = R * t.value / (v_eos - bmix) - amix / denom;
        let zc = p * v_eos / (R * t.value);
        let bi = b[i];
        let sum = self.mixing.aij_sum(&a, &b, z, i, t.value, P_REF);
        let (logterm, bd) = self.log_term(zc, t.value, bmix, p);
        let coef = self.attractive_coeff(bmix, t.value);
        let term = 2.0 * sum / amix - bi / bmix;
        let lnphi = bi / bmix * (zc - 1.0) - (zc - bd).ln() - coef * amix * term * logterm;
        Ok(lnphi)
    }

    fn molar_enthalpy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEnergy, ThermoError> {
        let v_eos = v.value + self.c_mix(z);
        let (amix, bmix, _, _) = self.mix_params(t.value, z);
        let (u, w) = (self.model.u(), self.model.w());
        let denom = v_eos * v_eos + u * bmix * v_eos + w * bmix * bmix;
        if denom <= 0.0 || v_eos <= bmix {
            return Err(ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::NumericalIssueReason::NonPhysical,
            )));
        }
        let p = R * t.value / (v_eos - bmix) - amix / denom;
        let zc = p * v_eos / (R * t.value);
        // Numerical d(amix)/dT.
        let dt = t.value.abs().max(1.0) * 1e-3;
        let amix_p = self.mix_params(t.value + dt, z).0;
        let amix_m = self.mix_params(t.value - dt, z).0;
        let da_dt = (amix_p - amix_m) / (2.0 * dt);
        let (logterm, _bd) = self.log_term(zc, t.value, bmix, p);
        // H^R = R T (Z - 1) + (T da/dT - a) / (2√2 b) · ln-term  (volume `b`).
        let coef_h = match self.model {
            CubicModel::PengRobinson => 1.0 / (2.0 * SQRT2 * bmix),
            CubicModel::SoaveRedlichKwong => 1.0 / bmix,
        };
        let h_res = R * t.value * (zc - 1.0) + (t.value * da_dt - amix) * coef_h * logterm;
        // Peneloux volume-translation correction to internal energy, hence H.
        let h = h_res - p * self.c_mix(z);
        Ok(MolarEnergy::new::<joule_per_mole>(h))
    }

    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarEntropy, ThermoError> {
        let v_eos = v.value + self.c_mix(z);
        let (amix, bmix, _, _) = self.mix_params(t.value, z);
        let (u, w) = (self.model.u(), self.model.w());
        let denom = v_eos * v_eos + u * bmix * v_eos + w * bmix * bmix;
        if denom <= 0.0 || v_eos <= bmix {
            return Err(ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::NumericalIssueReason::NonPhysical,
            )));
        }
        let p = R * t.value / (v_eos - bmix) - amix / denom;
        let zc = p * v_eos / (R * t.value);
        let dt = t.value.abs().max(1.0) * 1e-3;
        let amix_p = self.mix_params(t.value + dt, z).0;
        let amix_m = self.mix_params(t.value - dt, z).0;
        let da_dt = (amix_p - amix_m) / (2.0 * dt);
        let (logterm, bd) = self.log_term(zc, t.value, bmix, p);
        let coef = self.attractive_coeff(bmix, t.value);
        // S^R/R = ln(Z - B) + (da/dT) · coef · T · logterm - ln Z.
        let s_res = (zc - bd).ln() - zc.ln() + da_dt * coef * logterm;
        // Ideal-gas (reference) entropy: mixing + -ln(P/P_ref), with zero reference Cp.
        let mixing: f64 = -R
            * z.iter()
                .filter(|&&zi| zi > 0.0)
                .map(|&zi| zi * zi.ln())
                .sum::<f64>();
        let p_ref = P_REF;
        let s_ig = mixing - R * (p / p_ref).ln();
        Ok(molar_entropy(s_res * R + s_ig))
    }

    fn molar_isobaric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        let cp = self.cp_cv(t, v, z)?.0;
        Ok(MolarHeatCapacity::new::<joule_per_kelvin_mole>(cp))
    }

    fn molar_isochoric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        let cv = self.cp_cv(t, v, z)?.1;
        Ok(MolarHeatCapacity::new::<joule_per_kelvin_mole>(cv))
    }

    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<Velocity, ThermoError> {
        let (cp, cv) = self.cp_cv(t, v, z)?;
        if cv <= 0.0 {
            return Err(ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::NumericalIssueReason::SingularJacobian,
            )));
        }
        let dp_dv = self.dp_dv(t, v, z)?;
        let m_mass: f64 = self
            .molar_masses
            .iter()
            .zip(z.iter())
            .map(|(mi, zi)| mi * zi)
            .sum();
        if m_mass <= 0.0 {
            return Err(ThermoError::InvalidInput("non-positive molar mass"));
        }
        let a2 = (cp / cv) * (-v.value * v.value * dp_dv) / m_mass;
        if a2 <= 0.0 {
            return Err(ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::NumericalIssueReason::NonPhysical,
            )));
        }
        Ok(Velocity::new::<meter_per_second>(a2.sqrt()))
    }
}

impl CubicEos {
    /// `(c_p, c_v)` (J·mol⁻¹·K⁻¹) at `(T, v, z)` via finite differences of the
    /// residual internal energy and the thermodynamic identity
    /// `c_p − c_v = −T (∂P/∂T)_v² / (∂P/∂v)_T`.
    fn cp_cv(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<(f64, f64), ThermoError> {
        let dt = t.value.abs().max(1.0) * 1e-3;
        let u = |tt: f64| -> Result<f64, ThermoError> {
            let tv = Temperature::new::<kelvin>(tt);
            let hv = self.molar_enthalpy(tv, v, z)?;
            let (amix, bmix, _, _) = self.mix_params(tt, z);
            let ve = v.value + self.c_mix(z);
            let p = R * tt / (ve - bmix)
                - amix
                    / (ve * ve + self.model.u() * bmix * ve + self.model.w() * bmix * bmix);
            let zc2 = p * ve / (R * tt);
            // Residual internal energy u = H − R T (Z − 1) (the +R T is the
            // ideal-gas reference that cancels in derivatives but keeps u finite).
            Ok(hv.value - (zc2 - 1.0) * R * tt)
        };
        let u_p = u(t.value + dt)?;
        let u_m = u(t.value - dt)?;
        let cv = (u_p - u_m) / (2.0 * dt);
        let dp_dv = self.dp_dv(t, v, z)?;
        let dp_dt = self.dp_dt(t, v, z)?;
        if dp_dv.abs() < 1e-30 {
            return Err(ThermoError::Numerical(tpt_thermo_core::ConvergenceStatus::NumericalIssue(
                tpt_thermo_core::NumericalIssueReason::SingularJacobian,
            )));
        }
        let cp = cv - t.value * dp_dt * dp_dt / dp_dv;
        Ok((cp, cv))
    }
}
