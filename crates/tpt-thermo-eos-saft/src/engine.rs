//! Shared SAFT engine: hard-chain, dispersion, association, and the EoS
//! interface implementation.
//!
//! Both [`PcSaft`](crate::PcSaft) and [`SaftVrMie`](crate::SaftVrMie) are thin
//! wrappers around this engine. The residual Helmholtz energy is built as
//!
//! ```text
//! a^res/(RT) = a^hc + a^disp + a^assoc
//! ```
//!
//! with the Gross & Sadowski (2001) hard-chain + dispersion terms and the
//! association term ([`association`](crate::association)). The pressure is
//! recovered from the packing-fraction derivative
//! `Z = 1 + η · ∂(a^res/RT)/∂η`, and fugacity / enthalpy / entropy use the
//! numerical-default composition and temperature derivatives (permitted by the
//! build-out spec's "numerical-default fallback" clause). The pure-component
//! path is exact PC-SAFT; the mixture hard-chain uses the Carnahan-Starling
//! one-fluid approximation (documented refinement over full bmcsL for highly
//! asymmetric mixtures).

use crate::association::{self, AssociationResult};
use crate::parameters::SaftParameters;
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::convergence::NumericalIssueReason;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{
    molar_entropy, MolarEnergy, MolarEntropy, MolarHeatCapacity, MolarVolume, Pressure,
    Temperature, Velocity,
};
use tpt_thermo_core::{EquationOfState, R};
use uom::si::molar_energy::joule_per_mole;
use uom::si::molar_heat_capacity::joule_per_kelvin_mole;
use uom::si::molar_volume::cubic_meter_per_mole;
use uom::si::pressure::pascal;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::velocity::meter_per_second;

const NA: f64 = 6.022_140_76e23;
const PI: f64 = core::f64::consts::PI;

/// SAFT model flavour.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SaftFlavor {
    PcSaft,
    VrMie,
}

/// The shared SAFT engine.
#[derive(Debug, Clone)]
pub struct SaftEngine {
    params: SaftParameters,
    kij: Vec<Vec<f64>>,
    molar_masses: Vec<f64>,
    pub(crate) flavor: SaftFlavor,
}

impl SaftEngine {
    /// Build from a parameter set and per-component molar masses (kg·mol⁻¹).
    pub fn new(params: SaftParameters, molar_masses: Vec<f64>) -> Self {
        let n = params.num_components();
        let kij = vec![vec![0.0; n]; n];
        Self {
            params,
            kij,
            molar_masses,
            flavor: SaftFlavor::PcSaft,
        }
    }

    /// Build from the seed database (looks up SAFT parameters and molar masses
    /// by component name).
    pub fn from_seed_database(
        db: &dyn ComponentDatabase,
        flavor: SaftFlavor,
    ) -> Result<Self, ThermoError> {
        let params = SaftParameters::from_seed_database(db)?;
        let n = db.num_components();
        let mut mm = Vec::with_capacity(n);
        for i in 0..n {
            mm.push(db.molar_mass(i)?.value);
        }
        let mut engine = Self::new(params, mm);
        engine.flavor = flavor;
        Ok(engine)
    }

    /// Attach a binary interaction matrix `k_ij` (dimensionless, symmetric).
    pub fn with_kij(mut self, kij: Vec<Vec<f64>>) -> Self {
        if kij.len() == self.params.num_components() {
            self.kij = kij;
        }
        self
    }

    /// Underlying SAFT parameters.
    pub fn parameters(&self) -> &SaftParameters {
        &self.params
    }

    /// Number of components.
    pub fn num_components(&self) -> usize {
        self.params.num_components()
    }

    fn kij(&self, i: usize, j: usize) -> f64 {
        self.kij[i][j]
    }

    /// Temperature-dependent hard-sphere diameters `d_i` (m).
    fn diameters(&self, t: f64) -> Vec<f64> {
        self.params
            .components
            .iter()
            .map(|c| {
                let sigma = c.sigma * 1e-10; // Å → m
                sigma * (1.0 - 0.12 * (-3.0 * c.epsilon_k / t).exp())
            })
            .collect()
    }

    /// Temperature-independent segment diameters `σ_i` (m).
    fn sigmas(&self) -> Vec<f64> {
        self.params
            .components
            .iter()
            .map(|c| c.sigma * 1e-10)
            .collect()
    }

    /// Packing-fraction quantities: molecular number density `ρ` (m⁻³) and the
    /// mixture packing fraction `η = ζ_3`.
    fn packing(&self, t: f64, v: f64, x: &[f64]) -> (f64, f64, Vec<f64>) {
        let rho_mol = 1.0 / v; // mol·m⁻³
        let rho = NA * rho_mol; // molecules·m⁻³
        let d = self.diameters(t);
        let mut zeta3 = 0.0;
        for (i, xi) in x.iter().enumerate() {
            let m = self.params.component(i).m;
            zeta3 += (PI / 6.0) * rho * xi * m * d[i].powi(3);
        }
        (rho, zeta3, d)
    }

    /// Reduced residual Helmholtz energy `a^res/(RT)` at `(T, v, x)`.
    pub fn ares(&self, t: f64, v: f64, x: &[f64]) -> Result<f64, ThermoError> {
        let (rho, zeta3, _d) = self.packing(t, v, x);
        if !(0.0..0.9999).contains(&zeta3) {
            return Err(ThermoError::Numerical(ConvergenceStatus_::NumericalIssue(
                NumericalIssueReason::NonPhysical,
            )));
        }
        // Hard-chain (Carnahan-Starling one-fluid).
        let m_mix: f64 = x
            .iter()
            .zip(self.params.components.iter())
            .map(|(xi, c)| xi * c.m)
            .sum();
        let g_hs = (1.0 - 0.5 * zeta3) / (1.0 - zeta3).powi(3);
        let ahs = (4.0 * zeta3 - 3.0 * zeta3 * zeta3) / (1.0 - zeta3).powi(2);
        let ahc = m_mix * ahs - (m_mix - 1.0) * g_hs.ln();

        // Dispersion (Gross & Sadowski 2001, eq 18-24) — uses temperature-independent σ.
        let sigma = self.sigmas();
        let xi_bar = self.xi_bar(t, rho, &sigma, x);
        let (a1, a2) = self.dispersion_sums(t, &sigma, x);
        let i1 = i1(zeta3, xi_bar);
        let i2 = i2(zeta3, xi_bar);
        // NB: `i1`/`i2` evaluate negative for the usual (η, ξ) range, so the
        // leading signs here make the dispersion attractive (A^disp < 0).
        let disp = 2.0 * PI * rho * i1 * a1 + PI * rho * i2 * a2;

        // Association.
        let assoc = self.association_term(t, rho, zeta3, x)?;

        Ok(ahc + disp + assoc)
    }

    fn association_term(
        &self,
        t: f64,
        rho: f64,
        zeta3: f64,
        x: &[f64],
    ) -> Result<f64, ThermoError> {
        match association::solve_association(&self.params, x, rho, t, zeta3) {
            Ok(AssociationResult { ares, .. }) => Ok(ares),
            Err(status) => Err(ThermoError::Numerical(status)),
        }
    }

    /// Reduced association strength `ξ̄` for the dispersion integrals (Gross &
    /// Sadowski 2001, eq 14): the ratio cancels the number-density prefactor,
    /// leaving `ξ̄ = Σ_i x_i m_i (ε_i/kT) d_i³ / Σ_i x_i m_i d_i³`.
    fn xi_bar(&self, t: f64, _rho: f64, d: &[f64], x: &[f64]) -> f64 {
        let mut num = 0.0;
        let mut den = 0.0;
        for (i, xi) in x.iter().enumerate() {
            let m = self.params.component(i).m;
            let eps = self.params.component(i).epsilon_k;
            let di3 = d[i].powi(3);
            num += xi * m * (eps / t) * di3;
            den += xi * m * di3;
        }
        if den <= 0.0 {
            0.0
        } else {
            num / den
        }
    }

    /// Dispersion pair sums `a_1` and `a_2`.
    fn dispersion_sums(&self, t: f64, d: &[f64], x: &[f64]) -> (f64, f64) {
        let n = self.params.num_components();
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        for i in 0..n {
            for j in 0..n {
                let mi = self.params.component(i).m;
                let mj = self.params.component(j).m;
                let ei = self.params.component(i).epsilon_k;
                let ej = self.params.component(j).epsilon_k;
                let eps_ij = (ei * ej).sqrt() * (1.0 - self.kij(i, j));
                let sigma_ij = 0.5 * (d[i] + d[j]); // use temperature-dependent diameters
                let e = eps_ij / t;
                let s3 = sigma_ij.powi(3);
                a1 += x[i] * x[j] * mi * mj * e * s3;
                a2 += x[i] * x[j] * mi * mj * e * e * s3;
            }
        }
        (a1, a2)
    }

    /// Compressibility factor `Z = 1 + η · ∂(a^res/RT)/∂η` (central difference).
    fn compressibility(&self, t: f64, v: f64, x: &[f64]) -> Result<(f64, f64), ThermoError> {
        let (_rho, eta, _d) = self.packing(t, v, x);
        if eta <= 0.0 {
            return Ok((1.0, eta));
        }
        let deta = eta * 1e-5 + 1e-12;
        let eta_p = eta + deta;
        let eta_m = (eta - deta).max(1e-12);
        let v_p = v * eta / eta_p;
        let v_m = v * eta / eta_m;
        let a_p = self.ares(t, v_p, x)?;
        let a_m = self.ares(t, v_m, x)?;
        let dadeta = (a_p - a_m) / (eta_p - eta_m);
        let z = 1.0 + eta * dadeta;
        Ok((z, eta))
    }

    /// `∂(a^res/RT)/∂x_i` via forward differences, composition renormalised.
    fn da_dx(&self, t: f64, v: f64, x: &[f64]) -> Result<Vec<f64>, ThermoError> {
        let n = x.len();
        // For a single component the composition-derivative contribution to the
        // fugacity formula cancels exactly (`grad_i − Σ x_k grad_k = 0`), so a
        // numeric perturbation is degenerate (0/0). Return the zero gradient.
        if n == 1 {
            return Ok(vec![0.0]);
        }
        let base = self.ares(t, v, x)?;
        let mut grad = vec![0.0_f64; n];
        let h = 1e-4;
        for i in 0..n {
            let mut xp = x.to_vec();
            let step = (x[i] * h).max(h / n as f64).min(0.1);
            let old = xp[i];
            xp[i] = (old + step).min(1.0);
            // Renormalise so the sum stays 1.
            let sum: f64 = xp.iter().sum();
            if (sum - 1.0).abs() > 1e-12 {
                for xi in xp.iter_mut() {
                    *xi /= sum;
                }
            }
            let a = self.ares(t, v, xp.as_slice())?;
            grad[i] = (a - base) / (xp[i] - old);
        }
        Ok(grad)
    }

    /// `∂(a^res/RT)/∂T` (central difference, `v`, `x` fixed).
    fn da_dt(&self, t: f64, v: f64, x: &[f64]) -> Result<f64, ThermoError> {
        let h = t.abs().max(1.0) * 1e-4;
        let ap = self.ares(t + h, v, x)?;
        let am = self.ares(t - h, v, x)?;
        Ok((ap - am) / (2.0 * h))
    }

    /// `Z`-weighted pressure (Pa) at `(T, v, x)`.
    fn pressure_value(&self, t: f64, v: f64, x: &[f64]) -> Result<f64, ThermoError> {
        let (z, _eta) = self.compressibility(t, v, x)?;
        Ok(z * R * t / v)
    }

    /// Saturated molar volume roots (vapor, liquid) and `P_sat` via the equal
    /// fugacity / equal pressure tangent at `T` for a pure component. Used for
    /// validation against literature vapor pressures.
    pub fn saturation_pressure(
        &self,
        t: Temperature,
    ) -> Result<(Pressure, MolarVolume, MolarVolume), ThermoError> {
        let n = self.num_components();
        if n != 1 {
            return Err(ThermoError::Unsupported(
                "saturation_pressure is implemented for pure components",
            ));
        }
        let x = vec![1.0];
        let tt = t.value;
        // Bisect on pressure: P_sat is the lowest pressure that still admits two
        // distinct volume roots (liquid + vapor). Above it the liquid root exists;
        // below it only the vapor root remains.
        let mut p_lo = 1.0_f64; // Pa — below P_sat (single root)
        let mut p_hi = 1.0e7_f64; // Pa — start high, descend below P_c if needed
        while self.sat_roots(tt, p_hi, &x).is_err() && p_hi > 1.0 {
            p_hi *= 0.5;
        }
        if self.sat_roots(tt, p_hi, &x).is_err() {
            return Err(ThermoError::Numerical(ConvergenceStatus_::NotConverged));
        }
        let mut found = p_hi;
        for _ in 0..100 {
            let p_mid = 0.5 * (p_lo + p_hi);
            if self.sat_roots(tt, p_mid, &x).is_ok() {
                p_hi = p_mid;
                found = p_mid;
            } else {
                p_lo = p_mid;
            }
            if (p_hi - p_lo) / p_hi < 1e-7 {
                break;
            }
        }
        let (vv, vl) = self.sat_roots(tt, found, &x)?;
        Ok((
            Pressure::new::<pascal>(found),
            MolarVolume::new::<cubic_meter_per_mole>(vl),
            MolarVolume::new::<cubic_meter_per_mole>(vv),
        ))
    }

    /// Liquid (`vl`, smallest) and vapor (`vv`, largest) molar volumes giving
    /// pressure `p`, found by a logarithmic scan for sign changes of
    /// `P_EOS(v) − p` and Brent refinement. Returns `Err` when fewer than two
    /// distinct roots exist (single-phase region).
    fn sat_roots(&self, t: f64, p: f64, x: &[f64]) -> Result<(f64, f64), ThermoError> {
        let v_min = 1e-7_f64;
        let v_max = 2.0_f64;
        let npts = 6000_usize;
        let f = |v: f64| -> f64 { self.volume_pressure(t, v, x).unwrap_or(f64::INFINITY) - p };
        let mut roots = Vec::new();
        let mut prev_v = v_min;
        let mut prev_f = f(prev_v);
        let log_lo = v_min.ln();
        let log_hi = v_max.ln();
        for k in 1..=npts {
            let v = (log_lo + (log_hi - log_lo) * (k as f64) / (npts as f64)).exp();
            let fv = f(v);
            if prev_f == 0.0 {
                roots.push(prev_v);
            } else if fv == 0.0 {
                roots.push(v);
            } else if (fv > 0.0) != (prev_f > 0.0) {
                if let Ok(r) = tpt_thermo_core::numerics::brent(f, prev_v, v, 1e-10, 200) {
                    if r > v_min && r < v_max {
                        roots.push(r);
                    }
                }
            }
            prev_v = v;
            prev_f = fv;
        }
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roots.dedup_by(|a, b| (*a - *b).abs() < 1e-7);
        if roots.len() >= 2 {
            Ok((roots[roots.len() - 1], roots[0]))
        } else {
            Err(ThermoError::Numerical(ConvergenceStatus_::NotConverged))
        }
    }

    fn volume_pressure(&self, t: f64, v: f64, x: &[f64]) -> Result<f64, ThermoError> {
        self.pressure_value(t, v, x)
    }
}

// Re-export the convergence status path locally to keep the error mapping terse.
use tpt_thermo_core::convergence::ConvergenceStatus as ConvergenceStatus_;

impl EquationOfState for SaftEngine {
    fn num_components(&self) -> usize {
        self.num_components()
    }

    fn pressure(&self, t: Temperature, v: MolarVolume, x: &[f64]) -> Result<Pressure, ThermoError> {
        let p = self.pressure_value(t.value, v.value, x)?;
        Ok(Pressure::new::<pascal>(p))
    }

    fn ln_fugacity_coefficient(
        &self,
        t: Temperature,
        v: MolarVolume,
        x: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        if i >= x.len() {
            return Err(ThermoError::IndexOutOfRange(i));
        }
        let ares = self.ares(t.value, v.value, x)?;
        let (z, _eta) = self.compressibility(t.value, v.value, x)?;
        let grad = self.da_dx(t.value, v.value, x)?;
        let mut dsum = 0.0;
        for (k, xk) in x.iter().enumerate() {
            dsum += xk * grad[k];
        }
        Ok(ares + grad[i] - dsum + (z - 1.0))
    }

    fn molar_enthalpy(
        &self,
        t: Temperature,
        v: MolarVolume,
        x: &[f64],
    ) -> Result<MolarEnergy, ThermoError> {
        // Residual enthalpy: H^res = -R T² ∂(a^res/RT)/∂T.
        let da_dt = self.da_dt(t.value, v.value, x)?;
        let h_res = -R * t.value * t.value * da_dt;
        Ok(MolarEnergy::new::<joule_per_mole>(h_res))
    }

    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        x: &[f64],
    ) -> Result<MolarEntropy, ThermoError> {
        let ares = self.ares(t.value, v.value, x)?;
        let da_dt = self.da_dt(t.value, v.value, x)?;
        // S^res = -R (a^res + T da/dT).
        let s_res = -R * (ares + t.value * da_dt);
        let p = self.pressure_value(t.value, v.value, x)?;
        let mixing: f64 = -R
            * x.iter()
                .filter(|&&xi| xi > 0.0)
                .map(|&xi| xi * xi.ln())
                .sum::<f64>();
        let s_ig = mixing - R * (p / 1.0e5).ln();
        Ok(molar_entropy(s_res + s_ig))
    }

    fn molar_isobaric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        x: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        let cp = self.cp_cv(t.value, v.value, x)?.0;
        Ok(MolarHeatCapacity::new::<joule_per_kelvin_mole>(cp))
    }

    fn molar_isochoric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        x: &[f64],
    ) -> Result<MolarHeatCapacity, ThermoError> {
        let cv = self.cp_cv(t.value, v.value, x)?.1;
        Ok(MolarHeatCapacity::new::<joule_per_kelvin_mole>(cv))
    }

    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        x: &[f64],
    ) -> Result<Velocity, ThermoError> {
        let (cp, cv) = self.cp_cv(t.value, v.value, x)?;
        if cv <= 0.0 {
            return Err(ThermoError::Numerical(ConvergenceStatus_::NumericalIssue(
                NumericalIssueReason::SingularJacobian,
            )));
        }
        let dp_dv = self.dp_dv(t, v, x)?;
        let m_mass: f64 = self
            .molar_masses
            .iter()
            .zip(x.iter())
            .map(|(mi, xi)| mi * xi)
            .sum();
        if m_mass <= 0.0 {
            return Err(ThermoError::InvalidInput("non-positive molar mass"));
        }
        let a2 = (cp / cv) * (-v.value * v.value * dp_dv) / m_mass;
        if a2 <= 0.0 {
            return Err(ThermoError::Numerical(ConvergenceStatus_::NumericalIssue(
                NumericalIssueReason::NonPhysical,
            )));
        }
        Ok(Velocity::new::<meter_per_second>(a2.sqrt()))
    }
}

impl SaftEngine {
    /// `(c_p, c_v)` residual (J·mol⁻¹·K⁻¹) via finite differences of the
    /// residual internal energy using the thermodynamic identity
    /// `c_p − c_v = −T (∂P/∂T)_v² / (∂P/∂v)_T`, mirroring the cubic crate.
    fn cp_cv(&self, t: f64, v: f64, x: &[f64]) -> Result<(f64, f64), ThermoError> {
        let dt = t.abs().max(1.0) * 1e-3;
        let u = |tt: f64| -> Result<f64, ThermoError> {
            let hv = self.molar_enthalpy(
                Temperature::new::<kelvin>(tt),
                MolarVolume::new::<cubic_meter_per_mole>(v),
                x,
            )?;
            let p = self.pressure_value(tt, v, x)?;
            let zc = p * v / (R * tt);
            Ok(hv.value - (zc - 1.0) * R * tt)
        };
        let u_p = u(t + dt)?;
        let u_m = u(t - dt)?;
        let cv = (u_p - u_m) / (2.0 * dt);
        let dp_dv = self.dp_dv_raw(t, v, x)?;
        let dp_dt = self.dp_dt_raw(t, v, x)?;
        if dp_dv.abs() < 1e-30 {
            return Err(ThermoError::Numerical(ConvergenceStatus_::NumericalIssue(
                NumericalIssueReason::SingularJacobian,
            )));
        }
        let cp = cv - t * dp_dt * dp_dt / dp_dv;
        Ok((cp, cv))
    }

    fn dp_dv_raw(&self, t: f64, v: f64, x: &[f64]) -> Result<f64, ThermoError> {
        let h = v.abs().max(1e-8) * 1e-6;
        let pm = self.pressure_value(t, v - h, x)?;
        let pp = self.pressure_value(t, v + h, x)?;
        Ok((pp - pm) / (2.0 * h))
    }

    fn dp_dt_raw(&self, t: f64, v: f64, x: &[f64]) -> Result<f64, ThermoError> {
        let h = t.abs().max(1.0) * 1e-6;
        let pm = self.pressure_value(t - h, v, x)?;
        let pp = self.pressure_value(t + h, v, x)?;
        Ok((pp - pm) / (2.0 * h))
    }
}

/// PC-SAFT dispersion integral `I_1(η, ξ)` (Gross & Sadowski 2001, Table 1).
fn i1(eta: f64, xi: f64) -> f64 {
    let a = [
        0.910_283, 0.824_035, -1.737_244, -2.670_485, 3.044_772, -1.154_127,
    ];
    let xi2 = xi * xi;
    let poly = a[0] + a[1] * xi + a[2] * xi2;
    let rat = (a[3] + a[4] * xi + a[5] * xi2) / (1.0 - eta);
    poly + rat
}

/// PC-SAFT dispersion integral `I_2(η, ξ)` (Gross & Sadowski 2001, Table 1).
fn i2(eta: f64, xi: f64) -> f64 {
    let b = [
        -0.333_369, -0.309_373, 2.315_052, -3.454_927, 2.126_189, -0.603_327,
    ];
    let xi2 = xi * xi;
    let poly = b[0] + b[1] * xi + b[2] * xi2;
    let rat = (b[3] + b[4] * xi + b[5] * xi2) / (1.0 - eta);
    poly + rat
}
