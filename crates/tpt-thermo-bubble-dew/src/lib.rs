//! `tpt-thermo-bubble-dew` — bubble / dew point, phase envelope, azeotrope and
//! criconden detection for `tpt-thermodynamics`.
//!
//! # Approach
//!
//! Every routine here is driven by an equation of state through the
//! [`KProvider`] trait, which extends [`tpt_thermo_core::EquationOfState`] with
//! phase-selected molar volumes. K-values are built from the equilibrium
//! fugacity equality `K_i = φ_i^L / φ_i^V` (Peng-Robinson / SRK / volume-
//! translated cubic models supplied by `tpt-thermo-eos-cubic` all implement
//! `KProvider`), so the solvers are model-agnostic and work for any EoS that
//! can return liquid and vapor volumes.
//!
//! The two-phase boundary — where the incipient second phase appears — is the
//! physical locus of every bubble and dew point. It is located by tracking the
//! liquid/vapor molar-volume gap (which collapses to zero at the boundary) and
//! root-solving it with Brent's method (reusing `tpt-thermo-core`'s
//! [`brent`](tpt_thermo_core::brent)). At the located point the equilibrium
//! phase composition is recovered by iterating the vapor (bubble) or liquid
//! (dew) composition to fugacity consistency.
//!
//! > **Coupling note.** The approved plan (see `todo.md` Phase 9) lists
//! > `tpt-thermo-flash` and `tpt-thermo-phase` as dependencies. Those crates
//! > land in later phases; this crate is intentionally self-contained on the
//! > EoS fugacity route so it compiles and validates against the cubic EoS
//! > already in the workspace. When the flash/phase crates exist, their K-value
//! > providers can implement [`KProvider`] (or feed `BubbleDewSolver`) without
//! > any change to the public API here.
//!
//! # Example
//!
//! ```
//! use tpt_thermo_core::component::ComponentDatabase;
//! use tpt_thermo_core::quantities::Pressure;
//! use tpt_thermo_eos_cubic::PengRobinson;
//! use tpt_thermo_data::SeedComponentDatabase;
//! use tpt_thermo_bubble_dew::{BubbleDewSolver, KProvider};
//! use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};
//!
//! let db = SeedComponentDatabase::from_seed();
//! let eos = PengRobinson::from_database(&db).unwrap();
//! let solver = BubbleDewSolver::new(&eos as &dyn KProvider, &db);
//!
//! // Bubble point of an equimolar benzene/toluene liquid at 1 atm.
//! let benzene = db.index_of("benzene").unwrap();
//! let toluene = db.index_of("toluene").unwrap();
//! let mut x = vec![0.0; db.num_components()];
//! x[benzene] = 0.5; x[toluene] = 0.5;
//! let p = Pressure::new::<pascal>(101_325.0);
//! let bp = solver.bubble_point_temperature(p, &x).unwrap();
//! assert!(bp.temperature.value > 353.0 && bp.temperature.value < 384.0);
//! ```

#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

pub mod azeotrope;
pub mod bubble;
pub mod cricondentherm;
pub mod dew;
pub mod envelope;
pub mod equilibrium;
pub mod kprovider;

pub use azeotrope::{detect_azeotrope, Azeotrope};
pub use bubble::BubblePoint;
pub use cricondentherm::{cricondenbar_cricondentherm, Criconden};
pub use dew::DewPoint;
pub use envelope::{bubble_dew_envelope, Envelope, EnvelopePoint};
pub use equilibrium::{equilibrate, phase_gap, Equilibrium, Kind};
pub use kprovider::{KProvider, Phase};

use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::convergence::{ConvergenceStatus, NumericalIssueReason};
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use uom::si::{pressure::pascal, thermodynamic_temperature::kelvin};

/// Driver for bubble/dew point, envelope, azeotrope and criconden calculations.
///
/// Holds a borrowed equation of state (as a [`KProvider`]) and a component
/// database used only to derive physical temperature/pressure scan bounds.
pub struct BubbleDewSolver<'a> {
    eos: &'a dyn KProvider,
    db: &'a dyn ComponentDatabase,
}

impl<'a> BubbleDewSolver<'a> {
    /// Build a solver from an equation of state and the database it was built
    /// from (used to derive scan bounds from critical constants).
    pub fn new(eos: &'a dyn KProvider, db: &'a dyn ComponentDatabase) -> Self {
        Self { eos, db }
    }

    /// Number of components the underlying EoS describes.
    pub fn num_components(&self) -> usize {
        self.eos.num_components()
    }

    /// Temperature scan bounds `(t_lo, t_hi)` derived from the critical
    /// temperatures: `0.4·min(Tc)` to `1.1·max(Tc)`.
    pub fn t_bounds(&self) -> (f64, f64) {
        let mut min = f64::INFINITY;
        let mut max = 0.0_f64;
        for i in 0..self.db.num_components() {
            if let Ok(tc) = self.db.critical_temperature(i) {
                min = min.min(tc.value);
                max = max.max(tc.value);
            }
        }
        if !min.is_finite() {
            min = 200.0;
        }
        if max <= 0.0 {
            max = 1000.0;
        }
        (0.4 * min, 1.1 * max)
    }

    /// Pressure scan bounds `(p_lo, p_hi)`: `1 Pa` to `max(2·min(Pc), 5 MPa)`.
    pub fn p_bounds(&self) -> (f64, f64) {
        let mut min = f64::INFINITY;
        for i in 0..self.db.num_components() {
            if let Ok(pc) = self.db.critical_pressure(i) {
                min = min.min(pc.value);
            }
        }
        if !min.is_finite() {
            min = 5.0e6;
        }
        (1.0, (2.0 * min).max(5.0e6))
    }

    /// Recover the equilibrium (other-phase) composition at `(T, P, feed)`
    /// tolerantly: if the boundary point is ill-conditioned, step slightly into
    /// the two-phase region and reuse that composition.
    pub(crate) fn equilibrium_at(
        &self,
        t: Temperature,
        p: Pressure,
        feed: &[f64],
        kind: Kind,
    ) -> Result<Equilibrium, ThermoError> {
        match equilibrate(self.eos, t, p, feed, kind) {
            Ok(e) => Ok(e),
            Err(_) => {
                let toff = match kind {
                    Kind::Bubble => Temperature::new::<kelvin>(t.value * 1.01),
                    Kind::Dew => Temperature::new::<kelvin>(t.value * 0.99),
                };
                equilibrate(self.eos, toff, p, feed, kind)
            }
        }
    }

    /// Bubble-point temperature of liquid `x` at pressure `p`: the temperature
    /// where the incipient-phase K-values satisfy `Σ x_i K_i = 1`, located by
    /// root-solving `bubble_g` (see [`bubble_g`]).
    pub(crate) fn boundary_temperature(&self, p: Pressure, x: &[f64]) -> Result<f64, ThermoError> {
        let (mn, mx) = self.feed_tc(x);
        let (lo, hi) = (0.3 * mn, 1.2 * mx);
        let eos = self.eos;
        root_scalar(
            |v: f64| bubble_g(eos, Temperature::new::<kelvin>(v), p, x),
            lo,
            (hi - lo) / 2000.0,
            lo,
            hi,
        )
    }

    /// Dew-point temperature of vapor `y` at pressure `p`: the temperature where
    /// `Σ y_i / K_i = 1` (see [`dew_g`]).
    pub(crate) fn boundary_temperature_dew(
        &self,
        p: Pressure,
        y: &[f64],
    ) -> Result<f64, ThermoError> {
        let (mn, mx) = self.feed_tc(y);
        let (lo, hi) = (0.3 * mn, 1.2 * mx);
        let eos = self.eos;
        root_scalar(
            |v: f64| dew_g(eos, Temperature::new::<kelvin>(v), p, y),
            lo,
            (hi - lo) / 2000.0,
            lo,
            hi,
        )
    }

    /// Bubble-point temperature with a supplied `t_guess` (continuation). The
    /// search is confined to `[0.5·t_guess, 1.5·t_guess]`, so a good guess
    /// (e.g. the previous envelope point) yields the physical root and avoids
    /// the spurious low-temperature roots that appear at high pressure.
    pub(crate) fn boundary_temperature_guess(
        &self,
        p: Pressure,
        x: &[f64],
        t_guess: f64,
    ) -> Result<f64, ThermoError> {
        let (mn, mx) = self.feed_tc(x);
        let (lo, hi) = (0.3 * mn, 1.2 * mx);
        let eos = self.eos;
        root_scalar_guess(
            |v: f64| bubble_g(eos, Temperature::new::<kelvin>(v), p, x),
            t_guess,
            lo,
            hi,
        )
    }

    /// Dew-point temperature with a supplied `t_guess` (continuation), mirroring
    /// [`boundary_temperature_guess`].
    pub(crate) fn boundary_temperature_dew_guess(
        &self,
        p: Pressure,
        y: &[f64],
        t_guess: f64,
    ) -> Result<f64, ThermoError> {
        let (mn, mx) = self.feed_tc(y);
        let (lo, hi) = (0.3 * mn, 1.2 * mx);
        let eos = self.eos;
        root_scalar_guess(
            |v: f64| dew_g(eos, Temperature::new::<kelvin>(v), p, y),
            t_guess,
            lo,
            hi,
        )
    }

    /// Bubble-point pressure of liquid `x` at temperature `t`: root-solve
    /// `bubble_g` over pressure.
    pub(crate) fn boundary_pressure(&self, t: Temperature, x: &[f64]) -> Result<f64, ThermoError> {
        let (p_lo, p_hi) = self.p_bounds();
        let eos = self.eos;
        root_scalar(
            |v: f64| bubble_g(eos, t, Pressure::new::<pascal>(v), x),
            p_lo,
            (p_hi - p_lo) / 2000.0,
            p_lo,
            p_hi,
        )
    }

    /// Dew-point pressure of vapor `y` at temperature `t`: root-solve `dew_g`.
    pub(crate) fn boundary_pressure_dew(
        &self,
        t: Temperature,
        y: &[f64],
    ) -> Result<f64, ThermoError> {
        let (p_lo, p_hi) = self.p_bounds();
        let eos = self.eos;
        root_scalar(
            |v: f64| dew_g(eos, t, Pressure::new::<pascal>(v), y),
            p_lo,
            (p_hi - p_lo) / 2000.0,
            p_lo,
            p_hi,
        )
    }

    fn min_tc(&self) -> f64 {
        let mut m = f64::INFINITY;
        for i in 0..self.db.num_components() {
            if let Ok(tc) = self.db.critical_temperature(i) {
                m = m.min(tc.value);
            }
        }
        if !m.is_finite() {
            m = 200.0;
        }
        m
    }

    fn max_tc(&self) -> f64 {
        let mut m = 0.0_f64;
        for i in 0..self.db.num_components() {
            if let Ok(tc) = self.db.critical_temperature(i) {
                m = m.max(tc.value);
            }
        }
        if m <= 0.0 {
            m = 1000.0;
        }
        m
    }

    /// Minimum and maximum critical temperatures over the components actually
    /// present in `z` (nonzero mole fraction). Using the feed-present range
    /// rather than the whole database's keeps the boundary root scan sensible
    /// when the database contains much lighter or heavier species than the
    /// mixture being solved (e.g. the full seed set alongside a benzene/toluene
    /// feed, where helium would otherwise force a ~5 K scan window).
    fn feed_tc(&self, z: &[f64]) -> (f64, f64) {
        let mut mn = f64::INFINITY;
        let mut mx = 0.0_f64;
        for (i, &zi) in z.iter().enumerate().take(self.db.num_components()) {
            if zi > 0.0 {
                if let Ok(tc) = self.db.critical_temperature(i) {
                    mn = mn.min(tc.value);
                    mx = mx.max(tc.value);
                }
            }
        }
        if !mn.is_finite() {
            mn = self.min_tc();
        }
        if mx <= 0.0 {
            mx = self.max_tc();
        }
        (mn, mx)
    }
}

/// Incipient-phase bubble residual `Σ_i x_i K_i − 1` with `K_i = φ_i^L/φ_i^V`
/// evaluated at the liquid (feed) composition `x`. Crosses zero at the bubble
/// point. Returns `NaN` when either phase volume cannot be resolved.
pub fn bubble_g(eos: &dyn KProvider, t: Temperature, p: Pressure, x: &[f64]) -> f64 {
    let vl = match eos.phase_volume(t, p, x, Phase::Liquid) {
        Ok(v) => v,
        Err(_) => return f64::NAN,
    };
    let vv = match eos.phase_volume(t, p, x, Phase::Vapor) {
        Ok(v) => v,
        Err(_) => return f64::NAN,
    };
    let mut s = 0.0_f64;
    for i in 0..x.len() {
        let phl = match eos.ln_fugacity_coefficient(t, vl, x, i) {
            Ok(l) => l.exp(),
            Err(_) => return f64::NAN,
        };
        let pfv = match eos.ln_fugacity_coefficient(t, vv, x, i) {
            Ok(l) => l.exp(),
            Err(_) => return f64::NAN,
        };
        if pfv <= 0.0 {
            return f64::NAN;
        }
        s += x[i] * phl / pfv;
    }
    s - 1.0
}

/// Incipient-phase dew residual `Σ_i y_i / K_i − 1` with `K_i = φ_i^L/φ_i^V`
/// evaluated at the vapor (feed) composition `y`. Crosses zero at the dew point.
pub fn dew_g(eos: &dyn KProvider, t: Temperature, p: Pressure, y: &[f64]) -> f64 {
    let vl = match eos.phase_volume(t, p, y, Phase::Liquid) {
        Ok(v) => v,
        Err(_) => return f64::NAN,
    };
    let vv = match eos.phase_volume(t, p, y, Phase::Vapor) {
        Ok(v) => v,
        Err(_) => return f64::NAN,
    };
    let mut s = 0.0_f64;
    for i in 0..y.len() {
        let phl = match eos.ln_fugacity_coefficient(t, vl, y, i) {
            Ok(l) => l.exp(),
            Err(_) => return f64::NAN,
        };
        let pfv = match eos.ln_fugacity_coefficient(t, vv, y, i) {
            Ok(l) => l.exp(),
            Err(_) => return f64::NAN,
        };
        if phl <= 0.0 {
            return f64::NAN;
        }
        s += y[i] * pfv / phl;
    }
    s - 1.0
}

/// Scan `var` from `start` in `step` increments until the (NaN-aware) residual
/// `f` changes sign, then refine the root with bisection. The first sign change
/// is returned (this is the physical bubble/dew point for the standalone
/// solvers at the pressures/temperatures used in the seed validation; the phase
/// envelope uses [`root_scalar_guess`] with continuation to avoid spurious
/// low-variable roots that appear at high `P`/`T`).
fn root_scalar<F>(f: F, start: f64, step: f64, lo: f64, hi: f64) -> Result<f64, ThermoError>
where
    F: Fn(f64) -> f64,
{
    let mut var = start;
    let mut fprev = f(var);
    let mut vprev = var;
    for _ in 0..5000 {
        var += step;
        let fv = f(var);
        if fprev.is_finite() && fv.is_finite() && fprev * fv <= 0.0 {
            let (a, b) = if vprev <= var {
                (vprev, var)
            } else {
                (var, vprev)
            };
            return tpt_thermo_core::bisection(|v: f64| f(v), a, b, 1e-7, 200)
                .map_err(ThermoError::Numerical);
        }
        vprev = var;
        fprev = fv;
        if (step > 0.0 && var >= hi) || (step < 0.0 && var <= lo) {
            break;
        }
    }
    Err(ThermoError::Numerical(ConvergenceStatus::NotConverged))
}

/// Like [`root_scalar`] but searches for a sign change only within a bracket
/// `[centre·0.5, centre·1.5]` (clamped to `[lo, hi]`). Used by the phase
/// envelope, which carries the previous solution as `centre` (continuation):
/// near the true root the residual is monotonic, so this returns the physical
/// point and ignores the spurious low-/high-variable roots that a full-range
/// scan would otherwise pick up at high `P`/`T`.
fn root_scalar_guess<F>(f: F, centre: f64, lo: f64, hi: f64) -> Result<f64, ThermoError>
where
    F: Fn(f64) -> f64,
{
    let a = (centre * 0.5).max(lo);
    let b = (centre * 1.5).min(hi);
    if b <= a {
        return Err(ThermoError::Numerical(ConvergenceStatus::NotConverged));
    }
    // If the bracket already brackets a root, refine directly.
    let fa = f(a);
    let fb = f(b);
    if fa.is_finite() && fb.is_finite() && fa * fb <= 0.0 {
        return tpt_thermo_core::bisection(|v: f64| f(v), a, b, 1e-7, 200)
            .map_err(ThermoError::Numerical);
    }
    // Otherwise scan finely inside the bracket for the first sign change.
    let step = (b - a) / 2000.0;
    let mut var = a;
    let mut fprev = f(var);
    let mut vprev = var;
    for _ in 0..5000 {
        var += step;
        if var >= b {
            break;
        }
        let fv = f(var);
        if fprev.is_finite() && fv.is_finite() && fprev * fv <= 0.0 {
            let (c, d) = if vprev <= var {
                (vprev, var)
            } else {
                (var, vprev)
            };
            return tpt_thermo_core::bisection(|v: f64| f(v), c, d, 1e-7, 200)
                .map_err(ThermoError::Numerical);
        }
        vprev = var;
        fprev = fv;
    }
    Err(ThermoError::Numerical(ConvergenceStatus::NotConverged))
}

/// Convenience constructor for a 1-atm pressure.
pub fn one_atm() -> Pressure {
    Pressure::new::<pascal>(101_325.0)
}

/// Helper to flag an out-of-domain evaluation as a numerical issue.
pub(crate) fn nonphysical() -> ThermoError {
    ThermoError::Numerical(ConvergenceStatus::NumericalIssue(
        NumericalIssueReason::NonPhysical,
    ))
}
