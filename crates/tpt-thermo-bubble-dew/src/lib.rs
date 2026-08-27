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

pub use azeotrope::{Azeotrope, detect_azeotrope};
pub use bubble::BubblePoint;
pub use cricondentherm::{Criconden, cricondenbar_cricondentherm};
pub use dew::DewPoint;
pub use envelope::{Envelope, EnvelopePoint, bubble_dew_envelope};
pub use equilibrium::{Equilibrium, equilibrate, phase_gap, Kind};
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

    /// Locate the two-phase boundary temperature for a fixed pressure and feed
    /// composition, scanning in the direction `up` (`true` = rising temperature,
    /// used for bubble points; `false` = falling, used for dew points).
    pub(crate) fn boundary_temperature(
        &self,
        p: Pressure,
        z: &[f64],
        up: bool,
    ) -> Result<f64, ThermoError> {
        let (t_lo, t_hi) = self.t_bounds();
        let eos = self.eos;
        let is_tp = |v: f64| phase_gap(eos, Temperature::new::<kelvin>(v), p, z) > 0.0;
        let (a, b) = bracket_transition(is_tp, t_lo, t_hi, up)?;
        tpt_thermo_core::bisection(
            |v: f64| phase_gap(eos, Temperature::new::<kelvin>(v), p, z),
            a,
            b,
            1e-7,
            200,
        )
        .map_err(ThermoError::Numerical)
    }

    /// Locate the two-phase boundary pressure for a fixed temperature and feed
    /// composition, scanning in the direction `up` (`true` = rising pressure,
    /// used for dew points; `false` = falling, used for bubble points).
    pub(crate) fn boundary_pressure(
        &self,
        t: Temperature,
        z: &[f64],
        up: bool,
    ) -> Result<f64, ThermoError> {
        let (p_lo, p_hi) = self.p_bounds();
        let eos = self.eos;
        let is_tp = |v: f64| phase_gap(eos, t, Pressure::new::<pascal>(v), z) > 0.0;
        let (a, b) = bracket_transition(is_tp, p_lo, p_hi, up)?;
        tpt_thermo_core::bisection(
            |v: f64| phase_gap(eos, t, Pressure::new::<pascal>(v), z),
            a,
            b,
            1e-4,
            200,
        )
        .map_err(ThermoError::Numerical)
    }
}

/// Scan `var` from `start` (here `lo` or `hi`, chosen by `up`) in steps of
/// `(hi-lo)/2000` until `is_tp` flips, returning the two adjacent sample points
/// that bracket the two-phase boundary.
fn bracket_transition<F>(
    mut is_tp: F,
    lo: f64,
    hi: f64,
    up: bool,
) -> Result<(f64, f64), ThermoError>
where
    F: FnMut(f64) -> bool,
{
    let step = if up { (hi - lo) / 2000.0 } else { -(hi - lo) / 2000.0 };
    let mut var = if up { lo } else { hi };
    let mut prev = is_tp(var);
    let mut prev_var = var;
    for _ in 0..5000 {
        var += step;
        let tp = is_tp(var);
        if tp != prev {
            let (a, b) = if prev_var <= var {
                (prev_var, var)
            } else {
                (var, prev_var)
            };
            return Ok((a, b));
        }
        if (up && var >= hi) || (!up && var <= lo) {
            return Err(ThermoError::Numerical(ConvergenceStatus::NotConverged));
        }
        prev = tp;
        prev_var = var;
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
