//! Mixing-rule and coupling traits.
//!
//! * [`MixingRule`] combines pure-component parameters into mixture ones.
//! * [`ExcessGibbsModel`] and [`StabilityTest`] are **forward-declared** here so
//!   the cubic crate (Phase 4) and activity crate (Phase 5) can reference them
//!   as bounds, and the flash crate (Phase 7) / phase crate (Phase 8) can
//!   implement them — without creating a cyclic dependency between those crates.

use crate::error::ThermoError;
use crate::quantities::{MolarEnergy, Pressure, Temperature};
use alloc::vec::Vec;

/// Combines pure-component parameters into mixture parameters (e.g. van der
/// Waals one-fluid). `combine` returns the mixture `a` and `b` (energy and
/// co-volume) given pure values and the composition.
pub trait MixingRule {
    /// Combined attractive (`a`) and co-volume (`b`) parameters at `(T, z)`.
    fn combine(
        &self,
        t: Temperature,
        a_pure: &[f64],
        b_pure: &[f64],
        z: &[f64],
    ) -> Result<(f64, f64), ThermoError>;
}

/// An excess Gibbs free energy model `g^E(T, P, x)`.
///
/// Forward-declared in the core so cubic EoS mixing rules (Phase 4) can be
/// generic over it; implemented by the activity crate (Phase 5).
pub trait ExcessGibbsModel: Send + Sync {
    /// Number of components.
    fn num_components(&self) -> usize;

    /// Reduced excess Gibbs energy `g^E / (R T)` at `(T, P, x)`.
    fn reduced_excess_gibbs(
        &self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
    ) -> Result<f64, ThermoError>;

    /// Natural log of the activity coefficient of component `i` at `(T, P, x)`.
    fn ln_gamma(
        &self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError>;
}

/// Outcome of a phase-stability test.
#[derive(Debug, Clone, PartialEq)]
pub struct StabilityResult {
    /// Whether the phase is stable at the tested conditions.
    pub stable: bool,
    /// Tangent-plane-distance minimiser compositions (one per trial), when
    /// computed.
    pub trial_compositions: Vec<Vec<f64>>,
    /// True if any trial indicated an incipient second phase.
    pub found_second_phase: bool,
}

/// A phase-stability (tangent-plane-distance) test.
///
/// Forward-declared in the core so the flash crate (Phase 7) can require it as
/// a bound; implemented by the phase crate (Phase 8).
pub trait StabilityTest: Send + Sync {
    /// Test whether `composition` at `(T, P)` is stable.
    fn test(
        &self,
        t: Temperature,
        p: Pressure,
        composition: &[f64],
    ) -> Result<StabilityResult, ThermoError>;

    /// Excess molar Gibbs energy contribution used by the TPD (J·mol⁻¹).
    fn excess_gibbs(
        &self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
    ) -> Result<MolarEnergy, ThermoError>;
}
