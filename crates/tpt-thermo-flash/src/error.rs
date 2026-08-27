//! Errors and convergence reporting for the flash solvers.

use tpt_thermo_core::convergence::ConvergenceStatus;
use tpt_thermo_core::error::ThermoError;

/// Errors that can arise during a flash calculation.
#[derive(Debug, Clone, PartialEq)]
pub enum FlashError {
    /// An underlying thermodynamic evaluation failed.
    Thermo(ThermoError),
    /// The iteration did not meet tolerance within its budget.
    NotConverged(ConvergenceStatus),
    /// The feed composition was invalid (empty or not normalised).
    InvalidFeed,
}

impl core::fmt::Display for FlashError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            FlashError::Thermo(e) => write!(f, "flash thermo error: {e}"),
            FlashError::NotConverged(s) => write!(f, "flash did not converge: {s:?}"),
            FlashError::InvalidFeed => write!(f, "flash: invalid feed composition"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for FlashError {}

impl From<ThermoError> for FlashError {
    fn from(e: ThermoError) -> Self {
        FlashError::Thermo(e)
    }
}

impl From<ConvergenceStatus> for FlashError {
    fn from(s: ConvergenceStatus) -> Self {
        FlashError::NotConverged(s)
    }
}
