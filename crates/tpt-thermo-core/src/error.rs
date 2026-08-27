//! Error type for the core crate.

use crate::convergence::ConvergenceStatus;

/// Errors raised by core thermodynamic operations.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermoError {
    /// A numerical solver failed; carries its [`ConvergenceStatus`].
    Numerical(ConvergenceStatus),
    /// An input failed a physical/domain sanity check.
    InvalidInput(&'static str),
    /// A requested operation is not implemented for this model.
    Unsupported(&'static str),
    /// A component/phase index was out of range.
    IndexOutOfRange(usize),
    /// A database lookup failed for the given component/property.
    Database(&'static str),
}

impl core::fmt::Display for ThermoError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ThermoError::Numerical(s) => write!(f, "numerical failure: {s:?}"),
            ThermoError::InvalidInput(m) => write!(f, "invalid input: {m}"),
            ThermoError::Unsupported(m) => write!(f, "unsupported: {m}"),
            ThermoError::IndexOutOfRange(i) => write!(f, "index out of range: {i}"),
            ThermoError::Database(m) => write!(f, "database error: {m}"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for ThermoError {}

/// Convenience constructor for a numerical failure.
#[allow(dead_code)]
pub(crate) fn numerical(status: ConvergenceStatus) -> ThermoError {
    ThermoError::Numerical(status)
}

#[allow(dead_code)]
pub(crate) fn invalid_input(msg: &'static str) -> ThermoError {
    ThermoError::InvalidInput(msg)
}
