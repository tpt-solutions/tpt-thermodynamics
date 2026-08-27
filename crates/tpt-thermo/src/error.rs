//! Unified error type for the umbrella crate.
//!
//! Every constituent crate ultimately returns [`tpt_thermo_core::ThermoError`], so
//! the workspace-wide error is simply that type, re-exported here.

pub use tpt_thermo_core::ThermoError;
