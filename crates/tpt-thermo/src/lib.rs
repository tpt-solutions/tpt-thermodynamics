//! `tpt-thermo` — umbrella crate for the `tpt-thermodynamics` workspace.
//!
//! This crate re-exports [`core`](tpt_thermo_core) and
//! [`data`](tpt_thermo_data) unconditionally, and every other constituent crate
//! behind a flat, non-implied feature flag (per spec sec3): `cubic`, `activity`,
//! `saft`, `phase`, `bubble-dew`, `flash`, `transport`, `electrolyte`, `polymer`.
//!
//! Enabling a feature pulls in that crate's public API under the matching module
//! name (e.g. `tpt_thermo::flash::FlashCalculator` with `feature = "flash"`). A
//! convenience high-level API lives in [`api`].
//!
//! # Feature matrix
//!
//! | Feature      | Module re-export            | Crate                |
//! |--------------|----------------------------|----------------------|
//! | `cubic`      | `eos_cubic`                 | `tpt-thermo-eos-cubic`   |
//! | `activity`   | `eos_activity`              | `tpt-thermo-eos-activity`|
//! | `saft`       | `eos_saft`                  | `tpt-thermo-eos-saft`    |
//! | `phase`      | `phase`                    | `tpt-thermo-phase`       |
//! | `bubble-dew` | `bubble_dew`               | `tpt-thermo-bubble-dew`  |
//! | `flash`      | `flash`                    | `tpt-thermo-flash`       |
//! | `transport`  | `transport`                | `tpt-thermo-transport`   |
//! | `electrolyte`| `electrolyte`              | `tpt-thermo-electrolyte` |
//! | `polymer`    | `polymer`                  | `tpt-thermo-polymer`     |
//!
//! `cargo build --no-default-features` yields only `core` + `data`; each feature
//! builds standalone and `--all-features` builds them all.

pub use tpt_thermo_core as core;
pub use tpt_thermo_data as data;

#[cfg(feature = "bubble-dew")]
pub use tpt_thermo_bubble_dew as bubble_dew;
#[cfg(feature = "electrolyte")]
pub use tpt_thermo_electrolyte as electrolyte;
#[cfg(feature = "activity")]
pub use tpt_thermo_eos_activity as eos_activity;
#[cfg(feature = "cubic")]
pub use tpt_thermo_eos_cubic as eos_cubic;
#[cfg(feature = "saft")]
pub use tpt_thermo_eos_saft as eos_saft;
#[cfg(feature = "flash")]
pub use tpt_thermo_flash as flash;
#[cfg(feature = "phase")]
pub use tpt_thermo_phase as phase;
#[cfg(feature = "polymer")]
pub use tpt_thermo_polymer as polymer;
#[cfg(feature = "transport")]
pub use tpt_thermo_transport as transport;

pub mod api;
pub mod error;
