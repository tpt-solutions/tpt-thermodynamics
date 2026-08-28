//! High-level convenience API: re-exports of the most-used entry points, gated by
//! their feature so `tpt_thermo::api::FlashCalculator` only exists when `flash` is on.

#[cfg(feature = "flash")]
pub use tpt_thermo_flash::{
    flash_ph, flash_pt, flash_pu, flash_pv, flash_ts, flash_tv, FlashCalculator, FlashResult,
};

#[cfg(feature = "bubble-dew")]
pub use tpt_thermo_bubble_dew::{BubbleDewSolver, KProvider};

#[cfg(feature = "transport")]
pub use tpt_thermo_transport::{
    conductivity, diffusivity, mixing_rules, residual_entropy_scaling, viscosity,
};
