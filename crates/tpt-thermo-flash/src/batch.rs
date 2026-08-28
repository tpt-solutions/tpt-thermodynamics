//! Batch PT flash over a table of feed compositions.

use alloc::vec::Vec;
use tpt_thermo_core::eos::EquationOfState;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{Pressure, Temperature};

use crate::pt::{flash_pt_impl, FlashResult};

/// Run [`flash_pt`](crate::pt::flash_pt) over every feed in `feeds` at the same
/// `(T, P)`. The first feed's length determines the component count; all feeds must
/// match. Returns one [`FlashResult`] per feed (failures are returned as
/// [`ThermoError`]s immediately, aborting the batch).
pub fn flash_pt_batch<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn tpt_thermo_core::component::ComponentDatabase>,
    t: Temperature,
    p: Pressure,
    feeds: &[Vec<f64>],
) -> Result<Vec<FlashResult>, ThermoError> {
    let nc = eos.num_components();
    let mut out = alloc::vec::Vec::with_capacity(feeds.len());
    for z in feeds {
        if z.len() != nc {
            return Err(ThermoError::InvalidInput("feed length mismatch in batch"));
        }
        let r = flash_pt_impl(
            eos,
            db,
            nc,
            t,
            p,
            z,
            crate::pt::PT_MAX_ITER,
            crate::pt::PT_TOL,
        )
        .map_err(|e| match e {
            crate::FlashError::Thermo(te) => te,
            crate::FlashError::NotConverged(_) => ThermoError::Numerical(
                tpt_thermo_core::convergence::ConvergenceStatus::NotConverged,
            ),
            crate::FlashError::InvalidFeed => ThermoError::InvalidInput("invalid feed"),
        })?;
        out.push(r);
    }
    Ok(out)
}
