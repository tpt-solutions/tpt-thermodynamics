//! Batch PT flash over a table of feed compositions.

use alloc::vec::Vec;
use tpt_thermo_core::eos::EquationOfState;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{Pressure, Temperature};

use crate::pt::{flash_pt_impl, FlashResult};

/// Convert a [`crate::FlashError`] into a [`ThermoError`] for batch aggregation.
fn map_err(e: crate::FlashError) -> ThermoError {
    match e {
        crate::FlashError::Thermo(te) => te,
        crate::FlashError::NotConverged(_) => {
            ThermoError::Numerical(tpt_thermo_core::convergence::ConvergenceStatus::NotConverged)
        }
        crate::FlashError::InvalidFeed => ThermoError::InvalidInput("invalid feed"),
    }
}

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
        .map_err(map_err)?;
        out.push(r);
    }
    Ok(out)
}

/// Parallel variant of [`flash_pt_batch`] (requires the `std` feature). Feeds are
/// partitioned across the available hardware threads and solved concurrently;
/// results are returned in the original feed order. Each feed is an independent
/// non-linear solve, so this is the practical realisation of the deferred
/// explicit-SIMD/vectorised batch (the per-feed inner loop is itself iterative
/// and not directly SIMD-able).
#[cfg(feature = "std")]
pub fn flash_pt_batch_parallel<E: EquationOfState + ?Sized + Sync>(
    eos: &E,
    db: Option<&dyn tpt_thermo_core::component::ComponentDatabase>,
    t: Temperature,
    p: Pressure,
    feeds: &[Vec<f64>],
) -> Result<Vec<FlashResult>, ThermoError> {
    let nc = eos.num_components();
    for z in feeds {
        if z.len() != nc {
            return Err(ThermoError::InvalidInput("feed length mismatch in batch"));
        }
    }
    let n = feeds.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let n_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .max(1);
    let chunk_size = n.div_ceil(n_threads);
    type Slot = std::sync::Mutex<Option<Vec<Result<FlashResult, ThermoError>>>>;
    let slots: Vec<Slot> = (0..n_threads)
        .map(|_| std::sync::Mutex::new(None))
        .collect();
    std::thread::scope(|s| {
        for (ci, chunk) in feeds.chunks(chunk_size).enumerate() {
            let slot = &slots[ci];
            s.spawn(move || {
                let mut local = Vec::with_capacity(chunk.len());
                for z in chunk {
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
                    .map_err(map_err);
                    local.push(r);
                }
                *slot.lock().expect("flash batch slot poisoned") = Some(local);
            });
        }
    });
    let mut out = Vec::with_capacity(n);
    for slot in slots {
        if let Some(local) = slot.into_inner().expect("flash batch slot poisoned") {
            for r in local {
                out.push(r?);
            }
        }
    }
    Ok(out)
}
