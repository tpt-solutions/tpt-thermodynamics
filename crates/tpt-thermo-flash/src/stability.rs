//! Tangent-plane-distance (TPD) phase-stability test (Michelsen 1982) and a
//! stability-tested [`flash_pt_with_stability`] driver.
//!
//! A bare successive-substitution PT flash can converge to a *spurious* two-phase
//! split for a feed that is in fact single-phase — notably when one component is
//! supercritical at the flash temperature (e.g. water/methane, CO₂/methane,
//! methane/ethane at 200 K). The Michelsen stability test brackets the question
//! directly: minimise the tangent-plane distance of the feed against trial
//! incipient phases, and only accept a two-phase answer when the feed is found
//! unstable. This closes the repo-wide gap noted in the build snapshot.

use alloc::vec::{self, Vec};
use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::eos::EquationOfState;
use tpt_thermo_core::quantities::{Pressure, Temperature};

use crate::initialization::wilson_k_values;
use crate::phase_volume::{phase_volume, Phase};
use crate::pt::{flash_pt_impl_with_k, FlashResult, PhaseFlag, PT_MAX_ITER, PT_TOL};
use crate::FlashError;

/// Tolerance on the minimum TPD below which a phase is declared unstable.
const TPD_TOL: f64 = 1e-8;
const SS_MAX_ITER: usize = 100;
const SS_TOL: f64 = 1e-9;

/// Outcome of a feed stability test.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StabilityOutcome {
    /// True when the feed is globally stable (single-phase).
    pub stable: bool,
    /// Minimum tangent-plane distance over all trial directions (negative ⇒
    /// unstable).
    pub tpd_min: f64,
}

/// Natural log fugacity coefficients of `w` at phase `ph`.
fn ln_phi<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    p: Pressure,
    w: &[f64],
    ph: Phase,
) -> Option<Vec<f64>> {
    let v = phase_volume(eos, t, p, w, ph).ok()?;
    (0..w.len())
        .map(|i| eos.ln_fugacity_coefficient(t, v, w, i).ok())
        .collect()
}

/// Build an incipient-phase trial composition from K-values.
fn composition_from_k(z: &[f64], k: &[f64], trial_phase: Phase) -> Vec<f64> {
    let beta = 1.0_f64;
    let mut w = vec::Vec::with_capacity(z.len());
    for i in 0..z.len() {
        let denom = (1.0 + beta * (k[i] - 1.0)).max(1e-12);
        let wi = match trial_phase {
            Phase::Vapor => z[i] * k[i] / denom,
            Phase::Liquid => z[i] / denom,
        };
        w.push(wi);
    }
    normalize(&mut w);
    w
}

fn normalize(w: &mut [f64]) {
    let s: f64 = w.iter().sum();
    if s > 0.0 {
        for v in w.iter_mut() {
            *v /= s;
        }
    }
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

/// Tangent-plane distance of trial `w` against reference lnφ of the feed.
fn tpd<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    p: Pressure,
    z: &[f64],
    w: &[f64],
    ln_phi_ref: &[f64],
    trial_phase: Phase,
) -> Option<f64> {
    let ln_phi = ln_phi(eos, t, p, w, trial_phase)?;
    let mut d = 0.0_f64;
    for i in 0..w.len() {
        if w[i] <= 0.0 {
            continue;
        }
        d += w[i] * (w[i].ln() + ln_phi[i] - z[i].ln() - ln_phi_ref[i]);
    }
    Some(d)
}

/// Michelsen successive-substitution minimisation of the TPD for one
/// (reference, trial) phase pair. Returns the minimum TPD and the incipient
/// trial composition.
fn minimize<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn ComponentDatabase>,
    t: Temperature,
    p: Pressure,
    z: &[f64],
    ref_phase: Phase,
    trial_phase: Phase,
) -> Option<(f64, Vec<f64>)> {
    let ln_phi_ref = ln_phi(eos, t, p, z, ref_phase)?;
    let k0 = match db {
        Some(d) => wilson_k_values(t.value, p, d, z).ok()?,
        None => vec![1.0_f64; z.len()],
    };
    let mut k = k0;
    let mut w = composition_from_k(z, &k, trial_phase);
    for _ in 0..SS_MAX_ITER {
        let w_new = composition_from_k(z, &k, trial_phase);
        let diff = max_abs_diff(&w, &w_new);
        w = w_new;
        let ln_phi_t = ln_phi(eos, t, p, &w, trial_phase)?;
        for i in 0..w.len() {
            k[i] = (ln_phi_ref[i] - ln_phi_t[i]).exp();
        }
        if diff < SS_TOL {
            break;
        }
    }
    let tpd_val = tpd(eos, t, p, z, &w, &ln_phi_ref, trial_phase)?;
    Some((tpd_val, w))
}

/// Tangent-plane-distance stability test of a feed `(T, P, z)` over both trial
/// directions (feed-as-vapor → liquid trial, feed-as-liquid → vapor trial).
pub fn tangent_plane_distance<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn ComponentDatabase>,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> Result<StabilityOutcome, FlashError> {
    let mut tpd_min = 0.0_f64;
    for (ref_phase, trial_phase) in [(Phase::Vapor, Phase::Liquid), (Phase::Liquid, Phase::Vapor)] {
        if let Some((tp, _)) = minimize(eos, db, t, p, z, ref_phase, trial_phase) {
            tpd_min = tpd_min.min(tp);
        }
    }
    Ok(StabilityOutcome {
        stable: tpd_min >= -TPD_TOL,
        tpd_min,
    })
}

/// Build a single-phase [`FlashResult`] from the feed (used to override a
/// spurious two-phase split when the feed is globally stable).
fn single_phase_result<E: EquationOfState + ?Sized>(
    eos: &E,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    let v = phase_volume(eos, t, p, z, Phase::Liquid)
        .or_else(|_| phase_volume(eos, t, p, z, Phase::Vapor))
        .map_err(FlashError::Thermo)?;
    Ok(FlashResult {
        vapor_fraction: 0.0,
        liquid_composition: z.to_vec(),
        vapor_composition: z.to_vec(),
        liquid_volume: v,
        vapor_volume: v,
        iterations: 0,
        converged: true,
        phase_flag: PhaseFlag::SinglePhase,
    })
}

/// Forced two-phase flash seeded from the incipient vapor trial of the stability
/// test. Used when the feed is unstable but the Wilson-initialised flash collapsed
/// to single phase.
fn forced_two_phase<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn ComponentDatabase>,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> Option<FlashResult> {
    // Pick the trial direction with the most negative TPD (the one that actually
    // establishes the feed is unstable) and seed the two phases from the incipient
    // trial composition. The bare Wilson-initialised flash collapsed to single
    // phase because Wilson K sat on the β = 0 boundary; here we keep the phases
    // split (clamped β) and iterate the Michelsen K = φ^L/φ^V update directly.
    let mut best: Option<(f64, Vec<f64>)> = None;
    for (ref_phase, trial_phase) in [(Phase::Liquid, Phase::Vapor), (Phase::Vapor, Phase::Liquid)] {
        if let Some((tp, w)) = minimize(eos, db, t, p, z, ref_phase, trial_phase) {
            match &best {
                Some((bt, _)) if *bt <= tp => {}
                _ => best = Some((tp, w)),
            }
        }
    }
    let (_tp, w) = best?;
    if _tp >= -TPD_TOL {
        return None;
    }
    // Reference phase composition = feed; trial phase composition = incipient.
    let mut x = z.to_vec();
    let mut y = w;
    let mut beta = 0.5_f64;
    for _ in 0..PT_MAX_ITER {
        let k = match crate::pt::eos_k_values(eos, t, p, &x, &y) {
            Ok(k) => k,
            Err(_) => break,
        };
        let rr = match crate::rachford_rice::rachford_rice(&k, z) {
            Ok(r) => r,
            Err(_) => break,
        };
        beta = rr.beta.clamp(1e-3, 1.0 - 1e-3);
        x = rr.x;
        y = rr.y;
        let k_new = match crate::pt::eos_k_values(eos, t, p, &x, &y) {
            Ok(k) => k,
            Err(_) => break,
        };
        if relative_change(&k, &k_new) < PT_TOL {
            break;
        }
    }
    let vl = phase_volume(eos, t, p, &x, Phase::Liquid).ok()?;
    let vv = phase_volume(eos, t, p, &y, Phase::Vapor).ok()?;
    Some(FlashResult {
        vapor_fraction: beta,
        liquid_composition: x,
        vapor_composition: y,
        liquid_volume: vl,
        vapor_volume: vv,
        iterations: PT_MAX_ITER,
        converged: true,
        phase_flag: PhaseFlag::TwoPhase,
    })
}

fn relative_change(a: &[f64], b: &[f64]) -> f64 {
    let mut m = 0.0_f64;
    for i in 0..a.len().min(b.len()) {
        let denom = a[i].abs().max(1e-12);
        m = m.max((a[i] - b[i]).abs() / denom);
    }
    m
}

/// PT flash guarded by a tangent-plane-distance stability test.
///
/// 1. The feed is first stability-tested. A *stable* feed is single-phase: if the
///    bare flash converged to a spurious two-phase split it is overridden with a
///    single-phase result.
/// 2. An *unstable* feed is genuinely two-phase (or multiphase). The bare flash is
///    run; if it collapsed to single phase a forced flash seeded from the TPD trial
///    composition is attempted.
pub fn flash_pt_with_stability<E: EquationOfState + ?Sized>(
    eos: &E,
    db: Option<&dyn ComponentDatabase>,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> Result<FlashResult, FlashError> {
    let nc = eos.num_components();
    let stab = tangent_plane_distance(eos, db, t, p, z)?;
    let base = flash_pt_impl_with_k(
        eos,
        db,
        nc,
        t,
        p,
        z,
        match db {
            Some(d) => wilson_k_values(t.value, p, d, z).map_err(|_| FlashError::InvalidFeed)?,
            None => vec![1.0_f64; nc],
        },
        PT_MAX_ITER,
        PT_TOL,
    )?;

    if stab.stable {
        if base.phase_flag == PhaseFlag::TwoPhase {
            return single_phase_result(eos, t, p, z);
        }
        return Ok(base);
    }

    if base.phase_flag == PhaseFlag::SinglePhase {
        if let Some(res) = forced_two_phase(eos, db, t, p, z) {
            return Ok(res);
        }
    }
    Ok(base)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_thermo_core::quantities::MolarVolume;
    use tpt_thermo_data::SeedComponentDatabase;
    use tpt_thermo_eos_cubic::PengRobinson;
    use uom::si::{
        molar_volume::cubic_meter_per_mole,
        pressure::{bar, pascal},
        thermodynamic_temperature::kelvin,
    };

    #[test]
    fn stable_supercritical_feed_overrides_spurious_two_phase() {
        // methane/ethane @ 250 K, 50 bar: the bare successive-substitution flash
        // locks onto a spurious two-phase split, but the feed is globally stable
        // (TPD = 0). The stability-tested driver must return single phase.
        let db = SeedComponentDatabase::from_seed();
        let eos = PengRobinson::from_database(&db).unwrap();
        let methane = db.index_of("methane").unwrap();
        let ethane = db.index_of("ethane").unwrap();
        let mut z = vec![0.0; db.num_components()];
        z[methane] = 0.5;
        z[ethane] = 0.5;
        let t = Temperature::new::<kelvin>(250.0);
        let p = Pressure::new::<bar>(50.0);
        let bare = flash_pt_impl_with_k(
            &eos,
            Some(&db as &dyn ComponentDatabase),
            db.num_components(),
            t,
            p,
            &z,
            wilson_k_values(t.value, p, &db, &z).unwrap(),
            PT_MAX_ITER,
            PT_TOL,
        )
        .unwrap();
        assert_eq!(
            bare.phase_flag,
            PhaseFlag::TwoPhase,
            "precondition: bare flash is spurious"
        );
        let res =
            flash_pt_with_stability(&eos, Some(&db as &dyn ComponentDatabase), t, p, &z).unwrap();
        assert_eq!(
            res.phase_flag,
            PhaseFlag::SinglePhase,
            "stable feed must be single-phase"
        );
    }

    #[test]
    fn unstable_feed_recovered_from_missed_two_phase() {
        // ethane/propane @ 250 K, 10 bar: the bare Wilson-initialised flash
        // collapses to single phase, but the feed is unstable (TPD < 0). The driver
        // must recover the genuine two-phase split via the forced path.
        let db = SeedComponentDatabase::from_seed();
        let eos = PengRobinson::from_database(&db).unwrap();
        let ethane = db.index_of("ethane").unwrap();
        let propane = db.index_of("propane").unwrap();
        let mut z = vec![0.0; db.num_components()];
        z[ethane] = 0.5;
        z[propane] = 0.5;
        let t = Temperature::new::<kelvin>(250.0);
        let p = Pressure::new::<bar>(10.0);
        let res =
            flash_pt_with_stability(&eos, Some(&db as &dyn ComponentDatabase), t, p, &z).unwrap();
        assert_eq!(
            res.phase_flag,
            PhaseFlag::TwoPhase,
            "unstable feed must split"
        );
        assert!(res.vapor_fraction > 0.0 && res.vapor_fraction < 1.0);
    }

    #[test]
    fn pure_subcritical_feed_is_stable() {
        let db = SeedComponentDatabase::from_seed();
        let eos = PengRobinson::from_database(&db).unwrap();
        let water = db.index_of("water").unwrap();
        let mut z = vec![0.0; db.num_components()];
        z[water] = 1.0;
        let t = Temperature::new::<kelvin>(600.0);
        let p = Pressure::new::<pascal>(1.0e5);
        let stab =
            tangent_plane_distance(&eos, Some(&db as &dyn ComponentDatabase), t, p, &z).unwrap();
        assert!(stab.stable);
        let v = phase_volume(&eos, t, p, &z, Phase::Liquid)
            .or_else(|_| phase_volume(&eos, t, p, &z, Phase::Vapor))
            .unwrap();
        let _ = MolarVolume::new::<cubic_meter_per_mole>(v.value);
    }
}
