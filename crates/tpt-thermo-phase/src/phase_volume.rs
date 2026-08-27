//! Volume solvers: map `(T, P, z, phase)` to a molar volume of the requested
//! phase, used by the tangent-plane-distance machinery.
//!
//! For the in-repo cubic EoS we delegate to the engine's exact `solve_phase`;
//! for any other EoS a generic bracketing + Brent fallback inverts `pressure`.

use alloc::vec::Vec;
use tpt_thermo_core::brent;
use tpt_thermo_core::quantities::{MolarVolume, Pressure, Temperature};
use tpt_thermo_core::EquationOfState;
use tpt_thermo_eos_cubic::cubic_solver::Phase;
use tpt_thermo_eos_cubic::{CubicEos, PengRobinson, SoaveRedlichKwong, VolumeTranslated};
use uom::si::molar_volume::cubic_meter_per_mole;

/// Resolve the molar volume of a given phase at `(T, P, z)`.
pub trait PhaseVolume: Send + Sync {
    /// Molar volume of `phase` at `(T, P, z)`, or `None` if that phase has no
    /// real root (e.g. asking for a liquid root in the supercritical region).
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Option<MolarVolume>;
}

impl PhaseVolume for CubicEos {
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Option<MolarVolume> {
        self.solve_phase(t, p, z, phase).ok()
    }
}

impl PhaseVolume for PengRobinson {
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Option<MolarVolume> {
        self.engine().phase_volume(t, p, z, phase)
    }
}

impl PhaseVolume for SoaveRedlichKwong {
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Option<MolarVolume> {
        self.engine().phase_volume(t, p, z, phase)
    }
}

impl PhaseVolume for VolumeTranslated {
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Option<MolarVolume> {
        self.engine().phase_volume(t, p, z, phase)
    }
}

/// Generic fallback that inverts `pressure` for *any* EoS via bracketing + Brent.
///
/// Scans `pressure(T, v) − P` in log-volume, collects every sign change (a cubic
/// EoS yields two: liquid then vapor), and picks the root matching `phase`.
pub struct BrentPhaseVolume<'a> {
    eos: &'a dyn EquationOfState,
}

impl<'a> BrentPhaseVolume<'a> {
    /// Wrap any EoS so it satisfies [`PhaseVolume`].
    pub fn new(eos: &'a dyn EquationOfState) -> Self {
        Self { eos }
    }
}

impl PhaseVolume for BrentPhaseVolume<'_> {
    fn phase_volume(
        &self,
        t: Temperature,
        p: Pressure,
        z: &[f64],
        phase: Phase,
    ) -> Option<MolarVolume> {
        solve_generic(self.eos, t, p, z, phase).map(MolarVolume::new::<cubic_meter_per_mole>)
    }
}

fn solve_generic(
    eos: &dyn EquationOfState,
    t: Temperature,
    p: Pressure,
    z: &[f64],
    phase: Phase,
) -> Option<f64> {
    let v_min = 1e-7;
    let v_max = 1.0;
    let steps = 400usize;
    let press = |v: f64| -> f64 {
        eos.pressure(t, MolarVolume::new::<cubic_meter_per_mole>(v), z)
            .map(|pp| pp.value - p.value)
            .unwrap_or(f64::INFINITY)
    };
    let mut prev_v = v_min;
    let mut prev = press(v_min);
    let mut roots: Vec<f64> = Vec::new();
    for k in 1..=steps {
        let v = v_min * (v_max / v_min).powf(k as f64 / steps as f64);
        let fv = press(v);
        if prev.is_finite() && fv.is_finite() && prev * fv <= 0.0 {
            if let Ok(r) = brent(|x: f64| press(x), prev_v, v, 1e-9, 300) {
                roots.push(r);
            }
        }
        prev_v = v;
        prev = fv;
    }
    match phase {
        Phase::Liquid => roots.first().copied(),
        Phase::Vapor => roots.last().copied(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_thermo_core::component::ComponentDatabase;
    use tpt_thermo_data::SeedComponentDatabase;
    use tpt_thermo_eos_cubic::PengRobinson;
    use uom::si::pressure::pascal;
    use uom::si::thermodynamic_temperature::kelvin;

    #[test]
    fn cubic_liquid_and_vapor_roots_distinct() {
        let db = SeedComponentDatabase::from_seed();
        let eos = PengRobinson::from_database(&db).unwrap();
        let methane = db.index_of("methane").unwrap();
        let mut z = vec![0.0; db.num_components()];
        z[methane] = 1.0;
        let t = Temperature::new::<kelvin>(150.0);
        let p = Pressure::new::<pascal>(1.0e6);
        let v_l = eos
            .phase_volume(t, p, &z, Phase::Liquid)
            .expect("liquid root");
        let v_v = eos
            .phase_volume(t, p, &z, Phase::Vapor)
            .expect("vapor root");
        assert!(v_v > v_l, "vapor volume must exceed liquid volume");
    }

    #[test]
    fn brent_fallback_matches_cubic() {
        let db = SeedComponentDatabase::from_seed();
        let eos = PengRobinson::from_database(&db).unwrap();
        let methane = db.index_of("methane").unwrap();
        let mut z = vec![0.0; db.num_components()];
        z[methane] = 1.0;
        let t = Temperature::new::<kelvin>(150.0);
        let p = Pressure::new::<pascal>(1.0e6);
        let direct = eos.phase_volume(t, p, &z, Phase::Vapor).unwrap().value;
        let fb = BrentPhaseVolume::new(&eos)
            .phase_volume(t, p, &z, Phase::Vapor)
            .unwrap()
            .value;
        assert!((direct - fb).abs() / direct < 1e-6);
    }
}
