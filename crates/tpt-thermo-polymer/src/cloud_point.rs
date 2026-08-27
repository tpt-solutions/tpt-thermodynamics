//! Cloud-point (liquid–liquid phase split) detection for polymer solutions.
//!
//! Reuses the tangent-plane-distance (TPD) stability machinery in
//! [`tpt_thermo_phase`] to test whether a feed `(T, P, z)` is unstable to forming a
//! second *liquid* phase. A negative minimised TPD at the (Liquid, Liquid) stationarity
//! point indicates a cloud point (incipient L–L split).

use alloc::vec::Vec;
use tpt_thermo_core::quantities::{Pressure, Temperature};
use tpt_thermo_core::{ComponentDatabase, EquationOfState};
use tpt_thermo_eos_cubic::cubic_solver::Phase;
use tpt_thermo_phase::{BrentPhaseVolume, PhaseVolume, TangentPlaneDistance};

/// TPD tolerance below which a feed is declared unstable to L–L splitting.
const TPD_TOL: f64 = 1e-8;

/// Result of a cloud-point test.
#[derive(Debug, Clone)]
pub struct CloudPointResult {
    /// Whether an incipient second liquid phase was found.
    pub unstable: bool,
    /// Incipient (second-liquid) composition, if found.
    pub incipient_composition: Option<Vec<f64>>,
    /// Minimum tangent-plane distance (negative ⇒ unstable).
    pub tpd: f64,
}

/// Test `z` at `(T, P)` for a liquid–liquid cloud point using `eos`.
pub fn cloud_point<E>(
    eos: &E,
    db: &dyn ComponentDatabase,
    t: Temperature,
    p: Pressure,
    z: &[f64],
) -> CloudPointResult
where
    E: EquationOfState + Send + Sync + ?Sized,
{
    let volume = BrentPhaseVolume::new(eos);
    let calc = TangentPlaneDistance::new(eos, &volume, db, t, p, z.to_vec());
    match calc.minimize(Phase::Liquid, Phase::Liquid) {
        Some(sol) => CloudPointResult {
            unstable: sol.tpd < -TPD_TOL,
            incipient_composition: Some(sol.composition),
            tpd: sol.tpd,
        },
        None => CloudPointResult {
            unstable: false,
            incipient_composition: None,
            tpd: 0.0,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_thermo_core::quantities::MolarVolume;
    use tpt_thermo_data::SeedComponentDatabase;
    use tpt_thermo_eos_cubic::PengRobinson;
    use uom::si::{molar_volume::cubic_meter_per_mole, pressure::pascal, thermodynamic_temperature::kelvin};

    #[test]
    fn pure_feed_has_no_cloud_point() {
        // A single-component feed cannot undergo an L–L split.
        let db = SeedComponentDatabase::from_seed();
        let eos = PengRobinson::from_database(&db).unwrap();
        let methane = db.index_of("methane").unwrap();
        let mut z = vec![0.0; db.num_components()];
        z[methane] = 1.0;
        let res = cloud_point(&eos, &db, Temperature::new::<kelvin>(300.0), Pressure::new::<pascal>(1.0e5), &z);
        assert!(!res.unstable, "pure methane should be stable");
        let _ = MolarVolume::new::<cubic_meter_per_mole>(0.025);
    }
}
