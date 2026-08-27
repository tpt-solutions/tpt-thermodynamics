//! Shared helpers: per-component Lennard-Jones-style parameters derived from
//! the critical constants in a [`ComponentDatabase`].

use tpt_thermo_core::component::ComponentDatabase;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::MolarMass;

/// Lennard-Jones-like parameters used by the Chapman–Enskog / Chung
/// correlations, derived from critical constants.
#[derive(Debug, Clone, Copy)]
pub struct LjParams {
    /// Collision diameter `σ` (Å).
    pub sigma_a: f64,
    /// Energy parameter `ε/k` (K).
    pub eps_ok: f64,
    /// Molar mass (kg·mol⁻¹).
    pub molar_mass: f64,
    /// Acentric factor `ω`.
    pub omega: f64,
}

impl LjParams {
    /// Derive parameters for component `i` of `db`.
    ///
    /// `σ = 0.809386·V_c^{1/3}` with `V_c = Z_c R T_c / P_c` (`Z_c = 0.290`,
    /// `V_c` in cm³·mol⁻¹) and `ε/k = T_c / 1.2593`, the standard
    /// corresponding-states relations (Reid–Poling–Shields). `1.2593 T/T_c` is
    /// the reduced temperature fed to the collision-integral fits.
    pub fn from_database(db: &dyn ComponentDatabase, i: usize) -> Result<Self, ThermoError> {
        let tc = db.critical_temperature(i)?.value;
        let pc = db.critical_pressure(i)?.value;
        let omega = db.acentric_factor(i)?;
        let molar_mass: f64 = db.molar_mass(i)?.value;
        // V_c in cm³·mol⁻¹ (R = 83.14 cm³·bar·mol⁻¹·K⁻¹; Pc in bar).
        let pc_bar = pc / 1.0e5;
        let vc = 0.290 * 83.14 * tc / pc_bar.max(1.0);
        let sigma_a = 0.809386 * vc.powf(1.0 / 3.0);
        let eps_ok = tc / 1.2593;
        Ok(Self {
            sigma_a,
            eps_ok,
            molar_mass,
            omega,
        })
    }

    /// Reduced temperature `T* = 1.2593·T/T_c`.
    pub fn t_star(&self, t: f64) -> f64 {
        t / self.eps_ok
    }
}

/// Convenience: build a `Vec<LjParams>` for all components of `db`.
pub fn lj_params_for(db: &dyn ComponentDatabase) -> Result<Vec<LjParams>, ThermoError> {
    (0..db.num_components())
        .map(|i| LjParams::from_database(db, i))
        .collect()
}

/// Unit-safe accessor for a component's molar mass (kg·mol⁻¹) used by callers.
pub fn molar_mass_kg(db: &dyn ComponentDatabase, i: usize) -> Result<f64, ThermoError> {
    let _ = MolarMass::new::<uom::si::molar_mass::kilogram_per_mole>(0.0);
    Ok(db.molar_mass(i)?.value)
}
