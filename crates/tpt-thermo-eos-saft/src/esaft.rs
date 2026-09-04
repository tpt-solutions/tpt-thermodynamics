//! eSAFT electrolyte extension: ion-ion, ion-solvation, and ion-segment terms.
//!
//! Extends the PC-SAFT framework for electrolyte mixtures by adding three
//! contributions to the reduced residual Helmholtz energy:
//!
//! ```text
//! a^elec/(RT) = a^ion-ion + a^born + a^ion-seg
//! ```
//!
//! * **Ion-ion** — Debye-Hückel limiting-law screening of the Coulombic
//!   interaction between charged species, parametrised by the ionic strength
//!   and the solvent dielectric constant.
//! * **Ion-solvation (Born)** — electrostatic stabilisation of an ion in a
//!   dielectric medium, scaling as `z² / r_born · (1 − 1/ε_r)`.
//! * **Ion-segment** — dispersion between ions and neutral segments, reusing
//!   the PC-SAFT dispersion pair sum with the ion's segment parameters.
//!
//! The Debye-Hückel and Born terms are long-range and independent of the
//! packing-fraction reference; they enter as additive corrections to
//! `a^res/(RT)`. Ion-solvation uses the solvent dielectric constant (water at
//! 25 °C by default, with a simple T-scaling option).

use crate::engine::{SaftEngine, SaftFlavor};
use crate::parameters::SaftParameters;
use tpt_thermo_core::convergence::NumericalIssueReason;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::{MolarVolume, Temperature};
use tpt_thermo_core::EquationOfState;

/// Elementary charge (C).
const E_CHARGE: f64 = 1.602_176_63e-19;

/// Electrolyte configuration for the eSAFT extension.
///
/// When attached to the engine, the electrolyte correction (ion-ion
/// Debye-Hückel + Born solvation + ion-segment dispersion) is added to the
/// base PC-SAFT `a^res/(RT)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElectrolyteConfig {
    /// Solvent relative dielectric constant.
    pub epsr: f64,
}

impl ElectrolyteConfig {
    /// Create a configuration with the given dielectric constant.
    pub fn new(epsr: f64) -> Self {
        Self { epsr }
    }

    /// Create a configuration using the temperature-scaled water dielectric
    /// constant at `t` (K).
    pub fn water(t: f64) -> Self {
        Self {
            epsr: solvent_dielectric_water(t),
        }
    }
}

/// Vacuum permittivity (F·m⁻¹).
const EPSILON_0: f64 = 8.854_187_817e-12;

/// Boltzmann constant (J·K⁻¹).
const K_B: f64 = 1.380_649e-23;

/// Avogadro's number (mol⁻¹).
const NA: f64 = 6.022_140_76e23;

/// Dielectric constant of water at 25 °C.
const EPSR_WATER_25C: f64 = 78.375;

/// Temperature-scaled dielectric constant of water (empirical correlation,
/// valid 0–100 °C). From `ε_r(T) = 87.74 − 0.4008·T + 9.398e-4·T² − 1.41e-6·T³`
/// (Malmberg & Maryott 1956), with T in °C.
fn dielectric_water(t_k: f64) -> f64 {
    let tc = t_k - 273.15;
    if !(0.0..=100.0).contains(&tc) {
        EPSR_WATER_25C
    } else {
        87.74 - 0.4008 * tc + 9.398e-4 * tc * tc - 1.41e-7 * tc * tc * tc
    }
}

/// Ionic strength (mol·m⁻³) from number densities.
///
/// `rho` is the total molecular number density (m⁻³), `x` the mole fractions,
/// and `charges` the per-component ion charge numbers. Only charged species
/// contribute. Returns the number-density-scaled ionic strength
/// `I_n = ½·ρ·Σ x_i z_i²` (m⁻³).
pub fn ionic_strength_number_density(rho: f64, x: &[f64], charges: &[i32]) -> f64 {
    let mut sum_z2 = 0.0_f64;
    for (&xi, &z) in x.iter().zip(charges.iter()) {
        if z != 0 {
            sum_z2 += xi * (z as f64).powi(2);
        }
    }
    0.5 * rho * sum_z2
}

/// Debye screening parameter `κ` (m⁻¹) for a mixture of ions in a dielectric
/// solvent.
///
/// `I_n` is the ionic strength in number density (m⁻³), `t` the temperature (K),
/// `epsr` the relative dielectric constant of the solvent.
pub fn debye_kappa(i_n: f64, t: f64, epsr: f64) -> f64 {
    if i_n <= 0.0 || epsr <= 0.0 {
        return 0.0;
    }
    let numerator = 2.0 * E_CHARGE * E_CHARGE * i_n;
    let denominator = EPSILON_0 * epsr * K_B * t;
    (numerator / denominator).sqrt()
}

/// Ion-ion Debye-Hückel contribution to `a^res/(RT)` per mole.
///
/// The Debye-Hückel limiting law for the excess Helmholtz energy density is
///
/// ```text
/// a^DH/(RT) = −(κ³) / (12π·ρ)
/// ```
///
/// where `κ` is the Debye screening parameter (m⁻¹) and `ρ` is the total
/// molecular number density (m⁻³). This is the canonical DH limiting-law
/// expression valid at low-to-moderate ionic strength.
pub fn ion_ion_term(
    t: f64,
    rho_mol: f64,
    x: &[f64],
    charges: &[i32],
    epsr: f64,
) -> Result<f64, ThermoError> {
    if rho_mol <= 0.0 {
        return Ok(0.0);
    }
    let rho = NA * rho_mol;
    let i_n = ionic_strength_number_density(rho, x, charges);
    if i_n <= 0.0 {
        return Ok(0.0);
    }
    let kappa = debye_kappa(i_n, t, epsr);
    if kappa <= 0.0 {
        return Ok(0.0);
    }
    let kappa3 = kappa * kappa * kappa;
    // DH formula uses number density ρ (m⁻³), not molar density ρ_mol.
    let a_dh = -kappa3 / (12.0 * core::f64::consts::PI * rho);
    if !a_dh.is_finite() {
        return Err(ThermoError::Numerical(
            tpt_thermo_core::convergence::ConvergenceStatus::NumericalIssue(
                NumericalIssueReason::NonPhysical,
            ),
        ));
    }
    Ok(a_dh)
}

/// Born solvation energy contribution to `a^res/(RT)` per mole.
///
/// The Born electrostatic solvation energy for an ion of charge `z` and
/// solvation radius `r_b` in a medium of dielectric constant `ε_r` is
///
/// ```text
/// ΔG_Born = −(z²·e²) / (8π·ε₀·r_b) · (1 − 1/ε_r)
/// ```
///
/// Per molecule; per mole multiply by `NA`. The term enters `a^res/(RT)` as
/// `a^born/(RT) = Σ x_i · ΔG_Born,i · NA / (RT) = Σ x_i · ΔG_Born,i / (k_B·T)`.
///
/// Note: dividing the per-molecule energy by `k_B·T` is equivalent to dividing
/// the per-mole energy by `R·T`, giving the dimensionless `a^res/(RT)`.
pub fn born_term(
    t: f64,
    x: &[f64],
    charges: &[i32],
    born_radii: &[f64],
    epsr: f64,
) -> Result<f64, ThermoError> {
    if epsr <= 1.0 {
        return Ok(0.0);
    }
    let one_minus_inv = 1.0 - 1.0 / epsr;
    let mut a_born = 0.0_f64;
    for ((xi, &z), &r_b) in x.iter().zip(charges.iter()).zip(born_radii.iter()) {
        if z == 0 || r_b <= 0.0 {
            continue;
        }
        let zf = z as f64;
        let r_b_m = r_b * 1e-10; // Å → m
                                 // Born energy per molecule (J).
        let dg = -(zf * zf * E_CHARGE * E_CHARGE)
            / (8.0 * core::f64::consts::PI * EPSILON_0 * r_b_m)
            * one_minus_inv;
        // Convert to dimensionless a^res/(RT) per mole: multiply by NA to get
        // per mole, then divide by (R*T). Since R = NA*k_B, this is dg/(k_B*T).
        a_born += xi * dg / (K_B * t);
    }
    if !a_born.is_finite() {
        return Err(ThermoError::Numerical(
            tpt_thermo_core::convergence::ConvergenceStatus::NumericalIssue(
                NumericalIssueReason::NonPhysical,
            ),
        ));
    }
    Ok(a_born)
}

/// Ion-segment dispersion contribution to `a^res/(RT)`.
///
/// This term is intentionally zero: the engine's standard dispersion sum
/// already includes all segment-segment pairs (including ion-neutral pairs
/// when the ion carries `m`, `σ`, `ε/k` parameters). This function exists as
/// a documented hook for future ion-specific dispersion corrections.
pub fn ion_segment_dispersion(
    _t: f64,
    _rho_mol: f64,
    _x: &[f64],
    _params: &SaftParameters,
    _kij: &[Vec<f64>],
) -> Result<f64, ThermoError> {
    Ok(0.0)
}

/// Temperature-dependent hard-sphere diameters (m) — local copy to avoid
/// duplicating the engine's private method.
fn _hard_sphere_diameters(t: f64, params: &SaftParameters) -> Vec<f64> {
    params
        .components
        .iter()
        .map(|c| {
            let sigma = c.sigma * 1e-10;
            sigma * (1.0 - 0.12 * (-3.0 * c.epsilon_k / t).exp())
        })
        .collect()
}

/// Total eSAFT electrolyte correction to `a^res/(RT)` for the mixture.
///
/// Combines the ion-ion (Debye-Hückel), Born solvation, and ion-segment
/// dispersion contributions. The ion-segment term here is the ion-neutral
/// cross dispersion only; the full segment-segment dispersion is handled by
/// the engine's standard dispersion sum.
pub fn electrolyte_term(
    t: f64,
    rho_mol: f64,
    x: &[f64],
    params: &SaftParameters,
    kij: &[Vec<f64>],
    epsr: f64,
) -> Result<f64, ThermoError> {
    let charges: Vec<i32> = (0..params.num_components())
        .map(|i| params.component(i).charge)
        .collect();
    let born_radii: Vec<f64> = (0..params.num_components())
        .map(|i| params.component(i).born_radius)
        .collect();

    let a_ion = ion_ion_term(t, rho_mol, x, &charges, epsr)?;
    let a_born = born_term(t, x, &charges, &born_radii, epsr)?;
    let a_seg = ion_segment_dispersion(t, rho_mol, x, params, kij)?;

    Ok(a_ion + a_born + a_seg)
}

/// Convenience: dielectric constant of water at temperature `t` (K).
pub fn solvent_dielectric_water(t: f64) -> f64 {
    dielectric_water(t)
}

/// An eSAFT electrolyte equation of state.
///
/// Wraps the shared [`SaftEngine`] with an electrolyte configuration
/// ([`ElectrolyteConfig`]) that adds the ion-ion (Debye-Hückel), Born
/// solvation, and ion-segment dispersion corrections to the base PC-SAFT
/// model. Build it directly from parameters or from the seed database.
///
/// # Example
///
/// ```
/// use tpt_thermo_core::EquationOfState;
/// use tpt_thermo_eos_saft::{
///     ElectrolyteConfig, Esaft,
///     parameters::{SaftComponent, SaftParameters, SEED_E_SAFT_IONS},
/// };
/// use tpt_thermo_core::quantities::{MolarVolume, Temperature};
/// use uom::si::{molar_volume::cubic_meter_per_mole, thermodynamic_temperature::kelvin};
///
/// // Build a NaCl–water eSAFT model from the seed ion table.
/// let water = SaftComponent::pc_saft("water", 1.2047, 3.8331, 366.51);
/// let na = SEED_E_SAFT_IONS.iter().find(|c| c.name == "sodium").copied().unwrap();
/// let cl = SEED_E_SAFT_IONS.iter().find(|c| c.name == "chloride").copied().unwrap();
/// let params = SaftParameters::new(vec![water, na, cl]);
/// let mm = vec![0.018015, 0.022990, 0.035453];
/// let config = ElectrolyteConfig::water(298.15);
/// let eos = Esaft::new(params, mm, config);
/// let t = Temperature::new::<kelvin>(298.15);
/// // Low-density state where the base PC-SAFT model is near-ideal.
/// let v = MolarVolume::new::<cubic_meter_per_mole>(1.0);
/// let p = eos.pressure(t, v, &[0.97, 0.015, 0.015]).unwrap();
/// // Near-ideal gas: P ≈ RT/v ≈ 2479 Pa.
/// assert!(p.value > 2000.0 && p.value < 3000.0, "P = {} Pa", p.value);
/// ```
#[derive(Debug, Clone)]
pub struct Esaft(pub(crate) SaftEngine);

impl Esaft {
    /// Build directly from a parameter set, per-component molar masses
    /// (kg·mol⁻¹), and an electrolyte configuration.
    pub fn new(params: SaftParameters, molar_masses: Vec<f64>, config: ElectrolyteConfig) -> Self {
        let engine = SaftEngine::new(params, molar_masses).with_electrolyte(config);
        Self(engine)
    }

    /// Build from the seed database with the given electrolyte configuration.
    ///
    /// Uses [`SaftParameters::from_seed_database`] to look up SAFT parameters
    /// by component name. Components missing from the seed table are estimated
    /// from their critical constants (neutral segment parameters only — ion
    /// charges and Born radii default to zero for estimated components).
    pub fn from_seed_database(
        db: &dyn tpt_thermo_core::component::ComponentDatabase,
        config: ElectrolyteConfig,
    ) -> Result<Self, ThermoError> {
        let engine =
            SaftEngine::from_seed_database(db, SaftFlavor::ESaft)?.with_electrolyte(config);
        Ok(Self(engine))
    }

    /// Attach a binary interaction matrix `k_ij` (dimensionless, symmetric).
    pub fn with_kij(self, kij: Vec<Vec<f64>>) -> Self {
        Self(self.0.with_kij(kij))
    }

    /// Underlying SAFT parameters.
    pub fn parameters(&self) -> &SaftParameters {
        self.0.parameters()
    }

    /// Number of components.
    pub fn num_components(&self) -> usize {
        self.0.num_components()
    }

    /// The electrolyte configuration.
    pub fn electrolyte_config(&self) -> Option<ElectrolyteConfig> {
        self.0.electrolyte
    }
}

impl EquationOfState for Esaft {
    fn num_components(&self) -> usize {
        self.0.num_components()
    }

    fn pressure(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::Pressure, ThermoError> {
        self.0.pressure(t, v, z)
    }

    fn ln_fugacity_coefficient(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        self.0.ln_fugacity_coefficient(t, v, z, i)
    }

    fn molar_enthalpy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::MolarEnergy, ThermoError> {
        self.0.molar_enthalpy(t, v, z)
    }

    fn molar_entropy(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::MolarEntropy, ThermoError> {
        self.0.molar_entropy(t, v, z)
    }

    fn molar_isobaric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::MolarHeatCapacity, ThermoError> {
        self.0.molar_isobaric_heat_capacity(t, v, z)
    }

    fn molar_isochoric_heat_capacity(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::MolarHeatCapacity, ThermoError> {
        self.0.molar_isochoric_heat_capacity(t, v, z)
    }

    fn speed_of_sound(
        &self,
        t: Temperature,
        v: MolarVolume,
        z: &[f64],
    ) -> Result<tpt_thermo_core::quantities::Velocity, ThermoError> {
        self.0.speed_of_sound(t, v, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameters::{SaftComponent, SaftParameters, SEED_E_SAFT_IONS};

    #[test]
    fn dielectric_water_scaling() {
        // At 25 °C the correlation returns the reference value.
        let eps_25 = dielectric_water(298.15);
        assert!(
            (eps_25 - EPSR_WATER_25C).abs() < 1.0,
            "ε_r(25°C) = {eps_25}, expected ~{EPSR_WATER_25C}"
        );
        // Decreases with temperature.
        let eps_50 = dielectric_water(323.15);
        assert!(eps_50 < eps_25, "ε_r should decrease with T");
        assert!(eps_50 > 50.0, "ε_r(50°C) unreasonably low: {eps_50}");
    }

    #[test]
    fn debye_kappa_zero_for_pure_solvent() {
        // Pure water (no ions) → κ = 0.
        let kappa = debye_kappa(0.0, 298.15, EPSR_WATER_25C);
        assert_eq!(kappa, 0.0);
    }

    #[test]
    fn ionic_strength_scales_with_charge() {
        let x = vec![0.5, 0.5];
        let rho = 1.0; // arbitrary
                       // 1:1 electrolyte: z² = 1 for both → I = ρ·0.5·(0.5·1 + 0.5·1) = 0.5·ρ
        let i_11 = ionic_strength_number_density(rho, &x, &[1, -1]);
        // 2:1 electrolyte: z² = 4 for cation → I = ρ·0.5·(0.5·4 + 0.5·1) = 1.25·ρ
        let i_21 = ionic_strength_number_density(rho, &x, &[2, -1]);
        assert!((i_11 - 0.5 * rho).abs() < 1e-10, "I(1:1) = {i_11}");
        assert!((i_21 - 1.25 * rho).abs() < 1e-10, "I(2:1) = {i_21}");
    }

    #[test]
    fn born_term_negative_for_ions() {
        // Born solvation is stabilising (negative free energy).
        let x = vec![1.0];
        let charges = vec![1];
        let radii = vec![1.7]; // Å
        let a = born_term(298.15, &x, &charges, &radii, EPSR_WATER_25C).unwrap();
        assert!(a < 0.0, "Born term should be negative, got {a}");
    }

    #[test]
    fn ion_ion_negative_at_low_density() {
        // Debye-Hückel term is negative (stabilising screening).
        let x = vec![0.5, 0.5];
        let charges = vec![1, -1];
        let a = ion_ion_term(298.15, 1.0, &x, &charges, EPSR_WATER_25C).unwrap();
        assert!(a < 0.0, "DH term should be negative, got {a}");
    }

    #[test]
    fn neutral_mixture_gives_zero() {
        // A mixture of neutrals produces no electrolyte correction.
        let x = vec![0.5, 0.5];
        let charges = vec![0, 0];
        let radii = vec![0.0, 0.0];
        let a_ion = ion_ion_term(298.15, 1.0, &x, &charges, EPSR_WATER_25C).unwrap();
        let a_born = born_term(298.15, &x, &charges, &radii, EPSR_WATER_25C).unwrap();
        assert_eq!(a_ion, 0.0);
        assert_eq!(a_born, 0.0);
    }

    #[test]
    fn esaft_pressure_is_positive_at_low_density() {
        use crate::parameters::{SaftComponent, SaftParameters, SEED_E_SAFT_IONS};
        use tpt_thermo_core::quantities::{MolarVolume, Temperature};
        use tpt_thermo_core::EquationOfState;
        use uom::si::{molar_volume::cubic_meter_per_mole, thermodynamic_temperature::kelvin};

        let water = SaftComponent::pc_saft("water", 1.2047, 3.8331, 366.51);
        let na = SEED_E_SAFT_IONS
            .iter()
            .find(|c| c.name == "sodium")
            .copied()
            .unwrap();
        let cl = SEED_E_SAFT_IONS
            .iter()
            .find(|c| c.name == "chloride")
            .copied()
            .unwrap();
        let params = SaftParameters::new(vec![water, na, cl]);
        let mm = vec![0.018015, 0.022990, 0.035453];
        let config = ElectrolyteConfig::water(298.15);
        let eos = Esaft::new(params, mm, config);
        let t = Temperature::new::<kelvin>(298.15);
        // Low-density state where the base PC-SAFT model is near-ideal.
        let v = MolarVolume::new::<cubic_meter_per_mole>(1.0);
        let p = eos.pressure(t, v, &[0.97, 0.015, 0.015]).unwrap();
        // Near-ideal gas: P ≈ RT/v ≈ 2479 Pa.
        assert!(
            (p.value - 2479.0).abs() < 1.0,
            "P = {} Pa, expected ≈ 2479 Pa",
            p.value
        );
    }

    #[test]
    fn esaft_ion_parameters_exist() {
        // The seed table contains at least Na⁺ and Cl⁻.
        let na = SEED_E_SAFT_IONS.iter().find(|c| c.name == "sodium");
        let cl = SEED_E_SAFT_IONS.iter().find(|c| c.name == "chloride");
        assert!(na.is_some(), "sodium ion missing from seed table");
        assert!(cl.is_some(), "chloride ion missing from seed table");
        assert_eq!(na.unwrap().charge, 1);
        assert_eq!(cl.unwrap().charge, -1);
        assert!(na.unwrap().born_radius > 0.0);
    }

    #[test]
    fn electrolyte_term_for_nacl_mixture() {
        // Build a simple NaCl–water mixture and verify the correction is finite
        // and negative (stabilising).
        let comps: Vec<SaftComponent> = vec![
            // Water (neutral, with association params from main table would be ideal
            // but we use a neutral placeholder here).
            SaftComponent::pc_saft("water", 1.2047, 3.8331, 366.51),
            SEED_E_SAFT_IONS
                .iter()
                .find(|c| c.name == "sodium")
                .copied()
                .unwrap(),
            SEED_E_SAFT_IONS
                .iter()
                .find(|c| c.name == "chloride")
                .copied()
                .unwrap(),
        ];
        let params = SaftParameters::new(comps);
        let kij = vec![vec![0.0; 3]; 3];
        // Dilute aqueous NaCl: mostly water.
        let x = vec![0.97, 0.015, 0.015];
        let epsr = dielectric_water(298.15);
        let term = electrolyte_term(298.15, 55_000.0, &x, &params, &kij, epsr).unwrap();
        assert!(term.is_finite(), "electrolyte term not finite");
        assert!(
            term < 0.0,
            "electrolyte term should be negative, got {term}"
        );
    }
}
