//! Electrolyte model parameters.

/// Pitzer (1973) virial parameters for a single electrolyte (molality scale).
///
/// `beta0`, `beta1`, `beta2` are the second-virial terms, `cphi` the
/// third-virial term, and `alpha1`/`alpha2` the ion-size parameters for the
/// exponential term (commonly `alpha1 = 2.0`, `alpha2 = 0.0` for 1:1 salts).
#[derive(Debug, Clone, Copy)]
pub struct PitzerParams {
    /// Second-virial `β⁰`.
    pub beta0: f64,
    /// Second-virial `β¹`.
    pub beta1: f64,
    /// Second-virial `β²`.
    pub beta2: f64,
    /// Third-virial `C^φ`.
    pub cphi: f64,
    /// Ion-size parameter `α¹`.
    pub alpha1: f64,
    /// Ion-size parameter `α²`.
    pub alpha2: f64,
}

impl PitzerParams {
    /// Standard 1:1 parameters (NaCl-like defaults); caller overrides per salt.
    pub const fn one_to_one(beta0: f64, beta1: f64, cphi: f64) -> Self {
        Self {
            beta0,
            beta1,
            beta2: 0.0,
            cphi,
            alpha1: 2.0,
            alpha2: 0.0,
        }
    }
}

/// HKF (1981) standard partial-molal property parameters for a species.
///
/// `a1`–`a4` are the high-T coefficients, `c1`/`c2` the temperature/pressure
/// integrals, and `omega` the omega term driving the solvation contribution.
#[derive(Debug, Clone, Copy)]
pub struct HkfParams {
    /// `a₁` (J·bar⁻¹·mol⁻¹).
    pub a1: f64,
    /// `a₂` (J·mol⁻¹).
    pub a2: f64,
    /// `a₃` (J·K·mol⁻¹).
    pub a3: f64,
    /// `a₄` (J·K·mol⁻¹).
    pub a4: f64,
    /// `c₁` (J·K·mol⁻¹).
    pub c1: f64,
    /// `c₂` (J·mol⁻¹).
    pub c2: f64,
    /// `ω` (J·mol⁻¹).
    pub omega: f64,
    /// Reference molar volume `V°` at 298.15 K, 1 bar (cm³·mol⁻¹).
    pub v_ref: f64,
    /// Reference heat capacity `Cp°` at 298.15 K, 1 bar (J·mol⁻¹·K⁻¹).
    pub cp_ref: f64,
}
