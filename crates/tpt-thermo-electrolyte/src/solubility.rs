//! Gas solubility in electrolyte solutions: the Setschenow equation.

/// Setschenow salting-out/loading constant `k_s` (L·mol⁻¹) for a gas in a given
/// electrolyte. Returns the ratio of the gas solubility in the electrolyte to that in
/// pure water: `S/S₀ = exp(−k_s · c_salt)`.
///
/// * `ks` — Setschenow constant (L·mol⁻¹).
/// * `c_salt` — salt concentration (mol·L⁻¹).
pub fn setschenow_ratio(ks: f64, c_salt: f64) -> f64 {
    (-ks * c_salt).exp()
}

/// Solubility of a gas (mol·L⁻¹) in an electrolyte of salt concentration `c_salt`,
/// given the pure-water solubility `s0` and Setschenow constant `ks`.
pub fn gas_solubility_setschenow(s0: f64, ks: f64, c_salt: f64) -> f64 {
    s0 * setschenow_ratio(ks, c_salt)
}
