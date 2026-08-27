//! Ion association: Bjerrum criterion and mass-action association constants.

/// Bjerrum length `q` (m): the distance at which Coulombic and thermal energies are
/// equal, `q = |z_i z_j| e²/(4π ε₀ ε_r k T)`. Given the static relative permittivity
/// `eps_r` and temperature `t_kelvin`, returns `q` in metres.
pub fn bjerrum_length(eps_r: f64, t_kelvin: f64) -> f64 {
    // e²/(4π ε₀) = 2.307e-28 J·m; k = 1.380649e-23 J/K.
    let e2_4pi_eps0 = 2.307e-28;
    let k = 1.380_649e-23;
    (1.0 * 1.0) * e2_4pi_eps0 / (eps_r * k * t_kelvin)
}

/// Association (contact) distance `a` (m) conventionally taken as `q / 2` (or a fixed
/// ion-pair distance). Returns the Bjerrum criterion `q / a`.
pub fn bjerrum_criterion(eps_r: f64, t_kelvin: f64, a: f64) -> f64 {
    bjerrum_length(eps_r, t_kelvin) / a.max(1e-15)
}

/// Mass-action association constant `K_a` (m³·mol⁻¹) for a 1:1 ion pair via the Bjerrum
/// integral `K_a = 4π N_A ∫_a^q exp(u(r)) r² dr`, evaluated in closed form. `a` and
/// the returned constant are in metres / m³·mol⁻¹ respectively.
pub fn association_constant(eps_r: f64, t_kelvin: f64, a: f64) -> f64 {
    let q = bjerrum_length(eps_r, t_kelvin);
    let na = 6.022_140_76e23;
    let integral = if q > a {
        q - a + a * a / q - a // first-order Bjerrum approximation of the exponential integral
    } else {
        0.0
    };
    4.0 * core::f64::consts::PI * na * integral * 1.0e-3
}
