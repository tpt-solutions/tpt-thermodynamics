//! Compile-time typed thermodynamic quantity aliases.
//!
//! These re-export the `f64`-backed `uom` SI quantities used throughout the
//! crate and add a few convenience aliases the spec mandates
//! ([`Temperature`], [`EnergyPerMol`], [`MolarEntropy`]) that are absent from
//! the thin `tpt-math-units` 0.1.0 surface.

pub use uom::si::f64::{
    AmountOfSubstance, Area, DiffusionCoefficient, DynamicViscosity, Energy, HeatCapacity, Length,
    Mass, MolarEnergy, MolarHeatCapacity, MolarMass, MolarVolume, Pressure, Ratio,
    ThermalConductivity, ThermodynamicTemperature, Time, Velocity, Volume,
};

/// Absolute thermodynamic temperature (alias for `ThermodynamicTemperature`).
pub use uom::si::f64::ThermodynamicTemperature as Temperature;

/// Molar energy (alias whose name matches the `EquationOfState` trait
/// signature).
pub use uom::si::f64::MolarEnergy as EnergyPerMol;

/// Molar entropy (J·mol⁻¹·K⁻¹). `uom` 0.38 has no named `molar_entropy`
/// quantity, so it is built from the ISQ dimension (L²·M·T⁻²·Θ⁻¹·N⁻¹), which is
/// dimensionally identical to molar heat capacity.
pub type MolarEntropy = uom::si::Quantity<
    uom::si::ISQ<
        typenum::P2,
        typenum::P1,
        typenum::N2,
        typenum::Z0,
        typenum::N1,
        typenum::N1,
        typenum::Z0,
    >,
    uom::si::SI<f64>,
    f64,
>;

/// The universal gas constant as a molar heat-capacity quantity (J·mol⁻¹·K⁻¹).
pub fn gas_constant() -> MolarHeatCapacity {
    MolarHeatCapacity::new::<uom::si::molar_heat_capacity::joule_per_kelvin_mole>(crate::R)
}

/// Construct a [`MolarEntropy`] from a raw value in J·mol⁻¹·K⁻¹.
pub fn molar_entropy(value: f64) -> MolarEntropy {
    MolarEntropy::new::<uom::si::molar_heat_capacity::joule_per_kelvin_mole>(value)
}
