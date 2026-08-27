//! Composition representations: mole/mass fraction and molality newtypes, plus
//! a [`Composition`] helper that normalises and converts between bases.

use crate::error::ThermoError;
use crate::quantities::MolarMass;
use alloc::vec::Vec;

/// A mole fraction, `x ∈ [0, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MoleFraction(pub f64);

impl MoleFraction {
    /// Construct, validating `0 ≤ x ≤ 1`.
    pub fn new(x: f64) -> Result<Self, CompositionError> {
        if !x.is_finite() || !(0.0..=1.0).contains(&x) {
            return Err(CompositionError::OutOfRange);
        }
        Ok(Self(x))
    }

    /// Raw value.
    pub fn get(&self) -> f64 {
        self.0
    }
}

/// A mass fraction, `w ∈ [0, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassFraction(pub f64);

impl MassFraction {
    /// Construct, validating `0 ≤ w ≤ 1`.
    pub fn new(w: f64) -> Result<Self, CompositionError> {
        if !w.is_finite() || !(0.0..=1.0).contains(&w) {
            return Err(CompositionError::OutOfRange);
        }
        Ok(Self(w))
    }

    /// Raw value.
    pub fn get(&self) -> f64 {
        self.0
    }
}

/// A molality, moles of solute per kilogram of solvent (`m ≥ 0`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Molality(pub f64);

impl Molality {
    /// Construct, validating `m ≥ 0`.
    pub fn new(m: f64) -> Result<Self, CompositionError> {
        if !m.is_finite() || m < 0.0 {
            return Err(CompositionError::OutOfRange);
        }
        Ok(Self(m))
    }

    /// Raw value (mol·kg⁻¹).
    pub fn get(&self) -> f64 {
        self.0
    }
}

/// Errors raised while building or converting compositions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompositionError {
    /// A fraction was outside its valid range.
    OutOfRange,
    /// The composition was empty.
    Empty,
    /// The provided molar masses did not match the component count.
    LengthMismatch,
}

impl core::fmt::Display for CompositionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CompositionError::OutOfRange => write!(f, "composition value out of range"),
            CompositionError::Empty => write!(f, "empty composition"),
            CompositionError::LengthMismatch => {
                write!(f, "molar-mass / composition length mismatch")
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for CompositionError {}

/// A mixture composition with conversion utilities between mole, mass, and
/// molality bases. The canonical storage is normalised mole fractions.
#[derive(Debug, Clone, PartialEq)]
pub struct Composition {
    mole_fractions: Vec<f64>,
}

impl Composition {
    /// Build from raw fractions, normalising so they sum to 1.
    pub fn from_mole_fractions(x: Vec<f64>) -> Result<Self, CompositionError> {
        if x.is_empty() {
            return Err(CompositionError::Empty);
        }
        let sum: f64 = x.iter().sum();
        if !sum.is_finite() || sum <= 0.0 {
            return Err(CompositionError::OutOfRange);
        }
        let mole_fractions = x.iter().map(|v| v / sum).collect();
        Ok(Self { mole_fractions })
    }

    /// Number of components.
    pub fn len(&self) -> usize {
        self.mole_fractions.len()
    }

    /// True if there are no components.
    pub fn is_empty(&self) -> bool {
        self.mole_fractions.is_empty()
    }

    /// Normalised mole fractions.
    pub fn mole_fractions(&self) -> &[f64] {
        &self.mole_fractions
    }

    /// Mole fraction of component `i`.
    pub fn mole_fraction(&self, i: usize) -> Result<f64, ThermoError> {
        self.mole_fractions
            .get(i)
            .copied()
            .ok_or(ThermoError::IndexOutOfRange(i))
    }

    /// Convert to mass fractions given per-component molar masses (kg·mol⁻¹).
    pub fn mass_fractions(&self, molar_masses: &[MolarMass]) -> Result<Vec<f64>, CompositionError> {
        if molar_masses.len() != self.mole_fractions.len() {
            return Err(CompositionError::LengthMismatch);
        }
        let m: Vec<f64> = molar_masses.iter().map(|mm| mm.value).collect();
        let xm: Vec<f64> = self
            .mole_fractions
            .iter()
            .zip(m.iter())
            .map(|(x, mm)| x * mm)
            .collect();
        let total: f64 = xm.iter().sum();
        if total <= 0.0 {
            return Err(CompositionError::OutOfRange);
        }
        Ok(xm.iter().map(|v| v / total).collect())
    }

    /// Convert to molalities (mol·kg⁻¹) with respect to solvent component `0`.
    ///
    /// Returns `m_i = x_i / (x_0 · M_0)` for solutes; the solvent entry (index
    /// `0`) is `0`. Requires molar masses (kg·mol⁻¹).
    pub fn molalities(&self, molar_masses: &[MolarMass]) -> Result<Vec<f64>, CompositionError> {
        if molar_masses.len() != self.mole_fractions.len() {
            return Err(CompositionError::LengthMismatch);
        }
        let m0 = molar_masses[0].value;
        let x0 = self.mole_fractions[0];
        if m0 <= 0.0 || x0 <= 0.0 {
            return Err(CompositionError::OutOfRange);
        }
        Ok(self
            .mole_fractions
            .iter()
            .enumerate()
            .map(|(i, xi)| if i == 0 { 0.0 } else { xi / (x0 * m0) })
            .collect())
    }

    /// Build a composition from mass fractions and molar masses (kg·mol⁻¹).
    pub fn from_mass_fractions(
        w: Vec<f64>,
        molar_masses: &[MolarMass],
    ) -> Result<Self, CompositionError> {
        if w.is_empty() {
            return Err(CompositionError::Empty);
        }
        if molar_masses.len() != w.len() {
            return Err(CompositionError::LengthMismatch);
        }
        let m: Vec<f64> = molar_masses.iter().map(|mm| mm.value).collect();
        let wm: Vec<f64> = w.iter().zip(m.iter()).map(|(wi, mi)| wi / mi).collect();
        Self::from_mole_fractions(wm)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quantities::MolarMass;
    use uom::si::molar_mass::kilogram_per_mole;

    #[test]
    fn normalises_and_round_trips() {
        let c = Composition::from_mole_fractions(alloc::vec![1.0, 3.0]).unwrap();
        assert!((c.mole_fractions()[0] - 0.25).abs() < 1e-12);
        assert!((c.mole_fractions()[1] - 0.75).abs() < 1e-12);
        // Sum is 1 by construction.
        let s: f64 = c.mole_fractions().iter().sum();
        assert!((s - 1.0).abs() < 1e-12);
    }

    #[test]
    fn mole_to_mass_fraction() {
        // 50/50 mole methane (16) / water-like (18) → mass skewed to water.
        let c = Composition::from_mole_fractions(alloc::vec![0.5, 0.5]).unwrap();
        let mm = [
            MolarMass::new::<kilogram_per_mole>(0.016),
            MolarMass::new::<kilogram_per_mole>(0.018),
        ];
        let w = c.mass_fractions(&mm).unwrap();
        assert!((w[0] - 0.016 / 0.034).abs() < 1e-12);
        assert!((w[1] - 0.018 / 0.034).abs() < 1e-12);
        let s: f64 = w.iter().sum();
        assert!((s - 1.0).abs() < 1e-12);
    }

    #[test]
    fn mass_to_mole_round_trip() {
        let mm = [
            MolarMass::new::<kilogram_per_mole>(0.016),
            MolarMass::new::<kilogram_per_mole>(0.018),
        ];
        let original = alloc::vec![0.3_f64, 0.7];
        let c = Composition::from_mole_fractions(original.clone()).unwrap();
        let w = c.mass_fractions(&mm).unwrap();
        let back = Composition::from_mass_fractions(w, &mm).unwrap();
        for (a, b) in original.iter().zip(back.mole_fractions().iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn molality_vs_solvent() {
        let c = Composition::from_mole_fractions(alloc::vec![0.9, 0.1]).unwrap();
        let mm = [
            MolarMass::new::<kilogram_per_mole>(0.018),
            MolarMass::new::<kilogram_per_mole>(0.016),
        ];
        let m = c.molalities(&mm).unwrap();
        // m_0 (solvent) is 0; m_1 = x_1 / (x_0 M_0) = 0.1 / (0.9 * 0.018) ≈ 6.1728.
        assert!((m[0]).abs() < 1e-12);
        assert!((m[1] - 0.1 / (0.9 * 0.018)).abs() < 1e-9);
    }
}
