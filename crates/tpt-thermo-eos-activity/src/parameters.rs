//! Temperature-dependent interaction-parameter machinery shared by the activity
//! models.
//!
//! Three forms appear across the models in this crate:
//!
//! * [`TdParam`] — the general form `a + b/T + c·ln(T) + d·T + e/T²`, used by
//!   NRTL τ-parameters. The spec's "A + B/T + C·ln(T)" helper is the `(a, b, c)`
//!   subset; `d`/`e` are provided for completeness and default to zero.
//! * [`InteractionMatrix`] — a full `n×n` matrix of asymmetric binary
//!   energy-like parameters (units of K) used by Wilson and UNIQUAC in the
//!   reduced form `exp(-(a_ij − a_ii)/T)`.
//! * [`Volumes`] — component molar volumes (m³·mol⁻¹) required by Wilson's
//!   `Λ_ij = (V_j / V_i) · exp(-(a_ij − a_ii)/(R T))`.

use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::quantities::Temperature;

/// A temperature-dependent binary parameter:
/// `value(T) = a + b/T + c·ln(T) + d·T + e/T²`.
///
/// `T` is absolute temperature in kelvin. The default is the constant-zero
/// parameter; the spec's `A + B/T + C·ln(T)` three-term form is
/// [`TdParam::new`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TdParam {
    /// Constant term (K).
    pub a: f64,
    /// `1/T` term coefficient (K²).
    pub b: f64,
    /// `ln(T)` term coefficient (K).
    pub c: f64,
    /// Linear-in-`T` term coefficient (dimensionless).
    pub d: f64,
    /// `1/T²` term coefficient (K³).
    pub e: f64,
}

impl TdParam {
    /// The spec's three-term form `a + b/T + c·ln(T)`.
    pub const fn new(a: f64, b: f64, c: f64) -> Self {
        Self {
            a,
            b,
            c,
            d: 0.0,
            e: 0.0,
        }
    }

    /// Evaluate at absolute temperature `t` (kelvin).
    pub fn value(&self, t: Temperature) -> f64 {
        let tk = t.value;
        self.a + self.b / tk + self.c * libm::log(tk) + self.d * tk + self.e / (tk * tk)
    }
}

impl Default for TdParam {
    fn default() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }
}

/// A full `n×n` matrix of [`TdParam`]s (need not be symmetric).
#[derive(Debug, Clone)]
pub struct TdMatrix {
    n: usize,
    data: Vec<Vec<TdParam>>,
}

impl TdMatrix {
    /// A zero matrix of size `n×n`.
    pub fn zeros(n: usize) -> Self {
        Self {
            n,
            data: vec![vec![TdParam::default(); n]; n],
        }
    }

    /// Build from a full `n×n` matrix of parameters.
    pub fn from_full(data: Vec<Vec<TdParam>>) -> Result<Self, ThermoError> {
        let n = data.len();
        if n == 0 || data.iter().any(|row| row.len() != n) {
            return Err(ThermoError::InvalidInput("TdMatrix must be square and non-empty"));
        }
        Ok(Self { n, data })
    }

    /// Number of components.
    pub fn len(&self) -> usize {
        self.n
    }

    /// `true` if the matrix is empty (never the case for a valid model).
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Set entry `(i, j)`.
    pub fn set(&mut self, i: usize, j: usize, p: TdParam) {
        self.data[i][j] = p;
    }

    /// Get entry `(i, j)`.
    pub fn get(&self, i: usize, j: usize) -> TdParam {
        self.data[i][j]
    }

    /// Evaluate entry `(i, j)` at temperature `t`.
    pub fn value_at(&self, i: usize, j: usize, t: Temperature) -> f64 {
        self.data[i][j].value(t)
    }
}

/// A full `n×n` matrix of asymmetric energy-like binary parameters in units of
/// kelvin, used in the reduced form `exp(-(a_ij − a_ii)/T)`.
#[derive(Debug, Clone)]
pub struct InteractionMatrix {
    n: usize,
    data: Vec<Vec<f64>>,
}

impl InteractionMatrix {
    /// A zero matrix of size `n×n`.
    pub fn zeros(n: usize) -> Self {
        Self {
            n,
            data: vec![vec![0.0; n]; n],
        }
    }

    /// Build from a full `n×n` matrix.
    pub fn from_full(data: Vec<Vec<f64>>) -> Result<Self, ThermoError> {
        let n = data.len();
        if n == 0 || data.iter().any(|row| row.len() != n) {
            return Err(ThermoError::InvalidInput("InteractionMatrix must be square and non-empty"));
        }
        Ok(Self { n, data })
    }

    /// Number of components.
    pub fn len(&self) -> usize {
        self.n
    }

    /// `true` if the matrix is empty.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Set entry `(i, j)` (K).
    pub fn set(&mut self, i: usize, j: usize, a: f64) {
        self.data[i][j] = a;
    }

    /// Get entry `(i, j)` (K).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i][j]
    }
}

/// Per-component molar volumes (m³·mol⁻¹) used by Wilson's size-ratio term.
#[derive(Debug, Clone)]
pub struct Volumes(pub Vec<f64>);

impl Volumes {
    /// Build from a list of molar volumes (m³·mol⁻¹).
    pub fn new(v: Vec<f64>) -> Self {
        Self(v)
    }

    /// `V_i` for component `i` (m³·mol⁻¹).
    pub fn get(&self, i: usize) -> f64 {
        self.0[i]
    }

    /// Number of components.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// Sum of a slice of `f64`s (the activity models are `no_std` and cannot use
/// iterator `.sum::<f64>()` on a plain slice without `iter()` — but they can;
/// this is a small named helper for readability at call sites).
pub(crate) fn sum(xs: &[f64]) -> f64 {
    let mut s = 0.0;
    for &x in xs {
        s += x;
    }
    s
}

/// Require that a composition is non-empty and normalised, returning an error
/// otherwise. Used by every model's trait methods.
pub(crate) fn check_composition(x: &[f64]) -> Result<(), ThermoError> {
    if x.is_empty() {
        return Err(ThermoError::InvalidInput("empty composition"));
    }
    let total: f64 = sum(x);
    if (total - 1.0).abs() > 1e-6 {
        return Err(ThermoError::InvalidInput("composition does not sum to 1"));
    }
    Ok(())
}

/// Absolute temperature in kelvin as a plain `f64`.
pub(crate) fn tk(t: Temperature) -> f64 {
    t.value
}
