//! Real-root solver for the compressibility cubic `a z³ + b z² + c z + d = 0`
//! (Cardano / trigonometric method), used to select physically-meaningful
//! (vapor / liquid) roots of a cubic EoS, plus the [`Phase`] selection helper.

/// Which cubic model a [`crate::CubicEos`] uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubicModel {
    /// Peng-Robinson (1976).
    PengRobinson,
    /// Soave-Redlich-Kwong (1972).
    SoaveRedlichKwong,
}

impl CubicModel {
    /// `a`-factor in `a_i = a_factor·R²·Tc²/Pc`.
    pub fn a_factor(&self) -> f64 {
        match self {
            CubicModel::PengRobinson => 0.457235529,
            CubicModel::SoaveRedlichKwong => 0.427480233,
        }
    }

    /// `b`-factor in `b_i = b_factor·R·Tc/Pc`.
    pub fn b_factor(&self) -> f64 {
        match self {
            CubicModel::PengRobinson => 0.077796074,
            CubicModel::SoaveRedlichKwong => 0.086640350,
        }
    }

    /// Linear coefficient of `b` in the attractive denominator (`u` in
    /// `v² + u·b·v + w·b²`).
    pub fn u(&self) -> f64 {
        match self {
            CubicModel::PengRobinson => 2.0,
            CubicModel::SoaveRedlichKwong => 1.0,
        }
    }

    /// Quadratic coefficient of `b²` in the attractive denominator (`w`).
    pub fn w(&self) -> f64 {
        match self {
            CubicModel::PengRobinson => -1.0,
            CubicModel::SoaveRedlichKwong => 0.0,
        }
    }
}

/// The physically desired phase when selecting a root of the cubic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    /// Liquid (smallest positive root).
    Liquid,
    /// Vapor (largest root).
    Vapor,
}

/// Solve `a z³ + b z² + c z + d = 0` for its real roots (sorted ascending).
///
/// For the cubic EoS compressibility equation the coefficients are real; this
/// returns 1 root (supercritical / single phase) or 3 roots (two-phase region,
/// where the smallest is liquid, the largest is vapor).
pub fn cubic_real_roots(a: f64, b: f64, c: f64, d: f64) -> alloc::vec::Vec<f64> {
    if a.abs() < 1e-15 {
        // Degenerate to quadratic `b z² + c z + d = 0`.
        if b.abs() < 1e-15 {
            if c.abs() < 1e-15 {
                return alloc::vec::Vec::new();
            }
            return alloc::vec![-d / c];
        }
        let disc = c * c - 4.0 * b * d;
        if disc < 0.0 {
            return alloc::vec::Vec::new();
        }
        let sq = disc.sqrt();
        return alloc::vec![(-c + sq) / (2.0 * b), (-c - sq) / (2.0 * b)];
    }
    let b = b / a;
    let c = c / a;
    let d = d / a;
    // Depressed cubic: t³ + p t + q = 0 with z = t - b/3.
    let p = c - b * b / 3.0;
    let q = 2.0 * b * b * b / 27.0 - b * c / 3.0 + d;
    let disc = q * q / 4.0 + p * p * p / 27.0;

    let shift = -b / 3.0;
    let mut roots = alloc::vec::Vec::new();
    if disc > 0.0 {
        // One real root.
        let sqrt_disc = disc.sqrt();
        let u = (-q / 2.0 + sqrt_disc).cbrt();
        let v = (-q / 2.0 - sqrt_disc).cbrt();
        roots.push(u + v + shift);
    } else {
        // Three real roots.
        let r = (-p / 3.0).max(0.0).sqrt();
        let phi = (q / (2.0 * r * r * r).max(1e-300)).clamp(-1.0, 1.0).acos();
        for k in 0..3 {
            let theta = (phi + 2.0 * core::f64::consts::PI * k as f64) / 3.0;
            roots.push(2.0 * r * theta.cos() + shift);
        }
    }
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap());
    roots
}

/// Compressibility-factor real roots for a given model's cubic, given the
/// dimensionless `A = a·α·P/(R²T²)` and `B = b·P/(R T)`.
pub fn compressibility_roots(model: CubicModel, a: f64, b: f64) -> alloc::vec::Vec<f64> {
    match model {
        CubicModel::PengRobinson => {
            // Z³ + (B-1) Z² + (A - 2B - 3B²) Z + (B² + B³ - A B) = 0
            cubic_real_roots(
                1.0,
                b - 1.0,
                a - 2.0 * b - 3.0 * b * b,
                b * b + b * b * b - a * b,
            )
        }
        CubicModel::SoaveRedlichKwong => {
            // Z³ - Z² + (A - B - B²) Z - A B = 0
            cubic_real_roots(1.0, -1.0, a - b - b * b, -a * b)
        }
    }
}

/// Select the root corresponding to `phase` from the real roots returned by
/// [`compressibility_roots`].
///
/// With a single root (supercritical / single phase) it is returned unchanged;
/// with three roots the smallest is the liquid and the largest is the vapor.
pub fn select_root(roots: &[f64], phase: Phase) -> Option<f64> {
    if roots.is_empty() {
        return None;
    }
    let mut r = roots.to_vec();
    r.sort_by(|x, y| x.partial_cmp(y).unwrap());
    match phase {
        Phase::Liquid => Some(r[0]),
        Phase::Vapor => Some(*r.last().unwrap()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_cubic_three_distinct_roots() {
        // z^3 - z = 0  →  roots {-1, 0, 1}.
        let r = cubic_real_roots(1.0, 0.0, -1.0, 0.0);
        assert_eq!(r.len(), 3);
        assert!((r[0] + 1.0).abs() < 1e-9);
        assert!(r[1].abs() < 1e-9);
        assert!((r[2] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn roots_satisfy_pr_cubic() {
        let (a, b) = (0.6_f64, 0.08_f64);
        let roots = compressibility_roots(CubicModel::PengRobinson, a, b);
        assert!(!roots.is_empty());
        for z in &roots {
            // z^3 + (b-1) z^2 + (a - 2b - 3b^2) z + (b^2 + b^3 - a b) = 0.
            let f = z.powi(3)
                + (b - 1.0) * z.powi(2)
                + (a - 2.0 * b - 3.0 * b * b) * z
                + (b * b + b * b * b - a * b);
            assert!(f.abs() < 1e-9, "root {z} leaves residual {f}");
        }
    }

    #[test]
    fn roots_satisfy_srk_cubic() {
        let (a, b) = (0.6_f64, 0.08_f64);
        let roots = compressibility_roots(CubicModel::SoaveRedlichKwong, a, b);
        assert!(!roots.is_empty());
        for z in &roots {
            // z^3 - z^2 + (a - b - b^2) z - a b = 0.
            let f = z.powi(3) - z.powi(2) + (a - b - b * b) * z - a * b;
            assert!(f.abs() < 1e-9, "root {z} leaves residual {f}");
        }
    }

    #[test]
    fn select_root_picks_extremes() {
        let roots = vec![0.1, 0.5, 0.9];
        assert!((select_root(&roots, Phase::Liquid).unwrap() - 0.1).abs() < 1e-12);
        assert!((select_root(&roots, Phase::Vapor).unwrap() - 0.9).abs() < 1e-12);
        assert!(select_root(&[], Phase::Liquid).is_none());
    }
}
