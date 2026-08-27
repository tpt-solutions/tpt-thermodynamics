//! Minimal dense linear-algebra helpers used by the solvers in this crate.

use alloc::vec;
use alloc::vec::Vec;

/// Solve `A x = b` for a square `n×n` system via Gaussian elimination with
/// partial pivoting. Returns `None` if the matrix is (near) singular.
pub fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    if n == 0 || a.len() != n || a.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut m = vec![vec![0.0_f64; n + 1]; n];
    for i in 0..n {
        m[i][..n].copy_from_slice(&a[i]);
        m[i][n] = b[i];
    }
    for col in 0..n {
        let mut piv = col;
        let mut max = m[col][col].abs();
        for r in (col + 1)..n {
            let val = m[r][col].abs();
            if val > max {
                max = val;
                piv = r;
            }
        }
        if max < 1e-30 {
            return None;
        }
        m.swap(col, piv);
        let d = m[col][col];
        for j in col..=n {
            m[col][j] /= d;
        }
        for r in 0..n {
            if r != col {
                let f = m[r][col];
                if f != 0.0 {
                    for j in col..=n {
                        m[r][j] -= f * m[col][j];
                    }
                }
            }
        }
    }
    Some((0..n).map(|i| m[i][n]).collect())
}
