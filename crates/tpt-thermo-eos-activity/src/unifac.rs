//! UNIFAC — UNIQUalquac Functional-group Activity Coefficients (Fredenslund et
//! al., 1975), original and Dortmund-modified variants, with a **seed** group
//! table.
//!
//! Group-contribution activity coefficients:
//!
//! ```text
//! R_i = Σ_k ν_ik R_k,   Q_i = Σ_k ν_ik Q_k
//! X_k = (Σ_i x_i ν_ik) / (Σ_i x_i Σ_m ν_im)            (group mole fraction)
//! Θ_k = Q_k X_k / (Σ_m Q_m X_m)                        (group surface fraction)
//! ψ_mn = exp(−a_mn(T)/T)                               (a_mn in K)
//! Γ_k = Q_k[ 1 − ln(Σ_m Θ_m ψ_mk) − Σ_m Θ_m ψ_km / (Σ_n Θ_n ψ_nm) ]
//! ln γ_i^R = Σ_k ν_ik [ ln Γ_k − ln Γ_k^(i) ]         (Γ_k^(i): pure-i group fractions)
//! ln γ_i^C = 1 − V_i + ln V_i − 5 Q_i[ 1 − V_i/F_i + ln(V_i/F_i) ]
//!            V_i = R_i / Σ_j x_j R_j,  F_i = Q_i / Σ_j x_j Q_j
//! ln γ_i = ln γ_i^C + ln γ_i^R
//! ```
//!
//! The **Dortmund** modification (Gmehling) extends `a_mn(T)` with the full
//! `a + b·T + c·T² + d/T + e·ln T` temperature dependence and a different
//! parameter set; here the same machinery is reused with a
//! [`TdParam`](crate::parameters::TdParam) interaction matrix (constant `a` for
//! the Original variant). The shipped group table is a *seed* subset (see
//! `seed_group_table`); full group coverage is Deferred Scope.

use crate::parameters::{self, TdParam};
use alloc::vec::Vec;
use tpt_thermo_core::error::ThermoError;
use tpt_thermo_core::mixing::ExcessGibbsModel;
use tpt_thermo_core::quantities::{Pressure, Temperature};

/// Which UNIFAC parametrisation the interaction matrix represents.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnifacVariant {
    /// Original UNIFAC (constant `a_mn`).
    Original,
    /// Dortmund-modified UNIFAC (temperature-dependent `a_mn(T)`).
    Dortmund,
}

/// A UNIFAC group-parameter table: per-group `R_k`, `Q_k`, and the `a_mn(T)`
/// interaction matrix between groups.
#[derive(Debug, Clone)]
pub struct GroupTable {
    /// `R_k` volume parameter per group.
    pub r: Vec<f64>,
    /// `Q_k` surface-area parameter per group.
    pub q: Vec<f64>,
    /// `a_mn(T)` interaction parameters (K) between groups, full `n_g×n_g`.
    pub interaction: Vec<Vec<TdParam>>,
}

impl GroupTable {
    /// Build from `R_k`, `Q_k` and a full interaction matrix.
    pub fn new(
        r: Vec<f64>,
        q: Vec<f64>,
        interaction: Vec<Vec<TdParam>>,
    ) -> Result<Self, ThermoError> {
        let ng = r.len();
        if q.len() != ng || interaction.len() != ng || interaction.iter().any(|row| row.len() != ng)
        {
            return Err(ThermoError::InvalidInput(
                "UNIFAC group table dimension mismatch",
            ));
        }
        Ok(Self { r, q, interaction })
    }

    /// Number of groups.
    pub fn num_groups(&self) -> usize {
        self.r.len()
    }

    fn psi(&self, m: usize, n: usize, t: Temperature) -> f64 {
        let tk = parameters::tk(t);
        libm::exp(-self.interaction[m][n].value(t) / tk)
    }
}

/// A UNIFAC model for `n` components, each described by a list of
/// `(group_index, count)` contributions.
#[derive(Debug, Clone)]
pub struct UnifacModel {
    n: usize,
    /// `groups[i]` = `[(group_index, count), …]` for component `i`.
    groups: Vec<Vec<(usize, f64)>>,
    table: GroupTable,
    variant: UnifacVariant,
}

impl UnifacModel {
    /// Build for `n` components given per-component group counts and a table.
    pub fn new(
        n: usize,
        groups: Vec<Vec<(usize, f64)>>,
        table: GroupTable,
        variant: UnifacVariant,
    ) -> Result<Self, ThermoError> {
        if groups.len() != n {
            return Err(ThermoError::InvalidInput(
                "UNIFAC component/group count mismatch",
            ));
        }
        let ng = table.num_groups();
        for (i, g) in groups.iter().enumerate() {
            for &(gi, count) in g {
                if gi >= ng {
                    return Err(ThermoError::InvalidInput("UNIFAC group index out of range"));
                }
                if count < 0.0 {
                    return Err(ThermoError::InvalidInput("UNIFAC negative group count"));
                }
                let _ = i;
            }
        }
        Ok(Self {
            n,
            groups,
            table,
            variant,
        })
    }

    /// The chosen variant.
    pub fn variant(&self) -> UnifacVariant {
        self.variant
    }

    /// `R_i` for component `i`.
    pub fn r_i(&self, i: usize) -> f64 {
        self.groups[i]
            .iter()
            .map(|&(g, c)| c * self.table.r[g])
            .sum()
    }

    /// `Q_i` for component `i`.
    pub fn q_i(&self, i: usize) -> f64 {
        self.groups[i]
            .iter()
            .map(|&(g, c)| c * self.table.q[g])
            .sum()
    }

    fn group_fractions(&self, x: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let ng = self.table.num_groups();
        let mut num = vec![0.0f64; ng];
        let mut denom = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            let mut tot = 0.0;
            for &(g, c) in &self.groups[i] {
                num[g] += xi * c;
                tot += c;
            }
            denom += xi * tot;
        }
        let mut xg = vec![0.0f64; ng];
        for g in 0..ng {
            xg[g] = if denom > 0.0 { num[g] / denom } else { 0.0 };
        }
        let mut theta = vec![0.0f64; ng];
        let mut qs = 0.0;
        for g in 0..ng {
            theta[g] = self.table.q[g] * xg[g];
            qs += theta[g];
        }
        for g in 0..ng {
            theta[g] /= qs;
        }
        (xg, theta)
    }

    fn group_fractions_pure(&self, i: usize) -> (Vec<f64>, Vec<f64>) {
        let ng = self.table.num_groups();
        let mut xg = vec![0.0f64; ng];
        let mut tot = 0.0;
        for &(g, c) in &self.groups[i] {
            xg[g] = c;
            tot += c;
        }
        for g in 0..ng {
            xg[g] /= tot;
        }
        let mut theta = vec![0.0f64; ng];
        let mut qs = 0.0;
        for g in 0..ng {
            theta[g] = self.table.q[g] * xg[g];
            qs += theta[g];
        }
        for g in 0..ng {
            theta[g] /= qs;
        }
        (xg, theta)
    }

    fn gamma_k(&self, theta: &[f64], t: Temperature) -> Vec<f64> {
        let ng = theta.len();
        let mut gamma = vec![0.0f64; ng];
        for k in 0..ng {
            let mut s1 = 0.0;
            for m in 0..ng {
                s1 += theta[m] * self.table.psi(m, k, t);
            }
            let mut s2 = 0.0;
            for m in 0..ng {
                let mut denom = 0.0;
                for n in 0..ng {
                    denom += theta[n] * self.table.psi(n, m, t);
                }
                s2 += theta[m] * self.table.psi(k, m, t) / denom;
            }
            gamma[k] = self.table.q[k] * (1.0 - libm::log(s1) - s2);
        }
        gamma
    }

    /// Reduced excess Gibbs energy `g^E/(R T) = Σ_i x_i ln γ_i`.
    pub fn reduced_excess_gibbs_at(&self, t: Temperature, x: &[f64]) -> Result<f64, ThermoError> {
        parameters::check_composition(x)?;
        if x.len() != self.n {
            return Err(ThermoError::InvalidInput("composition length mismatch"));
        }
        let mut total = 0.0;
        for i in 0..self.n {
            total += x[i] * self.ln_gamma_at(t, x, i)?;
        }
        Ok(total)
    }

    /// Natural log of the activity coefficient of component `i`.
    pub fn ln_gamma_at(&self, t: Temperature, x: &[f64], i: usize) -> Result<f64, ThermoError> {
        parameters::check_composition(x)?;
        if x.len() != self.n || i >= self.n {
            return Err(ThermoError::InvalidInput("composition/ index mismatch"));
        }
        let r_i = self.r_i(i);
        let q_i = self.q_i(i);

        // Combinatorial term.
        let mut rs = 0.0;
        let mut qs = 0.0;
        for j in 0..self.n {
            rs += x[j] * self.r_i(j);
            qs += x[j] * self.q_i(j);
        }
        let v_i = r_i / rs;
        let f_i = q_i / qs;
        let comb =
            1.0 - v_i + libm::log(v_i) - 5.0 * q_i * (1.0 - v_i / f_i + libm::log(v_i / f_i));

        // Residual term.
        let (_xg, theta) = self.group_fractions(x);
        let gamma = self.gamma_k(&theta, t);
        let (_xp, theta_p) = self.group_fractions_pure(i);
        let gamma_p = self.gamma_k(&theta_p, t);
        let mut resid = 0.0;
        for &(g, c) in &self.groups[i] {
            resid += c * (gamma[g] - gamma_p[g]);
        }

        Ok(comb + resid)
    }

    /// Activity coefficient (not log) of component `i`.
    pub fn gamma_at(&self, t: Temperature, x: &[f64], i: usize) -> Result<f64, ThermoError> {
        Ok(libm::exp(self.ln_gamma_at(t, x, i)?))
    }
}

impl ExcessGibbsModel for UnifacModel {
    fn num_components(&self) -> usize {
        self.n
    }

    fn reduced_excess_gibbs(
        &self,
        t: Temperature,
        _p: Pressure,
        x: &[f64],
    ) -> Result<f64, ThermoError> {
        self.reduced_excess_gibbs_at(t, x)
    }

    fn ln_gamma(
        &self,
        t: Temperature,
        _p: Pressure,
        x: &[f64],
        i: usize,
    ) -> Result<f64, ThermoError> {
        self.ln_gamma_at(t, x, i)
    }
}

/// A minimal seed UNIFAC group table (groups with `R_k`/`Q_k` and a small set of
/// seeded main-group interactions). **Illustrative only** — full UNIFAC/Dortmund
/// group coverage is Deferred Scope.
pub fn seed_group_table() -> GroupTable {
    // Indices: 0 CH3, 1 CH2, 2 ACH, 3 OH, 4 H2O, 5 CH3OH, 6 COOH, 7 ACCH2.
    let r = vec![
        0.9011, 0.6744, 0.5313, 1.0000, 0.9200, 1.4311, 1.3013, 1.0396,
    ];
    let q = vec![
        0.8480, 0.5400, 0.4000, 1.2000, 1.4000, 1.4320, 1.2240, 0.6600,
    ];
    let mut inter = vec![vec![TdParam::default(); r.len()]; r.len()];
    // Helper to set a symmetric-ish pair (asymmetric a_mn, a_nm).
    let set = |inter: &mut Vec<Vec<TdParam>>, m: usize, n: usize, amn: f64, anm: f64| {
        inter[m][n] = TdParam::new(amn, 0.0, 0.0);
        inter[n][m] = TdParam::new(anm, 0.0, 0.0);
    };
    // Alkane (CH2) — water.
    set(&mut inter, 1, 4, 1318.7, 300.9);
    // Alkane (CH2) — OH.
    set(&mut inter, 1, 3, 986.5, 156.4);
    // Water — OH (alcohol).
    set(&mut inter, 4, 3, -229.1, 119.7);
    // Aromatic (ACH) — water.
    set(&mut inter, 2, 4, 657.0, 94.4);
    // Aromatic (ACH) — OH.
    set(&mut inter, 2, 3, 537.0, -184.0);
    // Aromatic (ACH) — alkane (CH2).
    set(&mut inter, 2, 1, -24.4, -110.0);
    // Methanol (CH3OH) — water.
    set(&mut inter, 5, 4, 292.8, 135.5);
    // Carboxyl (COOH) — water.
    set(&mut inter, 6, 4, 613.0, -144.0);
    // Carboxyl (COOH) — alkane (CH2).
    set(&mut inter, 6, 1, 339.5, 127.0);
    GroupTable::new(r, q, inter).unwrap()
}

/// Build a few predefined molecules' group lists for tests/examples.
pub fn seed_molecules() -> Vec<Vec<(usize, f64)>> {
    // [CH3, CH2, ACH, OH, H2O, CH3OH, COOH, ACCH2]
    vec![
        vec![(0, 2.0), (1, 4.0)],           // n-hexane
        vec![(0, 1.0), (1, 1.0), (3, 1.0)], // ethanol
        vec![(2, 6.0)],                     // benzene
        vec![(2, 5.0), (0, 1.0)],           // toluene
        vec![(4, 1.0)],                     // water
        vec![(5, 1.0)],                     // methanol
        vec![(0, 1.0), (6, 1.0)],           // acetic acid (CH3 + COOH)
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::thermodynamic_temperature::kelvin;

    #[test]
    fn pure_component_gives_zero_ln_gamma() {
        let table = seed_group_table();
        let mols = seed_molecules();
        let m = UnifacModel::new(7, mols, table, UnifacVariant::Original).unwrap();
        let t = Temperature::new::<kelvin>(298.15);
        for i in 0..7 {
            let mut x = vec![0.0f64; 7];
            x[i] = 1.0;
            let lg = m.ln_gamma_at(t, &x, i).unwrap();
            assert!(lg.abs() < 1e-9, "component {i} ln γ = {lg}");
        }
    }

    #[test]
    fn gibbs_duhem_consistency_ternary() {
        // Identity check on a 3-component mixture: g^E/(R T) == Σ_i x_i ln γ_i.
        let table = seed_group_table();
        let mols = seed_molecules();
        let m = UnifacModel::new(7, mols, table, UnifacVariant::Original).unwrap();
        let t = Temperature::new::<kelvin>(298.15);
        let ids = [0usize, 1, 4];
        let x = {
            let mut xx = vec![0.0f64; 7];
            xx[ids[0]] = 0.4;
            xx[ids[1]] = 0.3;
            xx[ids[2]] = 0.3;
            xx
        };
        let ge = m.reduced_excess_gibbs_at(t, &x).unwrap();
        let gd: f64 = ids
            .iter()
            .map(|&i| x[i] * m.ln_gamma_at(t, &x, i).unwrap())
            .sum();
        assert!(
            (ge - gd).abs() < 1e-9,
            "g^E/RT != Σ x lnγ: ge={ge}, gd={gd}"
        );
    }

    #[test]
    fn reduced_excess_gibbs_equals_sum_x_lngamma() {
        let table = seed_group_table();
        let mols = seed_molecules();
        let m = UnifacModel::new(7, mols, table, UnifacVariant::Original).unwrap();
        let t = Temperature::new::<kelvin>(333.15);
        let mut x = vec![0.0f64; 7];
        x[0] = 0.5;
        x[1] = 0.3;
        x[4] = 0.2;
        let ge = m.reduced_excess_gibbs_at(t, &x).unwrap();
        let sum = x[0] * m.ln_gamma_at(t, &x, 0).unwrap()
            + x[1] * m.ln_gamma_at(t, &x, 1).unwrap()
            + x[4] * m.ln_gamma_at(t, &x, 4).unwrap();
        assert!((ge - sum).abs() < 1e-12);
    }
}
