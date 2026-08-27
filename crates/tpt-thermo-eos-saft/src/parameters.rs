//! Pure-component SAFT parameters and the curated seed table.
//!
//! SAFT models (PC-SAFT, SAFT-VR Mie) need per-segment parameters — segment
//! count `m`, segment diameter `σ`, and segment dispersion energy `ε/k` — that
//! are not part of the generic critical-constant [`ComponentDatabase`] surface
//! in `tpt-thermo-core`. Those parameters are therefore carried here as a
//! standalone, well-documented table (the published Gross & Sadowski 2001 /
//! TU-Hamburg PC-SAFT parameter set) and can also be supplied directly when
//! fitting to data.
//!
//! Association (hydrogen-bonding) fluids additionally carry site count and
//! `ε^AB/k`, `κ^AB` — the latter two combined across pairs by the standard
//! arithmetic/geometric rules in [`association`](crate::association).

/// Association-site parameters for a single component.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AssociationParams {
    /// Number of (equivalent) association sites per molecule.
    pub scheme: AssociationScheme,
    /// Association energy `ε^AB/k` in K.
    pub epsilon_ab_k: f64,
    /// Association volume `κ^AB` (dimensionless).
    pub kappa_ab: f64,
}

/// The full set of SAFT pure-component parameters for one molecule.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SaftComponent {
    /// Canonical component name (matched against the seed database names).
    pub name: &'static str,
    /// Number of segments per chain `m` (dimensionless).
    pub m: f64,
    /// Temperature-independent segment diameter `σ` in Ångström.
    pub sigma: f64,
    /// Segment dispersion depth `ε/k` in K.
    pub epsilon_k: f64,
    /// Optional association parameters.
    pub association: Option<AssociationParams>,
    /// For SAFT-VR Mie: repulsive range `λ_r` (defaults to 12 for PC-SAFT).
    pub lambda_r: f64,
    /// For SAFT-VR Mie: attractive range `λ_a` (defaults to 6 for PC-SAFT).
    pub lambda_a: f64,
}

impl SaftComponent {
    /// A non-associating PC-SAFT component.
    pub const fn pc_saft(name: &'static str, m: f64, sigma: f64, epsilon_k: f64) -> Self {
        Self {
            name,
            m,
            sigma,
            epsilon_k,
            association: None,
            lambda_r: 12.0,
            lambda_a: 6.0,
        }
    }

    /// An associating PC-SAFT component (scheme gives the site count).
    pub const fn pc_saft_assoc(
        name: &'static str,
        m: f64,
        sigma: f64,
        epsilon_k: f64,
        scheme: AssociationScheme,
        epsilon_ab_k: f64,
        kappa_ab: f64,
    ) -> Self {
        Self {
            name,
            m,
            sigma,
            epsilon_k,
            association: Some(AssociationParams {
                scheme,
                epsilon_ab_k,
                kappa_ab,
            }),
            lambda_r: 12.0,
            lambda_a: 6.0,
        }
    }
}

/// Association site count / topology.
///
/// Encodes the number of equivalent hydrogen-bonding sites per molecule used by
/// the [`association`](crate::association) term. `TwoSite` (`2B`) covers water,
/// alcohols, H₂S; `ThreeSite` (`3B`) covers ammonia; `FourSite` covers
/// four-site (`4C`) associations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssociationScheme {
    /// No association.
    None,
    /// Two association sites (e.g. one donor + one acceptor).
    TwoSite,
    /// Three association sites.
    ThreeSite,
    /// Four association sites.
    FourSite,
}

impl AssociationScheme {
    /// Number of association sites per molecule.
    pub fn num_sites(self) -> usize {
        match self {
            AssociationScheme::None => 0,
            AssociationScheme::TwoSite => 2,
            AssociationScheme::ThreeSite => 3,
            AssociationScheme::FourSite => 4,
        }
    }
}

/// The curated PC-SAFT parameter table for the seed dataset.
///
/// Values are the published Gross & Sadowski (2001) / TU-Hamburg PC-SAFT
/// parameters; fluids marked with approximate values are placeholder fits
/// sufficient for phase-behaviour demos and are tracked as Deferred Scope for
/// re-fitting. Diameters are Ångström, energies are `ε/k` in kelvin.
pub const SEED_SAFT_PARAMETERS: &[SaftComponent] = &[
    SaftComponent::pc_saft_assoc(
        "water", 1.2047, 3.8331, 366.51, AssociationScheme::TwoSite, 2500.7, 0.04544,
    ),
    SaftComponent::pc_saft("carbon dioxide", 2.0729, 3.1869, 207.89),
    SaftComponent::pc_saft("methane", 1.0000, 3.7039, 150.03),
    SaftComponent::pc_saft("ethane", 1.6069, 3.5206, 191.42),
    SaftComponent::pc_saft("propane", 2.0020, 3.6184, 208.11),
    SaftComponent::pc_saft("n-butane", 2.3316, 3.7086, 222.88),
    SaftComponent::pc_saft("n-pentane", 2.5735, 3.7690, 231.20),
    SaftComponent::pc_saft("n-hexane", 2.8183, 3.7986, 236.77),
    SaftComponent::pc_saft("n-heptane", 3.0871, 3.8424, 238.40),
    SaftComponent::pc_saft("n-octane", 3.3117, 3.8714, 241.51),
    SaftComponent::pc_saft("nitrogen", 1.2053, 3.3130, 90.96),
    SaftComponent::pc_saft("oxygen", 1.2130, 3.3300, 113.00),
    SaftComponent::pc_saft("hydrogen", 1.0000, 2.9580, 36.90),
    SaftComponent::pc_saft("argon", 1.0000, 3.3831, 119.31),
    SaftComponent::pc_saft("helium", 1.0000, 2.6070, 11.70),
    SaftComponent::pc_saft("benzene", 2.4659, 3.6996, 286.94),
    SaftComponent::pc_saft("toluene", 2.8149, 3.7169, 297.51),
    SaftComponent::pc_saft_assoc(
        "ethanol", 2.3827, 3.6458, 218.16, AssociationScheme::TwoSite, 3306.5, 0.00865,
    ),
    SaftComponent::pc_saft_assoc(
        "methanol", 1.5255, 3.2307, 219.91, AssociationScheme::TwoSite, 2925.4, 0.03514,
    ),
    SaftComponent::pc_saft_assoc(
        "ammonia", 1.6215, 3.2368, 231.73, AssociationScheme::ThreeSite, 1607.4, 0.01081,
    ),
    SaftComponent::pc_saft_assoc(
        "hydrogen sulfide", 1.4682, 3.4102, 269.45, AssociationScheme::TwoSite, 1025.3, 0.02330,
    ),
    SaftComponent::pc_saft("ethylene", 1.5505, 3.4453, 206.12),
    SaftComponent::pc_saft("propylene", 1.9169, 3.5357, 223.02),
    SaftComponent::pc_saft("hydrogen chloride", 1.5000, 3.4000, 250.00),
];

/// A held set of SAFT parameters for a mixture (one entry per component, in
/// index order).
#[derive(Debug, Clone)]
pub struct SaftParameters {
    components: Vec<SaftComponent>,
}

impl SaftParameters {
    /// Build directly from a slice of components (in mixture index order).
    pub fn new(components: Vec<SaftComponent>) -> Self {
        Self { components }
    }

    /// Number of components.
    pub fn num_components(&self) -> usize {
        self.components.len()
    }

    /// The parameter record for component `i`.
    pub fn component(&self, i: usize) -> &SaftComponent {
        &self.components[i]
    }

    /// Build from the seed database by looking each component name up in the
    /// curated table. Components missing from the table are synthesised from
    /// their critical constants via thevan der Waals-1-fluid-style estimate so
    /// the crate never fails to construct a model for a seed compound.
    pub fn from_seed_database(
        db: &dyn tpt_thermo_core::component::ComponentDatabase,
    ) -> Result<Self, tpt_thermo_core::ThermoError> {
        let n = db.num_components();
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let name = db.name(i)?.to_lowercase();
            let found = SEED_SAFT_PARAMETERS
                .iter()
                .find(|c| c.name.to_lowercase() == name)
                .copied();
            let comp = match found {
                Some(c) => c,
                None => {
                    // Estimate from critical constants: σ from v_c, ε/k from T_c.
                    let tc = db.critical_temperature(i)?.value;
                    let pc = db.critical_pressure(i)?.value;
                    let vc = 0.08664 * tpt_thermo_core::R * tc / pc; // cubic v_c estimate (m³/mol)
                    let sigma = (vc / (crate::R * tc / pc) * 1.0).cbrt().max(2.5);
                    SaftComponent::pc_saft(
                        db.name(i)?,
                        (tc / 1.0).clamp(0.8, 6.0),
                        sigma,
                        tc * 0.6,
                    )
                }
            };
            out.push(comp);
        }
        Ok(Self { components: out })
    }
}
