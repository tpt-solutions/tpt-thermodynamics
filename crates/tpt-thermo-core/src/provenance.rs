//! Provenance for parameters and dataset values: where a number came from, so
//! every seeded/curated quantity can be traced.

use alloc::string::String;

/// A calendar date without a timezone (kept `no_std`-friendly: no `chrono`
/// dependency in the core crate).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SourceDate {
    /// Year (CE).
    pub year: u16,
    /// Month, 1-12.
    pub month: u8,
    /// Day of month, 1-31.
    pub day: u8,
}

/// Where a parameter or value originated.
#[derive(Debug, Clone, PartialEq)]
pub struct ParameterSource {
    /// Human-readable citation (e.g. "DIPPR 801", "Poling et al. 2001").
    pub description: String,
    /// Optional URL to the source.
    pub url: Option<String>,
    /// Optional retrieval date.
    pub retrieved: Option<SourceDate>,
}

impl ParameterSource {
    /// Construct with just a description.
    pub fn new(description: &str) -> Self {
        Self {
            description: String::from(description),
            url: None,
            retrieved: None,
        }
    }

    /// Attach a URL.
    pub fn with_url(mut self, url: &str) -> Self {
        self.url = Some(String::from(url));
        self
    }

    /// Attach a retrieval date.
    pub fn with_date(mut self, date: SourceDate) -> Self {
        self.retrieved = Some(date);
        self
    }
}

/// Provenance wrapper attaching a [`ParameterSource`] to a value.
#[derive(Debug, Clone, PartialEq)]
pub struct Provenance<T> {
    /// The value.
    pub value: T,
    /// Its source.
    pub source: ParameterSource,
    /// Optional free-text notes.
    pub notes: Option<String>,
}

impl<T> Provenance<T> {
    /// Wrap `value` with `source`.
    pub fn new(value: T, source: ParameterSource) -> Self {
        Self {
            value,
            source,
            notes: None,
        }
    }

    /// Attach a note.
    pub fn with_notes(mut self, notes: &str) -> Self {
        self.notes = Some(String::from(notes));
        self
    }
}

/// A binary interaction parameter `k_ij` (or any BIP-shaped value) with
/// provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct BipParameter {
    /// First component index.
    pub i: usize,
    /// Second component index.
    pub j: usize,
    /// The dimensionless interaction value.
    pub value: f64,
    /// Provenance of the value.
    pub source: ParameterSource,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn provenance_round_trip() {
        let src = ParameterSource::new("DIPPR 801")
            .with_url("https://example.org/dippr")
            .with_date(SourceDate {
                year: 2024,
                month: 1,
                day: 15,
            });
        let p = Provenance::new(0.95_f64, src.clone()).with_notes("estimated");
        assert_eq!(p.value, 0.95);
        assert_eq!(p.source.description, "DIPPR 801");
        assert_eq!(p.source.retrieved.unwrap().year, 2024);
    }

    #[test]
    fn bip_stores_indices() {
        let bip = BipParameter {
            i: 0,
            j: 1,
            value: 0.01,
            source: ParameterSource::new("fit"),
        };
        assert_eq!((bip.i, bip.j, bip.value), (0, 1, 0.01));
    }
}
