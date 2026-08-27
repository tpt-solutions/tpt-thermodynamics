//! The embedded curated seed dataset.

/// The curated seed component dataset, embedded at compile time from
/// `data/seed.toml`.
pub const SEED_TOML: &str = include_str!("../data/seed.toml");
