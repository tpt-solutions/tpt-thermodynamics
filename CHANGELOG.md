# Changelog

All notable changes to this workspace are documented here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/) once crates are published.

> Per the project decision recorded in `todo.md`, **no crates in this repo are
> published to crates.io as part of the build-out** — versioning is tracked
> locally until the user chooses to release.

## [Unreleased]

### Added
- Phase 0 repository bootstrap:
  - Workspace `Cargo.toml` (`resolver = "2"`, `[workspace.package]`,
    `[workspace.dependencies]` seeded with upstream version-string deps,
    `[workspace.lints]`, `[profile.release]`).
  - `rust-toolchain.toml` (stable + rustfmt/clippy + `wasm32-unknown-unknown`),
    `rustfmt.toml`, `deny.toml`, `.cargo/config.toml` (`xtask` alias),
    `.gitignore`.
  - `xtask/` task runner (`fmt`, `clippy`, `test`, `deny`, `wasm`,
    `check`, `new-crate`, `all`).
  - `justfile` mirroring the `cargo xtask` commands.
  - `LICENSE-MIT`, `LICENSE-APACHE` (dual `MIT OR Apache-2.0`).
  - `README.md`, `CONTRIBUTING.md`, `SECURITY.md`, this `CHANGELOG.md`.
  - `.github/workflows/ci.yml` (self-contained fmt / clippy / test+doc / wasm /
    deny / msrv) and `.github/dependabot.yml`.
  - `examples/` crate scaffold (populated per phase).
