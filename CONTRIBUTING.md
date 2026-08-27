# Contributing

Thanks for your interest in `tpt-thermodynamics`. This document covers the
workflow conventions used across the `tpt-*` family of repositories.

## Getting started

1. Install the pinned toolchain via `rustup` (see `rust-toolchain.toml`):

   ```bash
   rustup show   # picks up rust-toolchain.toml automatically
   ```

2. Build and run the workspace task runner:

   ```bash
   cargo xtask all   # fmt, clippy, test, deny, wasm, check
   ```

   `just` recipes mirror these (`just --list`).

## Conventions

- **Workspace lints are authoritative.** Every crate inherits
  `[lints] workspace = true`. `unsafe_code` is forbidden workspace-wide;
  `clippy::all` warns. Keep `clippy --all-targets --all-features -- -D warnings`
  clean.
- **Formatting:** `cargo fmt --all` before every commit (CI enforces
  `--check`).
- **Licensing:** every crate is `MIT OR Apache-2.0`. Header new files with the
  SPDX-style dual license noted in `README.md`. Do not introduce GPL/AGPL or
  other copyleft dependencies (see `deny.toml`).
- **Dependencies:** upstream `tpt-*` crates are crates.io version strings, never
  path/git deps across repo boundaries. New in-workspace crates are added to
  `[workspace.dependencies]` and registered in the `members` array (use
  `cargo xtask new-crate <name>` to scaffold).
- **`no_std`:** only `tpt-thermo-core` is `no_std` + `alloc`. `Vec`-returning
  APIs there are gated behind `alloc`. Other crates may assume `std`.
- **Tests & docs:** every public item needs a doctest or unit test and rustdoc.
  Validation targets are against the curated seed dataset described in
  `todo.md` "Known Deferred Scope".

## Dependency hygiene

- Run `cargo xtask deny` (cargo-deny) locally; CI runs the same.
- Dependabot opens grouped weekly patch/minor PRs; security advisories surface
  via cargo-deny in CI.

## Pull requests

- Keep changes focused; reference the relevant `todo.md` phase.
- Ensure `cargo xtask all` is green locally before opening a PR.
- CI must pass (fmt / clippy / test+doc / wasm / deny / msrv).

## Reporting issues

See `SECURITY.md` for vulnerability reporting. For general bugs, open an issue
with a minimal reproducer and the `rust-toolchain` version.
