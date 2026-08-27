# tpt-thermodynamics — task recipes.
#
# All recipes shell out to `cargo xtask` (see xtask/ and .cargo/config.toml),
# so the logic lives in one place. Requires https://github.com/casey/just.

default:
    @just --list

# Format the whole workspace
fmt:
    cargo xtask fmt

# Clippy with -D warnings, all features
clippy:
    cargo xtask clippy

# Tests + doctests, all features
test:
    cargo xtask test

# cargo-deny (advisories / bans / licenses / sources)
deny:
    cargo xtask deny

# Cross-check non-xtask crates for wasm32-unknown-unknown
wasm:
    cargo xtask wasm

# Fast compile check, all features
check:
    cargo xtask check

# Scaffold a new tpt-thermo-* crate
new-crate name:
    cargo xtask new-crate {{name}}

# Everything above, in order
all:
    cargo xtask all
