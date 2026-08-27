# Security Policy

## Supported versions

This repository is under active development (Phases 0-13 per `todo.md`). Only
the `master` branch is supported. Security fixes land on `master` and are
consumed downstream via crates.io releases managed by the maintainers.

## Reporting a vulnerability

**Please do not open a public GitHub issue for security vulnerabilities.**

Report privately via GitHub's security advisory feature
(`Security` → `Report a vulnerability` on the repository) or by contacting the
maintainers at the address published in the repository's `CODEOWNERS` /
maintainer metadata.

Include:

- A description of the vulnerability and its impact.
- A minimal reproducer or proof-of-concept.
- Affected crate(s) and version(s) / commit(s).
- If known, a suggested mitigation.

You can expect an acknowledgement within **5 business days**. We will keep you
informed of progress and will credit you in the advisory unless you request
anonymity.

## Scope notes

This is a numerical/scientific computing library. Typical risk classes of
interest:

- Numerical instability or panics that can be triggered by untrusted inputs
  passed to EoS / flash / phase-stability APIs.
- Unsound `unsafe` code (note: `unsafe_code` is forbidden workspace-wide).
- Denial-of-service via unbounded iteration / non-terminating solvers without
  the documented `ConvergenceStatus` / iteration-limit guards.
- Cryptographic or supply-chain issues in the build (see `deny.toml`).

Out of scope: physical/modelling accuracy of thermodynamic predictions (track
those as normal issues), and downstream application misuse.
