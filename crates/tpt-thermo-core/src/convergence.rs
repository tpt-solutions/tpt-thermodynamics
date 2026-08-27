//! Convergence and numerical-issue reporting shared by every solver in the
//! workspace.

/// Outcome of an iterative solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvergenceStatus {
    /// The iteration met its tolerance.
    Converged,
    /// The iteration left the basin of attraction / failed to tighten.
    Diverged(DivergenceReason),
    /// The iteration ran out of budget without meeting tolerance.
    NotConverged,
    /// A numerical precondition failed (singular Jacobian, out-of-domain
    /// evaluation, …) before divergence per se.
    NumericalIssue(NumericalIssueReason),
}

/// Why an iterative solve diverged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DivergenceReason {
    /// Exceeded the maximum iteration count.
    MaxIterations,
    /// Successive iterates oscillated without damping.
    Oscillation,
    /// A non-finite (NaN/inf) value was produced.
    NonFinite,
    /// The damped step collapsed to (near) zero without converging.
    StepTooSmall,
}

/// A numerical precondition that failed during a solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NumericalIssueReason {
    /// A linear solve encountered a singular (or near-singular) Jacobian.
    SingularJacobian,
    /// An evaluation was requested outside the model's valid domain.
    OutOfDomain,
    /// A pressure/volume/composition became non-physical (e.g. negative).
    NonPhysical,
    /// A required parameter was missing.
    MissingParameter,
}
