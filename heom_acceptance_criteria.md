# HEOM Acceptance Criteria

This document pre-registers acceptance criteria for HEOM validation in the Quantubulin v3.5.0 release.

## Scope

The HEOM layer is used as a nonperturbative validation ledger for short-window convergence and model-comparison diagnostics. It does not replace staged long-window full-system production trajectories.

## Acceptance criteria

1. Observable-level hierarchy-depth contraction must be monotonic or asymptotically decreasing across the tested NC ladder.
2. Nk refinement must be subdominant relative to NC refinement in the tested regime.
3. Full-system 1JFF two-depth entries are treated as ratio calibrations, not independently identified curve fits.
4. HEOM--Redfield discrepancies are interpreted as model-level differences only when larger than projected truncation uncertainty.
5. Long-window NC=7 production trajectories remain staged extensions unless explicitly generated and audited.
6. Bayesian HEOM contraction summaries must use log-scale positive-jump observables and hierarchical shrinkage; raw-level weakly identified MCMC fits are not accepted as inferential evidence.

## Current release status

Quantubulin v3.5.0 satisfies the short-window validation-ledger criteria and reports long-window full-system production trajectories as staged computational closure.
