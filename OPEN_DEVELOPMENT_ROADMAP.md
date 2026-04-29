# Open Development Roadmap (JOSS readiness)

This roadmap is the public maintenance plan used to sustain JOSS-ready quality over time.

## Objectives

1. Keep reproducibility workflows stable and reviewer-friendly.
2. Maintain transparent open development signals (issues, PRs, releases).
3. Improve documentation and evidence quality for software-paper claims.

## Monthly cycle

### Week 1 - Maintenance triage
- Review and triage all open issues.
- Label and prioritize bug, docs, CI, and reproducibility tasks.
- Select a bounded sprint scope.

### Week 2 - Implementation
- Deliver high-priority fixes and documentation updates.
- Ensure every behavior change has test or verification evidence.

### Week 3 - Reproducibility check
- Run canonical fast reproduction pipeline.
- Verify expected artifacts and basic numerical sanity.
- Update README/paper consistency if paths or workflows changed.

### Week 4 - Release and reporting
- Merge pending PRs.
- Publish release notes for changes shipped.
- Record status in the maintenance log and define next cycle goals.

## Quarterly checkpoints

- Reassess paper "State of the field" and "Research impact" text.
- Audit CI coverage and runtime budget.
- Review governance docs (`CONTRIBUTING`, `CODE_OF_CONDUCT`, `SECURITY`, `SUPPORT`).

## Definition of done (cycle)

- At least one merged maintenance PR with evidence.
- No unresolved critical reproducibility regressions.
- Documentation remains aligned with executable entrypoints.
- Maintenance status logged in the consolidated workspace log.
