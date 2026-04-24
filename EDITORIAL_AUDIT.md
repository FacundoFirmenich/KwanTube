# Editorial Audit Ledger

This document tracks the editorial consistency and metadata alignment of the Quantubulin repository following the v3.5.0 hardening pass.

## 1. Bibliographic Architecture

- **Manuscript Source**: `manuscript/quantubulin_prxlife_v1.tex` utilizes an embedded bibliography (`\bibitem{...}`) to ensure independence from external library fluctuations.
- **Software Metadata**: `paper.md` (JOSS-style metadata) and `README.md` utilize `paper.bib` for citation management.

The separation between the conceptual manuscript and the software documentation is explicit and intentional.

## 2. Consistency Status (v3.5.0)

- **Metadata Alignment**: `pyproject.toml`, `CITATION.cff`, and `.zenodo.json` have been synchronized with correct author names (Firmenich et al.) and institutional affiliations (CEDESUR & UNAJ).
- **ORCID Verification**: Primary author ORCID (0009-0002-6578-3811) has been verified across all metadata paths.
- **HEOM Guardrails**: The Bayesian HEOM v2 layer is explicitly documented as a validation summary tool, not a dynamic solver, ensuring scientific transparency.

## 3. Release Readiness

- **Version Alignment**: All artifacts are fixed to version **3.5.0**.
- **DOI Provisioning**: Placeholder DOIs in `CITATION.cff` and `.zenodo.json` must be updated upon Zenodo record activation.
- **Runbook Execution**: A final end-to-end execution of `RELEASE_RUNBOOK.md` is mandatory prior to the archival tag.

## 4. Pending Items

- [ ] Update Zenodo DOI in `CITATION.cff`.
- [ ] Final audit of `paper.bib` against the latest PRX Life citation requirements.
- [ ] Verify that `LIVING_SI.md` cryptographic hash matches the pre-registered `heom_acceptance_criteria.md`.

*End of Audit Ledger.*
