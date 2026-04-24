# Editorial Audit (paper.md / paper.bib / paper.tex)

This file tracks editorial consistency status after the final v3.5.0 hardening pass.

## 1) Bibliography architecture clarified

- `paper/paper.tex` uses an **embedded bibliography** (`\bibitem{...}` block in-file).
- `git_repo/paper.bib` is used by `git_repo/paper.md` (software-paper metadata path).

Therefore, `paper/paper.tex` is not blocked by `paper.bib` key coverage, because it does not import it.

## 2) Current consistency status

- `paper.md` placeholder metadata fields were resolved (authors, affiliations, ORCID, acknowledgements text).
- `paper.md` version mention aligned to **v3.5.0**.
- Bayesian HEOM v2 guardrail is documented in repository/manuscript scope as:
  - validation layer for existing convergence evidence,
  - not a replacement for new full-system HEOM runs.

## 3) Metadata consistency

Aligned:
- `pyproject.toml` version -> **3.5.0**
- `CITATION.cff` version -> **3.5.0**

Pending before archival DOI freeze:
- Replace `CITATION.cff` placeholder DOI (`10.5281/zenodo.XXXXX`) with final Zenodo record.

## 4) Release recommendation

Before final archival tag:
1. Mint/fetch Zenodo DOI and update `CITATION.cff`.
2. Run release runbook (`RELEASE_RUNBOOK.md`) end-to-end.
3. Keep bibliography split explicit: `paper.tex` (embedded), `paper.md` (`paper.bib`).
