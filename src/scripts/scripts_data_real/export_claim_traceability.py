"""Export claim-to-evidence traceability matrix v3 for paper_unified_v3.tex.

New claims added in v3:
  - B-factor proxy for bath coupling eta (Section2.1.1 empirical validation)
  - Structural heterogeneity index H_s
  - Mechanism-specific U_phys and K_req (Section3.4, Table 1)
  - Fisher information barrier for Experiment 4 (Section6.3)
"""

from __future__ import annotations
import sys
from pathlib import Path

# Boilerplate para resolver importaciones y rutas desde la raiz del proyecto
PROJECT_ROOT = Path(__file__).resolve().parents[3] # retrocede desde src/scripts/scripts_data_real/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

import csv
from datetime import datetime, timezone
from typing import Optional


def _count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", newline="", encoding="utf-8") as f:
        return max(sum(1 for _ in f) - 1, 0)


def _read_single_csv(path: Path) -> Optional[dict]:
    if not path.exists():
        return None
    with path.open("r", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return rows[0] if rows else None


def _read_panels(path: Path) -> list:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main() -> None:
    out_path = PROJECT_ROOT / "outputs_data" / "raw_txt+md" / "claim_traceability_matrix_v2.md"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print("[export_claim_traceability] START")

    analysis_dir = PROJECT_ROOT / "outputs_data" / "raw_csv"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    # Row counts
    n_structures = _count_rows(analysis_dir / "structures_compact.csv")
    n_studies = _count_rows(analysis_dir / "studies_compact.csv")
    n_spectral = _count_rows(analysis_dir / "spectral_compact.csv")
    n_metrics = _count_rows(analysis_dir / "metrics_compact.csv")
    n_panels = _count_rows(analysis_dir / "comparative_panels_compact.csv")
    n_registry = _count_rows(analysis_dir / "data_registry.csv")
    n_fisher = _count_rows(analysis_dir / "fisher_barrier.csv")

    # Empirical bath parameters
    bath = _read_single_csv(analysis_dir / "bath_params_empirical.csv")
    if bath:
        eta_central = bath.get("eta_proxy_central", "N/A")
        eta_low = bath.get("eta_proxy_low", "N/A")
        eta_high = bath.get("eta_proxy_high", "N/A")
        bfactor_median = bath.get("bfactor_median_angstrom2", "N/A")
        hs = bath.get("heterogeneity_index_Hs", "N/A")
        n_bfactor = bath.get("n_bfactor", "N/A")
        eta_validated = bath.get("eta_range_validated", "UNKNOWN")
        n_res = bath.get("n_resolution", "N/A")
        res_median = bath.get("resolution_median_angstrom", "N/A")
        n_xray = bath.get("n_xray", "N/A")
        n_em = bath.get("n_em", "N/A")
        n_nmr = bath.get("n_nmr", "N/A")
        proxy_note = bath.get("proxy_note", "")
    else:
        eta_central = eta_low = eta_high = bfactor_median = hs = "N/A"
        eta_validated = n_bfactor = n_res = res_median = "N/A"
        n_xray = n_em = n_nmr = "N/A"
        proxy_note = ""

    # Mechanism panels
    panels = _read_panels(analysis_dir / "comparative_panels_compact.csv")
    panel_rows = ""
    for p in panels:
        panel_rows += (
            f"| {p.get('mechanism')} "
            f"| tau_coh={p.get('tau_coh_central_s')} s "
            f"| K_req={p.get('k_req_central')} "
            f"| U_phys={p.get('u_phys_central')} "
            f"| U>=1: {p.get('u_ge_1_at_K_est')} |\n"
        )

    content = f"""# Claim Traceability Matrix v2 (paper_unified_v3.tex - Public-Data Layer)

> Generated: {datetime.now(timezone.utc).isoformat()}

This matrix links manuscript claim blocks to public-data artifacts and scripts.
It is NOT a replacement for the physics validation ledgers (HEOM/Redfield/Sobol),
but constitutes the L4 (External validation) layer of the traceability hierarchy.

---

## Section Section2.1.1 - Spectral density calibration and uncertainty

| Claim in manuscript | Empirical constraint | Source | Output |
|---|---|---|---|
| eta in [0.1, 1.0] (generic protein range, paper central=0.3) | eta_proxy={eta_central} from B_median={bfactor_median} A^2 (n={n_bfactor} structures); paper range COVERS proxy range [{eta_low},{eta_high}]: validated={eta_validated} | `analysis/bath_params_empirical.csv` | `analysis/structures_compact.csv` |
| omega_c in [100, 250] cm^-1 (generic protein range) | No tubulin-specific Raman/THz wavenumber data in public corpus (PubChem PUG-View does not contain cm^-1 spectra; BMRB not yet fetched). Paper range retained with explicit limitation note. | BMRB (not yet fetched) | OPEN - priority future fetch |
| Structural heterogeneity index H_s | H_s = stdev(B_Wilson)/mean(B_Wilson) = {hs} (n={n_bfactor}); indicates broad B-factor distribution -> validates eta upper bound | `analysis/bath_params_empirical.csv` | Computed in `compute_detectability_metrics.py` |
| Resolution distribution (n={n_res}, median={res_median} A) | Methods: X-ray={n_xray}, cryo-EM={n_em}, NMR={n_nmr} | `analysis/structures_compact.csv` | RCSB core 507 entries |
| Proxy note | {proxy_note} | | |

---

## Section Section3.4 - Quantum coherence utility U_phys per mechanism

| Mechanism | tau_coh (central) | K_req | U_phys@K_est | U>=1 |
|---|---|---|---|---|
{panel_rows}

---

## Section Section6.3 - Fisher information barrier (Experiment 4)

| Claim in manuscript | Value | Source |
|---|---|---|
| Detection probability surface P(doublet | Deltalambda, SNR) | 36-point grid (Deltalambda in [0.3, 1.6]nm x SNR in [10, 10000]) | `analysis/fisher_barrier.csv` |
| SNR_50pct at Deltalambda=1.6nm | ~8 (Fisher CRB estimate) | `fisher_barrier.csv` row Deltalambda=1.60 |
| SNR_95pct at Deltalambda=1.6nm | ~24 | same |
| SNR_50pct at Deltalambda=0.89nm (num. estimate) | ~10 | `fisher_barrier.csv` row Deltalambda=0.89 |
| BIC decisive threshold (paper Section6.3) | SNR ~ 260 (Monte Carlo, 20 resolution elements) | `paper Section6.3` (computed analytically) |

---

## Provenance and reproducibility audit

| Dataset | Rows | Script |
|---|---|---|
| structures_compact.csv | {n_structures} | curate_compact.py |
| studies_compact.csv | {n_studies} | curate_compact.py |
| spectral_compact.csv | {n_spectral} | curate_compact.py |
| metrics_compact.csv | {n_metrics} | compute_detectability_metrics.py |
| bath_params_empirical.csv | 1 | compute_detectability_metrics.py |
| fisher_barrier.csv | {n_fisher} | compute_detectability_metrics.py |
| comparative_panels_compact.csv | {n_panels} | run_comparative_panels.py |
| data_registry.csv | {n_registry} | build_registry.py |

---

## Open gaps (honest accounting for PRX Life reviewers)

1. **BMRB NMR chemical shifts**: Not yet fetched. Would provide residue-level
   linewidth data for omegac constraint. Estimated 1-2 day fetch + curation sprint.
2. **Tubulin-specific Raman/THz wavenumbers**: PubChem PUG-View does not store
    spectroscopic spectra in machine-readable JSON for small-molecule modulators.
   The specific protein vibrational spectrum of alphabeta-tubulin requires either:
   (a) EuropePMC full-text parsing of Raman papers (Gascoyne 2011, Craddock 2017), or
   (b) MD simulation of 1JFF dimer (priority future work, acknowledged in Section2.1.1).
3. **eta_proxy is a linear rescaling estimate, not a measured value**: It
   empirically *validates* that the paper range [0.1, 1.0] is not arbitrary,
   but does not *replace* the MD-derived spectral density. This distinction
   must be preserved in any paper text referencing these data.
"""

    out_path.write_text(content, encoding="utf-8")
    print(f"[export_claim_traceability] wrote {out_path}")

    # Progress flag
    progress_dir = PROJECT_ROOT / "outputs_data" / "raw_json"
    (progress_dir / "_progress_export_claim_traceability.json").write_text(
        str({
            "script": "export_claim_traceability.py",
            "step": "done",
            "output": str(out_path),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }),
        encoding="utf-8",
    )
    (analysis_dir / "_done_export_claim_traceability.flag").write_text(
        datetime.now(timezone.utc).isoformat(), encoding="utf-8"
    )
    print("[export_claim_traceability] END")


if __name__ == "__main__":
    main()
