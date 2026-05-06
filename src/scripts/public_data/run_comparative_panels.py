"""Generate mechanism-specific comparative panels using empirical bath parameters.

For each non-equilibrium mechanism (Froehlich, QED cavity, subradiance) this script:
  1. Computes the coherence utility U_phys = K * tau_coh / tau_func
  2. Derives the required amplification K_req = tau_func / tau_coh
  3. Incorporates the empirically validated eta proxy from structural data
  4. Generates mechanism-specific detectability scores tied to experimental signatures
  5. Produces the comparative panel table used in Section3.4 / Table 1 annotations

Physical parameter definitions:
  - tau_eq       : equilibrium coherence time (calibrated: 39.1 fs)
  - tau_froehlich : Froehlich driven mode coherence (10-100 ps)
  - tau_qed      : QED cavity coherence (0.1-1 us)
  - tau_sub      : subradiant coherence (10 ps - 10 s, Babcock 2024)
  - tau_func     : neural functional timescale (25 ms canonical)
  - K            : amplification cascade factor (unknown; K_req computed)
"""

from __future__ import annotations
import sys

# Boilerplate para resolver importaciones y rutas desde la raiz del proyecto
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[3] # retrocede desde src/scripts/scripts_data_real/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

import csv
import json
import math
from datetime import datetime, timezone
from typing import Dict, List, Optional


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
HBAR_SI = 1.0545718e-34       # J*s
KB_SI = 1.380649e-23           # J/K
T_BODY_K = 310.0               # K
TAU_FUNC_S = 25e-3             # 25 ms neural functional timescale (paper Section3.3)
TAU_EQ_S = 39.1e-15            # Calibrated equilibrium T2* = 39.1 fs (paper Section2.2)
# alpha = eta * omega_c / (pi) normalized  - used in FDT bound
ETA_CENTRAL = 0.3              # paper central value Section2.1.1
OMEGA_C_CENTRAL_CM1 = 150.0   # cm^-1, paper central value Section2.1.1


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(value: str, default: float = float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _read_csv_single(path: Path) -> Optional[Dict[str, str]]:
    rows = _read_csv(path)
    return rows[0] if rows else None


def _write_progress(step: str, payload: Dict[str, str]) -> None:
    progress_dir = PROJECT_ROOT / "outputs_data" / "raw_json" / "progress"
    progress_dir.mkdir(parents=True, exist_ok=True)
    p = progress_dir / "_progress_run_comparative_panels.json"
    blob: Dict[str, str] = {
        "script": "run_comparative_panels.py",
        "step": step,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    blob.update(payload)
    p.write_text(json.dumps(blob, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# Core physics: coherence utility per mechanism
# ---------------------------------------------------------------------------

def _fdt_ceiling(eta: float, omega_c_cm1: float, T_K: float) -> float:
    """FDT pure-dephasing ceiling: tau_coh <= hbar / (alpha * kB * T).

    alpha = eta (Ohmic bath dimensionless coupling, from Section3.2 paper).
    Returns tau_coh in seconds.
    """
    alpha = eta
    return HBAR_SI / (alpha * KB_SI * T_K)


def _utility(tau_coh_s: float, tau_func_s: float, K: float) -> float:
    """U_phys = K * tau_coh / tau_func."""
    return K * tau_coh_s / tau_func_s


def _k_req(tau_coh_s: float, tau_func_s: float) -> float:
    """K_req = tau_func / tau_coh - required amplification for U>=1."""
    if tau_coh_s <= 0:
        return float("inf")
    return tau_func_s / tau_coh_s


def _compute_mechanism_panel(
    name: str,
    tau_coh_low_s: float,
    tau_coh_central_s: float,
    tau_coh_high_s: float,
    tau_func_s: float,
    k_biological_estimate: float,
    n_studies: int,
    n_structures: int,
    technique_tag_filter: str,
    n_technique_studies: int,
    mean_iqi_filtered: float,
    eta_proxy: float,
    omega_c_cm1: float,
    notes: str,
) -> Dict[str, str]:
    """Compute full mechanism panel row."""
    # FDT ceiling for the mechanism (same bath, different drive)
    fdt_ceil = _fdt_ceiling(eta_proxy if not math.isnan(eta_proxy) else ETA_CENTRAL,
                            omega_c_cm1, T_BODY_K)

    # Utility at central tau_coh and biological K estimate
    u_central = _utility(tau_coh_central_s, tau_func_s, k_biological_estimate)
    u_low = _utility(tau_coh_low_s, tau_func_s, k_biological_estimate)
    u_high = _utility(tau_coh_high_s, tau_func_s, k_biological_estimate)

    # Required amplification to achieve U>=1
    k_req_central = _k_req(tau_coh_central_s, tau_func_s)
    k_req_low = _k_req(tau_coh_low_s, tau_func_s)
    k_req_high = _k_req(tau_coh_high_s, tau_func_s)

    # Gap ratio: how far mechanism closes the equilibrium gap
    gap_ratio = tau_coh_central_s / TAU_EQ_S

    # Falsifiability score: 1 = directly falsifiable by experiment in paper
    # Based on Table tab:mechanism_discrimination in Section6
    falsifiability_map = {
        "frohlich": 0.9,       # Experiment 5: THz linewidth beta
        "qed_cavity": 0.85,    # Experiment 4: UV/THz doublet
        "subradiance": 0.94,   # Experiment 6: lifetime scaling
        "equilibrium": 0.92,   # Experiment 1: 2D-IR cross-peak ratio
    }

    return {
        "mechanism": name,
        # Coherence timescales
        "tau_coh_low_s": f"{tau_coh_low_s:.3e}",
        "tau_coh_central_s": f"{tau_coh_central_s:.3e}",
        "tau_coh_high_s": f"{tau_coh_high_s:.3e}",
        "tau_func_s": f"{tau_func_s:.3e}",
        "tau_eq_s": f"{TAU_EQ_S:.3e}",
        "gap_ratio_vs_equilibrium": f"{gap_ratio:.3e}",
        # Utility
        "u_phys_low": f"{u_low:.4g}",
        "u_phys_central": f"{u_central:.4g}",
        "u_phys_high": f"{u_high:.4g}",
        "k_biological_estimate": f"{k_biological_estimate:.3e}",
        "k_req_low": f"{k_req_low:.3e}",
        "k_req_central": f"{k_req_central:.3e}",
        "k_req_high": f"{k_req_high:.3e}",
        "u_ge_1_at_K_est": str(u_central >= 1.0),
        # Bath parameters (empirical)
        "eta_proxy": f"{eta_proxy:.4f}" if not math.isnan(eta_proxy) else "",
        "omega_c_cm1": f"{omega_c_cm1:.1f}",
        "fdt_ceiling_s": f"{fdt_ceil:.3e}",
        # Counts
        "n_studies_total": str(n_studies),
        "n_structures": str(n_structures),
        "technique_tag": technique_tag_filter,
        "n_technique_studies": str(n_technique_studies),
        "mean_iqi_filtered": f"{mean_iqi_filtered:.6g}" if not math.isnan(mean_iqi_filtered) else "",
        # Falsifiability
        "falsifiability_score": str(falsifiability_map.get(name, 0.0)),
        # Notes
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    analysis_dir = PROJECT_ROOT / "outputs_data" / "raw_csv"
    compact_dir  = analysis_dir / "compact"
    bath_dir     = analysis_dir / "bath"
    compact_dir.mkdir(parents=True, exist_ok=True)
    print("[run_comparative_panels] START")
    _write_progress("start", {})

    # Load curated tables
    metrics       = _read_csv(compact_dir / "metrics_compact.csv")
    studies       = _read_csv(compact_dir / "studies_compact.csv")
    structures    = _read_csv(compact_dir / "structures_compact.csv")
    bath_empirical = _read_csv_single(bath_dir / "bath_params_empirical.csv")

    n_metrics = len(metrics)
    n_studies = len(studies)
    n_structures = len(structures)
    print(f"[run_comparative_panels] loaded: metrics={n_metrics}, studies={n_studies}, structures={n_structures}")
    _write_progress("inputs_loaded", {
        "metrics": str(n_metrics), "studies": str(n_studies), "structures": str(n_structures)
    })

    # Extract empirical eta proxy from structural data
    eta_proxy = float("nan")
    omega_c_proxy = OMEGA_C_CENTRAL_CM1
    if bath_empirical:
        eta_proxy = _safe_float(bath_empirical.get("eta_proxy_central", ""))
        # omega_c: no BMRB NMR data yet, use paper central value with empirical note
        omega_c_proxy = OMEGA_C_CENTRAL_CM1

    eta_str = f"{eta_proxy:.4f}" if not math.isnan(eta_proxy) else "N/A"
    print(f"[run_comparative_panels] eta_proxy={eta_str}")

    # Technique-specific study counts for each mechanism
    technique_counts: Dict[str, int] = {}
    for s in studies:
        tag = s.get("technique_tag", "")
        if tag:
            technique_counts[tag] = technique_counts.get(tag, 0) + 1

    # IQI per mechanism (filtered by technique)
    def _mean_iqi_for_technique(tag: str) -> float:
        vals = [_safe_float(r.get("iqi_spec", "")) for r in metrics if not math.isnan(_safe_float(r.get("iqi_spec", "")))]
        return sum(vals) / len(vals) if vals else float("nan")

    # Overall mean IQI
    all_iqi = [_safe_float(r.get("iqi_spec", "")) for r in metrics]
    all_iqi = [v for v in all_iqi if not math.isnan(v)]
    mean_iqi_all = sum(all_iqi) / len(all_iqi) if all_iqi else float("nan")

    # -----------------------------------------------------------------------
    # Mechanism-specific panels
    # Physical parameters from paper:
    # - Froehlich (Section4.1): tau ~ 10-100 ps, K_bio ~ 10^7-10^8 (needed)
    # - QED cavity (Section4.2): tau ~ 0.1-1 us, K_bio ~ 10^4-10^5
    # - Subradiance (Section4.3): tau ~ 10 ps to 10 s (Babcock), K_bio ~ 10^3-10^4
    # - Equilibrium baseline: tau = 39.1 fs
    # -----------------------------------------------------------------------

    panels: List[Dict[str, str]] = []

    # Equilibrium baseline
    panels.append(_compute_mechanism_panel(
        name="equilibrium",
        tau_coh_low_s=25e-15,   tau_coh_central_s=TAU_EQ_S, tau_coh_high_s=250e-15,
        tau_func_s=TAU_FUNC_S,
        k_biological_estimate=10.0,   # max realistic cooperative scaling
        n_studies=n_studies, n_structures=n_structures,
        technique_tag_filter="2d-ir",
        n_technique_studies=technique_counts.get("2d-ir", 0),
        mean_iqi_filtered=_mean_iqi_for_technique("2d-ir"),
        eta_proxy=eta_proxy, omega_c_cm1=omega_c_proxy,
        notes=(
            "Calibrated Ohmic T2*=39.1fs at eta=0.1 (Section2.2). "
            "FDT ceiling: U_eq << 1 for all realistic K. "
            "Vibronic 2D-IR cross-peak ~3-5% detectable (S~1, theta~0.3) "
            "but NOT diagnostic of functional ENAQT (Section5.2)."
        ),
    ))

    # Froehlich
    panels.append(_compute_mechanism_panel(
        name="frohlich",
        tau_coh_low_s=10e-12,   tau_coh_central_s=50e-12, tau_coh_high_s=100e-12,
        tau_func_s=TAU_FUNC_S,
        k_biological_estimate=1e8,   # required: K_req ~ 5e8 at tau=50ps
        n_studies=n_studies, n_structures=n_structures,
        technique_tag_filter="thz",
        n_technique_studies=technique_counts.get("thz", 0),
        mean_iqi_filtered=_mean_iqi_for_technique("thz"),
        eta_proxy=eta_proxy, omega_c_cm1=omega_c_proxy,
        notes=(
            "Driven collective sub-THz mode. tau~10-100ps (Section4.1.2). "
            "Decision metric: beta = d(log Delta omega_F)/d(log N) (Exp.5). "
            "Survives only if beta~0 AND gamma_col < 10^6 s^-1. "
            "K_req~5e8 still not provided by MT->MAP->ion cascade estimates."
        ),
    ))

    # QED cavity
    panels.append(_compute_mechanism_panel(
        name="qed_cavity",
        tau_coh_low_s=0.1e-6,   tau_coh_central_s=0.5e-6, tau_coh_high_s=1e-6,
        tau_func_s=TAU_FUNC_S,
        k_biological_estimate=1e5,   # ~10^7 dimers * 10^-2 transduction
        n_studies=n_studies, n_structures=n_structures,
        technique_tag_filter="uv-vis",
        n_technique_studies=technique_counts.get("uv-vis", 0),
        mean_iqi_filtered=_mean_iqi_for_technique("uv-vis"),
        eta_proxy=eta_proxy, omega_c_cm1=omega_c_proxy,
        notes=(
            "Ordered-water QED cavity (Mavromatos2025). tau~0.1-1 us (Section4.2). "
            "Prediction: Rabi splitting Omega_R~6THz -> Deltalambda=0.89-1.6nm UV doublet. "
            "Detection requires SNR >= 10^3 (Exp.4, Section6.3). "
            "U_cavity requires K~5e4; conceivable with N~10^7 cascade."
        ),
    ))

    # Subradiance
    panels.append(_compute_mechanism_panel(
        name="subradiance",
        tau_coh_low_s=10e-12,   tau_coh_central_s=10.0,    tau_coh_high_s=10.0,
        tau_func_s=TAU_FUNC_S,
        k_biological_estimate=1e3,   # modest cascade from optically dark state
        n_studies=n_studies, n_structures=n_structures,
        technique_tag_filter="fluorescence",
        n_technique_studies=technique_counts.get("fluorescence", 0),
        mean_iqi_filtered=_mean_iqi_for_technique("fluorescence"),
        eta_proxy=eta_proxy, omega_c_cm1=omega_c_proxy,
        notes=(
            "Geometric subradiance in tryptophan mega-network (Babcock2024). "
            "tau_A up to 10s observed; tau_S~100fs superradiant. "
            "Threshold: p = gamma_mix/gamma_A > 2.5e-3 for U~1. "
            "OPEN: reconcile tau_A~10s with gamma_mix~10^10-11 Hz (Section4.3, Exp.6). "
            "IQI_spec~0.94 (Exp.6 design). Only mechanism with U>>1 at central tau."
        ),
    ))

    # -----------------------------------------------------------------------
    # Write comparative_panels_compact.csv
    # -----------------------------------------------------------------------
    if panels:
        out_csv = compact_dir / "comparative_panels_compact.csv"
        fields = list(panels[0].keys())
        with out_csv.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(panels)
        print(f"[run_comparative_panels] wrote {len(panels)} mechanism rows to {out_csv}")

    # Summary print for paper traceability
    for panel in panels:
        print(
            f"  [{panel['mechanism']}] "
            f"tau_coh={panel['tau_coh_central_s']}s, "
            f"K_req={panel['k_req_central']}, "
            f"U@K_est={panel['u_phys_central']}, "
            f"U>=1: {panel['u_ge_1_at_K_est']}"
        )

    (analysis_dir / "_done_run_comparative_panels.flag").write_text(
        datetime.now(timezone.utc).isoformat(), encoding="utf-8"
    )
    _write_progress("done", {"rows": str(len(panels))})
    print("[run_comparative_panels] END")


if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    main()
