#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_heom_production_figure.py
Dedicated extraction, metric computation, and publication-ready figure generation
for the completed 30 ps HEOM production trajectory (1JFF, NC=7, Nk=1).

Outputs:
  - Terminal report with full metric disclosure
  - Audit trail appended to outputs_data/raw_txt+md/execution_memory.log.txt
  - heom_production_comparison.pdf / .png in outputs_data/figures_final/

This script operates independently of the global validation pipeline and does not
modify LIVING_SI.md or validation_report.json.
"""

import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime, timezone

# =============================================================================
# PATH RESOLUTION & AUDIT INFRASTRUCTURE
# =============================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent.parent # Adjusted to src/scripts/analysis/ -> KwanTube/
NPZ_PATH = PROJECT_ROOT / "outputs_data" / "raw_npz" / "master_results.npz"
FIG_DIR = PROJECT_ROOT / "outputs_data" / "figures_final"
LOG_PATH = PROJECT_ROOT / "outputs_data" / "raw_txt+md" / "execution_memory.log.txt"

FIG_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH.parent.mkdir(parents=True, exist_ok=True)

def audit_log(msg: str) -> None:
    """Appends a timestamped audit entry to execution_memory.log.txt"""
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")
    line = f"[{ts}] [extract_heom_production_figure.py] {msg}\n"
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(line)

# =============================================================================
# EXTRACTION & METRIC COMPUTATION
# =============================================================================
def run_extraction_and_plot() -> int:
    audit_log("[RUN_AUDIT] START script=extract_heom_production_figure.py")
    t0 = time.time()

    if not NPZ_PATH.exists():
        err_msg = f"[FAIL] Master dataset not found at {NPZ_PATH}. Run assemble_master_results.py first."
        print(err_msg)
        audit_log(err_msg)
        return 1

    data = np.load(NPZ_PATH)
    t_fs = data["tlist"]          # shape (3001,)
    pops = data["populations"]    # shape (8, 3001)

    # Initial excitation site: B:103 corresponds to internal index 5
    init_idx = 5
    P_init = pops[init_idx, :]

    # Time-resolved metrics
    purity = np.sum(pops**2, axis=0)
    ipr = 1.0 / np.maximum(np.sum(pops**2, axis=0), 1e-15)

    # Interpolate at 500 fs for direct Redfield comparison
    P_init_500fs = float(np.interp(500.0, t_fs, P_init))
    purity_30ps = float(purity[-1])
    ipr_30ps = float(ipr[-1])

    # Redfield baseline (validated ledger value)
    redfield_P_init_500fs = 0.6248
    discrepancy_pct = float(abs(P_init_500fs - redfield_P_init_500fs) / P_init_500fs * 100)

    # =====================================================================
    # TERMINAL REPORT
    # =====================================================================
    print("=" * 72)
    print("HEOM PRODUCTION TRAJECTORY EXTRACTION & METRIC DISCLOSURE")
    print("=" * 72)
    print(f"Dataset      : {NPZ_PATH.name}")
    print(f"Time window  : 0 → {t_fs[-1]:.0f} fs ({len(t_fs)} points)")
    print(f"Init site    : B:103 (internal index {init_idx})")
    print("-" * 72)
    print(f"P_init(500 fs) [HEOM]      : {P_init_500fs:.4f}")
    print(f"P_init(500 fs) [Redfield]  : {redfield_P_init_500fs:.4f}")
    print(f"Relative discrepancy       : {discrepancy_pct:.2f}%")
    print(f"Purity(30 ps)              : {purity_30ps:.3f}")
    print(f"IPR(30 ps)                 : {ipr_30ps:.3f}")
    print(f"Regime diagnosis           : {'NON-EQUILIBRIUM TRANSIENT' if purity_30ps < 0.25 else 'NEAR-THERMAL'}")
    print("=" * 72)

    # =====================================================================
    # FIGURE GENERATION (ACADEMIC STYLING)
    # =====================================================================
    plt.rcParams.update({
        "font.size": 10,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.linewidth": 0.8,
        "grid.linestyle": ":",
        "grid.alpha": 0.45,
        "figure.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05
    })

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3))

    # Panel A: Initial site population decay vs Redfield baseline
    ax1 = axes[0]
    ax1.plot(t_fs, P_init, color="#1f77b4", linewidth=1.6, label="HEOM ($N_C=7$)")
    ax1.axhline(y=redfield_P_init_500fs, color="#d62728", linestyle="--", linewidth=1.2,
                label="Redfield baseline (500 fs)")
    ax1.set_xlabel("Time (fs)")
    ax1.set_ylabel("Initial Site Population $P_{\\mathrm{init}}(t)$")
    ax1.set_xlim(0, 30000)
    ax1.set_ylim(0, 1.05)
    ax1.grid(True)
    ax1.legend(frameon=False, loc="upper right", fontsize=9)

    # Panel B: Purity and Inverse Participation Ratio
    ax2 = axes[1]
    ax2.plot(t_fs, purity, color="#2ca02c", linewidth=1.6, label="Purity $\\mathcal{P}(t)$")
    ax2.plot(t_fs, ipr, color="#9467bd", linestyle="--", linewidth=1.6, label="IPR$(t)$")
    ax2.set_xlabel("Time (fs)")
    ax2.set_ylabel("Metric")
    ax2.set_xlim(0, 30000)
    ax2.set_ylim(0, max(purity.max(), ipr.max()) * 1.1)
    ax2.grid(True)
    ax2.legend(frameon=False, loc="upper right", fontsize=9)

    plt.tight_layout(pad=1.2)

    pdf_path = FIG_DIR / "heom_production_comparison.pdf"
    png_path = FIG_DIR / "heom_production_comparison.png"
    fig.savefig(pdf_path, transparent=False)
    fig.savefig(png_path, transparent=False)
    plt.close(fig)

    elapsed = time.time() - t0
    success_msg = (f"[SUCCESS] Figures saved: {pdf_path.name}, {png_path.name} | "
                   f"Wall-time: {elapsed:.2f}s")
    print(success_msg)
    audit_log(success_msg)
    audit_log("[RUN_AUDIT] END status=ok")
    return 0

# =============================================================================
# ENTRY POINT
# =============================================================================
if __name__ == "__main__":
    sys.exit(run_extraction_and_plot())
