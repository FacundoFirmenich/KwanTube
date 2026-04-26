"""
sbc_report.py -- Simulation-Based Calibration Report Generator.
Validates the Nested Sampling engine and exports rank histogram (PNG/PDF)
and a clean JSON validation report.
STANDARD CONFIGURATION: n_sim=1000, confidence=0.95
"""
import numpy as np
import json
import sys
from pathlib import Path

# Boilerplate to resolve imports from package root
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_dist

from qmc_mt.nested_sampling import nested_sample
from qmc_mt.sbc import simulation_based_calibration
from qmc_mt.meta import _norm_ppf

# ========== STANDARD CONFIGURATION ==========
SIGMA_DATA = 1.0
PRIOR_SD   = 2.0
N_SIM      = 1000      # ← ESTÁNDAR: 1000 simulaciones
L_POST     = 99        # ← 99 muestras posteriores
N_BINS     = 20        # ← 20 bins (dof=19)
SEED       = 42        # ← Reproducibilidad
CONFIDENCE = 0.95      # ← Nivel de confianza 95%
FIG_DIR    = PROJECT_ROOT / "figures_final"
# ============================================

def prior_sampler(rng: np.random.Generator) -> float:
    return float(rng.normal(0.0, PRIOR_SD))

def data_sampler(theta: float, rng: np.random.Generator) -> np.ndarray:
    return rng.normal(theta, SIGMA_DATA, size=1)

def posterior_sampler_ns(y: np.ndarray, L: int, rng: np.random.Generator) -> np.ndarray:
    """Uses Nested Sampling to draw L samples from the posterior."""
    val_y = float(y[0])
    s = SIGMA_DATA
    
    def loglike(theta: np.ndarray) -> float:
        t = float(theta[0])
        return -0.5 * ((val_y - t) / s) ** 2 - 0.5 * np.log(2 * np.pi * s**2)
    
    def prior_transform(u: np.ndarray) -> np.ndarray:
        return np.array([PRIOR_SD * _norm_ppf(float(u[0]))])

    res = nested_sample(loglike, prior_transform, ndim=1, n_live=100, 
                       seed=rng.integers(0, int(1e6)))
    idx = rng.choice(len(res.samples), size=L, p=res.weights)
    return res.samples[idx].ravel()

def run_sbc_report(n_sim: int = N_SIM) -> dict:
    print(f"Starting SBC for Nested Sampling Engine")
    print(f"  Configuration: n_sim={n_sim}, L={L_POST}, n_bins={N_BINS}")
    print(f"  Confidence level: {CONFIDENCE*100:.0f}%")
    print(f"  Expected chi2 range: [{chi2_dist.ppf((1-CONFIDENCE)/2, N_BINS-1):.2f}, "
          f"{chi2_dist.ppf((1+CONFIDENCE)/2, N_BINS-1):.2f}]")
    
    res = simulation_based_calibration(
        prior_sampler=prior_sampler,
        data_sampler=data_sampler,
        posterior_sampler=posterior_sampler_ns,
        n_sim=n_sim, L=L_POST, n_bins=N_BINS, seed=SEED
    )
    
    # Expected range for confidence interval
    chi2_lower = chi2_dist.ppf((1-CONFIDENCE)/2, N_BINS-1)
    chi2_upper = chi2_dist.ppf((1+CONFIDENCE)/2, N_BINS-1)
    
    # ---- PLOT GENERATION ----
    plt.figure(figsize=(10, 6))
    plt.bar(res.bin_edges[:-1], res.bin_counts, width=np.diff(res.bin_edges),
            align='edge', color='skyblue', edgecolor='navy', alpha=0.7,
            label=f'Observed (chi2={res.chi2:.2f}, p={res.p_value:.4f})')
    plt.axhline(n_sim / N_BINS, color='red', linestyle='--', linewidth=2, 
                label=f'Ideal Uniformity (n_sim/n_bins = {n_sim/N_BINS:.1f})')
    
    # Confidence band
    expected_per_bin = n_sim / N_BINS
    std_per_bin = np.sqrt(expected_per_bin)
    plt.axhspan(expected_per_bin - 1.96*std_per_bin, 
                expected_per_bin + 1.96*std_per_bin,
                alpha=0.2, color='green', label='95% Poisson CI per bin')
    
    plt.xlabel("Rank", fontsize=12, fontweight='bold')
    plt.ylabel("Frequency", fontsize=12, fontweight='bold')
    plt.title(f"SBC Rank Histogram: Nested Sampling Engine\n"
              f"p-value = {res.p_value:.4f}, chi2 = {res.chi2:.2f} (dof={N_BINS-1})\n"
              f"Expected chi2 range ({CONFIDENCE*100:.0f}%): [{chi2_lower:.2f}, {chi2_upper:.2f}]",
              fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout()

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        out_path = FIG_DIR / f"sbc_calibration_ns.{ext}"
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"  [OK] Saved {out_path}")
    plt.close()

    # ---- JSON REPORT (Clean keys, no trailing spaces) ----
    in_expected_range = bool(chi2_lower <= res.chi2 <= chi2_upper)
    
    report = {
        "engine": "nested_sampling.py",
        "n_simulations": int(res.n_sim),
        "n_post_samples_per_sim": int(res.L),
        "n_bins": int(N_BINS),
        "confidence_level": float(CONFIDENCE),
        "chi2_stat": float(round(res.chi2, 4)),
        "chi2_expected_range": [float(round(chi2_lower, 4)), float(round(chi2_upper, 4))],
        "chi2_in_expected_range": in_expected_range,
        "p_value": float(round(res.p_value, 6)),
        "is_calibrated_at_0.05": bool(0.05 < res.p_value < 0.95),
        "interpretation": _get_interpretation(res.chi2, res.p_value, chi2_lower, chi2_upper),
        "variability_check": {
            "expected_chi2_mean": int(N_BINS - 1),
            "observed_chi2": float(res.chi2),
            "deviation_from_expected": float(round((res.chi2 - (N_BINS-1)) / (N_BINS-1) * 100, 2)),
            "status": "NORMAL" if in_expected_range else "UNUSUAL"
        }
    }

    report_path = PROJECT_ROOT / "sbc_report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(f"  [OK] Saved {report_path}")
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"SBC CALIBRATION SUMMARY")
    print(f"{'='*60}")
    print(f"chi2 statistic:      {res.chi2:.4f}")
    print(f"Expected range:      [{chi2_lower:.2f}, {chi2_upper:.2f}]")
    print(f"p-value:             {res.p_value:.4f}")
    print(f"In expected range:   {in_expected_range}")
    print(f"Interpretation:      {_get_interpretation(res.chi2, res.p_value, chi2_lower, chi2_upper)}")
    print(f"{'='*60}\n")
    
    return report

def _get_interpretation(chi2: float, p_value: float, lower: float, upper: float) -> str:
    if not (lower <= chi2 <= upper):
        if chi2 < lower:
            return "UNUSUALLY_LOW_VARIANCE - bins too uniform (possible overfitting or bug)"
        else:
            return "HIGH_VARIANCE - possible miscalibration"
    elif p_value < 0.05:
        return "MISALIGNED - reject calibration at 5% level"
    elif p_value > 0.95:
        return "UNUSUALLY_GOOD - variance lower than expected (check for issues)"
    else:
        return "WELL_CALIBRATED - consistent with uniform ranks"

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate SBC report for KwanTube")
    parser.add_argument("--n-sim", type=int, default=N_SIM, 
                       help=f"Number of SBC simulations (default: {N_SIM})")
    args = parser.parse_args()
    run_sbc_report(n_sim=args.n_sim)