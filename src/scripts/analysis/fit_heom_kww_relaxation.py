#!/usr/bin/env python3
"""Fit stretched-exponential relaxation diagnostics to HEOM trajectories.

This script addresses the reviewer request for a Kohlrausch--Williams--Watts
(KWW) relaxation analysis without inventing observables that are not present in
the archived data. It supports two evidence layers:

1. The 30 ps production master dataset (``master_results.npz``), which stores
   site populations. From this file the script fits population-diagonal
   observables such as ``sum_i p_i(t)^2`` and the initially excited population.
2. Optional HEOM pickle dumps containing full density matrices (``rho_t``). If
   available, the script additionally fits quantum purity and L1 coherence.

Outputs are written as JSON, CSV, and publication-ready figures. The JSON report
explicitly labels each observable as either population-diagonal or full-density,
so downstream manuscript text can avoid overclaiming.

Usage
-----
    # From the project root (KwanTube/):
    python paper/fit_heom_kww_relaxation.py

    # With explicit paths:
    python paper/fit_heom_kww_relaxation.py \\
        --project-root /path/to/KwanTube \\
        --master-npz outputs_data/raw_npz/master_results.npz \\
        --density-pkl outputs_data/raw_pkl/heom_1JFF.pkl \\
        --bootstrap 500

Interpretation note
-------------------
The KWW fit on a 30 ps window captures the early-time non-Markovian transient.
A stretched exponent beta < 1 indicates deviation from simple exponential
kinetics but does NOT by itself prove a glassy bath; it is presented as a
lower-bound descriptor of non-Markovian relaxation complexity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import pickle
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FitResult:
    """Container for a single KWW fit result."""

    metric: str
    data_layer: str
    n_points: int
    fit_start_fs: float
    fit_end_fs: float
    y_inf: float
    amplitude: float
    tau_fs: float
    beta: float
    rmse: float
    r2: float
    aic: float
    bic: float
    stderr: dict
    bootstrap_ci95: dict
    status: str
    warning: str


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def resolve_project_root(explicit: Optional[str] = None) -> Path:
    """Resolve the KwanTube project root.

    Args:
        explicit: Optional explicit root path supplied by the caller.

    Returns:
        Absolute path to the project root.
    """
    if explicit:
        return Path(explicit).expanduser().resolve()
    here = Path(__file__).resolve()
    candidates = [
        here.parents[3],
        Path.cwd().resolve(),
        Path.cwd().resolve() / "KwanTube",
    ]
    for candidate in candidates:
        if (candidate / "outputs_data").exists() and (candidate / "src").exists():
            return candidate.resolve()
    return here.parents[3].resolve()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Compute a SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_log(project_root: Path, message: str) -> None:
    """Append a timestamped execution line to the canonical memory log."""
    log_path = (
        project_root
        / "outputs_data"
        / "raw_txt+md"
        / "logs"
        / "execution_memory.log.txt"
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(
            f"[{timestamp}] [fit_heom_kww_relaxation.py] {message}\n"
        )


# ---------------------------------------------------------------------------
# KWW model and fitting primitives
# ---------------------------------------------------------------------------

def kww_model(
    t_fs: np.ndarray,
    y_inf: float,
    amplitude: float,
    tau_fs: float,
    beta: float,
) -> np.ndarray:
    """Kohlrausch--Williams--Watts relaxation model.

    y(t) = y_inf + amplitude * exp[-(t/tau)^beta].

    The model is sign-flexible through ``amplitude`` so it can fit both
    decays and rises toward an asymptote.
    """
    safe_t = np.maximum(np.asarray(t_fs, dtype=float), 0.0)
    safe_tau = max(float(tau_fs), np.finfo(float).eps)
    return y_inf + amplitude * np.exp(-np.power(safe_t / safe_tau, beta))


def finite_window(
    t_fs: np.ndarray,
    y: np.ndarray,
    fit_start_fs: float,
    fit_end_fs: Optional[float],
) -> tuple:
    """Select a finite, non-NaN fitting window.

    Returns:
        (t_relative, y_selected) where t_relative starts at 0.
    """
    mask = np.isfinite(t_fs) & np.isfinite(y) & (t_fs >= fit_start_fs)
    if fit_end_fs is not None:
        mask &= t_fs <= fit_end_fs
    t_sel = np.asarray(t_fs[mask], dtype=float)
    y_sel = np.asarray(y[mask], dtype=float)
    if t_sel.size < 8:
        raise ValueError(
            "KWW fitting requires at least 8 finite points in the selected window."
        )
    return t_sel - t_sel[0], y_sel


def initial_guess_and_bounds(
    t_rel: np.ndarray, y: np.ndarray
) -> tuple:
    """Create robust initial guesses and parameter bounds for KWW fitting."""
    y_min = float(np.min(y))
    y_max = float(np.max(y))
    y_range = max(y_max - y_min, 1e-12)
    tail_n = max(3, int(0.1 * y.size))
    y_inf0 = float(np.mean(y[-tail_n:]))
    amp0 = float(y[0] - y_inf0)
    tau0 = max(float((t_rel[-1] - t_rel[0]) / 3.0), 1.0)
    beta0 = 0.7

    lower = [y_min - 2.0 * y_range, -3.0 * y_range, 1e-9, 0.05]
    upper = [
        y_max + 2.0 * y_range,
        3.0 * y_range,
        max(t_rel[-1] * 100.0, 10.0),
        3.0,
    ]
    return [y_inf0, amp0, tau0, beta0], (lower, upper)


def information_criteria(
    y: np.ndarray, residuals: np.ndarray, n_params: int
) -> tuple:
    """Compute Gaussian AIC and BIC from residuals."""
    n = y.size
    rss = float(np.sum(residuals**2))
    sigma2 = max(rss / n, np.finfo(float).tiny)
    log_likelihood = -0.5 * n * (math.log(2.0 * math.pi * sigma2) + 1.0)
    aic = 2 * n_params - 2 * log_likelihood
    bic = n_params * math.log(n) - 2 * log_likelihood
    return float(aic), float(bic)


def bootstrap_kww(
    t_rel: np.ndarray,
    y_fit: np.ndarray,
    y_hat: np.ndarray,
    popt: np.ndarray,
    bounds: tuple,
    n_bootstrap: int,
    seed: int,
) -> dict:
    """Residual-bootstrap confidence intervals for KWW parameters.

    The returned intervals are exploratory because HEOM trajectories are
    autocorrelated time series. They are useful for manuscript-scale
    uncertainty disclosure, but should not be described as
    independent-sample confidence intervals.
    """
    if n_bootstrap <= 0:
        return {}
    rng = np.random.default_rng(seed)
    residuals = y_fit - y_hat
    estimates = []
    for _ in range(n_bootstrap):
        sampled = rng.choice(residuals, size=residuals.size, replace=True)
        y_star = y_hat + sampled
        try:
            p_star, _ = curve_fit(
                kww_model,
                t_rel,
                y_star,
                p0=popt,
                bounds=bounds,
                maxfev=50000,
            )
            estimates.append(p_star)
        except Exception:
            continue
    if len(estimates) < max(10, int(0.1 * n_bootstrap)):
        return {}
    arr = np.vstack(estimates)
    names = ["y_inf", "amplitude", "tau_fs", "beta"]
    return {
        name: [float(v) for v in np.percentile(arr[:, idx], [2.5, 97.5])]
        for idx, name in enumerate(names)
    }


def fit_single_metric(
    metric: str,
    data_layer: str,
    t_fs: np.ndarray,
    y: np.ndarray,
    fit_start_fs: float,
    fit_end_fs: Optional[float],
    n_bootstrap: int,
    seed: int,
) -> FitResult:
    """Fit a single time-resolved observable to the KWW model."""
    try:
        t_rel, y_fit = finite_window(t_fs, y, fit_start_fs, fit_end_fs)
        p0, bounds = initial_guess_and_bounds(t_rel, y_fit)
        popt, pcov = curve_fit(
            kww_model, t_rel, y_fit, p0=p0, bounds=bounds, maxfev=100000
        )
        y_hat = kww_model(t_rel, *popt)
        residuals = y_fit - y_hat
        ss_res = float(np.sum(residuals**2))
        ss_tot = float(np.sum((y_fit - np.mean(y_fit)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        rmse = float(np.sqrt(np.mean(residuals**2)))
        aic, bic = information_criteria(y_fit, residuals, n_params=4)
        stderr_values = (
            np.sqrt(np.clip(np.diag(pcov), 0.0, np.inf))
            if pcov.size
            else np.full(4, np.nan)
        )
        stderr = {
            name: float(stderr_values[idx])
            for idx, name in enumerate(["y_inf", "amplitude", "tau_fs", "beta"])
        }
        # Locate the actual start time in original axis
        t_orig_start = float(
            t_fs[np.isfinite(t_fs) & (t_fs >= fit_start_fs)][0]
        )
        ci95 = bootstrap_kww(
            t_rel, y_fit, y_hat, popt, bounds,
            n_bootstrap=n_bootstrap, seed=seed,
        )
        return FitResult(
            metric=metric,
            data_layer=data_layer,
            n_points=int(y_fit.size),
            fit_start_fs=t_orig_start,
            fit_end_fs=float(t_rel[-1] + t_orig_start),
            y_inf=float(popt[0]),
            amplitude=float(popt[1]),
            tau_fs=float(popt[2]),
            beta=float(popt[3]),
            rmse=rmse,
            r2=float(r2),
            aic=aic,
            bic=bic,
            stderr=stderr,
            bootstrap_ci95=ci95,
            status="ok",
            warning=(
                "Residual bootstrap intervals are exploratory for "
                "autocorrelated HEOM time series."
            ),
        )
    except Exception as exc:
        return FitResult(
            metric=metric,
            data_layer=data_layer,
            n_points=0,
            fit_start_fs=float(fit_start_fs),
            fit_end_fs=float(fit_end_fs) if fit_end_fs is not None else float("nan"),
            y_inf=float("nan"),
            amplitude=float("nan"),
            tau_fs=float("nan"),
            beta=float("nan"),
            rmse=float("nan"),
            r2=float("nan"),
            aic=float("nan"),
            bic=float("nan"),
            stderr={},
            bootstrap_ci95={},
            status="failed",
            warning=str(exc),
        )


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_population_master(
    path: Path, init_site_index: int
) -> tuple:
    """Load the 30 ps population master dataset and compute observables.

    Returns:
        (t_fs, metrics_dict, metadata_dict)
    """
    data = np.load(path, allow_pickle=False)
    t_fs = np.asarray(data["tlist"], dtype=float)
    populations = np.asarray(data["populations"], dtype=float)
    if populations.ndim != 2:
        raise ValueError(
            f"Expected populations with shape (n_sites, n_times), "
            f"got {populations.shape}."
        )
    # Handle both (n_sites, n_times) and (n_times, n_sites)
    if populations.shape[1] != t_fs.size and populations.shape[0] == t_fs.size:
        populations = populations.T
    if populations.shape[1] != t_fs.size:
        raise ValueError("Population/time axes are incompatible.")

    p_sum = np.sum(populations, axis=0)
    diagonal_purity = np.sum(populations**2, axis=0)
    participation_ratio = 1.0 / np.maximum(diagonal_purity, 1e-15)
    clipped = np.clip(populations, 1e-15, None)
    shannon_entropy = -np.sum(clipped * np.log(clipped), axis=0)
    safe_idx = min(int(init_site_index), populations.shape[0] - 1)
    metrics = {
        "population_diagonal_purity": diagonal_purity,
        "participation_ratio": participation_ratio,
        "initial_site_population": populations[safe_idx, :],
        "population_shannon_entropy": shannon_entropy,
        "population_trace": p_sum,
    }
    metadata = {
        "source": str(path),
        "sha256": sha256_file(path),
        "n_sites": int(populations.shape[0]),
        "n_times": int(populations.shape[1]),
        "final_time_fs": float(t_fs[-1]),
        "data_layer": "population_diagonal_production_30ps",
        "init_site_index": int(init_site_index),
        "trace_min": float(np.min(p_sum)),
        "trace_max": float(np.max(p_sum)),
    }
    return t_fs, metrics, metadata


def qobj_to_array(rho: object) -> np.ndarray:
    """Convert a QuTiP Qobj-like density matrix to a dense complex ndarray."""
    if hasattr(rho, "full"):
        return np.asarray(rho.full(), dtype=complex)
    return np.asarray(rho, dtype=complex)


def load_density_pickle(
    path: Path, init_site_index: int
) -> tuple:
    """Load a HEOM pickle containing full density matrices and compute observables.

    Returns:
        (t_fs, metrics_dict, metadata_dict)
    """
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    t_fs = np.asarray(payload["tlist"], dtype=float)
    rho_t = payload["rho_t"]

    quantum_purity: list = []
    l1_coherence: list = []
    diagonal_purity: list = []
    initial_population: list = []

    for rho in rho_t:
        arr = qobj_to_array(rho)
        diagonal = np.real(np.diag(arr))
        quantum_purity.append(float(np.real(np.trace(arr @ arr))))
        l1_coherence.append(
            float(np.sum(np.abs(arr)) - np.sum(np.abs(np.diag(arr))))
        )
        diagonal_purity.append(float(np.sum(diagonal**2)))
        safe_idx = min(int(init_site_index), diagonal.size - 1)
        initial_population.append(float(diagonal[safe_idx]))

    metrics = {
        "quantum_purity": np.asarray(quantum_purity, dtype=float),
        "l1_coherence": np.asarray(l1_coherence, dtype=float),
        "density_diagonal_purity": np.asarray(diagonal_purity, dtype=float),
        "density_initial_site_population": np.asarray(initial_population, dtype=float),
    }
    metadata = {
        "source": str(path),
        "sha256": sha256_file(path),
        "n_times": int(t_fs.size),
        "final_time_fs": float(t_fs[-1]),
        "data_layer": "full_density_pickle",
        "init_site_index": int(init_site_index),
    }
    return t_fs, metrics, metadata


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_timeseries_csv(path: Path, series: list) -> None:
    """Write all computed observables in long CSV format."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["dataset", "data_layer", "metric", "time_fs", "value"])
        for dataset, layer, t_fs, metrics in series:
            for metric_name, values in metrics.items():
                for time_value, metric_value in zip(t_fs, values):
                    writer.writerow([
                        dataset,
                        layer,
                        metric_name,
                        f"{float(time_value):.12g}",
                        f"{float(metric_value):.12g}",
                    ])


def result_to_dict(result: FitResult) -> dict:
    """Convert a dataclass fit result into a JSON-serializable dictionary."""
    return {
        "metric": result.metric,
        "data_layer": result.data_layer,
        "n_points": result.n_points,
        "fit_start_fs": result.fit_start_fs,
        "fit_end_fs": result.fit_end_fs,
        "parameters": {
            "y_inf": result.y_inf,
            "amplitude": result.amplitude,
            "tau_fs": result.tau_fs,
            "tau_ps": (
                result.tau_fs / 1000.0
                if math.isfinite(result.tau_fs)
                else float("nan")
            ),
            "beta": result.beta,
        },
        "stderr": result.stderr,
        "bootstrap_ci95": result.bootstrap_ci95,
        "fit_quality": {
            "rmse": result.rmse,
            "r2": result.r2,
            "aic": result.aic,
            "bic": result.bic,
        },
        "status": result.status,
        "warning": result.warning,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_fits(
    fig_prefix: Path,
    fit_results: list,
    series_lookup: dict,
) -> None:
    """Generate a compact multi-panel figure for successful KWW fits."""
    import matplotlib.pyplot as plt  # lazy import — optional dependency

    successful = [res for res in fit_results if res.status == "ok"]
    if not successful:
        return
    fig_prefix.parent.mkdir(parents=True, exist_ok=True)
    n = len(successful)
    n_cols = 2
    n_rows = int(math.ceil(n / n_cols))
    plt.rcParams.update({
        "font.size": 9,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "figure.dpi": 300,
        "axes.linewidth": 0.8,
    })
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(10.5, 3.6 * n_rows), squeeze=False
    )
    for ax, result in zip(axes.ravel(), successful):
        key = (result.data_layer, result.metric)
        t_fs, values = series_lookup[key]
        ax.plot(t_fs / 1000.0, values, color="#1f77b4", lw=1.4, label="HEOM data")
        try:
            t_fit, y_sel = finite_window(
                t_fs, values, result.fit_start_fs, result.fit_end_fs
            )
            y_hat = kww_model(
                t_fit,
                result.y_inf,
                result.amplitude,
                result.tau_fs,
                result.beta,
            )
            ax.plot(
                (t_fit + result.fit_start_fs) / 1000.0,
                y_hat,
                color="#d62728",
                lw=1.6,
                label="KWW fit",
            )
        except Exception:
            pass
        beta_str = f"{result.beta:.3f}"
        tau_str = f"{result.tau_fs / 1000.0:.3g}"
        r2_str = f"{result.r2:.3f}"
        ax.set_title(
            f"{result.metric}\n"
            r"$\beta$=" + beta_str + r", $\tau$=" + tau_str + r" ps, $R^2$=" + r2_str
        )
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel("Value")
        ax.grid(True, alpha=0.35, linestyle=":")
        ax.legend(frameon=False, fontsize=8)
    for ax in axes.ravel()[len(successful):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(str(fig_prefix) + ".png", bbox_inches="tight")
    fig.savefig(str(fig_prefix) + ".pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description=(
            "Fit KWW relaxation models to HEOM trajectory diagnostics. "
            "Reads master_results.npz (populations) and optionally a "
            "full-density HEOM pickle."
        )
    )
    parser.add_argument(
        "--project-root", default=None,
        help="Optional KwanTube project root override.",
    )
    parser.add_argument(
        "--master-npz",
        default="outputs_data/raw_npz/master_results.npz",
        help="Population production NPZ path relative to project root.",
    )
    parser.add_argument(
        "--density-pkl",
        default="outputs_data/raw_pkl/heom_1JFF.pkl",
        help="Optional full-density HEOM pickle path relative to project root.",
    )
    parser.add_argument(
        "--skip-density-pkl", action="store_true",
        help="Skip full-density pickle analysis even if the file exists.",
    )
    parser.add_argument(
        "--init-site-index", type=int, default=5,
        help="Initially excited site index used in the 1JFF production run.",
    )
    parser.add_argument(
        "--fit-start-fs", type=float, default=0.0,
        help="Start time for fitting window in femtoseconds.",
    )
    parser.add_argument(
        "--fit-end-fs", type=float, default=None,
        help="End time for fitting window in femtoseconds (default: full range).",
    )
    parser.add_argument(
        "--bootstrap", type=int, default=200,
        help="Residual-bootstrap replicates for exploratory parameter intervals.",
    )
    parser.add_argument(
        "--seed", type=int, default=12345,
        help="Random seed for bootstrap reproducibility.",
    )
    parser.add_argument(
        "--output-json",
        default="outputs_data/raw_json/metrics/heom_kww_relaxation_fit.json",
        help="Output JSON path relative to project root.",
    )
    parser.add_argument(
        "--output-csv",
        default="outputs_data/raw_csv/heom_kww_relaxation_timeseries.csv",
        help="Output CSV path relative to project root.",
    )
    parser.add_argument(
        "--fig-prefix",
        default="outputs_data/figures_final/heom_kww_relaxation",
        help="Figure prefix relative to project root (no extension).",
    )
    args = parser.parse_args()

    project_root = resolve_project_root(args.project_root)
    audit_log(project_root, "[RUN_AUDIT] START script=fit_heom_kww_relaxation.py")
    start = time.time()

    # ------------------------------------------------------------------
    # Load datasets
    # ------------------------------------------------------------------
    datasets = []  # (name, layer, t_fs, metrics, metadata)

    master_path = project_root / args.master_npz
    if not master_path.exists():
        raise FileNotFoundError(
            f"Production master dataset not found: {master_path}"
        )
    t_master, metrics_master, meta_master = load_population_master(
        master_path, init_site_index=args.init_site_index
    )
    datasets.append((
        "master_results",
        str(meta_master["data_layer"]),
        t_master,
        metrics_master,
        meta_master,
    ))
    print(
        f"[INFO] master_results.npz loaded: "
        f"{meta_master['n_sites']} sites x {meta_master['n_times']} points, "
        f"final_t={meta_master['final_time_fs']:.0f} fs"
    )

    density_path = project_root / args.density_pkl
    density_warning = ""
    if not args.skip_density_pkl and density_path.exists():
        try:
            t_density, metrics_density, meta_density = load_density_pickle(
                density_path, init_site_index=args.init_site_index
            )
            datasets.append((
                "heom_density_pickle",
                str(meta_density["data_layer"]),
                t_density,
                metrics_density,
                meta_density,
            ))
            print(
                f"[INFO] heom_1JFF.pkl loaded: "
                f"{meta_density['n_times']} time points, "
                f"final_t={meta_density['final_time_fs']:.0f} fs"
            )
        except Exception as exc:
            density_warning = (
                f"Full-density pickle analysis skipped after load failure: {exc}"
            )
            print(f"[WARN] {density_warning}")
    elif not args.skip_density_pkl:
        density_warning = f"Full-density pickle not found: {density_path}"
        print(f"[WARN] {density_warning}")

    # ------------------------------------------------------------------
    # Fit KWW to each observable
    # ------------------------------------------------------------------
    fit_results = []
    series_lookup = {}  # (layer, metric) -> (t_fs, values)

    for _, layer, t_fs, metrics, _ in datasets:
        for metric_name, values in metrics.items():
            if metric_name == "population_trace":
                # Constant by conservation; skip
                continue
            result = fit_single_metric(
                metric=metric_name,
                data_layer=layer,
                t_fs=t_fs,
                y=values,
                fit_start_fs=args.fit_start_fs,
                fit_end_fs=args.fit_end_fs,
                n_bootstrap=args.bootstrap,
                seed=args.seed,
            )
            fit_results.append(result)
            series_lookup[(layer, metric_name)] = (t_fs, values)

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    output_csv = project_root / args.output_csv
    write_timeseries_csv(
        output_csv,
        [(name, layer, t_fs, metrics) for name, layer, t_fs, metrics, _ in datasets],
    )

    fig_prefix = project_root / args.fig_prefix
    try:
        plot_fits(fig_prefix, fit_results, series_lookup)
    except ImportError:
        print("[WARN] matplotlib not available; figures skipped.")

    report = {
        "script": "fit_heom_kww_relaxation.py",
        "version": "1.1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(project_root),
        "important_interpretation_note": (
            "The 30 ps master_results.npz layer contains populations only; "
            "its fitted purity is population-diagonal purity, not full "
            "density-matrix purity. Full quantum purity is reported only for "
            "pickle datasets containing rho_t. A stretched exponent beta < 1 "
            "indicates deviation from simple Markovian kinetics and is "
            "presented as a lower-bound descriptor of non-Markovian relaxation "
            "complexity, not as evidence of a glassy bath."
        ),
        "density_pickle_warning": density_warning,
        "datasets": [metadata for _, _, _, _, metadata in datasets],
        "fits": [result_to_dict(result) for result in fit_results],
        "outputs": {
            "json": str(project_root / args.output_json),
            "csv": str(output_csv),
            "figure_png": str(fig_prefix) + ".png",
            "figure_pdf": str(fig_prefix) + ".pdf",
        },
        "elapsed_seconds": time.time() - start,
    }
    output_json = project_root / args.output_json
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------
    print("=" * 78)
    print("HEOM KWW RELAXATION FIT COMPLETE")
    print("=" * 78)
    for result in fit_results:
        if result.status == "ok":
            print(
                f"  [{result.data_layer[:28]:28s}] "
                f"{result.metric[:32]:32s} "
                f"beta={result.beta:.4f}  "
                f"tau={result.tau_fs / 1000.0:.4g} ps  "
                f"R2={result.r2:.4f}"
            )
        else:
            print(
                f"  [{result.data_layer[:28]:28s}] "
                f"{result.metric[:32]:32s} FAILED: {result.warning}"
            )
    if density_warning:
        print(f"[WARN] {density_warning}")
    print(f"JSON : {output_json}")
    print(f"CSV  : {output_csv}")
    print(f"FIG  : {str(fig_prefix)}.png")
    audit_log(
        project_root,
        f"[RUN_AUDIT] END status=ok output={output_json.name} "
        f"elapsed={time.time() - start:.1f}s",
    )
    return 0


if __name__ == "__main__":
    from pathlib import Path as _P
    import sys as _sys
    for _parent in _P(__file__).resolve().parents:
        if (_parent / "qmc_mt" / "run_audit.py").exists():
            _sys.path.insert(0, str(_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    _sys.exit(main())
