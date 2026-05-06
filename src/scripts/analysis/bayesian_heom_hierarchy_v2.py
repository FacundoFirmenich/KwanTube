#!/usr/bin/env python3
from __future__ import annotations
"""
Bayesian HEOM hierarchy v2: stable contraction model for small-N HEOM validation data.
"""

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False

EPS = 1e-15
LOG2PI = math.log(2.0 * math.pi)

# Boilerplate para resolver rutas Tier-0
PROJECT_ROOT = Path(__file__).resolve().parents[3] # retrocede desde src/scripts/analysis/ a la raiz

def q(x: np.ndarray, p: float) -> float:
    return float(np.quantile(np.asarray(x, dtype=float), p))

def normal_logpdf_vec(x: np.ndarray, mu: float, sd: float) -> np.ndarray:
    if sd <= 0:
        return np.full_like(x, -np.inf, dtype=float)
    z = (x - mu) / sd
    return -0.5 * (LOG2PI + 2.0 * math.log(sd) + z * z)

@dataclass
class Row:
    observable: str
    group: str
    time_fs: float
    nc: int
    value: float
    sigma_obs: float
    reference: Optional[float]

@dataclass
class Group:
    name: str
    observable: str
    time_fs: float
    nc: np.ndarray
    value: np.ndarray
    sigma_obs: np.ndarray
    reference: Optional[float]
    n0: int
    is_jump: bool

def load_groups(path: str) -> List[Group]:
    rows: List[Row] = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ref_raw = row.get("reference", "")
            ref = None if ref_raw in (None, "") else float(ref_raw)
            obs = row["observable"].strip()
            gname = (row.get("group") or f"{obs}@{row['time_fs']}").strip()
            rows.append(Row(
                observable=obs,
                group=gname,
                time_fs=float(row["time_fs"]),
                nc=int(float(row["nc"])),
                value=float(row["value"]),
                sigma_obs=float(row.get("sigma_obs", 0.0) or 0.0),
                reference=ref,
            ))
    by: Dict[str, List[Row]] = {}
    for r in rows:
        by.setdefault(r.group, []).append(r)
    groups: List[Group] = []
    for name, vals in by.items():
        vals = sorted(vals, key=lambda z: z.nc)
        groups.append(Group(
            name=name,
            observable=vals[0].observable,
            time_fs=vals[0].time_fs,
            nc=np.array([v.nc for v in vals], dtype=int),
            value=np.array([v.value for v in vals], dtype=float),
            sigma_obs=np.array([v.sigma_obs for v in vals], dtype=float),
            reference=vals[0].reference,
            n0=min(v.nc for v in vals),
            is_jump=vals[0].observable.startswith("jump_"),
        ))
    return sorted(groups, key=lambda g: g.name)

def weighted_line_fit(x, y, sigma):
    X = np.column_stack([np.ones_like(x, dtype=float), x.astype(float)])
    w = 1.0 / np.maximum(sigma, 1e-12) ** 2
    XtW = X.T * w
    A = XtW @ X
    B = XtW @ y
    rank = int(np.linalg.matrix_rank(A))
    coef = np.linalg.lstsq(A, B, rcond=None)[0]
    pred = X @ coef
    resid = y - pred
    if len(y) > 2:
        rss = float(np.sum(w * resid * resid))
        dof = max(len(y) - 2, 1)
        resid_sd = math.sqrt(max(rss / dof, 0.0) / max(float(np.mean(w)), 1e-12))
    else:
        resid_sd = float("nan")
    return float(coef[0]), float(coef[1]), resid_sd, rank

def fit_jump_group_mc(g, draws, rng, min_log_sigma, max_log_sigma):
    if np.any(g.value <= 0):
        raise ValueError(f"Jump group {g.name} contains non-positive values; cannot log-fit.")
    x = g.nc.astype(float) - float(g.n0)
    logy = np.log(g.value)
    siglog = np.clip(g.sigma_obs / np.maximum(np.abs(g.value), EPS), min_log_sigma, max_log_sigma)
    a_hat, b_hat, resid_sd, rank = weighted_line_fit(x, logy, siglog)
    a_draw, b_draw = [], []
    attempts, max_attempts = 0, draws * 30
    while len(b_draw) < draws and attempts < max_attempts:
        attempts += 1
        ystar = rng.normal(logy, siglog)
        a, b, _, _ = weighted_line_fit(x, ystar, siglog)
        if b < 0.0 and np.isfinite(a) and np.isfinite(b):
            a_draw.append(a); b_draw.append(b)
    if len(b_draw) < draws:
        idx = rng.integers(0, len(b_draw), size=draws)
        a_draw, b_draw = np.asarray(a_draw)[idx], np.asarray(b_draw)[idx]
    return {
        "name": g.name, "observable": g.observable, "time_fs": g.time_fs, "n_points": len(g.nc),
        "n0": g.n0, "nc_min": int(np.min(g.nc)), "nc_max": int(np.max(g.nc)),
        "a_hat": a_hat, "b_hat": b_hat, "resid_sd_hat": resid_sd, "rank": rank,
        "logA_draw": np.asarray(a_draw), "logr_draw": np.asarray(b_draw),
        "r_draw": np.exp(np.asarray(b_draw)), "beta_draw": -np.asarray(b_draw),
        "siglog_median": float(np.median(siglog)), "siglog_max": float(np.max(siglog)),
    }

def hierarchy_grid(group_fits, draws, rng):
    m = np.array([float(np.mean(f["logr_draw"])) for f in group_fits])
    se = np.array([float(np.std(f["logr_draw"], ddof=1)) for f in group_fits])
    se = np.maximum(se, 1e-4)
    mu_grid = np.linspace(-1.8, -0.05, 500)
    tau_grid = np.linspace(0.001, 1.0, 500)
    MU, TAU = np.meshgrid(mu_grid, tau_grid, indexing="ij")
    logp = normal_logpdf_vec(MU, math.log(0.5), 0.70) - 0.5 * (TAU / 0.50) ** 2
    for mi, sei in zip(m, se):
        sd = np.sqrt(TAU * TAU + sei * sei)
        logp += -0.5 * (LOG2PI + 2.0 * np.log(sd) + ((mi - MU) / sd) ** 2)
    logp -= np.nanmax(logp)
    w = np.exp(logp); w /= np.sum(w)
    idx = rng.choice(w.size, size=draws, replace=True, p=w.ravel())
    mu_draw, tau_draw = MU.ravel()[idx], TAU.ravel()[idx]
    b_group = np.zeros((len(group_fits), draws), dtype=float)
    for gi, (mi, sei) in enumerate(zip(m, se)):
        prec = 1.0 / (sei * sei) + 1.0 / np.maximum(tau_draw, 1e-6) ** 2
        var = 1.0 / prec
        mean = var * (mi / (sei * sei) + mu_draw / np.maximum(tau_draw, 1e-6) ** 2)
        b_group[gi] = np.minimum(rng.normal(mean, np.sqrt(var)), -1e-8)
    return {"mu_logr": mu_draw, "global_r": np.exp(mu_draw), "tau_logr": tau_draw, "group_logr": b_group, "group_r": np.exp(b_group), "m_logr": m, "se_logr": se, "global_beta": -mu_draw, "group_beta": -b_group}

def level_reference_checks(groups, draws, rng):
    rows = []
    for g in groups:
        if g.is_jump: continue
        vals = np.vstack([rng.normal(y, max(float(sy), 1e-12), size=draws) for y, sy in zip(g.value, g.sigma_obs)])
        ref = float(g.reference) if g.reference is not None else float("nan")
        delta_last = vals[-1] - ref if np.isfinite(ref) else np.full(draws, np.nan)
        rows.append({"group": g.name, "observable": g.observable, "time_fs": g.time_fs, "n_points": len(g.nc), "reference": ref, "last_nc": int(g.nc[-1]), "last_value": float(g.value[-1]), "delta_last_mean": float(np.mean(delta_last)), "delta_last_q50": q(delta_last, 0.5)})
    return rows

def _summary_stats(values: np.ndarray) -> Tuple[float, float, float, float]:
    arr = np.asarray(values, dtype=float)
    return float(np.mean(arr)), q(arr, 0.025), q(arr, 0.50), q(arr, 0.975)


def _write_posterior_plots(fig_dir: Path, hier: Dict[str, np.ndarray]) -> None:
    """Write posterior diagnostic plots in PNG and PDF formats."""
    if not HAS_MPL:
        print("WARNING: matplotlib unavailable; posterior plots were not generated.")
        return
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    plot_specs = [
        ("global_r", "Global contraction r"),
        ("global_beta", "Global decay beta"),
        ("tau_logr", "Between-group heterogeneity tau"),
    ]
    for ax, (key, title) in zip(axes, plot_specs):
        vals = np.asarray(hier[key], dtype=float)
        ax.hist(vals, bins=50, color="#3B82F6", alpha=0.78, edgecolor="white")
        lo, med, hi = q(vals, 0.025), q(vals, 0.50), q(vals, 0.975)
        ax.axvline(med, color="#111827", lw=1.5, label="median")
        ax.axvspan(lo, hi, color="#F59E0B", alpha=0.22, label="95% CI")
        ax.set_title(title)
        ax.set_ylabel("posterior draws")
        ax.grid(alpha=0.25)
    axes[0].legend(loc="best", fontsize=8)
    fig.suptitle("Bayesian HEOM hierarchy v2 posterior diagnostics", fontsize=13)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out_path = fig_dir / f"posterior_plots_v2.{ext}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Posterior plot written to: {out_path.resolve()}")
    plt.close(fig)


def write_outputs(csv_dir, txt_dir, fig_dir, groups, jump_fits, hier, level_rows, thresholds, draws):
    csv_dir.mkdir(parents=True, exist_ok=True)
    txt_dir.mkdir(parents=True, exist_ok=True)
    with (csv_dir / "group_loglinear_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["group", "observable", "time_fs", "logr_mean", "logr_q025", "logr_q50", "logr_q975", "r_mean", "r_q025", "r_q50", "r_q975", "beta_mean", "beta_q025", "beta_q50", "beta_q975"])
        for fit in jump_fits:
            logr = _summary_stats(fit["logr_draw"])
            rvals = _summary_stats(fit["r_draw"])
            beta = _summary_stats(fit["beta_draw"])
            w.writerow([fit["name"], fit["observable"], fit["time_fs"], *logr, *rvals, *beta])
    with (csv_dir / "hierarchy_global_contraction.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "mean", "q025", "q50", "q975"])
        for key in ("global_r", "global_beta", "tau_logr"):
            w.writerow([key, *_summary_stats(hier[key])])
    with (csv_dir / "level_reference_checks.csv").open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["group", "observable", "time_fs", "n_points", "reference", "last_nc", "last_value", "delta_last_mean", "delta_last_q50"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(level_rows)
    with (csv_dir / "extrapolated_jumps.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["group", "observable", "nc", "pred_mean"])
        for fit in jump_fits:
            a, b, n0 = fit["logA_draw"], fit["logr_draw"], int(fit["n0"])
            for nc in range(int(fit["nc_min"]), int(fit["nc_max"]) + 4):
                pred = np.exp(a + b * (nc - n0))
                w.writerow([fit["name"], fit["observable"], nc, float(np.mean(pred))])
    diagnostics = txt_dir / "diagnostics_v2.txt"
    diagnostics.write_text(
        "Bayesian HEOM hierarchy v2 diagnostics\n"
        f"draws={draws}\n"
        f"n_groups={len(groups)}\n"
        f"n_jump_groups={len(jump_fits)}\n"
        f"thresholds={thresholds}\n"
        f"csv_dir={csv_dir.resolve()}\n"
        f"fig_dir={fig_dir.resolve()}\n",
        encoding="utf-8",
    )
    _write_posterior_plots(fig_dir, hier)
    print(f"CSV outputs written to: {csv_dir.resolve()}")
    print(f"Diagnostics written to: {diagnostics.resolve()}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", default=None)
    ap.add_argument("--csv-dir", default=str(PROJECT_ROOT / "outputs_data" / "raw_csv" / "heom_+bayesian_analysis"))
    ap.add_argument("--txt-dir", default=str(PROJECT_ROOT / "outputs_data" / "raw_txt+md" / "reports"))
    ap.add_argument("--fig-dir", default=str(PROJECT_ROOT / "outputs_data" / "figures_final"))
    ap.add_argument("--output-dir", default=None, help="Deprecated; ignored. Outputs are written to standard outputs_data subdirectories.")
    ap.add_argument("--draws", type=int, default=10000)
    args = ap.parse_args()
    rng = np.random.default_rng(20260424)
    input_csv = Path(args.input_csv) if args.input_csv else _default_input_csv()
    if not input_csv.exists():
        print(f"ERROR: Input CSV {input_csv} not found.")
        return
    groups = load_groups(str(input_csv))
    jump_groups = [g for g in groups if g.is_jump]
    jump_fits = [fit_jump_group_mc(g, args.draws, rng, 0.03, 0.75) for g in jump_groups]
    hier = hierarchy_grid(jump_fits, args.draws, rng)
    level_rows = level_reference_checks(groups, args.draws, rng)
    write_outputs(Path(args.csv_dir), Path(args.txt_dir), Path(args.fig_dir), groups, jump_fits, hier, level_rows, [0.01, 0.005], args.draws)


def _default_input_csv() -> Path:
    """Return the canonical HEOM Bayes input CSV."""
    return PROJECT_ROOT / "outputs_data" / "raw_csv" / "heom_+bayesian_analysis" / "heom_bayes_input_current.csv"

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
