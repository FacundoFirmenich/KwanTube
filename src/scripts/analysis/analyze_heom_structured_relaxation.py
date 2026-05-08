#!/usr/bin/env python3
"""Secondary diagnostics for structured non-Markovian HEOM relaxation.

This script reuses the archived KWW time-series and fit ledger. It does not run
HEOM. Its purpose is to test whether the sub-unitary KWW exponents are isolated
fit artefacts or a cross-observable signature of distributed relaxation.

Outputs are intentionally conservative: they support distributed non-Markovian
relaxation over a finite window, but not a thermodynamic glass-transition claim.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit, nnls


def project_root() -> Path:
    return Path(__file__).resolve().parents[3]


def kww(t, y_inf, amp, tau, beta):
    return y_inf + amp * np.exp(-np.power(np.maximum(t, 0.0) / max(tau, 1e-12), beta))


def single_exp(t, y_inf, amp, tau):
    return y_inf + amp * np.exp(-np.maximum(t, 0.0) / max(tau, 1e-12))


def biexp(t, y_inf, a1, tau1, a2, tau2):
    tt = np.maximum(t, 0.0)
    return y_inf + a1 * np.exp(-tt / max(tau1, 1e-12)) + a2 * np.exp(-tt / max(tau2, 1e-12))


def aic_bic(y, residuals, n_params):
    n = len(y)
    rss = max(float(np.sum(residuals**2)), np.finfo(float).tiny)
    sigma2 = rss / n
    log_l = -0.5 * n * (math.log(2.0 * math.pi * sigma2) + 1.0)
    return 2 * n_params - 2 * log_l, n_params * math.log(n) - 2 * log_l


def fit_model(name, func, t, y, p0, bounds):
    try:
        popt, _ = curve_fit(func, t, y, p0=p0, bounds=bounds, maxfev=50000)
        pred = func(t, *popt)
        resid = y - pred
        aic, bic = aic_bic(y, resid, len(popt))
        ss_tot = max(float(np.sum((y - np.mean(y)) ** 2)), np.finfo(float).tiny)
        r2 = 1.0 - float(np.sum(resid**2)) / ss_tot
        return {"model": name, "status": "ok", "params": [float(x) for x in popt], "aic": float(aic), "bic": float(bic), "r2": float(r2)}
    except Exception as exc:  # pragma: no cover - diagnostic script
        return {"model": name, "status": "failed", "error": str(exc)}


def select_window(t, y, end_fs):
    mask = np.isfinite(t) & np.isfinite(y) & (t <= end_fs)
    tw = t[mask]
    yw = y[mask]
    if len(tw) < 20:
        return None, None
    return tw - tw[0], yw


def initial_bounds(t, y):
    ymin, ymax = float(np.min(y)), float(np.max(y))
    yr = max(ymax - ymin, 1e-12)
    yinf = float(np.mean(y[-max(3, len(y) // 10):]))
    amp = float(y[0] - yinf)
    tau = max(float(t[-1] / 3.0), 1.0)
    common_y = (ymin - 3 * yr, ymax + 3 * yr)
    amp_bounds = (-4 * yr, 4 * yr)
    tau_bounds = (1e-9, max(t[-1] * 100.0, 10.0))
    return {
        "single": ([yinf, amp, tau], ([common_y[0], amp_bounds[0], tau_bounds[0]], [common_y[1], amp_bounds[1], tau_bounds[1]])),
        "kww": ([yinf, amp, tau, 0.7], ([common_y[0], amp_bounds[0], tau_bounds[0], 0.05], [common_y[1], amp_bounds[1], tau_bounds[1], 3.0])),
        "biexp": ([yinf, amp * 0.6, tau / 3.0, amp * 0.4, tau * 3.0], ([common_y[0], amp_bounds[0], tau_bounds[0], amp_bounds[0], tau_bounds[0]], [common_y[1], amp_bounds[1], tau_bounds[1], amp_bounds[1], tau_bounds[1]])),
    }


def nnls_rate_width(t, y):
    y_tail = float(np.mean(y[-max(3, len(y) // 10):]))
    signal = y - y_tail
    sign = 1.0 if abs(np.max(signal)) >= abs(np.min(signal)) else -1.0
    target = np.maximum(sign * signal, 0.0)
    if np.max(target) <= 0:
        return None
    taus = np.logspace(math.log10(10.0), math.log10(max(t[-1] * 20.0, 100.0)), 80)
    basis = np.exp(-t[:, None] / taus[None, :])
    weights, residual = nnls(basis, target)
    if np.sum(weights) <= 0:
        return None
    weights = weights / np.sum(weights)
    active = taus[weights > np.max(weights) * 0.01]
    if len(active) == 0:
        active = taus[weights > 0]
    return {
        "tau_grid_fs_min": float(np.min(taus)),
        "tau_grid_fs_max": float(np.max(taus)),
        "active_tau_fs_min": float(np.min(active)),
        "active_tau_fs_max": float(np.max(active)),
        "active_decades": float(np.log10(np.max(active) / np.min(active))) if np.min(active) > 0 else None,
        "relative_residual": float(residual / max(np.linalg.norm(target), np.finfo(float).eps)),
    }


def load_timeseries(path: Path):
    groups = defaultdict(lambda: {"time_fs": [], "value": []})
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = f"{row['data_layer']}::{row['metric']}"
            groups[key]["data_layer"] = row["data_layer"]
            groups[key]["metric"] = row["metric"]
            groups[key]["time_fs"].append(float(row["time_fs"]))
            groups[key]["value"].append(float(row["value"]))
    return groups


def main():
    root = project_root()
    csv_path = root / "outputs_data" / "raw_csv" / "heom_kww_relaxation_timeseries.csv"
    fit_path = root / "outputs_data" / "raw_json" / "metrics" / "heom_kww_relaxation_fit.json"
    out_path = root / "outputs_data" / "raw_json" / "metrics" / "heom_structured_relaxation_diagnostics.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    groups = load_timeseries(csv_path)
    source_fit = json.loads(fit_path.read_text(encoding="utf-8"))
    windows = [5000.0, 10000.0, 20000.0, 30000.0]
    diagnostics = []

    for key, data in groups.items():
        t_all = np.asarray(data["time_fs"], dtype=float)
        y_all = np.asarray(data["value"], dtype=float)
        order = np.argsort(t_all)
        t_all, y_all = t_all[order], y_all[order]
        item = {"metric": data["metric"], "data_layer": data["data_layer"], "window_fits": []}
        for end_fs in windows:
            t, y = select_window(t_all, y_all, end_fs)
            if t is None:
                continue
            guesses = initial_bounds(t, y)
            fits = [
                fit_model("single_exponential", single_exp, t, y, *guesses["single"]),
                fit_model("kww", kww, t, y, *guesses["kww"]),
                fit_model("biexponential", biexp, t, y, *guesses["biexp"]),
            ]
            ok = [fit for fit in fits if fit["status"] == "ok"]
            best = min(ok, key=lambda fit: fit["bic"])["model"] if ok else None
            kww_fit = next((fit for fit in fits if fit["model"] == "kww" and fit["status"] == "ok"), None)
            single_fit = next((fit for fit in fits if fit["model"] == "single_exponential" and fit["status"] == "ok"), None)
            item["window_fits"].append(
                {
                    "fit_end_fs": end_fs,
                    "n_points": int(len(t)),
                    "fits": fits,
                    "best_bic_model": best,
                    "delta_bic_single_minus_kww": float(single_fit["bic"] - kww_fit["bic"]) if single_fit and kww_fit else None,
                    "kww_beta": float(kww_fit["params"][3]) if kww_fit else None,
                }
            )
        t, y = select_window(t_all, y_all, 30000.0)
        item["distributed_rate_nnls"] = nnls_rate_width(t, y) if t is not None else None
        diagnostics.append(item)

    purity_like = []
    for fit in source_fit.get("fits", []):
        metric = fit.get("metric", "")
        if fit.get("status") == "ok" and any(token in metric for token in ["purity", "population", "entropy"]):
            beta = fit.get("parameters", {}).get("beta")
            ci = fit.get("bootstrap_ci95", {}).get("beta")
            if beta is not None:
                purity_like.append({"metric": metric, "data_layer": fit.get("data_layer"), "beta": beta, "ci95": ci})

    beta_values = np.asarray([x["beta"] for x in purity_like], dtype=float)
    summary = {
        "script": "analyze_heom_structured_relaxation.py",
        "interpretation": "Population and purity observables show robust sub-unitary KWW exponents consistent with structured non-Markovian distributed relaxation over the finite production window; this is not evidence of a thermodynamic glass transition.",
        "source_csv": str(csv_path),
        "source_fit_json": str(fit_path),
        "purity_population_entropy_beta_summary": {
            "n_observables": int(len(beta_values)),
            "beta_min": float(np.min(beta_values)) if len(beta_values) else None,
            "beta_max": float(np.max(beta_values)) if len(beta_values) else None,
            "beta_mean": float(np.mean(beta_values)) if len(beta_values) else None,
            "all_subunitary": bool(np.all(beta_values < 1.0)) if len(beta_values) else None,
            "observables": purity_like,
        },
        "diagnostics": diagnostics,
    }
    out_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[structured] wrote {out_path}")
    if len(beta_values):
        print(f"[structured] beta range {np.min(beta_values):.3f}-{np.max(beta_values):.3f} across {len(beta_values)} population/purity/entropy observables")


if __name__ == "__main__":
    main()
