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
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/scripts/ a la raiz

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

def write_outputs(out, groups, jump_fits, hier, level_rows, thresholds, draws):
    out.mkdir(parents=True, exist_ok=True)
    with (out / "group_loglinear_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["group", "observable", "time_fs", "logr_mean", "r_mean", "beta_mean"])
        for fit in jump_fits:
            w.writerow([fit["name"], fit["observable"], fit["time_fs"], float(np.mean(fit["logr_draw"])), float(np.mean(fit["r_draw"])), float(np.mean(fit["beta_draw"]))])
    with (out / "extrapolated_jumps.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["group", "observable", "nc", "pred_mean"])
        for fit in jump_fits:
            a, b, n0 = fit["logA_draw"], fit["logr_draw"], int(fit["n0"])
            for nc in range(int(fit["nc_min"]), int(fit["nc_max"]) + 4):
                pred = np.exp(a + b * (nc - n0))
                w.writerow([fit["name"], fit["observable"], nc, float(np.mean(pred))])
    print(f"Outputs written to: {out.resolve()}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", default=str(PROJECT_ROOT / "outputs_data" / "raw_csv" / "heom_bayes_input_current.csv"))
    ap.add_argument("--output-dir", default="outputs_data/heom_bayes_out_v2")
    ap.add_argument("--draws", type=int, default=10000)
    args = ap.parse_args()
    rng = np.random.default_rng(20260424)
    if not Path(args.input_csv).exists():
        print(f"ERROR: Input CSV {args.input_csv} not found.")
        return
    groups = load_groups(args.input_csv)
    jump_groups = [g for g in groups if g.is_jump]
    jump_fits = [fit_jump_group_mc(g, args.draws, rng, 0.03, 0.75) for g in jump_groups]
    hier = hierarchy_grid(jump_fits, args.draws, rng)
    level_rows = level_reference_checks(groups, args.draws, rng)
    out = PROJECT_ROOT / args.output_dir
    write_outputs(out, groups, jump_fits, hier, level_rows, [0.01, 0.005], args.draws)

if __name__ == "__main__":
    main()
