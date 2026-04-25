#!/usr/bin/env python3
from __future__ import annotations

"""
Bayesian HEOM hierarchy v2: stable contraction model for small-N HEOM validation data.

Why v2 exists
-------------
The previous raw-level nonlinear hierarchy tried to estimate theta, alpha, beta and sigma
inside each group, even for two-point groups. That posterior is weakly identified and can
produce very poor R-hat. This version uses the structure of the HEOM validation problem:

1) Jump/convergence observables are positive and should contract toward zero:
       y_{g,n} > 0,      y_{g,n} ~= A_g * r_g^(n-n0_g),     0 < r_g < 1
   so we fit them on the log scale.

2) The population-level Richardson check is not pooled with jump magnitudes. It is treated
   as an independent reference/asymptote check.

3) The hierarchy is applied to contraction slopes log(r_g), not to raw observables with
   incomparable scales. Group-level log(r_g) estimates are combined with a Bayesian
   random-effects meta-model on a deterministic grid, avoiding fragile MCMC.

Dependencies: Python stdlib + numpy + matplotlib optional.

Input CSV columns:
    observable,time_fs,nc,value,group,sigma_obs,reference

Outputs:
    group_loglinear_summary.csv
    hierarchy_global_contraction.csv
    hierarchy_group_shrinkage.csv
    extrapolated_jumps.csv
    level_reference_checks.csv
    diagnostics_v2.txt
    posterior_plots_v2.png   (if matplotlib is available)
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


def q(x: np.ndarray, p: float) -> float:
    return float(np.quantile(np.asarray(x, dtype=float), p))


def normal_logpdf_vec(x: np.ndarray, mu: float, sd: float) -> np.ndarray:
    if sd <= 0:
        return np.full_like(x, -np.inf, dtype=float)
    z = (x - mu) / sd
    return -0.5 * (LOG2PI + 2.0 * math.log(sd) + z * z)


def normal_logpdf_scalar(x: float, mu: float, sd: float) -> float:
    if sd <= 0:
        return -np.inf
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


def weighted_line_fit(x: np.ndarray, y: np.ndarray, sigma: np.ndarray) -> Tuple[float, float, float, int]:
    """Weighted linear regression y = a + b x. Returns a,b,resid_sd,rank."""
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
        # Weighted residual scatter on log scale. Conservative floor is applied elsewhere.
        rss = float(np.sum(w * resid * resid))
        dof = max(len(y) - 2, 1)
        resid_sd = math.sqrt(max(rss / dof, 0.0) / max(float(np.mean(w)), 1e-12))
    else:
        resid_sd = float("nan")
    return float(coef[0]), float(coef[1]), resid_sd, rank


def log_sigma_from_value(y: np.ndarray, sy: np.ndarray, min_log_sigma: float, max_log_sigma: float) -> np.ndarray:
    # Delta-method: sd(log y) ~= sd(y)/y. Add a floor because sigma_obs is often a numerical
    # proxy rather than a calibrated measurement error.
    raw = sy / np.maximum(np.abs(y), EPS)
    return np.clip(raw, min_log_sigma, max_log_sigma)


def fit_jump_group_mc(
    g: Group,
    draws: int,
    rng: np.random.Generator,
    min_log_sigma: float,
    max_log_sigma: float,
) -> Dict[str, np.ndarray | float | int | str]:
    if np.any(g.value <= 0):
        raise ValueError(f"Jump group {g.name} contains non-positive values; cannot log-fit.")
    x = g.nc.astype(float) - float(g.n0)
    logy = np.log(g.value)
    siglog = log_sigma_from_value(g.value, g.sigma_obs, min_log_sigma, max_log_sigma)
    a_hat, b_hat, resid_sd, rank = weighted_line_fit(x, logy, siglog)

    # Monte Carlo posterior approximation: perturb log observations by their calibrated log error,
    # fit line, retain convergent fits b<0. This is intentionally conservative for tiny groups.
    a_draw = []
    b_draw = []
    attempts = 0
    max_attempts = draws * 30
    while len(b_draw) < draws and attempts < max_attempts:
        attempts += 1
        ystar = rng.normal(logy, siglog)
        a, b, _, _ = weighted_line_fit(x, ystar, siglog)
        if b < 0.0 and np.isfinite(a) and np.isfinite(b):
            a_draw.append(a)
            b_draw.append(b)
    if len(b_draw) < max(100, draws // 20):
        raise RuntimeError(f"Could not obtain enough convergent MC fits for {g.name}.")
    if len(b_draw) < draws:
        # Resample with replacement to requested draw count.
        idx = rng.integers(0, len(b_draw), size=draws)
        a_draw = np.asarray(a_draw)[idx]
        b_draw = np.asarray(b_draw)[idx]
    else:
        a_draw = np.asarray(a_draw)
        b_draw = np.asarray(b_draw)

    r_draw = np.exp(b_draw)
    beta_draw = -b_draw
    return {
        "name": g.name,
        "observable": g.observable,
        "time_fs": g.time_fs,
        "n_points": len(g.nc),
        "n0": g.n0,
        "nc_min": int(np.min(g.nc)),
        "nc_max": int(np.max(g.nc)),
        "a_hat": a_hat,
        "b_hat": b_hat,
        "resid_sd_hat": resid_sd,
        "rank": rank,
        "logA_draw": a_draw,
        "logr_draw": b_draw,
        "r_draw": r_draw,
        "beta_draw": beta_draw,
        "siglog_median": float(np.median(siglog)),
        "siglog_max": float(np.max(siglog)),
    }


def hierarchy_grid(
    group_fits: List[Dict[str, np.ndarray | float | int | str]],
    draws: int,
    rng: np.random.Generator,
) -> Dict[str, np.ndarray]:
    # Data are group estimates of b_g = log r_g < 0.
    m = np.array([float(np.mean(f["logr_draw"])) for f in group_fits])
    se = np.array([float(np.std(f["logr_draw"], ddof=1)) for f in group_fits])
    se = np.maximum(se, 1e-4)

    # Grid in mu_b and tau_b. Prior: mu centered at log(0.5), broad; tau half-normal scale 0.5.
    mu_grid = np.linspace(-1.8, -0.05, 500)       # r ~ 0.17 to 0.95
    tau_grid = np.linspace(0.001, 1.0, 500)
    MU, TAU = np.meshgrid(mu_grid, tau_grid, indexing="ij")
    logp = np.zeros_like(MU)
    # Prior on mu: centered around geometric contraction r=0.5 but broad.
    logp += normal_logpdf_vec(MU, math.log(0.5), 0.70)
    # Half-normal prior on tau, plus constant irrelevant.
    logp += -0.5 * (TAU / 0.50) ** 2
    for mi, sei in zip(m, se):
        sd = np.sqrt(TAU * TAU + sei * sei)
        logp += -0.5 * (LOG2PI + 2.0 * np.log(sd) + ((mi - MU) / sd) ** 2)
    logp -= np.nanmax(logp)
    w = np.exp(logp)
    w /= np.sum(w)
    flat_w = w.ravel()
    idx = rng.choice(flat_w.size, size=draws, replace=True, p=flat_w)
    mu_draw = MU.ravel()[idx]
    tau_draw = TAU.ravel()[idx]

    # Conditional posterior draws for each latent group slope b_g.
    b_group = np.zeros((len(group_fits), draws), dtype=float)
    for gi, (mi, sei) in enumerate(zip(m, se)):
        prec = 1.0 / (sei * sei) + 1.0 / np.maximum(tau_draw, 1e-6) ** 2
        var = 1.0 / prec
        mean = var * (mi / (sei * sei) + mu_draw / np.maximum(tau_draw, 1e-6) ** 2)
        bg = rng.normal(mean, np.sqrt(var))
        # enforce convergence support; rare if priors/data are reasonable
        bg = np.minimum(bg, -1e-8)
        b_group[gi] = bg

    return {
        "mu_logr": mu_draw,
        "global_r": np.exp(mu_draw),
        "global_beta": -mu_draw,
        "tau_logr": tau_draw,
        "group_logr": b_group,
        "group_r": np.exp(b_group),
        "group_beta": -b_group,
        "m_logr": m,
        "se_logr": se,
    }


def level_reference_checks(groups: List[Group], draws: int, rng: np.random.Generator) -> List[Dict[str, float | str | int]]:
    rows = []
    for g in groups:
        if g.is_jump:
            continue
        # For level observables, do not infer contraction from 2 points. Check agreement with reference/asymptote.
        vals = []
        for y, sy in zip(g.value, g.sigma_obs):
            sy_eff = max(float(sy), 1e-12)
            vals.append(rng.normal(y, sy_eff, size=draws))
        vals = np.vstack(vals)
        last = vals[-1]
        mean_level = np.mean(vals, axis=0)
        ref = float(g.reference) if g.reference is not None else float("nan")
        delta_last = last - ref if np.isfinite(ref) else np.full(draws, np.nan)
        delta_mean = mean_level - ref if np.isfinite(ref) else np.full(draws, np.nan)
        rows.append({
            "group": g.name,
            "observable": g.observable,
            "time_fs": g.time_fs,
            "n_points": len(g.nc),
            "reference": ref,
            "last_nc": int(g.nc[-1]),
            "last_value": float(g.value[-1]),
            "delta_last_mean": float(np.mean(delta_last)),
            "delta_last_q025": q(delta_last, 0.025),
            "delta_last_q50": q(delta_last, 0.5),
            "delta_last_q975": q(delta_last, 0.975),
            "p_abs_delta_last_lt_1e_minus_4": float(np.mean(np.abs(delta_last) < 1e-4)),
            "mean_value": float(np.mean(g.value)),
            "delta_mean_mean": float(np.mean(delta_mean)),
            "delta_mean_q025": q(delta_mean, 0.025),
            "delta_mean_q50": q(delta_mean, 0.5),
            "delta_mean_q975": q(delta_mean, 0.975),
        })
    return rows


def write_outputs(
    out: Path,
    groups: List[Group],
    jump_fits: List[Dict[str, np.ndarray | float | int | str]],
    hier: Dict[str, np.ndarray],
    level_rows: List[Dict[str, float | str | int]],
    thresholds: List[float],
    draws: int,
    min_log_sigma: float,
    max_log_sigma: float,
):
    out.mkdir(parents=True, exist_ok=True)

    with (out / "group_loglinear_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "group", "observable", "time_fs", "n_points", "nc_min", "nc_max", "n0",
            "logr_mean", "logr_q025", "logr_q50", "logr_q975",
            "r_mean", "r_q025", "r_q50", "r_q975",
            "beta_mean", "beta_q025", "beta_q50", "beta_q975",
            "logA_mean", "logA_q50", "resid_sd_hat", "siglog_median", "siglog_max",
            "warning"
        ])
        for fit in jump_fits:
            b = np.asarray(fit["logr_draw"])
            r = np.asarray(fit["r_draw"])
            be = np.asarray(fit["beta_draw"])
            a = np.asarray(fit["logA_draw"])
            n_points = int(fit["n_points"])
            warning = "two_point_slope_only_no_internal_residual" if n_points < 3 else ""
            w.writerow([
                fit["name"], fit["observable"], fit["time_fs"], n_points, fit["nc_min"], fit["nc_max"], fit["n0"],
                float(np.mean(b)), q(b,0.025), q(b,0.5), q(b,0.975),
                float(np.mean(r)), q(r,0.025), q(r,0.5), q(r,0.975),
                float(np.mean(be)), q(be,0.025), q(be,0.5), q(be,0.975),
                float(np.mean(a)), q(a,0.5), fit["resid_sd_hat"], fit["siglog_median"], fit["siglog_max"],
                warning,
            ])

    with (out / "hierarchy_global_contraction.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "mean", "q025", "q50", "q975"])
        for name in ["mu_logr", "global_r", "global_beta", "tau_logr"]:
            arr = np.asarray(hier[name])
            w.writerow([name, float(np.mean(arr)), q(arr,0.025), q(arr,0.5), q(arr,0.975)])

    with (out / "hierarchy_group_shrinkage.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "group", "observable", "n_points", "raw_logr_mean", "raw_logr_se",
            "shrunk_logr_mean", "shrunk_logr_q025", "shrunk_logr_q50", "shrunk_logr_q975",
            "shrunk_r_mean", "shrunk_r_q025", "shrunk_r_q50", "shrunk_r_q975",
            "shrunk_beta_mean", "shrunk_beta_q025", "shrunk_beta_q50", "shrunk_beta_q975",
            "warning"
        ])
        for gi, fit in enumerate(jump_fits):
            bg = hier["group_logr"][gi]
            rr = hier["group_r"][gi]
            be = hier["group_beta"][gi]
            warning = "two_point_shrunk_by_global_hierarchy" if int(fit["n_points"]) < 3 else ""
            w.writerow([
                fit["name"], fit["observable"], fit["n_points"], hier["m_logr"][gi], hier["se_logr"][gi],
                float(np.mean(bg)), q(bg,0.025), q(bg,0.5), q(bg,0.975),
                float(np.mean(rr)), q(rr,0.025), q(rr,0.5), q(rr,0.975),
                float(np.mean(be)), q(be,0.025), q(be,0.5), q(be,0.975),
                warning,
            ])

    with (out / "extrapolated_jumps.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        header = ["group", "observable", "nc", "pred_mean", "pred_q025", "pred_q50", "pred_q975"]
        for t in thresholds:
            header.append(f"p_pred_lt_{t:g}")
        w.writerow(header)
        for fit in jump_fits:
            a = np.asarray(fit["logA_draw"])
            b = np.asarray(fit["logr_draw"])
            n0 = int(fit["n0"])
            start = int(fit["nc_min"])
            stop = int(fit["nc_max"]) + 3
            for nc in range(start, stop + 1):
                pred = np.exp(a + b * (nc - n0))
                row = [fit["name"], fit["observable"], nc, float(np.mean(pred)), q(pred,0.025), q(pred,0.5), q(pred,0.975)]
                row.extend([float(np.mean(pred < t)) for t in thresholds])
                w.writerow(row)

    if level_rows:
        keys = list(level_rows[0].keys())
        with (out / "level_reference_checks.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(level_rows)

    with (out / "diagnostics_v2.txt").open("w", encoding="utf-8") as f:
        f.write("Bayesian HEOM hierarchy v2 diagnostics\n")
        f.write("=====================================\n")
        f.write(f"jump_groups = {len(jump_fits)}\n")
        f.write(f"level_reference_groups = {len(level_rows)}\n")
        f.write(f"posterior_draws = {draws}\n")
        f.write(f"min_log_sigma = {min_log_sigma}\n")
        f.write(f"max_log_sigma = {max_log_sigma}\n")
        f.write("\nMethod:\n")
        f.write("- Jump observables are modeled on log scale as log(y_gn)=log(A_g)+(n-n0)log(r_g).\n")
        f.write("- theta=0 is imposed for jump magnitudes; no free asymptotic floor is estimated from small-N data.\n")
        f.write("- Group slopes log(r_g) are pooled by a Bayesian random-effects grid. No fragile MCMC is used.\n")
        f.write("- Population/asymptote level checks are not pooled with jump magnitudes.\n")
        f.write("\nWarnings:\n")
        for fit in jump_fits:
            if int(fit["n_points"]) < 3:
                f.write(f"- {fit['name']}: only {fit['n_points']} points; slope is ratio-calibration, not independent curve-shape evidence.\n")
        global_r = np.asarray(hier["global_r"])
        f.write("\nGlobal contraction ratio r = exp(mu_logr):\n")
        f.write(f"mean={float(np.mean(global_r)):.6g}, q025={q(global_r,0.025):.6g}, q50={q(global_r,0.5):.6g}, q975={q(global_r,0.975):.6g}\n")
        f.write("\nInterpretation:\n")
        f.write("This output is designed for small HEOM validation ledgers. It is not a substitute for new HEOM runs; it compresses existing convergence evidence into a stable posterior over contraction ratios.\n")

    if HAS_MPL:
        n = len(jump_fits)
        fig, axes = plt.subplots(nrows=n, ncols=1, figsize=(7.5, max(3.0, 2.8*n)), squeeze=False)
        for gi, fit in enumerate(jump_fits):
            ax = axes[gi, 0]
            g = next(gr for gr in groups if gr.name == fit["name"])
            x = g.nc
            ax.scatter(x, g.value, label="observed")
            a = np.asarray(fit["logA_draw"])
            b = np.asarray(fit["logr_draw"])
            grid = np.arange(int(np.min(g.nc)), int(np.max(g.nc))+4)
            preds = np.exp(a[:,None] + b[:,None]*(grid[None,:]-int(fit["n0"])))
            ax.plot(grid, np.median(preds, axis=0), label="median fit")
            ax.fill_between(grid, np.quantile(preds,0.025,axis=0), np.quantile(preds,0.975,axis=0), alpha=0.25)
            ax.set_yscale("log")
            ax.set_title(str(fit["name"]))
            ax.set_xlabel("NC index / jump target depth")
            ax.set_ylabel("jump magnitude")
            ax.legend(fontsize=8)
        plt.tight_layout()
        fig.savefig(out / "posterior_plots_v2.png", dpi=180)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Stable Bayesian hierarchy for HEOM contraction ledgers")
    ap.add_argument("input_csv")
    ap.add_argument("--output-dir", default="heom_bayes_v2_out")
    ap.add_argument("--draws", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260424)
    ap.add_argument("--min-log-sigma", type=float, default=0.03, help="floor for sd(log y), avoids overconfidence in pseudo-errors")
    ap.add_argument("--max-log-sigma", type=float, default=0.75, help="cap for sd(log y), avoids unstable near-zero weights")
    ap.add_argument("--thresholds", default="0.01,0.005,0.001,0.0005", help="comma-separated jump thresholds for extrapolation probabilities")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    groups = load_groups(args.input_csv)
    jump_groups = [g for g in groups if g.is_jump]
    if not jump_groups:
        raise SystemExit("No jump_* groups found. v2 needs positive contraction/jump observables.")
    jump_fits = [fit_jump_group_mc(g, args.draws, rng, args.min_log_sigma, args.max_log_sigma) for g in jump_groups]
    hier = hierarchy_grid(jump_fits, args.draws, rng)
    level_rows = level_reference_checks(groups, args.draws, rng)
    thresholds = [float(x) for x in args.thresholds.split(",") if x.strip()]
    out = Path(args.output_dir)
    write_outputs(out, groups, jump_fits, hier, level_rows, thresholds, args.draws, args.min_log_sigma, args.max_log_sigma)
    print(f"Completed HEOM Bayesian hierarchy v2: {len(jump_fits)} jump groups, {len(level_rows)} level-reference groups")
    print(f"Outputs written to: {out.resolve()}")


if __name__ == "__main__":
    main()
