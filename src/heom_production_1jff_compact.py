#!/usr/bin/env python3
"""
heom_production_1jff_compact.py

Compact production harness for 1JFF HEOM runs.

Goals:
1. Keep the scientifically relevant logic of the original production script.
2. Reduce I/O overhead by checkpointing only the resumable ADO state.
3. Save observables window-by-window instead of rewriting one giant trajectory file.
4. Make path resolution more robust on Windows and mixed working directories.

This script is intentionally conservative: it does not change the bath model,
Hamiltonian loading, or physical constants. It only changes orchestration and
persistence strategy.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import qutip as qt
from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver


CM_TO_RADFS = 2 * np.pi * 2.9979e-5
LAM_RADFS = 35.0 * CM_TO_RADFS
GAM_RADFS = 53.0 * CM_TO_RADFS
T_RADFS = 300.0 * 0.69503 * CM_TO_RADFS


def resolve_project_root(explicit_root: str | None = None) -> Path:
    """Resolve project root robustly across direct runs and copied scripts."""
    candidates: list[Path] = []

    if explicit_root:
        candidates.append(Path(explicit_root).expanduser().resolve())

    here = Path(__file__).resolve()
    candidates.extend(
        [
            here.parent,
            here.parent.parent,
            Path.cwd().resolve(),
            Path.cwd().resolve().parent,
        ]
    )

    seen = set()
    uniq: list[Path] = []
    for c in candidates:
        if c in seen:
            continue
        uniq.append(c)
        seen.add(c)

    for root in uniq:
        if (root / "H_1JFF.npz").exists():
            return root
        if (root / "git_repo").exists() and (root / "H_1JFF.npz").exists():
            return root

    # Last-resort heuristic: climb upwards from cwd looking for H_1JFF.npz.
    for base in [Path.cwd().resolve(), here.parent]:
        current = base
        for _ in range(5):
            if (current / "H_1JFF.npz").exists():
                return current
            current = current.parent

    raise FileNotFoundError(
        "Could not resolve project root containing H_1JFF.npz. "
        "Pass --project-root explicitly."
    )


def load_1jff(project_root: Path) -> tuple[qt.Qobj, list[str]]:
    npz_path = project_root / "H_1JFF.npz"
    if not npz_path.exists():
        raise FileNotFoundError(f"Missing Hamiltonian file: {npz_path}")
    data = np.load(npz_path, allow_pickle=True)
    H_cm = data["H_cm1"]
    labels = list(data["labels"])
    H_rad = (H_cm - np.mean(np.diag(H_cm)) * np.eye(H_cm.shape[0])) * CM_TO_RADFS
    return qt.Qobj(H_rad), labels


def build_solver(H_S: qt.Qobj, labels: list[str], nk: int, nc: int, store_ados: bool) -> HEOMSolver:
    N = len(labels)
    baths = []
    L_total = qt.liouvillian(H_S)
    for i in range(N):
        Q = qt.basis(N, i) * qt.basis(N, i).dag()
        bath = DrudeLorentzPadeBath(Q, lam=LAM_RADFS, gamma=GAM_RADFS, T=T_RADFS, Nk=nk)
        baths.append(bath)
        _, L_term = bath.terminator()
        L_total += L_term

    solver = HEOMSolver(
        L_total,
        baths,
        max_depth=nc,
        options={
            "nsteps": 100_000,
            "store_ados": store_ados,
            # Observables are what we care about in production stitching.
            "store_states": False,
            "progress_bar": False,
        },
    )
    return solver


def save_checkpoint(path: Path, payload: dict) -> None:
    with gzip.open(path, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_checkpoint(path: Path) -> dict:
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def window_file_name(index: int) -> str:
    return f"window_{index:03d}.npz"


def append_log(log_file: Path, msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compact 1JFF HEOM production harness")
    p.add_argument("--project-root", type=str, default=None)
    p.add_argument("--nc", type=int, default=7)
    p.add_argument("--nk", type=int, default=1)
    p.add_argument("--target-ps", type=float, default=30.0)
    p.add_argument("--window-ps", type=float, default=5.0)
    p.add_argument("--dt-fs", type=float, default=10.0)
    p.add_argument("--start-site-label", type=str, default="B:103")
    p.add_argument("--run-tag", type=str, default="nc7_nk1_compact")
    p.add_argument(
        "--store-ados",
        action="store_true",
        help="Store ADO trajectory in solver internals. Keep enabled if your QuTiP build requires it for final_ados_state.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    project_root = resolve_project_root(args.project_root)
    out_dir = project_root / f"heom_prod_{args.run_tag}"
    windows_dir = out_dir / "windows"
    windows_dir.mkdir(parents=True, exist_ok=True)

    ckpt_file = out_dir / "checkpoint.pkl.gz"
    summary_file = out_dir / "run_summary.json"
    log_file = out_dir / "production_progress.log"

    append_log(log_file, f"=== HEOM PRODUCTION START ({args.run_tag}) ===")
    append_log(log_file, f"project_root = {project_root}")
    append_log(log_file, f"out_dir = {out_dir}")
    append_log(log_file, f"config: NC={args.nc}, Nk={args.nk}, target_ps={args.target_ps}, window_ps={args.window_ps}, dt_fs={args.dt_fs}")

    H_S, labels = load_1jff(project_root)
    solver = build_solver(H_S, labels, nk=args.nk, nc=args.nc, store_ados=args.store_ados)
    N = len(labels)
    e_ops = [qt.basis(N, i) * qt.basis(N, i).dag() for i in range(N)]

    if ckpt_file.exists():
        ckpt = load_checkpoint(ckpt_file)
        current_state = ckpt["ados_state"]
        t_start_fs = float(ckpt["t_fs"])
        window_index = int(ckpt["window_index"])
        append_log(log_file, f"Resuming from checkpoint at t={t_start_fs:.1f} fs, window_index={window_index}")
    else:
        try:
            i0 = labels.index(args.start_site_label)
        except ValueError as exc:
            raise ValueError(f"Start-site label {args.start_site_label!r} not found in labels: {labels}") from exc
        rho0 = qt.basis(N, i0) * qt.basis(N, i0).dag()
        current_state = rho0
        t_start_fs = 0.0
        window_index = 0
        append_log(log_file, f"No checkpoint found. Starting from {args.start_site_label}")

    target_fs = args.target_ps * 1000.0
    window_fs = args.window_ps * 1000.0

    runtime_windows: list[float] = []
    finished = False

    while t_start_fs < target_fs - 1e-12:
        w_start = t_start_fs
        w_end = min(target_fs, w_start + window_fs)
        tlist = np.arange(w_start, w_end + args.dt_fs, args.dt_fs, dtype=float)
        if len(tlist) < 2:
            tlist = np.array([w_start, w_end], dtype=float)

        append_log(log_file, f"Window {window_index+1}: {w_start:.1f} -> {w_end:.1f} fs")
        t0 = time.time()
        result = solver.run(current_state, tlist, e_ops=e_ops)
        wall = time.time() - t0
        runtime_windows.append(wall)

        expect = np.asarray(result.expect, dtype=float)
        np.savez_compressed(
            windows_dir / window_file_name(window_index),
            tlist=tlist,
            expect=expect,
            labels=np.asarray(labels, dtype=object),
            nc=np.array(args.nc),
            nk=np.array(args.nk),
        )

        if not hasattr(result, "final_ados_state"):
            raise RuntimeError(
                "Result has no final_ados_state. Re-run with --store-ados if your QuTiP build requires it."
            )

        current_state = result.final_ados_state
        t_start_fs = w_end
        window_index += 1

        ms_per_fs = wall / max(w_end - w_start, 1e-12) * 1000.0
        append_log(log_file, f"Window done in {wall/60.0:.2f} min | speed={ms_per_fs:.1f} ms/fs")

        save_checkpoint(
            ckpt_file,
            {
                "t_fs": t_start_fs,
                "window_index": window_index,
                "ados_state": current_state,
                "labels": labels,
                "nc": args.nc,
                "nk": args.nk,
            },
        )
        append_log(log_file, f"Checkpoint saved at t={t_start_fs:.1f} fs")

    finished = True

    avg_wall = float(np.mean(runtime_windows)) if runtime_windows else 0.0
    avg_ms_per_fs = avg_wall / (window_fs if window_fs > 0 else 1.0) * 1000.0
    summary = {
        "finished": finished,
        "project_root": str(project_root),
        "run_tag": args.run_tag,
        "nc": args.nc,
        "nk": args.nk,
        "target_ps": args.target_ps,
        "window_ps": args.window_ps,
        "dt_fs": args.dt_fs,
        "num_windows": window_index,
        "avg_window_wall_s": avg_wall,
        "avg_ms_per_fs": avg_ms_per_fs,
        "output_dir": str(out_dir),
        "windows_dir": str(windows_dir),
        "checkpoint": str(ckpt_file),
    }
    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    append_log(log_file, "=== PRODUCTION FINISHED ===")
    append_log(log_file, f"Summary saved to {summary_file}")


if __name__ == "__main__":
    main()
