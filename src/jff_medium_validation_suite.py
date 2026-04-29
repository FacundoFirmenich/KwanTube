#!/usr/bin/env python3
"""
jff_medium_validation_suite.py
NC-SPECIFIC VALIDATION SUITE for the full 1JFF HEOM model.
"""

from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import qutip as qt
from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[1] # retrocede desde src/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

# ------------------------- physical parameters -------------------------
LAM_CM = 35.0
GAM_CM = 53.0
T_K = 300.0
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
T_RADFS = T_K * 0.69503 * CM_TO_RADFS
DEFAULT_THRESHOLD = 1e-3

def find_project_root(explicit_root: str | None) -> Path:
    if explicit_root:
        return Path(explicit_root).resolve()
    # Prioridad: Estructura Tier-0
    if (PROJECT_ROOT / "outputs_data" / "raw_npz" / "H_1JFF.npz").exists():
        return PROJECT_ROOT
    if (PROJECT_ROOT / "H_1JFF.npz").exists():
        return PROJECT_ROOT
    raise FileNotFoundError("Could not locate H_1JFF.npz. Pass --project-root explicitly.")

def load_hamiltonian(project_root: Path) -> Tuple[qt.Qobj, str, List[str]]:
    npz_path = project_root / "outputs_data" / "raw_npz" / "H_1JFF.npz"
    if not npz_path.exists():
        npz_path = project_root / "H_1JFF.npz"
    if not npz_path.exists():
        raise FileNotFoundError(f"Missing Hamiltonian file: {npz_path}")
    data = np.load(npz_path, allow_pickle=True)
    H_cm = data["H_cm1"]
    labels = list(data.get("labels", []))
    H = (H_cm - np.mean(np.diag(H_cm)) * np.eye(H_cm.shape[0])) * CM_TO_RADFS
    return qt.Qobj(H), str(npz_path), labels

def site_projectors(n_sites: int) -> List[qt.Qobj]:
    return [qt.basis(n_sites, n) * qt.basis(n_sites, n).dag() for n in range(n_sites)]

def initial_state(n_sites: int, site: int) -> qt.Qobj:
    psi = qt.basis(n_sites, site)
    return psi * psi.dag()

def run_heom(H_S: qt.Qobj, coupling_ops: List[qt.Qobj], NC: int, Nk: int, tlist: np.ndarray, rho0: qt.Qobj, label: str) -> Tuple[List[qt.Qobj], float]:
    lam_rad = LAM_CM * CM_TO_RADFS
    gam_rad = GAM_CM * CM_TO_RADFS
    baths = [DrudeLorentzPadeBath(Q=Q, lam=lam_rad, gamma=gam_rad, T=T_RADFS, Nk=Nk) for Q in coupling_ops]
    solver = HEOMSolver(H_S, baths, max_depth=NC, options={"nsteps": 100_000, "store_states": True, "progress_bar": False})
    print(f"[run {label}] NC={NC}, Nk={Nk}, tmax={tlist[-1]:.1f} fs", flush=True)
    t0 = time.time()
    result = solver.run(rho0, tlist)
    wall = time.time() - t0
    print(f"           finished in {wall:.1f} s", flush=True)
    return result.states, wall

def state_metrics(tlist: np.ndarray, states_a: List[qt.Qobj], states_b: List[qt.Qobj]) -> Dict[str, any]:
    dpop_matrix = []
    dcoh = []
    dfrob = []
    for a, b in zip(states_a, states_b):
        A, B = a.full(), b.full()
        D = A - B
        dpop_matrix.append(np.abs(np.real(np.diag(D))))
        off = D - np.diag(np.diag(D))
        dcoh.append(float(np.max(np.abs(off))))
        dfrob.append(float(np.linalg.norm(D, ord="fro")))
    
    dpop_matrix = np.array(dpop_matrix)
    max_pop_idx = np.unravel_index(np.argmax(dpop_matrix), dpop_matrix.shape)
    
    return {
        "dPop": np.max(dpop_matrix, axis=1),
        "dCoh": np.array(dcoh),
        "dFrob": np.array(dfrob),
        "report": {
            "t_dPop_max": float(tlist[max_pop_idx[0]]),
            "site_dPop_max": int(max_pop_idx[1]),
            "t_dCoh_max": float(tlist[np.argmax(dcoh)]),
            "t_dFrob_max": float(tlist[np.argmax(dfrob)]),
            "final_dPop": float(dpop_matrix[-1].max()),
            "final_dCoh": float(dcoh[-1]),
            "final_dFrob": float(dfrob[-1])
        }
    }

def window_summary(tlist: np.ndarray, metrics: Dict[str, any], windows_fs: List[int]) -> Dict[str, Dict[str, float]]:
    out = {}
    for w in windows_fs:
        mask = tlist <= float(w) + 1e-12
        out[f"0-{int(w)}fs"] = {
            "dPop_max": float(np.max(metrics["dPop"][mask])),
            "dCoh_max": float(np.max(metrics["dCoh"][mask])),
            "dFrob_max": float(np.max(metrics["dFrob"][mask])),
        }
    return out

def build_markdown_summary(payload: Dict) -> str:
    cfg = payload["config"]
    lines = ["# 1JFF NC Validation Summary\n", "## Configuration"]
    lines.append(f"- h_source: {cfg['h_source']}")
    lines.append(f"- Level: NC={cfg['NC']}, Nk={cfg['Nk']}")
    lines.append(f"- Init site: {cfg['init_site_index']} ({cfg['init_site_label']})")
    lines.append(f"- Horizon: {cfg['tmax_fs']} fs\n")
    lines.append("## Runtimes")
    for name, r in payload["runs"].items():
        lines.append(f"- {name}: NC={r['NC']}, Nk={r['Nk']}, wall={r['wall_s']:.1f} s")
    nc = payload["comparisons"]["NC_stress"]
    lines.append("\n## NC-stress summary")
    lines.append(f"- Max diffs: dPop={nc['full']['dPop_max']:.2e}, dCoh={nc['full']['dCoh_max']:.2e}, dFrob={nc['full']['dFrob_max']:.2e}")
    lines.append(f"- t_dPop_max: {nc['report']['t_dPop_max']} fs (site {nc['report']['site_dPop_max']})")
    lines.append(f"- Convergence ratio r = {payload['comparisons']['ratio_r']:.3f}\n")
    lines.append("| Window | dPop max | dCoh max | dFrob max |")
    lines.append("|---|---:|---:|---:|")
    for window, vals in nc["windows"].items():
        lines.append(f"| {window} | {vals['dPop_max']:.2e} | {vals['dCoh_max']:.2e} | {vals['dFrob_max']:.2e} |")
    return "\n".join(lines) + "\n"

def main():
    parser = argparse.ArgumentParser(description="NC-specific 1JFF HEOM validation suite")
    parser.add_argument("--project-root", type=str, default=None)
    parser.add_argument("--nc", type=int, default=5, help="Working NC level")
    parser.add_argument("--nk", type=int, default=1, help="Fixed Nk level")
    parser.add_argument("--site", type=int, default=5, help="Initial site (default 5 for B:103)")
    parser.add_argument("--tmax", type=float, default=200.0)
    parser.add_argument("--sample", type=float, default=5.0)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--out-json", type=str, default=str(PROJECT_ROOT / "outputs_data" / "raw_json" / "jff_medium_validation_suite.json"))
    parser.add_argument("--out-md", type=str, default=str(PROJECT_ROOT / "outputs_data" / "raw_txt+md" / "jff_medium_validation_suite.md"))
    args = parser.parse_args()

    project_root = find_project_root(args.project_root)
    H_S, h_source, labels = load_hamiltonian(project_root)
    N = H_S.shape[0]
    init_label = labels[args.site] if args.site < len(labels) else "N/A"
    tlist = np.linspace(0.0, args.tmax, int(args.tmax/args.sample) + 1)
    rho0 = initial_state(N, args.site)
    S_ops = site_projectors(N)

    runs = {}
    states_lo, wall_lo = run_heom(H_S, S_ops, args.nc - 1, args.nk, tlist, rho0, "NC-1")
    runs["NC_minus_1"] = {"NC": args.nc - 1, "Nk": args.nk, "wall_s": wall_lo}
    states_mid, wall_mid = run_heom(H_S, S_ops, args.nc, args.nk, tlist, rho0, "NC")
    runs["baseline"] = {"NC": args.nc, "Nk": args.nk, "wall_s": wall_mid}
    states_hi, wall_hi = run_heom(H_S, S_ops, args.nc + 1, args.nk, tlist, rho0, "NC+1")
    runs["NC_stress"] = {"NC": args.nc + 1, "Nk": args.nk, "wall_s": wall_hi}

    m_lo_mid = state_metrics(tlist, states_lo, states_mid)
    m_mid_hi = state_metrics(tlist, states_mid, states_hi)
    d_lo_mid = float(np.max(m_lo_mid["dFrob"]))
    d_mid_hi = float(np.max(m_mid_hi["dFrob"]))
    ratio_r = d_mid_hi / d_lo_mid if d_lo_mid > 0 else float("nan")

    payload = {
        "config": {"h_source": h_source, "NC": args.nc, "Nk": args.nk, "init_site_index": args.site, "init_site_label": init_label, "tmax_fs": args.tmax},
        "runs": runs,
        "comparisons": {
            "ratio_r": ratio_r,
            "NC_stress": {
                "full": {"dPop_max": float(np.max(m_mid_hi["dPop"])), "dCoh_max": float(np.max(m_mid_hi["dCoh"])), "dFrob_max": d_mid_hi},
                "report": m_mid_hi["report"],
                "windows": window_summary(tlist, m_mid_hi, [40, 80, 120, 160, 200])
            }
        }
    }

    p_json = Path(args.out_json)
    p_json.parent.mkdir(parents=True, exist_ok=True)
    p_json.write_text(json.dumps(payload, indent=2))
    
    p_md = Path(args.out_md)
    p_md.parent.mkdir(parents=True, exist_ok=True)
    p_md.write_text(build_markdown_summary(payload))
    
    p_txt = p_md.with_suffix(".txt")
    p_txt.write_text(build_markdown_summary(payload))

    print(f"\nNC Validation Complete. Ratio r = {ratio_r:.3f}")
    print(f"Results saved in outputs_data/")

if __name__ == "__main__":
    main()
