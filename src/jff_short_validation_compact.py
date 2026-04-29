"""
jff_short_validation_compact.py
Validacion compacta de (NC*, Nk*) para 1JFF.
"""

from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
import qutip as qt
from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver

# ----------------- parametros fisicos -----------------
LAM_CM = 35.0
GAM_CM = 53.0
T_K = 300.0
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
T_RADFS = T_K * 0.69503 * CM_TO_RADFS

DEFAULT_TMAX_FS = 120.0
DEFAULT_SAMPLE_FS = 5.0
DEFAULT_THRESHOLD = 1e-3

def find_project_root() -> Path:
    """Encuentra la raiz del proyecto de forma robusta."""
    here = Path(__file__).resolve()
    # La raiz es biofisicaquantiqaCLINE (2 niveles arriba de src/)
    root = here.parents[1]
    if (root / "outputs_data").exists():
        return root
    return root

PROJECT_ROOT = find_project_root()

# Import path robusto
src_path = PROJECT_ROOT / "src"
if src_path.exists() and str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

# Sincronizado con la arquitectura Tier-0
PADE_REPORT_CANDIDATES = [
    PROJECT_ROOT / "outputs_data" / "raw_json" / "pade_convergence_report.json",
    PROJECT_ROOT / "pade_convergence_report.json",
]

def load_best_pair(user_nc: int | None, user_nk: int | None) -> tuple[int, int, str]:
    if user_nc is not None and user_nk is not None:
        return user_nc, user_nk, "CLI"
    for p in PADE_REPORT_CANDIDATES:
        if p.exists():
            with open(p, "r", encoding="utf-8") as f:
                rep = json.load(f)
            v = rep.get("verdict", {})
            if isinstance(v, dict) and "NC" in v and "Nk" in v:
                return int(v["NC"]), int(v["Nk"]), str(p)
    return 5, 1, "fallback_default"

def H_1JFF_from_file() -> qt.Qobj:
    npz_path = PROJECT_ROOT / "outputs_data" / "raw_npz" / "H_1JFF.npz"
    if not npz_path.exists():
        npz_path = PROJECT_ROOT / "H_1JFF.npz"
        
    if not npz_path.exists():
        raise FileNotFoundError(f"No se encontro H_1JFF.npz en {npz_path}")
        
    d = np.load(npz_path, allow_pickle=True)
    H_cm = d["H_cm1"]
    H = (H_cm - np.mean(np.diag(H_cm)) * np.eye(H_cm.shape[0])) * CM_TO_RADFS
    return qt.Qobj(H)

def site_projectors(n: int) -> list[qt.Qobj]:
    return [qt.basis(n, i) * qt.basis(n, i).dag() for i in range(n)]

def coherence_op(n: int, i: int = 0, j: int = 1) -> qt.Qobj:
    return qt.basis(n, i) * qt.basis(n, j).dag()

def initial_state(n: int, site: int = 0) -> qt.Qobj:
    psi = qt.basis(n, site)
    return psi * psi.dag()

def run_heom(H_S, coupling_ops, NC, Nk, t_max_fs, sample_fs, rho0, label="") -> tuple[np.ndarray, list[qt.Qobj], float]:
    lam_rad = LAM_CM * CM_TO_RADFS
    gam_rad = GAM_CM * CM_TO_RADFS
    baths = [DrudeLorentzPadeBath(Q=Q, lam=lam_rad, gamma=gam_rad, T=T_RADFS, Nk=Nk) for Q in coupling_ops]
    solver = HEOMSolver(H_S, baths, max_depth=NC, options={"nsteps": 100_000, "store_states": True, "progress_bar": False})
    tlist = np.arange(0.0, t_max_fs + sample_fs, sample_fs)
    print(f"[run {label}] NC={NC}, Nk={Nk}, t_max={t_max_fs} fs, sample={sample_fs} fs ...", flush=True)
    t0 = time.time()
    result = solver.run(rho0, tlist)
    wall = time.time() - t0
    print(f"          wall = {wall:.1f} s, n_times={len(tlist)}")
    return tlist, result.states, wall

def max_frobenius_diff(states_A, states_B):
    diffs = [np.linalg.norm((a - b).full(), ord="fro") for a, b in zip(states_A, states_B)]
    return float(np.max(diffs)), int(np.argmax(diffs))

def max_population_diff(states_A, states_B, proj_ops):
    max_diff = -1.0
    argmax_t = 0
    for ti, (a, b) in enumerate(zip(states_A, states_B)):
        pa = np.array([(a * P).tr().real for P in proj_ops])
        pb = np.array([(b * P).tr().real for P in proj_ops])
        d = np.max(np.abs(pa - pb))
        if d > max_diff:
            max_diff = float(d)
            argmax_t = ti
    return max_diff, argmax_t

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nc", type=int, default=None)
    ap.add_argument("--nk", type=int, default=None)
    ap.add_argument("--tmax", type=float, default=DEFAULT_TMAX_FS)
    ap.add_argument("--sample", type=float, default=DEFAULT_SAMPLE_FS)
    ap.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    ap.add_argument("--with-nk-stress", action="store_true")
    ap.add_argument("--site", type=int, default=0)
    args = ap.parse_args()

    NC_star, Nk_star, source = load_best_pair(args.nc, args.nk)
    print("=" * 72)
    print("  VALIDACION COMPACTA 1JFF")
    print("=" * 72)
    print(f"PROJECT_ROOT: {PROJECT_ROOT}")
    print(f"Par base (NC*, Nk*): ({NC_star}, {Nk_star})   [fuente={source}]")

    H_S = H_1JFF_from_file()
    n = H_S.shape[0]
    S_ops = site_projectors(n)
    coh = coherence_op(n, 0, 1)
    rho0 = initial_state(n, site=args.site)

    t_b, st_b, w_b = run_heom(H_S, S_ops, NC_star, Nk_star, args.tmax, args.sample, rho0, "baseline")
    t_n, st_n, w_n = run_heom(H_S, S_ops, NC_star + 1, Nk_star, args.tmax, args.sample, rho0, "NC-stress")
    
    dF_nc, tiF_nc = max_frobenius_diff(st_b, st_n)
    dP_nc, tiP_nc = max_population_diff(st_b, st_n, S_ops)
    
    print("\n" + "-" * 72)
    print("NC-STRESS")
    print(f"  max dFrob = {dF_nc:.2e} @ t={t_b[tiF_nc]:.1f} fs   {'OK' if dF_nc < args.threshold else 'FAIL'}")
    print(f"  max dPop  = {dP_nc:.2e} @ t={t_b[tiP_nc]:.1f} fs")

    result = {
        "baseline": {"NC": NC_star, "Nk": Nk_star, "tmax_fs": args.tmax, "sample_fs": args.sample, "wall_s": w_b},
        "nc_stress": {"NC": NC_star + 1, "Nk": Nk_star, "wall_s": w_n, "dFrob_max": dF_nc, "dPop_max": dP_nc},
    }

    if args.with_nk_stress:
        t_k, st_k, w_k = run_heom(H_S, S_ops, NC_star, Nk_star + 1, args.tmax, args.sample, rho0, "Nk-stress")
        dF_nk, _ = max_frobenius_diff(st_b, st_k)
        result["nk_stress"] = {"NC": NC_star, "Nk": Nk_star + 1, "wall_s": w_k, "dFrob_max": dF_nk}

    out_path = PROJECT_ROOT / "outputs_data" / "raw_json" / "jff_short_validation_compact_report.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"\nReporte guardado en: {out_path}")

if __name__ == "__main__":
    main()
