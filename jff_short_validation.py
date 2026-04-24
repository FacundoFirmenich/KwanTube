"""
jff_short_validation.py

Valida que el par (NC*, Nk*) identificado en 6DPU sigue siendo
suficiente para el sistema completo 1JFF (8 sitios).
"""

import numpy as np
import qutip as qt
import sys
from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver
import pickle, json, time
from pathlib import Path

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parent.parent # La raíz es biofisicaquantiqaCLINE
if str(PROJECT_ROOT / "git_repo" / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "git_repo" / "src"))

# ----------------- parametros fisicos -----------------
LAM_CM   = 35.0
GAM_CM   = 53.0
T_K      = 300.0
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
T_RADFS  = T_K * 0.69503 * CM_TO_RADFS

T_MAX_VAL_FS = 200.0
DT_FS        = 1.0
THRESHOLD    = 1e-3

PADE_REPORT = PROJECT_ROOT / "git_repo" / "pade_convergence_report.json"

def H_1JFF_from_file():
    npz_path = PROJECT_ROOT / "H_1JFF.npz"
    if not npz_path.exists():
        print(f"ERROR: No se encuentra {npz_path}")
        sys.exit(1)
    d = np.load(npz_path, allow_pickle=True)
    H_cm = d["H_cm1"]
    H = (H_cm - np.mean(np.diag(H_cm)) * np.eye(H_cm.shape[0])) * CM_TO_RADFS
    return qt.Qobj(H)

def site_projectors(N):
    return [qt.basis(N, n) * qt.basis(N, n).dag() for n in range(N)]

def initial_state(N, site=0):
    psi = qt.basis(N, site)
    return psi * psi.dag()

def run_heom(H_S, coupling_ops, NC, Nk, t_max_fs, dt_fs, rho0, label=""):
    # Convertir lambdas a rad/fs para consistencia con H
    lam_rad = LAM_CM * CM_TO_RADFS
    gam_rad = GAM_CM * CM_TO_RADFS
    baths = [
        DrudeLorentzPadeBath(Q=Q, lam=lam_rad, gamma=gam_rad, T=T_RADFS, Nk=Nk)
        for Q in coupling_ops
    ]
    solver = HEOMSolver(H_S, baths, max_depth=NC,
                        options={"nsteps": 100_000,
                                 "store_states": True,
                                 "progress_bar": False})
    tlist = np.arange(0.0, t_max_fs + dt_fs, dt_fs)
    print(f"[run {label}] NC={NC}, Nk={Nk}, t_max={t_max_fs} fs ...", flush=True)
    t0 = time.time()
    result = solver.run(rho0, tlist)
    wall = time.time() - t0
    print(f"          wall = {wall:.1f} s")
    return tlist, result.states, wall

def max_frobenius_diff(states_A, states_B):
    diffs = [np.linalg.norm((a - b).full(), ord='fro')
             for a, b in zip(states_A, states_B)]
    return float(np.max(diffs)), int(np.argmax(diffs))

def main():
    if not PADE_REPORT.exists():
        print("ERROR: correr primero heom_pade_convergence.py")
        return
    with open(PADE_REPORT) as f:
        rep = json.load(f)
    v = rep.get("verdict")
    if not isinstance(v, dict):
        print(f"ERROR: verdict no es dict: {v}")
        return
    NC_star, Nk_star = int(v["NC"]), int(v["Nk"])
    
    H_S   = H_1JFF_from_file()
    N     = H_S.shape[0]
    S_ops = site_projectors(N)
    rho0  = initial_state(N, site=0)
    
    # 1. baseline
    t_b, st_b, w_b = run_heom(H_S, S_ops, NC_star, Nk_star, T_MAX_VAL_FS, DT_FS, rho0, "baseline")
    # 2. NC-stress
    t_n, st_n, w_n = run_heom(H_S, S_ops, NC_star+1, Nk_star, T_MAX_VAL_FS, DT_FS, rho0, "NC-stress")
    # 3. Nk-stress
    t_k, st_k, w_k = run_heom(H_S, S_ops, NC_star, Nk_star+1, T_MAX_VAL_FS, DT_FS, rho0, "Nk-stress")
    
    d_nc, ti_nc = max_frobenius_diff(st_b, st_n)
    d_nk, ti_nk = max_frobenius_diff(st_b, st_k)
    
    print("\n" + "="*60)
    print("  VALIDACION 1JFF (ventana corta)")
    print(f"  NC {NC_star}->{NC_star+1} dmax={d_nc:.2e} {'OK' if d_nc<THRESHOLD else 'FAIL'}")
    print(f"  Nk {Nk_star}->{Nk_star+1} dmax={d_nk:.2e} {'OK' if d_nk<THRESHOLD else 'FAIL'}")
    
    est_prod = w_b * (30000.0 / T_MAX_VAL_FS)
    print(f"\n  [cost] est production 30ps: {est_prod/60:.1f} min")

if __name__ == "__main__":
    main()
