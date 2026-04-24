"""
jff_calibration_check.py -- Calibration tests on 1JFF (8 sites)

1. Smoke test: NC=5, 500fs.
2. Conv test: NC=6, 1ps.
3. Compare ratio r(5->6) with fragment baseline (0.42-0.45).
"""

import numpy as np
import qutip as qt
import sys
from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver
import time, pickle, json
from pathlib import Path

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parent.parent # La raíz es biofisicaquantiqaCLINE
if str(PROJECT_ROOT / "git_repo" / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "git_repo" / "src"))

# --- Constants ---
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
LAM_RADFS = 35.0 * CM_TO_RADFS
GAM_RADFS = 53.0 * CM_TO_RADFS
T_RADFS   = 300.0 * 0.69503 * CM_TO_RADFS

# --- Load 1JFF ---
def load_1jff():
    npz_path = PROJECT_ROOT / "H_1JFF.npz"
    if not npz_path.exists():
        print(f"ERROR: No se encuentra {npz_path}")
        sys.exit(1)
    d = np.load(npz_path, allow_pickle=True)
    H_cm = d["H_cm1"]
    labels = list(d["labels"])
    # Center H
    H_rad = (H_cm - np.mean(np.diag(H_cm)) * np.eye(H_cm.shape[0])) * CM_TO_RADFS
    return qt.Qobj(H_rad), labels

def run_test(H_S, labels, NC, Nk, t_max_fs, dt_fs=1.0):
    N = len(labels)
    baths = []
    L_total = qt.liouvillian(H_S)
    
    # Start at B:103
    i0 = labels.index("B:103")
    rho0 = qt.basis(N, i0) * qt.basis(N, i0).dag()
    
    for i in range(N):
        Q = qt.basis(N, i) * qt.basis(N, i).dag()
        bath = DrudeLorentzPadeBath(Q, lam=LAM_RADFS, gamma=GAM_RADFS, T=T_RADFS, Nk=Nk)
        baths.append(bath)
        _, L_term = bath.terminator()
        L_total += L_term
        
    solver = HEOMSolver(L_total, baths, max_depth=NC, 
                        options={"nsteps": 100_000, "store_states": True})
    
    tlist = np.arange(0, t_max_fs + dt_fs, dt_fs)
    print(f"  [setup] NC={NC}, Nk={Nk}, ADOs={len(solver.ados.labels)}")
    
    t0 = time.time()
    result = solver.run(rho0, tlist)
    wall = time.time() - t0
    return tlist, result.states, wall

def get_pop(states):
    N = states[0].shape[0]
    return np.array([[s.full()[i,i].real for s in states] for i in range(N)])

def main():
    print("="*60)
    print("  1JFF CALIBRATION CHECK (Smoke + Convergence)")
    print("="*60)
    
    H_S, labels = load_1jff()
    
    # 1. Smoke test (NC=5, 500fs)
    print("\n[1] Smoke Test: NC=5, Nk=1, 500fs...")
    tlist5, states5, wall5 = run_test(H_S, labels, NC=5, Nk=1, t_max_fs=500.0)
    print(f"    Wall time: {wall5:.1f}s  ({wall5/500*1000:.1f} ms/fs)")
    
    # 2. Convergence test (NC=6, 500fs)
    print("\n[2] Conv Test: NC=6, Nk=1, 500fs...")
    tlist6, states6, wall6 = run_test(H_S, labels, NC=6, Nk=1, t_max_fs=500.0)
    print(f"    Wall time: {wall6:.1f}s  ({wall6/500*1000:.1f} ms/fs)")
    
    # 3. Ratio Calculation
    # Need NC=5 for 1ps too
    print("\n[3] Calculating ratio r(5->6) at 500fs...")
    # populations at 500fs
    idx500 = 500 # assuming dt=1.0
    P5 = get_pop(states5)
    P6 = get_pop(states6[:501])
    d56 = np.max(np.abs(P6 - P5))
    
    # We need NC=4 for 500fs to get the previous jump d45
    print("    Running NC=4 for d45 baseline...")
    _, states4, _ = run_test(H_S, labels, NC=4, Nk=1, t_max_fs=500.0)
    P4 = get_pop(states4)
    d45 = np.max(np.abs(P5 - P4))
    
    ratio = d56 / d45 if d45 > 0 else 0
    print(f"\n    Jump d(4->5): {d45:.2e}")
    print(f"    Jump d(5->6): {d56:.2e}")
    print(f"    Ratio r(1JFF): {ratio:.3f}  (Fragment was ~0.45)")
    
    # 4. Projections
    print("\n[4] Projections for 1JFF (t=30ps):")
    # NC=7 cost estimate
    # NC=7, L=16 -> 245,157 ADOs
    # ratio of ADOs NC7/NC6 = 245157 / 74613 (approx)
    # L=16, NC=6 -> (16+6)!/(16!6!) = 22!/(16!6!) = 74613
    ado_ratio = 245157 / 74613
    time_per_fs_7 = (wall6/1000) * ado_ratio
    total_time_30ps = time_per_fs_7 * 30000
    print(f"    Est. wall-time NC=7 (30ps): {total_time_30ps/3600:.1f} hours")
    
    # Truncation error projection
    eps7 = (d56 * ratio) / (1 - ratio) if ratio < 1 else 999
    
    # Save data for analysis
    out_path = PROJECT_ROOT / "git_repo" / "jff_calib_data.npz"
    np.savez(out_path, 
             tlist=tlist5, P5=P5, P6=P6, 
             ratio=ratio, d56=d56, eps7=eps7)
    print(f"\n[5] Data saved to {out_path.name}")

if __name__ == "__main__":
    main()
