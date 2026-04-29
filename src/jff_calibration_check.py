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

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[1] # retrocede desde src/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

# --- Constants ---
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
LAM_RADFS = 35.0 * CM_TO_RADFS
GAM_RADFS = 53.0 * CM_TO_RADFS
T_RADFS   = 300.0 * 0.69503 * CM_TO_RADFS

# --- Load 1JFF ---
def load_1jff():
    # Sincronizado con outputs_data/raw_npz/
    npz_path = PROJECT_ROOT / "outputs_data" / "raw_npz" / "H_1JFF.npz"
    if not npz_path.exists():
        # Fallback para la raiz
        npz_path = PROJECT_ROOT / "H_1JFF.npz"
        
    if not npz_path.exists():
        print(f"ERROR: No se encuentra H_1JFF.npz en {npz_path}")
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
    
    # Start at B:103
    try:
        i0 = labels.index("B:103")
    except ValueError:
        i0 = 0 # Fallback si los labels no coinciden
        
    rho0 = qt.basis(N, i0) * qt.basis(N, i0).dag()
    
    for i in range(N):
        Q = qt.basis(N, i) * qt.basis(N, i).dag()
        bath = DrudeLorentzPadeBath(Q, lam=LAM_RADFS, gamma=GAM_RADFS, T=T_RADFS, Nk=Nk)
        baths.append(bath)
        
    solver = HEOMSolver(H_S, baths, max_depth=NC, 
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
    print("\n[3] Calculating ratio r(5->6) at 500fs...")
    P5 = get_pop(states5)
    P6 = get_pop(states6[:501])
    d56 = np.max(np.abs(P6 - P5))
    
    # NC=4 for baseline
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
    ado_ratio = 245157 / 74613 # Ratio aproximado NC7/NC6
    time_per_fs_7 = (wall6/500) * ado_ratio
    total_time_30ps = time_per_fs_7 * 30000
    print(f"    Est. wall-time NC=7 (30ps): {total_time_30ps/3600:.1f} hours")
    
    # Truncation error projection
    eps7 = (d56 * ratio) / (1 - ratio) if ratio < 1 else 999
    
    # Save data
    out_path = PROJECT_ROOT / "outputs_data" / "raw_npz" / "jff_calib_data.npz"
    np.savez(out_path, 
             tlist=tlist5, P5=P5, P6=P6, 
             ratio=ratio, d56=d56, eps7=eps7)
    print(f"\n[5] Data saved to {out_path}")

if __name__ == "__main__":
    main()
