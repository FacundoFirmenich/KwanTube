"""
heom_production_1jff.py -- Production run with Windowed Checkpointing

Simulates 1JFF (8 sites) at NC=7, Nk=1 (Pade).
Saves full ADO stack every 2 ps to allow resumption.
"""

import numpy as np
import qutip as qt
import sys
from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver
import time, pickle, os
from pathlib import Path

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[1] # retrocede desde src/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

# --- Configuration ---
NC = 7
NK = 1
T_TARGET_PS = 10.0
WINDOW_PS   = 2.0
DT_FS       = 5.0
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
LAM_RADFS = 35.0 * CM_TO_RADFS
GAM_RADFS = 53.0 * CM_TO_RADFS
T_RADFS   = 300.0 * 0.69503 * CM_TO_RADFS

# Files - sincronizados con la nueva arquitectura
CKPT_FILE = PROJECT_ROOT / "outputs_data" / "raw_pkl" / "heom_1jff_prod_checkpoint.pkl"
TRAJ_FILE = PROJECT_ROOT / "outputs_data" / "raw_pkl" / "heom_1jff_prod_trajectory.pkl"
LOG_FILE  = PROJECT_ROOT / "outputs_data" / "raw_txt+md" / "production_progress.log"

def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}\n"
    print(line, end="")
    # Aseguramos que el directorio del log exista
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG_FILE, "a") as f:
        f.write(line)

def load_1jff():
    npz_path = PROJECT_ROOT / "outputs_data" / "raw_npz" / "H_1JFF.npz"
    if not npz_path.exists():
        npz_path = PROJECT_ROOT / "H_1JFF.npz"
        
    if not npz_path.exists():
        log(f"ERROR: No se encuentra H_1JFF.npz en {npz_path}")
        sys.exit(1)
        
    d = np.load(npz_path, allow_pickle=True)
    H_cm = d["H_cm1"]
    labels = list(d["labels"])
    H_rad = (H_cm - np.mean(np.diag(H_cm)) * np.eye(H_cm.shape[0])) * CM_TO_RADFS
    return qt.Qobj(H_rad), labels

def setup_solver(H_S, labels):
    N = len(labels)
    baths = []
    for i in range(N):
        Q = qt.basis(N, i) * qt.basis(N, i).dag()
        bath = DrudeLorentzPadeBath(Q, lam=LAM_RADFS, gamma=GAM_RADFS, T=T_RADFS, Nk=NK)
        baths.append(bath)
    solver = HEOMSolver(H_S, baths, max_depth=NC, 
                        options={"nsteps": 100_000, "store_ados": True})
    return solver

def main():
    log("=== HEOM PRODUCTION START (1JFF) ===")
    H_S, labels = load_1jff()
    solver = setup_solver(H_S, labels)
    N = len(labels)
    e_ops = [qt.basis(N, i) * qt.basis(N, i).dag() for i in range(N)]
    
    # Check for existing checkpoint
    if os.path.exists(CKPT_FILE):
        log(f"Resuming from checkpoint: {CKPT_FILE}")
        with open(CKPT_FILE, "rb") as f:
            ckpt = pickle.load(f)
        current_state = ckpt["ados_state"]
        t_start = ckpt["t_fs"]
        accumulated_expect = ckpt["expect"]
        accumulated_tlist  = ckpt["tlist"]
    else:
        log("No checkpoint found. Starting from B:103.")
        try:
            i0 = labels.index("B:103")
        except ValueError:
            i0 = 0
        rho0 = qt.basis(N, i0) * qt.basis(N, i0).dag()
        current_state = rho0
        t_start = 0.0
        accumulated_expect = [[] for _ in range(N)]
        accumulated_tlist  = []

    total_windows = int(np.ceil((T_TARGET_PS * 1000 - t_start) / (WINDOW_PS * 1000)))
    
    for w in range(total_windows):
        w_start = t_start
        w_end   = min(T_TARGET_PS * 1000, w_start + WINDOW_PS * 1000)
        if w_end <= w_start: break
        
        tlist = np.arange(w_start, w_end + DT_FS, DT_FS)
        log(f"Propagating window {w+1}/{total_windows}: {w_start:.0f} -> {w_end:.0f} fs")
        
        t0 = time.time()
        result = solver.run(current_state, tlist, e_ops=e_ops)
        wall = time.time() - t0
        
        for i in range(N):
            if len(accumulated_tlist) > 0:
                accumulated_expect[i].extend(result.expect[i][1:])
            else:
                accumulated_expect[i].extend(result.expect[i])
        
        if len(accumulated_tlist) > 0:
            accumulated_tlist.extend(tlist[1:])
        else:
            accumulated_tlist.extend(tlist)
            
        current_state = result.final_ados_state
        t_start = w_end
        
        log(f"Window done in {wall/60:.1f} min. Speed: {wall/(w_end-w_start)*1000:.1f} ms/fs")
        
        # Save Checkpoint
        CKPT_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(CKPT_FILE, "wb") as f:
            pickle.dump({
                "t_fs": t_start,
                "ados_state": current_state,
                "expect": accumulated_expect,
                "tlist": accumulated_tlist
            }, f)
            
        # Save partial trajectory
        with open(TRAJ_FILE, "wb") as f:
            pickle.dump({
                "tlist": accumulated_tlist,
                "expect": accumulated_expect,
                "labels": labels
            }, f)
        log(f"Checkpoint saved at t={t_start} fs")

    log("=== PRODUCTION FINISHED ===")

if __name__ == "__main__":
    main()
