#!/usr/bin/env python3
"""
heom_pade_convergence.py -- THE TRUTH VERSION
Refined convergence check with proper site selection and audit trail.
"""
import numpy as np
import qutip as qt
import sys
import pickle, time, json
from pathlib import Path

# Configuracion de rutas Tier-0
PROJECT_ROOT = Path(__file__).resolve().parents[1] # retrocede desde src/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

try:
    from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver
except ImportError:
    from qutip.nonmarkov.heom import HEOMSolver, DrudeLorentzBath as DrudeLorentzPadeBath

from qmc_mt.heom_benchmark import select_fragment

# --- Constantes Fisicas (Standard) ---
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
LAM_CM, GAM_CM, T_K = 35.0, 53.0, 300.0
KT_RADFS = T_K * 0.69503 * CM_TO_RADFS
LAM_RADFS = LAM_CM * CM_TO_RADFS
GAM_RADFS = GAM_CM * CM_TO_RADFS

def run_converged_check():
    # 1. Cargar el Hamiltoniano real de 6DPU - Sincronizado con raw_npz
    npz_path = PROJECT_ROOT / "outputs_data" / "raw_npz" / "H_6DPU.npz"
    if not npz_path.exists():
        npz_path = PROJECT_ROOT / "H_6DPU.npz"
        
    if not npz_path.exists():
        print(f"Error: No se encuentra H_6DPU.npz en {npz_path}")
        return
        
    d = np.load(npz_path, allow_pickle=True)
    H_cm = d["H_cm1"]; labels = list(d["labels"])
    
    # 2. SELECCIÓN CRÍTICA: Empezar en A:346 (el "hotspot" de acoplamiento)
    H_sub, lab_sub, i0_sub, _ = select_fragment(H_cm, labels, "A:346", n_sites=4)
    N = H_sub.shape[0]
    
    # Centrar Hamiltoniano
    H_rad = (H_sub - np.mean(np.diag(H_sub)) * np.eye(N)) * CM_TO_RADFS
    H_S = qt.Qobj(H_rad)
    rho0 = qt.basis(N, i0_sub) * qt.basis(N, i0_sub).dag()
    
    print(f"--- Iniciando Auditoria de Convergencia ---")
    print(f"Sitio inicial: {lab_sub[i0_sub]} | Fragmento: {lab_sub}")
    
    results_summary = {}

    def execute_run(nc, nk, label):
        print(f"\n> Corriendo {label} (NC={nc}, Nk={nk})...")
        baths = []
        for i in range(N):
            Q = qt.basis(N, i) * qt.basis(N, i).dag()
            bath = DrudeLorentzPadeBath(Q, lam=LAM_RADFS, gamma=GAM_RADFS, T=KT_RADFS, Nk=nk)
            baths.append(bath)
            
        solver = HEOMSolver(H_S, baths, max_depth=nc, 
                            options={"nsteps": 10000, "store_states": True})
        
        tlist = np.linspace(0, 500, 51)
        t0 = time.time()
        result = solver.run(rho0, tlist)
        wall = time.time() - t0
        
        final_rho = result.states[-1]
        diag = np.real(np.diag(final_rho.full()))
        tr = final_rho.tr().real
        # Coherence metric
        coh = np.sum(np.abs(final_rho.full())) - np.sum(np.abs(np.diag(final_rho.full())))

        print(f"  Completado en {wall:.2f}s")
        print(f"  Final Diag: {diag}")

        stats = {
            "NC": nc, "Nk": nk,
            "diag_final": diag.tolist(),
            "trace_final": float(tr),
            "coh_final": float(coh),
            "wall": float(wall),
            "rho_t": result.states,
            "tlist": tlist.tolist()
        }
        results_summary[f"NC{nc}_Nk{nk}"] = stats
        
        # Guardar PKL en raw_pkl
        fname = f"pade_refined_ckpt_NC{nc}_Nk{nk}.pkl"
        out_pkl = PROJECT_ROOT / "outputs_data" / "raw_pkl" / fname
        out_pkl.parent.mkdir(parents=True, exist_ok=True)
        with open(out_pkl, "wb") as f:
            pickle.dump({
                "NC": nc, "Nk": nk, "rho_t": result.states, "tlist": tlist, 
                "wall": wall, "labels": lab_sub, "stats": stats
            }, f)
        return stats

    # Sweep NC
    for nc in [3, 5, 7]:
        execute_run(nc, 1, f"Sweep_NC{nc}")
    
    # Validation Nk=2
    s1 = results_summary.get("NC5_Nk1")
    s2 = execute_run(5, 2, "Validation_Nk2")

    if s1 and s2:
        r1 = s1["rho_t"]; r2 = s2["rho_t"]; t = s1["tlist"]
        diffs = [ (r1[i] - r2[i]).norm() for i in range(len(t)) ]
        max_diff = np.max(diffs)
        
        results_summary["comparison_NC5"] = {
            "max_abs_diff": float(max_diff),
            "time_of_max_abs_diff": float(t[np.argmax(diffs)])
        }

    # Save summary JSON - raw_json
    json_path = PROJECT_ROOT / "outputs_data" / "raw_json" / "heom_convergence_summary.json"
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(json_path, "w") as f:
        json_payload = {k: {sk: sv for sk, sv in v.items() if sk != "rho_t"} for k, v in results_summary.items()}
        json.dump(json_payload, f, indent=2)

    # Save summary TXT - raw_txt+md
    txt_path = PROJECT_ROOT / "outputs_data" / "raw_txt+md" / "heom_convergence_report.txt"
    txt_path.parent.mkdir(parents=True, exist_ok=True)
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("=== HEOM CONVERGENCE AUDIT REPORT ===\n")
        f.write(f"Timestamp: {time.ctime()}\n\n")
        for k, v in results_summary.items():
            if k.startswith("NC"):
                f.write(f"RUN {k}:\n")
                f.write(f"  Final Diag: {v['diag_final']}\n")
                f.write(f"  Wallclock: {v['wall']:.2f}s\n\n")

    print(f"\nResultados guardados en outputs_data/")

if __name__ == "__main__":
    run_converged_check()
