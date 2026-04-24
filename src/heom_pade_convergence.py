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

# Configuración de rutas
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT / "git_repo" / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "git_repo" / "src"))

try:
    from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver
except ImportError:
    from qutip.nonmarkov.heom import HEOMSolver, DrudeLorentzBath as DrudeLorentzPadeBath

from qmc_mt.heom_benchmark import select_fragment

# --- Constantes Físicas (Standard) ---
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
LAM_CM, GAM_CM, T_K = 35.0, 53.0, 300.0
KT_RADFS = T_K * 0.69503 * CM_TO_RADFS
LAM_RADFS = LAM_CM * CM_TO_RADFS
GAM_RADFS = GAM_CM * CM_TO_RADFS

def run_converged_check():
    # 1. Cargar el Hamiltoniano real de 6DPU
    npz_path = PROJECT_ROOT / "H_6DPU.npz"
    if not npz_path.exists():
        print(f"Error: {npz_path} no encontrado.")
        return
        
    d = np.load(npz_path, allow_pickle=True)
    H_cm = d["H_cm1"]; labels = list(d["labels"])
    
    # 2. SELECCIÓN CRÍTICA: Empezar en A:346 (el "hotspot" de acoplamiento)
    H_sub, lab_sub, i0_sub, _ = select_fragment(H_cm, labels, "A:346", n_sites=4)
    N = H_sub.shape[0]
    
    # Centrar Hamiltoniano (Shift diagonal a cero)
    H_rad = (H_sub - np.mean(np.diag(H_sub)) * np.eye(N)) * CM_TO_RADFS
    H_S = qt.Qobj(H_rad)
    rho0 = qt.basis(N, i0_sub) * qt.basis(N, i0_sub).dag()
    
    print(f"--- Iniciando Auditoría de Convergencia ---")
    print(f"Sitio inicial: {lab_sub[i0_sub]} | Fragmento: {lab_sub}")
    print(f"Energías relativas (rad/fs): {np.diag(H_rad)}")
    
    results_summary = {}

    def execute_run(nc, nk, label):
        print(f"\n> Corriendo {label} (NC={nc}, Nk={nk})...")
        baths = []
        L_total = qt.liouvillian(H_S)
        
        for i in range(N):
            Q = qt.basis(N, i) * qt.basis(N, i).dag()
            bath = DrudeLorentzPadeBath(Q, lam=LAM_RADFS, gamma=GAM_RADFS, T=KT_RADFS, Nk=nk)
            baths.append(bath)
            try:
                _, L_term = bath.terminator()
                L_total += L_term
            except: pass
            
        solver = HEOMSolver(L_total, baths, max_depth=nc, 
                            options={"nsteps": 10000, "store_states": True})
        
        tlist = np.linspace(0, 500, 51)
        t0 = time.time()
        result = solver.run(rho0, tlist)
        wall = time.time() - t0
        
        final_rho = result.states[-1]
        diag = np.real(np.diag(final_rho.full()))
        tr = final_rho.tr().real
        herm = (final_rho - final_rho.dag()).norm()
        min_eig = min(final_rho.eigenenergies())
        # Coherence: sum of absolute off-diagonals
        coh = np.sum(np.abs(final_rho.full())) - np.sum(np.abs(np.diag(final_rho.full())))

        print(f"  Completado en {wall:.2f}s")
        print(f"  Final Diag: {diag}")
        print(f"  Trace: {tr:.6f} | Herm: {herm:.2e} | MinEig: {min_eig:.2e}")

        stats = {
            "NC": nc, "Nk": nk,
            "init_site_index": int(i0_sub),
            "init_site_label": lab_sub[i0_sub],
            "fragment_labels": lab_sub,
            "site0_final": float(diag[0]),
            "diag_final": diag.tolist(),
            "trace_final": float(tr),
            "hermiticity_error_final": float(herm),
            "min_eig_final": float(min_eig),
            "coh_final": float(coh),
            "wall": float(wall),
            "rho_t": result.states,
            "tlist": tlist.tolist()
        }
        results_summary[f"NC{nc}_Nk{nk}"] = stats
        
        fname = f"pade_refined_ckpt_NC{nc}_Nk{nk}.pkl"
        with open(PROJECT_ROOT / fname, "wb") as f:
            pickle.dump({
                "NC": nc, "Nk": nk, "rho_t": result.states, "tlist": tlist, 
                "wall": wall, "labels": lab_sub, "H_S": H_S, "stats": stats
            }, f)
        return stats

    # Runs
    for nc in [3, 5, 7]:
        execute_run(nc, 1, f"Sweep_NC{nc}")
    
    # Validation Nk=2
    s1 = results_summary.get("NC5_Nk1")
    s2 = execute_run(5, 2, "Validation_Nk2")

    if s1 and s2:
        r1 = s1["rho_t"]
        r2 = s2["rho_t"]
        t = s1["tlist"]
        
        diffs = [ (r1[i] - r2[i]).norm() for i in range(len(t)) ]
        max_diff = np.max(diffs)
        idx_max = np.argmax(diffs)
        
        pop_diffs = [ np.abs(np.real(np.diag(r1[i].full() - r2[i].full()))) for i in range(len(t)) ]
        pop_diffs = np.array(pop_diffs)
        max_pop_idx = np.unravel_index(np.argmax(pop_diffs), pop_diffs.shape)
        
        print("\n" + "="*40)
        print("COMPARISON NC=5 (Nk=1 vs Nk=2)")
        print("="*40)
        print(f"Max Abs Diff (Frobenius): {max_diff:.2e} at {t[idx_max]} fs")
        print(f"Max Pop Diff:             {np.max(pop_diffs):.2e} at {t[max_pop_idx[0]]} fs (site {max_pop_idx[1]})")
        
        results_summary["comparison_NC5"] = {
            "max_abs_diff": float(max_diff),
            "time_of_max_abs_diff": float(t[idx_max]),
            "site_of_max_pop_diff": int(max_pop_idx[1]),
            "max_pop_diff": float(np.max(pop_diffs))
        }

    # Save summary JSON
    json_path = PROJECT_ROOT / "heom_convergence_summary.json"
    with open(json_path, "w") as f:
        json_payload = {k: {sk: sv for sk, sv in v.items() if sk != "rho_t"} for k, v in results_summary.items()}
        json.dump(json_payload, f, indent=2)

    # Save summary TXT (Audit log)
    txt_path = PROJECT_ROOT / "heom_convergence_report.txt"
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("=== HEOM CONVERGENCE AUDIT REPORT ===\n")
        f.write(f"Timestamp: {time.ctime()}\n")
        f.write(f"Fragment: {results_summary.get('NC3_Nk1', {}).get('fragment_labels', [])}\n\n")
        for k, v in results_summary.items():
            if k.startswith("NC"):
                f.write(f"RUN {k}:\n")
                f.write(f"  Final Diag: {v['diag_final']}\n")
                f.write(f"  Trace: {v['trace_final']:.6f} | Herm: {v['hermiticity_error_final']:.2e} | MinEig: {v['min_eig_final']:.2e}\n")
                f.write(f"  Wallclock: {v['wall']:.2f}s\n\n")
        
        comp = results_summary.get("comparison_NC5")
        if comp:
            f.write("=== COMPARISON NC=5 (Nk=1 vs Nk=2) ===\n")
            f.write(f"Max Abs Diff (Frobenius): {comp['max_abs_diff']:.2e} at {comp['time_of_max_abs_diff']} fs\n")
            f.write(f"Max Pop Diff:             {comp['max_pop_diff']:.2e} at {comp['time_of_max_abs_diff']} fs (site {comp['site_of_max_pop_diff']})\n")

    print(f"\nResultados guardados en:\n  - {json_path.name}\n  - {txt_path.name}")

if __name__ == "__main__":
    run_converged_check()



