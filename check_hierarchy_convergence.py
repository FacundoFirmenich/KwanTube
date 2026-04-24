import numpy as np
import qutip as qt
import sys
from pathlib import Path
import time

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[1] # La raíz es biofisicaquantiqaCLINE
if str(PROJECT_ROOT / "git_repo" / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "git_repo" / "src"))

from qmc_mt.heom_benchmark import run_heom, select_fragment, analyze

def main():
    print("=== Coherent-Window Hierarchy Convergence Check (NC=3 vs NC=4) ===")
    
    # Intentar localizar el archivo .npz en la raíz verdadera del proyecto
    npz_path = PROJECT_ROOT / "H_6DPU.npz"
    if not npz_path.exists():
        # Fallback a la ubicación secundaria en git_repo
        npz_path = PROJECT_ROOT / "git_repo" / "H_6DPU.npz"
    
    if not npz_path.exists():
        print(f"ERROR: No se encuentra {npz_path}. Verifica la ubicación de los Hamiltonianos.")
        sys.exit(1)

    d = np.load(npz_path, allow_pickle=True)
    H_cm = d["H_cm1"]; labels = list(d["labels"])
    H_sub, lab_sub, i0_sub, idx = select_fragment(H_cm, labels, "A:346", n_sites=4)
    
    # Parámetros para test de convergencia coherente
    T_MAX_TEST = 500.0  # fs
    DT_TEST = 1.0       # fs (dt fino para captar oscilaciones)
    
    import qmc_mt.heom_benchmark as hb
    hb.T_MAX = T_MAX_TEST
    hb.DT = DT_TEST

    # Run NC=3
    print(f"\n[RUNNING NC=3, t_max={T_MAX_TEST}fs, dt={DT_TEST}fs]")
    tlist, _, _, _, _, rho_t_3, _, _ = run_heom(H_sub, lab_sub, i0_sub, "6DPU_NC3", nc_override=3)
    
    # Run NC=4
    print(f"\n[RUNNING NC=4, t_max={T_MAX_TEST}fs, dt={DT_TEST}fs]")
    tlist, _, _, _, _, rho_t_4, _, _ = run_heom(H_sub, lab_sub, i0_sub, "6DPU_NC4", nc_override=4)
    
    # Métrica: Norma de la diferencia integrada
    diffs = [ (rho_t_4[i] - rho_t_3[i]).norm() for i in range(len(tlist)) ]
    max_diff = np.max(diffs)
    mean_diff = np.mean(diffs)
    
    print("\n" + "="*45)
    print(f"Max difference (NC3 vs NC4):  {max_diff:.2e}")
    print(f"Mean difference (NC3 vs NC4): {mean_diff:.2e}")
    print(f"Threshold (Scientific):             1.00e-03")
    
    if max_diff < 1e-3:
        print("VERDICT: NC=3 is CONVERGED for coherent dynamics.")
    else:
        print(f"VERDICT: NC=3 is NOT converged. Error {max_diff:.2e} > 1e-3.")
    print("="*45)

if __name__ == "__main__":
    main()
