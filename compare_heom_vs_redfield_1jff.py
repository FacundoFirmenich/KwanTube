import numpy as np
import pickle
import sys
import json
from pathlib import Path

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[1] # La raíz es biofisicaquantiqaCLINE
if str(PROJECT_ROOT / "git_repo" / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "git_repo" / "src"))

def main():
    print("="*60)
    print("  DIAGNOSTICO CIENTIFICO: 1JFF (HEOM vs Redfield)")
    print("="*60)

    # 1. Load Redfield
    red_path = PROJECT_ROOT / "figures_final" / "redfield_1JFF.npz"
    if not red_path.exists():
        print(f"ERROR: No se encuentra {red_path}. Ejecuta primero redfield_tubulin.py.")
        return

    red = np.load(red_path, allow_pickle=True)
    t_red = red["t_fs"]
    P_site_red = red["P_site"] # (nt, N)

    # 2. Load HEOM NC=6 (from calibration data)
    calib_path = PROJECT_ROOT / "git_repo" / "jff_calib_data.npz"
    if not calib_path.exists():
        print(f"ERROR: No se encuentra {calib_path}.")
        return

    calib = np.load(calib_path, allow_pickle=True)
    t_heom = calib["tlist"]
    P6 = calib["P6"] # (N, nt)
    P6 = P6.T # (nt, N)
    
    # Compare at ~500fs using robust indexing
    idx_500_red = np.argmin(np.abs(t_red - 500.0))
    idx_500_heom = np.argmin(np.abs(t_heom - 500.0))
    
    p_red_500 = P_site_red[idx_500_red]
    p_heom_500 = P6[idx_500_heom]
    
    # B:103 es el sitio inicial en 1JFF (índice 5)
    site_idx = 5
    
    diff = p_heom_500 - p_red_500
    max_diff = np.max(np.abs(diff))
    rel_diff = max_diff / np.max(p_heom_500)
    
    print(f"Poblaciones a ~500fs (Sitio B:103 init):")
    print(f"  Redfield:  {p_red_500[site_idx]:.4f} (Sitio {site_idx})")
    print(f"  HEOM NC=6: {p_heom_500[site_idx]:.4f} (Sitio {site_idx})")
    print(f"  Discrepancia Absoluta Max: {max_diff:.2e}")
    print(f"  Discrepancia Relativa Max: {rel_diff*100:.2f}%")
    
    print(f"\nRatio r(1JFF): {calib['ratio']:.3f}")
    print(f"Error Truncacion NC=7: {calib['eps7']*100:.2f}%")
    
    # Save a report for the SI
    report = {
        "max_redfield_deviation": float(max_diff),
        "heom_ratio": float(calib['ratio']),
        "truncation_error_nc7": float(calib['eps7']),
        "p_red_500": p_red_500.tolist(),
        "p_heom_500": p_heom_500.tolist()
    }
    
    report_path = PROJECT_ROOT / "heom_vs_redfield_report.json"
    with report_path.open("w") as f:
        json.dump(report, f, indent=2)
    print(f"\nReporte guardado en: {report_path.name}")

if __name__ == "__main__":
    main()
