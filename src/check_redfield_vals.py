import numpy as np
import sys
from pathlib import Path

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parent # Está en git_repo/
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

# Localizar el archivo de resultados de Redfield
npz_path = PROJECT_ROOT / "figures_final" / "redfield_1JFF.npz"

if not npz_path.exists():
    print(f"ERROR: No se encuentra {npz_path}. Ejecuta primero redfield_tubulin.py.")
    sys.exit(1)

d = np.load(npz_path, allow_pickle=True)
t = d["t_fs"]
P_site = d["P_site"]
labels = list(d["labels"])

print(f"=== Auditoría de Resultados Redfield (1JFF) ===")
print(f"Redfield P0(500fs) = {P_site[500, 0]:.4f}")

try:
    i0 = labels.index('B:103')
    print(f"Initial site B:103 index: {i0}")
    print(f"P_init(500fs) = {P_site[500, i0]:.4f}")
except ValueError:
    print("WARNING: 'B:103' no encontrado en las etiquetas.")
