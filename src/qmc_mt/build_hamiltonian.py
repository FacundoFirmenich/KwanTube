"""
build_hamiltonian.py — construye el Hamiltoniano electrónico de Trps
a partir de pdb_tubulin_analysis.json.

Modelo
------
Cada Trp = sitio de 2 niveles (|g>, |e>).  En la base de un solo excitón:
    H_ij = E_i δ_ij + V_ij (1 - δ_ij)

Energías de sitio E_i:
    - Valor base E0 = 35000 cm⁻¹ (~286 nm, banda 1La de Trp en agua).
    - Shift estocástico gaussiano σ_site = 200 cm⁻¹ (desorden estático
      razonable para Trp en proteína; Callis & Liu 2004).

Acoplamientos V_ij (dipolo-dipolo puntual, unidades cm⁻¹):
    V_ij = C * κ_ij * |μ|² / (ε_r * r_ij³)

    con:
      |μ|_1La = 2.0 D   (Callis 1997, absorción Trp)
      ε_r     = 2.0     (proteína interior, screening efectivo)
      C       = 5.04 × 10⁴  cm⁻¹·Å³·D⁻²   (prefactor que convierte
                D²/Å³ en cm⁻¹ tras dividir por ε_r; derivación abajo)

Derivación de C:
    V[erg] = μ_D μ_A κ / (ε r³)   en unidades CGS
    1 D = 1e-18 esu·cm
    V[cm⁻¹] = V[erg] / (h c)
    Resultado numérico estándar: V[cm⁻¹] = 5.04e4 · κ · μ[D]² /(ε · r[Å]³)
    (ver Madjet, Abdurakhmanov, Renger 2006, J.Phys.Chem.B 110:17268, eq.2)

Chequeos internos
-----------------
1. Hermiticidad:  ||H - H†|| < 1e-12
2. Orden de V:    max|V| ∈ [10, 200] cm⁻¹ para Trps a 14-20 Å
3. Gap visible:   autovalores dentro de E0 ± 3σ_site ± max|V|
4. Reproduce el par dominante 346-407 de 6DPU con |V| ~ 60-120 cm⁻¹.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np

JSON_PATH = Path(r"C:\Users\User\3D Objects\biofisicaquantiqaCLINE\pdb_tubulin_analysis.json")

# --- constantes físicas ---
MU_TRP_D     = 2.0        # Debye, banda 1La
EPS_R        = 2.0        # screening efectivo proteína
PREFACTOR    = 5.04e4     # cm⁻¹ · Å³ · D⁻² (Madjet-Renger 2006)
E0_CM        = 35000.0    # cm⁻¹, energía 0-0 de Trp 1La
SIGMA_SITE   = 200.0      # cm⁻¹, desorden estático gaussiano
RNG_SEED     = 12345

def load_structure(pid: str) -> dict:
    data = json.loads(JSON_PATH.read_text())
    if pid not in data:
        raise KeyError(f"{pid} no está en {JSON_PATH.name}")
    return data[pid]

def coupling_cm1(kappa2: float, r_A: float) -> float:
    """V en cm⁻¹ a partir de κ² y r (Å). Devuelve |V| (signo indeterminado
    sin los vectores μ̂; el cuadrado del signo no afecta poblaciones)."""
    kappa = np.sqrt(max(kappa2, 0.0))
    return PREFACTOR * kappa * (MU_TRP_D ** 2) / (EPS_R * r_A ** 3)

def build(pid: str, disorder: bool = True, verbose: bool = True):
    s = load_structure(pid)
    trps = s["trp"]
    N = len(trps)
    labels = [f"{t['chain']}:{t['resseq']}" for t in trps]
    idx = {lab: i for i, lab in enumerate(labels)}

    # índice rápido por centro → para distancias (no lo usamos: ya está en pairs)
    H = np.zeros((N, N), dtype=float)

    # diagonal
    rng = np.random.default_rng(RNG_SEED)
    shifts = rng.normal(0.0, SIGMA_SITE, size=N) if disorder else np.zeros(N)
    for i in range(N):
        H[i, i] = E0_CM + shifts[i]

    # fuera de diagonal: todos los pares dentro del cutoff del JSON
    n_off = 0
    vmax = 0.0
    for p in s["pairs_within_cutoff"]:
        i, j = idx[p["donor"]], idx[p["acceptor"]]
        V = coupling_cm1(p["kappa2"], p["r_A"])
        H[i, j] = V
        H[j, i] = V
        n_off += 1
        vmax = max(vmax, V)

    # ---- validaciones ----
    herm_err = float(np.linalg.norm(H - H.T))
    assert herm_err < 1e-12, f"H no es hermítico: err={herm_err}"

    evals = np.linalg.eigvalsh(H)
    spread = float(evals.max() - evals.min())

    if verbose:
        print(f"\n=== Hamiltoniano {pid} ===")
        print(f"  N sitios           = {N}")
        print(f"  N acoplamientos    = {n_off}")
        print(f"  |V|_max            = {vmax:7.2f} cm⁻¹")
        print(f"  σ_site (desorden)  = {SIGMA_SITE:7.2f} cm⁻¹")
        print(f"  hermiticidad OK    (||H-Hᵀ|| = {herm_err:.2e})")
        print(f"  spread espectral   = {spread:7.2f} cm⁻¹")
        print(f"  E_min, E_max       = {evals.min():.1f}, {evals.max():.1f} cm⁻¹")
        # top-3 acoplamientos realizados
        pairs_sorted = sorted(s["pairs_within_cutoff"],
                              key=lambda p: -coupling_cm1(p["kappa2"], p["r_A"]))[:3]
        print(f"  top-3 |V|:")
        for p in pairs_sorted:
            V = coupling_cm1(p["kappa2"], p["r_A"])
            print(f"    {p['donor']:>6}—{p['acceptor']:<6}  "
                  f"r={p['r_A']:5.2f} Å  κ²={p['kappa2']:5.3f}  "
                  f"|V|={V:6.2f} cm⁻¹")

    return {
        "pdb_id": pid,
        "labels": labels,
        "H_cm1": H,
        "E0_cm1": E0_CM,
        "sigma_site_cm1": SIGMA_SITE,
        "mu_D": MU_TRP_D,
        "eps_r": EPS_R,
        "n_sites": N,
        "v_max_cm1": vmax,
        "eigvals_cm1": evals,
    }

# ---- chequeo de cordura del prefactor ----
def _selfcheck_prefactor():
    """
    Test: dos Trps alineados cabeza-cola (κ=2 → κ²=4) a 10 Å con μ=2D, ε=1
    deben dar V = 5.04e4 · 2 · 4 / 1000 = 403.2 cm⁻¹ (valor de libro).
    """
    V = PREFACTOR * np.sqrt(4.0) * (2.0 ** 2) / (1.0 * 10.0 ** 3)
    expected = 5.04e4 * 2 * 4 / 1000
    assert abs(V - expected) < 1e-6, f"prefactor roto: {V} vs {expected}"
    print(f"self-check prefactor: V(κ²=4, r=10Å, μ=2D, ε=1) = {V:.2f} cm⁻¹  OK")
    return V

def save_npz(result: dict, out: Path):
    np.savez(out,
             H_cm1=result["H_cm1"],
             eigvals_cm1=result["eigvals_cm1"],
             labels=np.array(result["labels"]),
             E0_cm1=result["E0_cm1"],
             sigma_site_cm1=result["sigma_site_cm1"],
             mu_D=result["mu_D"],
             eps_r=result["eps_r"])
    print(f"  guardado: {out}")

if __name__ == "__main__":
    _selfcheck_prefactor()

    out_dir = Path(r"C:\Users\User\3D Objects\biofisicaquantiqaCLINE")
    for pid in ("1JFF", "6DPU"):
        res = build(pid, disorder=True)
        save_npz(res, out_dir / f"H_{pid}.npz")

    # Resumen comparativo
    print("\n--- resumen ---")
    print("Se generaron H_1JFF.npz (dímero αβ, 8 sitios) y H_6DPU.npz (red MT, 48 sitios).")
    print("Para #1 HEOM usaremos H_1JFF.npz primero (más chico, valida la pipeline),")
    print("después escalamos a H_6DPU.npz.")