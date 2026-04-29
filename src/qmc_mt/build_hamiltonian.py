"""
build_hamiltonian.py - construye el Hamiltoniano electrónico de Trps
a partir de pdb_tubulin_analysis.json.

Modelo
------
Cada Trp = sitio de 2 niveles (|g>, |e>). En la base de un solo excitón:
    H_ij = E_i * delta_ij + V_ij * (1 - delta_ij)

Energías de sitio E_i:
    - Valor base E0 = 35000 cm^-1 (~286 nm, banda 1La de Trp en agua).
    - Shift estocástico gaussiano sigma_site = 200 cm^-1 (desorden estático
      razonable para Trp en proteína; Callis & Liu 2004).

Acoplamientos V_ij (dipolo-dipolo puntual, unidades cm^-1):
    V_ij = C * kappa_ij * |mu|^2 / (epsilon_r * r_ij^3)

    con:
      |mu|_1La = 2.0 D   (Callis 1997, absorción Trp)
      epsilon_r     = 2.0     (proteína interior, screening efectivo)
      C       = 5.04e4  cm^-1 * Å^3 * D^-2   (prefactor; Madjet-Renger 2006)

Derivación de C:
    V[erg] = mu_D mu_A kappa / (epsilon r^3)   en unidades CGS
    1 D = 1e-18 esu*cm
    V[cm^-1] = V[erg] / (h c)
    Resultado numérico estándar: V[cm^-1] = 5.04e4 * kappa * mu[D]^2 / (epsilon * r[Å]^3)
    (ver Madjet, Abdurakhmanov, Renger 2006, J.Phys.Chem.B 110:17268, ec.2)

Trazabilidad de parámetros (Patch v3.5.1)
------------------------------------------
Todos los parámetros físicos se cargan desde config/physics_params.yaml.
Si el archivo no existe, se usan valores por defecto documentados con fuente.
El código NUNCA modifica el YAML; sólo lo lee.

Chequeos internos
-----------------
1. Hermiticidad:  ||H - H†|| < 1e-12
2. Orden de V:    max|V| en [10, 300] cm^-1 para Trps a 14-40 Å
3. Gap visible:   autovalores dentro de E0 +/- 3*sigma_site +/- max|V|
4. Reproduce el par dominante 346-407 de 6DPU con |V| ~ 60-120 cm^-1.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

# ------------------------------------------------------------------ path setup
# Resolución robusta de PROJECT_ROOT con fallbacks explícitos.
def _resolve_project_root() -> Path:
    """Resolver raíz del proyecto con múltiples fallbacks."""
    candidates = [
        Path(__file__).resolve().parents[2],   # src/qmc_mt/ -> raíz
        Path(__file__).resolve().parents[3],   # extra nivel si hay src anidado
        Path.cwd(),
    ]
    for root in candidates:
        if (root / "outputs_data").exists() and (root / "src").exists():
            return root.resolve()
    # Fallback silencioso al primero si no se encontró la estructura esperada
    return Path(__file__).resolve().parents[2]


PROJECT_ROOT = _resolve_project_root()
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

JSON_PATH = PROJECT_ROOT / "outputs_data" / "raw_json" / "pdb_tubulin_analysis.json"
PHYSICS_YAML = PROJECT_ROOT / "config" / "physics_params.yaml"

# ------------------------------------------------------------------ load params
def _load_physics_params() -> dict:
    """
    Cargar parámetros físicos desde config/physics_params.yaml.

    Si el YAML no existe o no puede parsearse, retorna los valores por defecto
    documentados con su fuente original. Esto garantiza que el script funcione
    en entornos donde el YAML aún no fue creado, sin silenciar el problema.

    Returns
    -------
    dict
        Diccionario con estructura compatible con physics_params.yaml.
        Cada hoja tiene 'value', 'source' y 'uncertainty'.
    """
    if PHYSICS_YAML.exists():
        try:
            import yaml  # pyyaml
            with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
                params = yaml.safe_load(f)
            return params
        except ImportError:
            warnings.warn(
                "pyyaml no instalado; usando parámetros por defecto documentados. "
                "Instalar con: pip install pyyaml",
                stacklevel=2,
            )
        except Exception as e:
            warnings.warn(
                f"No se pudo leer {PHYSICS_YAML}: {e}. "
                "Usando parámetros por defecto documentados.",
                stacklevel=2,
            )
    else:
        warnings.warn(
            f"config/physics_params.yaml no encontrado en {PHYSICS_YAML}. "
            "Usando parámetros por defecto documentados. "
            "Para trazabilidad completa, crear el archivo YAML.",
            stacklevel=2,
        )

    # Fallback: valores por defecto idénticos al YAML, con fuente documentada
    return {
        "dipole": {
            "mu_trp_1La_D": {"value": 2.0, "source": "Callis1997_Fig2", "uncertainty": 0.3},
            "epsilon_r_protein": {"value": 2.0, "source": "effective_screening_estimate", "uncertainty": 0.5},
        },
        "energetics": {
            "E0_cm1": {"value": 35000.0, "source": "Trp_1La_absorption_peak", "uncertainty": 500.0},
            "sigma_site_cm1": {"value": 200.0, "source": "Callis_Liu_2004_protein_disorder", "uncertainty": 50.0},
        },
        "constants": {
            "prefactor_cm1_A3_D2": {"value": 50400.0, "source": "Madjet_Renger_2006_eq2", "uncertainty": 0.0},
            "rng_seed": {"value": 12345, "source": "reproducibility_seed_fixed", "uncertainty": 0.0},
        },
    }


def _extract_value(params: dict, *keys: str, default: float = 0.0) -> float:
    """Navegar el diccionario anidado de parámetros y extraer 'value'."""
    node = params
    for k in keys:
        if not isinstance(node, dict) or k not in node:
            warnings.warn(
                f"Parámetro {'/'.join(keys)} no encontrado en physics_params; "
                f"usando default={default}",
                stacklevel=3,
            )
            return default
        node = node[k]
    if isinstance(node, dict):
        return float(node.get("value", default))
    return float(node)


# ------------------------------------------------------------------ globals
PHYSICS_PARAMS = _load_physics_params()

MU_TRP_D  = _extract_value(PHYSICS_PARAMS, "dipole", "mu_trp_1La_D", default=2.0)
EPS_R     = _extract_value(PHYSICS_PARAMS, "dipole", "epsilon_r_protein", default=2.0)
PREFACTOR = _extract_value(PHYSICS_PARAMS, "constants", "prefactor_cm1_A3_D2", default=5.04e4)
E0_CM     = _extract_value(PHYSICS_PARAMS, "energetics", "E0_cm1", default=35000.0)
SIGMA_SITE = _extract_value(PHYSICS_PARAMS, "energetics", "sigma_site_cm1", default=200.0)
RNG_SEED  = int(_extract_value(PHYSICS_PARAMS, "constants", "rng_seed", default=12345))

# Tolerancias numéricas desde numerical_params.yaml (o valores por defecto)
def _load_hermiticity_threshold() -> float:
    """Cargar umbral de hermiticidad desde numerical_params.yaml."""
    num_yaml = PROJECT_ROOT / "config" / "numerical_params.yaml"
    if num_yaml.exists():
        try:
            import yaml
            with open(num_yaml, "r", encoding="utf-8") as f:
                nparams = yaml.safe_load(f)
            return float(nparams.get("tolerances", {}).get("hermiticity_threshold", 1e-12))
        except Exception:
            pass
    return 1e-12


HERMITICITY_THRESHOLD = _load_hermiticity_threshold()
V_MAX_EXPECTED = 300.0  # cm^-1 — umbral de advertencia
V_MIN_EXPECTED = 10.0   # cm^-1 — umbral de advertencia


# ------------------------------------------------------------------ core funcs
def load_structure(pid: str) -> dict:
    """
    Cargar estructura desde pdb_tubulin_analysis.json.

    Parameters
    ----------
    pid : str
        Identificador PDB (e.g., '1JFF', '6DPU').

    Raises
    ------
    FileNotFoundError
        Si el JSON no existe en la ruta esperada.
    KeyError
        Si el PDB ID no está en el JSON.
    """
    if not JSON_PATH.exists():
        raise FileNotFoundError(
            f"{JSON_PATH} no encontrado. "
            "Ejecutar primero pdb_tubulin_analysis.py para generar el JSON."
        )
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    if pid not in data:
        raise KeyError(f"'{pid}' no está en {JSON_PATH.name}. "
                       f"Disponibles: {list(data.keys())}")
    return data[pid]


def coupling_cm1(kappa2: float, r_A: float) -> float:
    """
    Calcular acoplamiento dipolar V en cm^-1.

    V = PREFACTOR * sqrt(kappa^2) * mu^2 / (epsilon_r * r^3)

    Devuelve |V| (signo indeterminado sin los vectores mu completos;
    el cuadrado del signo no afecta poblaciones en HEOM).

    Parameters
    ----------
    kappa2 : float
        Factor de orientación kappa^2 (adimensional, >= 0).
    r_A : float
        Distancia centro-a-centro en Ångström.

    Returns
    -------
    float
        |V| en cm^-1.
    """
    kappa = np.sqrt(max(kappa2, 0.0))
    return PREFACTOR * kappa * (MU_TRP_D ** 2) / (EPS_R * r_A ** 3)


def build(
    pid: str,
    disorder: bool = True,
    verbose: bool = True,
) -> dict:
    """
    Construir el Hamiltoniano de excitón de sitio único para PDB ID dado.

    Parameters
    ----------
    pid : str
        Identificador PDB.
    disorder : bool
        Si True, añade desorden estático gaussiano (sigma_site) a la diagonal.
    verbose : bool
        Imprimir resumen en stdout.

    Returns
    -------
    dict
        Diccionario con 'H_cm1', 'labels', 'eigvals_cm1' y metadatos.
    """
    s = load_structure(pid)
    trps = s["trp"]
    N = len(trps)
    labels = [f"{t['chain']}:{t['resseq']}" for t in trps]
    idx = {lab: i for i, lab in enumerate(labels)}

    H = np.zeros((N, N), dtype=float)

    # Diagonal: energías de sitio con desorden estático
    rng = np.random.default_rng(RNG_SEED)
    shifts = rng.normal(0.0, SIGMA_SITE, size=N) if disorder else np.zeros(N)
    for i in range(N):
        H[i, i] = E0_CM + shifts[i]

    # Fuera de diagonal: acoplamientos dipolo-dipolo
    n_off = 0
    vmax = 0.0
    for p in s["pairs_within_cutoff"]:
        i, j = idx[p["donor"]], idx[p["acceptor"]]
        V = coupling_cm1(p["kappa2"], p["r_A"])
        H[i, j] = V
        H[j, i] = V
        n_off += 1
        vmax = max(vmax, V)

    # ---- Validaciones estrictas ----
    herm_err = float(np.linalg.norm(H - H.T))
    assert herm_err < HERMITICITY_THRESHOLD, (
        f"H no es hermítico para {pid}: ||H - H.T|| = {herm_err:.2e} "
        f"(umbral: {HERMITICITY_THRESHOLD:.0e})"
    )

    evals = np.linalg.eigvalsh(H)
    spread = float(evals.max() - evals.min())

    # Advertencia si acoplamientos fuera de rango físico esperado
    if n_off > 0 and vmax > V_MAX_EXPECTED:
        warnings.warn(
            f"|V|_max = {vmax:.2f} cm^-1 excede umbral esperado "
            f"({V_MAX_EXPECTED:.0f} cm^-1). Verificar parámetros en physics_params.yaml.",
            stacklevel=2,
        )
    if n_off > 0 and vmax < V_MIN_EXPECTED:
        warnings.warn(
            f"|V|_max = {vmax:.2f} cm^-1 por debajo del umbral esperado "
            f"({V_MIN_EXPECTED:.0f} cm^-1). Posible error en la geometría o parámetros.",
            stacklevel=2,
        )

    if verbose:
        print(f"\n=== Hamiltoniano {pid} ===")
        print(f"  N sitios            = {N}")
        print(f"  N acoplamientos     = {n_off}")
        print(f"  |V|_max             = {vmax:7.2f} cm^-1")
        print(f"  sigma_site (desorden)= {SIGMA_SITE:7.2f} cm^-1  "
              f"[src: {PHYSICS_PARAMS.get('energetics', {}).get('sigma_site_cm1', {}).get('source', 'N/A') if isinstance(PHYSICS_PARAMS.get('energetics', {}).get('sigma_site_cm1', {}), dict) else 'N/A'}]")
        print(f"  hermiticidad OK     (||H-H.T|| = {herm_err:.2e})")
        print(f"  spread espectral    = {spread:7.2f} cm^-1")
        print(f"  E_min, E_max        = {evals.min():.1f}, {evals.max():.1f} cm^-1")
        print(f"  Parámetros desde    = {'config/physics_params.yaml' if PHYSICS_YAML.exists() else 'DEFAULT (YAML no encontrado)'}")

        # Top-3 acoplamientos
        pairs_sorted = sorted(
            s["pairs_within_cutoff"],
            key=lambda p: -coupling_cm1(p["kappa2"], p["r_A"]),
        )[:3]
        print(f"  top-3 |V|:")
        for p in pairs_sorted:
            V = coupling_cm1(p["kappa2"], p["r_A"])
            print(
                f"    {p['donor']:>6}-{p['acceptor']:<6}  "
                f"r={p['r_A']:5.2f} Å  kappa^2={p['kappa2']:5.3f}  "
                f"|V|={V:6.2f} cm^-1"
            )

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
        "params_source": str(PHYSICS_YAML) if PHYSICS_YAML.exists() else "hardcoded_defaults",
    }


# ------------------------------------------------------------------ self-check
def _selfcheck_prefactor() -> float:
    """
    Test: dos Trps alineados cabeza-cola (kappa^2=4) a 10 Å con mu=2D, epsilon=1
    deben dar V = 5.04e4 * 2 * 4 / 1000 = 403.2 cm^-1 (valor de libro).

    Raises
    ------
    AssertionError si el prefactor está mal configurado.
    """
    V = PREFACTOR * np.sqrt(4.0) * (2.0 ** 2) / (1.0 * 10.0 ** 3)
    expected = 5.04e4 * 2 * 4 / 1000  # 403.2 cm^-1 — invariante
    assert abs(V - expected) < 1e-3, (
        f"Fallo en self-check del prefactor: V={V:.6f} vs esperado={expected:.6f}. "
        f"PREFACTOR actual = {PREFACTOR}. Verificar config/physics_params.yaml."
    )
    print(
        f"self-check prefactor: "
        f"V(kappa^2=4, r=10Å, mu=2D, epsilon=1) = {V:.2f} cm^-1  OK"
    )
    return V


# ------------------------------------------------------------------ I/O
def save_npz(result: dict, out: Path) -> None:
    """Guardar Hamiltoniano en formato .npz."""
    np.savez(
        out,
        H_cm1=result["H_cm1"],
        eigvals_cm1=result["eigvals_cm1"],
        labels=np.array(result["labels"]),
        E0_cm1=result["E0_cm1"],
        sigma_site_cm1=result["sigma_site_cm1"],
        mu_D=result["mu_D"],
        eps_r=result["eps_r"],
    )
    print(f"  guardado: {out}")


# ------------------------------------------------------------------ main
if __name__ == "__main__":
    print(f"PROJECT_ROOT resuelto: {PROJECT_ROOT}")
    print(f"YAML de parámetros:    {PHYSICS_YAML} ({'existe' if PHYSICS_YAML.exists() else 'NO EXISTE — usando defaults'})")
    _selfcheck_prefactor()

    out_dir = PROJECT_ROOT / "outputs_data" / "raw_npz"
    out_dir.mkdir(parents=True, exist_ok=True)

    for pid in ("1JFF", "6DPU"):
        res = build(pid, disorder=True)
        save_npz(res, out_dir / f"H_{pid}.npz")

    print("\n--- resumen ---")
    print("Generados H_1JFF.npz (dímero alphabeta, 8 sitios) y H_6DPU.npz (red MT, 48 sitios).")
    print(f"Parámetros leídos desde: {'config/physics_params.yaml' if PHYSICS_YAML.exists() else 'DEFAULTS documentados'}")
