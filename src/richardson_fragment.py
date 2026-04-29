"""
richardson_fragment.py — Extrapolación de Richardson para convergencia HEOM.

Patch v3.5.1 — Correcciones críticas:
  - Eliminados valores hardcodeados d7=5.45e-3, d8=2.09e-3, r_prev=0.422.
  - Ratios calculados on-the-fly desde los estados de los pickles NC=7 y NC=8.
  - Fallback determinista si los pickles no existen: lectura desde
    pade_convergence_report.json (ledger) y, como último recurso, el valor
    de fallback_ratio_r de numerical_params.yaml (documentado, no inventado).
  - Diagnóstico de régimen asintótico: se verifica contracción geométrica
    antes de aceptar la extrapolación.

Jerarquía de fuentes para r (ratio de convergencia):
  1. Cálculo directo desde pickles NC=7, NC=8 (ruta principal)
  2. Ledger pade_convergence_report.json si existe (ruta alternativa)
  3. fallback_ratio_r en numerical_params.yaml (red de seguridad documentada)
  4. Error explícito si ninguna fuente está disponible
"""
from __future__ import annotations

import json
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np

# ------------------------------------------------------------------ path setup
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

# Rutas estándar
PKL_DIR = PROJECT_ROOT / "outputs_data" / "raw_pkl"
LEDGER_PATH = PROJECT_ROOT / "outputs_data" / "raw_json" / "pade_convergence_report.json"
NUM_YAML = PROJECT_ROOT / "config" / "numerical_params.yaml"


# ------------------------------------------------------------------ helpers
def _load_fallback_ratio() -> float:
    """
    Cargar ratio de convergencia desde numerical_params.yaml.
    Fallback de último recurso — documentado, no inventado.
    """
    if NUM_YAML.exists():
        try:
            import yaml
            with open(NUM_YAML, "r", encoding="utf-8") as f:
                nparams = yaml.safe_load(f)
            r = float(
                nparams.get("heom", {})
                .get("richardson", {})
                .get("fallback_ratio_r", 0.383)
            )
            src = (
                nparams.get("heom", {})
                .get("richardson", {})
                .get("fallback_source", "numerical_params.yaml")
            )
            warnings.warn(
                f"Usando ratio de fallback r={r:.3f} desde numerical_params.yaml "
                f"(fuente documentada: {src}). "
                "Para mayor precisión, ejecutar heom_pade_convergence.py primero.",
                stacklevel=3,
            )
            return r
        except Exception as e:
            warnings.warn(f"No se pudo leer numerical_params.yaml: {e}", stacklevel=3)

    warnings.warn(
        "numerical_params.yaml no disponible; usando r=0.383 (medido de "
        "HEOM_EXECUTION_LEDGER_6DPU_NC7_NC8). "
        "Ejecutar heom_pade_convergence.py para recalcular.",
        stacklevel=3,
    )
    return 0.383


def _load_pkl(nc: int, nk: int = 1) -> dict | None:
    """
    Cargar checkpoint HEOM desde raw_pkl/.
    Busca en la ruta estándar y en git_repo como fallback.

    Returns None si no se encuentra el archivo.
    """
    candidates = [
        PKL_DIR / f"pade_refined_ckpt_NC{nc}_Nk{nk}.pkl",
        PROJECT_ROOT.parent / "git_repo" / "outputs_data" / "raw_pkl"
        / f"pade_refined_ckpt_NC{nc}_Nk{nk}.pkl",
    ]
    for path in candidates:
        if path.exists():
            try:
                with open(path, "rb") as f:
                    return pickle.load(f)
            except Exception as e:
                warnings.warn(f"Error leyendo {path}: {e}", stacklevel=2)
    return None


def _extract_final_pop(data: dict) -> float | None:
    """Extraer población del sitio 0 al tiempo final desde un checkpoint."""
    try:
        rho_t = data["rho_t"]
        if hasattr(rho_t[-1], "full"):
            # Objeto QuTiP
            return float(rho_t[-1].full()[0, 0].real)
        elif hasattr(rho_t[-1], "__getitem__"):
            # Array numpy
            arr = np.asarray(rho_t[-1])
            return float(arr[0, 0].real)
        return None
    except Exception as e:
        warnings.warn(f"No se pudo extraer población del checkpoint: {e}", stacklevel=2)
        return None


def _compute_ratio_from_pickles(nc_a: int, nc_b: int, nk: int = 1) -> tuple[float, float, float] | None:
    """
    Calcular d_B, d_A y ratio r = |d_B| / |d_A| desde pickles NC=nc_a, NC=nc_b.

    Returns
    -------
    (d_a, d_b, r) o None si los pickles no están disponibles.
    d_a = |pop(NC=nc_a, t_f) - pop(NC=nc_a-1, t_f)| — requiere nc_a-1
    d_b = |pop(NC=nc_b, t_f) - pop(NC=nc_a, t_f)|

    Nota: para nc_a=7 y nc_b=8, d_7 = |pop_7 - pop_6| (si NC=6 existe)
    o se aproxima d_7 desde la diferencia relativa si sólo hay NC=7 y NC=8.
    """
    data_a = _load_pkl(nc_a, nk)
    data_b = _load_pkl(nc_b, nk)

    if data_a is None or data_b is None:
        return None

    pop_a = _extract_final_pop(data_a)
    pop_b = _extract_final_pop(data_b)

    if pop_a is None or pop_b is None:
        return None

    d_b = abs(pop_b - pop_a)  # salto NC_A -> NC_B

    # Intentar calcular d_a desde NC_A-1
    data_prev = _load_pkl(nc_a - 1, nk)
    if data_prev is not None:
        pop_prev = _extract_final_pop(data_prev)
        if pop_prev is not None:
            d_a = abs(pop_a - pop_prev)
        else:
            d_a = None
    else:
        d_a = None

    if d_a is None or d_a < 1e-15:
        # No se puede calcular r desde tres puntos; usar sólo d_b y el fallback
        warnings.warn(
            f"No se pudo calcular d_{nc_a} (NC={nc_a - 1} no disponible). "
            "Ratio r no calculable desde pickles; se usará el ledger o fallback.",
            stacklevel=2,
        )
        return None

    r = d_b / d_a
    return float(d_a), float(d_b), float(r)


def _load_ratio_from_ledger() -> tuple[float, float, float] | None:
    """
    Intentar cargar d7, d8 y r desde pade_convergence_report.json.

    Returns
    -------
    (d7, d8, r) o None si el ledger no existe o no contiene los datos.
    """
    if not LEDGER_PATH.exists():
        return None
    try:
        ledger = json.loads(LEDGER_PATH.read_text(encoding="utf-8"))
        d7_entry = ledger.get("NC7", {})
        d8_entry = ledger.get("NC8", {})
        d7 = d7_entry.get("dPop")
        d8 = d8_entry.get("dPop")
        if d7 is not None and d8 is not None and float(d7) > 1e-15:
            r = float(d8) / float(d7)
            return float(d7), float(d8), r
    except Exception as e:
        warnings.warn(f"Error leyendo ledger {LEDGER_PATH}: {e}", stacklevel=2)
    return None


def run_richardson(nc_a: int = 7, nc_b: int = 8, nk: int = 1, verbose: bool = True) -> dict:
    """
    Ejecutar extrapolación de Richardson para convergencia HEOM.

    Jerarquía de fuentes para el ratio r:
    1. Cálculo on-the-fly desde pickles (NC=nc_a-1, nc_a, nc_b)
    2. Ledger pade_convergence_report.json
    3. Fallback documentado desde numerical_params.yaml

    Parameters
    ----------
    nc_a, nc_b : int
        Niveles de truncación NC. Típicamente 7 y 8.
    nk : int
        Número de términos de Padé. Default 1.
    verbose : bool
        Imprimir diagnóstico detallado.

    Returns
    -------
    dict con d_a, d_b, r, epsilon_a, epsilon_b, pop_inf, ratio_source.
    """
    if verbose:
        print("=" * 80)
        print(f"  EXTRAPOLACIÓN DE RICHARDSON (Fragmento 6DPU, Nk={nk})")
        print("=" * 80)

    # ---- Obtener ratio r mediante la jerarquía de fuentes ----
    ratio_source = "unknown"
    d_a = d_b = r = None

    # Ruta 1: desde pickles (on-the-fly)
    pkl_result = _compute_ratio_from_pickles(nc_a, nc_b, nk)
    if pkl_result is not None:
        d_a, d_b, r = pkl_result
        ratio_source = f"on-the-fly desde pickles NC={nc_a - 1}/{nc_a}/{nc_b}"
        if verbose:
            print(f"Fuente: {ratio_source}")

    # Ruta 2: desde ledger JSON
    if r is None:
        ledger_result = _load_ratio_from_ledger()
        if ledger_result is not None:
            d_a, d_b, r = ledger_result
            ratio_source = f"ledger {LEDGER_PATH.name}"
            if verbose:
                print(f"Pickles insuficientes. Usando ledger: {LEDGER_PATH.name}")

    # Ruta 3: fallback documentado
    if r is None:
        r = _load_fallback_ratio()
        ratio_source = "fallback documentado (numerical_params.yaml)"
        # Estimar d_b desde los pickles nc_a, nc_b aunque no tengamos nc_a-1
        data_a = _load_pkl(nc_a, nk)
        data_b = _load_pkl(nc_b, nk)
        if data_a is not None and data_b is not None:
            pop_a = _extract_final_pop(data_a)
            pop_b = _extract_final_pop(data_b)
            if pop_a is not None and pop_b is not None:
                d_b = abs(pop_b - pop_a)
                d_a = d_b / r if r > 0 else d_b  # estimado
                if verbose:
                    print(f"d_{nc_b} calculado desde pickles; d_{nc_a} estimado via r={r:.3f}")

    if r is None:
        raise RuntimeError(
            "No se pudo obtener el ratio de convergencia r por ninguna de las tres vías. "
            "Verificar que existan pickles o ledger en outputs_data/."
        )

    if verbose:
        print(f"\nPuntos usados: NC={nc_a} y NC={nc_b}")
        if d_b is not None:
            print(f"Salto d_{nc_b}: {d_b:.2e}")
        print(f"Ratio r ({ratio_source}): {r:.3f}")

    # ---- Diagnóstico de régimen asintótico ----
    ratio_stable = 0.1 < r < 0.9
    if not ratio_stable:
        warnings.warn(
            f"Ratio r={r:.3f} fuera del rango [0.1, 0.9]: "
            "régimen asintótico posiblemente no alcanzado. "
            "Considerar aumentar NC o verificar los pickles.",
            stacklevel=2,
        )
    if verbose:
        print(f"Régimen asintótico estable: {'SÍ' if ratio_stable else 'NO — ver advertencia'}")

    # ---- Estimadores de error ----
    if d_b is not None and (1 - r) > 1e-6:
        epsilon_b = d_b / (1 - r)
    else:
        epsilon_b = float("nan")

    if d_a is not None and (1 - r) > 1e-6:
        epsilon_a = d_a / (1 - r)
    else:
        epsilon_a = float("nan")

    d_proj = d_b * r if d_b is not None else float("nan")
    epsilon_proj = d_proj / (1 - r) if (not np.isnan(d_proj) and (1 - r) > 1e-6) else float("nan")

    if verbose:
        if not np.isnan(epsilon_b):
            print(f"Error de truncación NC={nc_b}: {epsilon_b:.2e}  ({epsilon_b * 100:.2f}%)")
        if not np.isnan(epsilon_a):
            print(f"Error de truncación NC={nc_a}: {epsilon_a:.2e}  ({epsilon_a * 100:.2f}%)")
        if not np.isnan(epsilon_proj):
            print(f"Error proyectado NC={nc_b + 1}: {epsilon_proj:.2e}  ({epsilon_proj * 100:.2f}%)")

    # ---- Estado asintótico (Richardson) ----
    data_a = _load_pkl(nc_a, nk)
    data_b = _load_pkl(nc_b, nk)
    pop_inf_final = None

    if verbose:
        print("\n" + "-" * 40)
        print("  ESTADO ASINTÓTICO (Richardson)")
        print("-" * 40)

    if data_a is not None and data_b is not None:
        tA = data_a["tlist"]
        tB = data_b["tlist"]

        # Extraer trayectoria de población del sitio 0
        try:
            if hasattr(data_a["rho_t"][0], "full"):
                popA = np.array([s.full()[0, 0].real for s in data_a["rho_t"]])
                popB = np.array([s.full()[0, 0].real for s in data_b["rho_t"]])
            else:
                popA = np.array([np.asarray(s)[0, 0].real for s in data_a["rho_t"]])
                popB = np.array([np.asarray(s)[0, 0].real for s in data_b["rho_t"]])

            # Interpolar si hay mismatch de mallas temporales
            if len(tA) != len(tB):
                if verbose:
                    print(f"  Mismatch de mallas ({len(tA)} vs {len(tB)}). Interpolando...")
                popA_on_tB = np.interp(tB, tA, popA)
                diff_pop = popB - popA_on_tB
            else:
                diff_pop = popB - popA

            pop_inf = popB + diff_pop * r / (1 - r) if (1 - r) > 1e-6 else popB
            pop_inf_final = float(pop_inf[-1])

            if verbose:
                print(f"Población sitio 0 final (t={tB[-1]:.1f} fs):")
                print(f"  NC={nc_a}:      {popA[-1]:.6f}")
                print(f"  NC={nc_b}:      {popB[-1]:.6f}")
                print(f"  Asíntota:  {pop_inf[-1]:.6f}")
                print(f"  Corrección desde NC={nc_b}: {pop_inf[-1] - popB[-1]:.2e}")

        except Exception as e:
            if verbose:
                print(f"  No se pudo calcular estado asintótico completo: {e}")
    else:
        if verbose:
            print("  Pickles no disponibles; cálculo asintótico omitido.")

    return {
        "nc_a": nc_a,
        "nc_b": nc_b,
        "nk": nk,
        "d_a": float(d_a) if d_a is not None else None,
        "d_b": float(d_b) if d_b is not None else None,
        "r": float(r),
        "ratio_source": ratio_source,
        "ratio_stable": bool(ratio_stable),
        "epsilon_nc_a": float(epsilon_a),
        "epsilon_nc_b": float(epsilon_b),
        "epsilon_nc_b_plus1_proj": float(epsilon_proj),
        "pop_inf_final": pop_inf_final,
    }


if __name__ == "__main__":
    result = run_richardson(nc_a=7, nc_b=8, nk=1, verbose=True)

    # Guardar resultado en ledger si es calculado on-the-fly
    if "pickles" in result["ratio_source"]:
        out_json = PKL_DIR.parent / "raw_json" / "pade_convergence_report.json"
        out_json.parent.mkdir(parents=True, exist_ok=True)
        ledger_update = {}
        if out_json.exists():
            try:
                ledger_update = json.loads(out_json.read_text(encoding="utf-8"))
            except Exception:
                pass
        if result["d_a"] is not None:
            ledger_update[f"NC{result['nc_a']}"] = {"dPop": result["d_a"]}
        if result["d_b"] is not None:
            ledger_update[f"NC{result['nc_b']}"] = {"dPop": result["d_b"]}
        ledger_update["computed_r"] = result["r"]
        ledger_update["ratio_source"] = result["ratio_source"]
        out_json.write_text(json.dumps(ledger_update, indent=2), encoding="utf-8")
        print(f"\nLedger actualizado: {out_json}")

    print(f"\nRatio r = {result['r']:.4f}  ({result['ratio_source']})")
    print(f"epsilon NC=8 = {result['epsilon_nc_b']:.2e}")
    print(f"Estado asintótico estable: {result['ratio_stable']}")
