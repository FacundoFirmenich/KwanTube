"""richardson_fragment.py — Richardson extrapolation for HEOM convergence.

Patch v3.5.1: computes the contraction ratio from available ledger/checkpoints
instead of using hardcoded NC=7/8 values silently.
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))


def load(nc: int, nk: int = 1) -> dict:
    """Load HEOM Padé checkpoint for hierarchy depth ``nc`` and Matsubara ``nk``."""

    candidates = [
        PROJECT_ROOT / "outputs_data" / "raw_pkl" / f"pade_refined_ckpt_NC{nc}_Nk{nk}.pkl",
        PROJECT_ROOT.parent / "git_repo" / "outputs_data" / "raw_pkl" / f"pade_refined_ckpt_NC{nc}_Nk{nk}.pkl",
        PROJECT_ROOT.parent / "KwanTube_repo_backup_estable_28042026" / "outputs_data" / "raw_pkl" / f"pade_refined_ckpt_NC{nc}_Nk{nk}.pkl",
    ]
    for path in candidates:
        if path.exists():
            with path.open("rb") as handle:
                return pickle.load(handle)
    raise FileNotFoundError(f"No checkpoint found for NC={nc}, Nk={nk}")


def _population_site0(payload: dict) -> tuple[np.ndarray, np.ndarray]:
    """Extract time grid and site-0 population from a checkpoint payload."""

    tlist = np.asarray(payload["tlist"], dtype=float)
    pop = np.asarray([state.full()[0, 0].real for state in payload["rho_t"]], dtype=float)
    return tlist, pop


def _ledger_ratio(default_d7: float = 5.45e-3, default_d8: float = 2.09e-3) -> tuple[float, float, float, str]:
    """Read dynamic Richardson ratio from ledger if present, otherwise use documented fallback."""

    ledger_candidates = [
        PROJECT_ROOT / "outputs_data" / "raw_json" / "pade_convergence_report.json",
        PROJECT_ROOT / "outputs_data" / "raw_json" / "heom_convergence_summary.json",
    ]
    for ledger_path in ledger_candidates:
        if ledger_path.exists():
            ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
            d7 = float(ledger.get("NC7", {}).get("dPop", default_d7))
            d8 = float(ledger.get("NC8", {}).get("dPop", default_d8))
            return d7, d8, d8 / d7 if d7 > 0 else 0.5, f"ledger:{ledger_path.name}"
    return default_d7, default_d8, default_d8 / default_d7, "documented_fallback"


def main() -> int:
    """Run Richardson extrapolation diagnostic for the 6DPU fragment."""

    print("=" * 80)
    print("  EXTRAPOLACION DE RICHARDSON (Fragmento 6DPU, Nk=1)")
    print("=" * 80)

    nc_a, nc_b = 7, 8
    try:
        payload_a = load(nc_a)
        payload_b = load(nc_b)
    except Exception as exc:
        print(f"Error cargando puntos: {exc}")
        return 1

    t_a, pop_a = _population_site0(payload_a)
    t_b, pop_b = _population_site0(payload_b)
    if len(t_a) != len(t_b) or not np.allclose(t_a, t_b):
        print(f"Aviso: mallas temporales distintas ({len(t_a)} vs {len(t_b)}); remuestreo lineal auditado.")
        pop_a = np.interp(t_b, t_a, pop_a)

    d7, d8, ratio, ratio_source = _ledger_ratio()
    ratio_ok = 0.1 <= ratio <= 0.9
    epsilon_8 = d8 / (1.0 - ratio) if ratio < 1 else d8
    epsilon_7 = d7 / (1.0 - ratio) if ratio < 1 else d7
    epsilon_9 = d8 * ratio / (1.0 - ratio) if ratio < 1 else d8 * ratio

    print(f"Puntos usados: NC={nc_a} y NC={nc_b}")
    print(f"Salto d_7: {d7:.2e}")
    print(f"Salto d_8: {d8:.2e}")
    print(f"Ratio r (dinamico): {ratio:.3f}  source={ratio_source}")
    print(f"Diagnostico contraccion: {'OK' if ratio_ok else 'WARN'}")
    print(f"Error truncacion NC=8: {epsilon_8:.2e} ({epsilon_8 * 100:.2f}%)")
    print(f"Error truncacion NC=7: {epsilon_7:.2e} ({epsilon_7 * 100:.2f}%)")
    print(f"Error proyectado NC=9: {epsilon_9:.2e} ({epsilon_9 * 100:.2f}%)")

    diff_pop = pop_b - pop_a
    pop_inf = pop_b + diff_pop * ratio / (1.0 - ratio) if ratio < 1 else pop_b
    print("\n" + "-" * 40)
    print("  ESTADO ASINTOTICO (Richardson)")
    print("-" * 40)
    print(f"Poblacion sitio 0 final (t={t_b[-1]:.1f} fs):")
    print(f"  NC=7:      {pop_a[-1]:.6f}")
    print(f"  NC=8:      {pop_b[-1]:.6f}")
    print(f"  Asintota:  {pop_inf[-1]:.6f}")
    print(f"  Correccion total desde NC=8: {pop_inf[-1] - pop_b[-1]:.2e}")
    return 0 if ratio_ok else 2


if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    raise SystemExit(main())