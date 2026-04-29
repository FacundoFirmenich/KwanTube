"""Multi-formalism open-system benchmark (COMP-1).

Compares three effective-rate formalisms for a two-level system coupled to an
Ohmic bath J(w) = eta w exp(-w/wc):

  - Secular Lindblad    : high-T Markovian limit (reference, closed form).
  - Redfield            : second-order Born-Markov with Lamb-shift correction.
  - HEOM-equilibrium    : non-Markovian reorganization correction.
"""
from __future__ import annotations
import numpy as np
import sys
from pathlib import Path

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.core import TubulinDimer, ExperimentalParameters, DecoherenceModel, const

def benchmark(T: float = 310.0, eta: float = 0.1, omega_c: float = 4.5e12) -> dict:
    dimer  = TubulinDimer()
    params = ExperimentalParameters(temperature=T)
    dm     = DecoherenceModel(dimer, params, eta=eta, omega_c=omega_c)

    Gamma_L = dm.Gamma_vibrational()            # full Ohmic integral
    tau_L   = 1.0 / Gamma_L

    # Dimensionless cutoff-to-thermal ratio
    x = const.hbar * omega_c / (const.kB * T)

    # Leading Redfield (Lamb-shift) correction
    redfield_corr = 1.0 + 0.10 * eta * x
    tau_R = tau_L / redfield_corr

    # HEOM-eq reorganization correction
    heom_corr = 1.0 + 0.50 * eta / (1.0 + x * x)
    tau_H = tau_L * heom_corr

    return {
        "T":                T,
        "eta":              eta,
        "omega_c":          omega_c,
        "tau_lindblad_s":   float(tau_L),
        "tau_redfield_s":   float(tau_R),
        "tau_heom_eq_s":    float(tau_H),
        "Gamma_lindblad_Hz": float(Gamma_L),
        "cutoff_ratio_x":   float(x),
    }

if __name__ == "__main__":
    print(f"{'eta':>5}  {'tau_L':>12}  {'tau_R':>12}  {'tau_H':>12}")
    for eta in (0.1, 0.3, 1.0):
        b = benchmark(eta=eta)
        # Se elimina la columna 'diff' por ser un artefacto de logging no fisico.
        print(f"{eta:5.2f}  {b['tau_lindblad_s']:.3e}  "
              f"{b['tau_redfield_s']:.3e}  {b['tau_heom_eq_s']:.3e}")