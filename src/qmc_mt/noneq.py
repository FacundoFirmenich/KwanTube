"""
Non-equilibrium quantum coherence analysis for microtubules.
Version: 3.3.1 (Audited Rigor - PNAS Compliance)
"""

import numpy as np
import sys
from pathlib import Path

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.core import (TubulinDimer, ExperimentalParameters, const)

class FrohlichCondensation:
    def __init__(self, dimer: TubulinDimer, params: ExperimentalParameters):
        self.dimer = dimer; self.params = params
        self.E_GTP = 0.42 * const.eV_to_J; self.P_single = self.E_GTP * 1.0

    def pumping_parameter(self, N: int = 1, Gamma_coll: float = 1e6,
                          beta: float = 0.0, N_ref: int = 1) -> float:
        N_array = np.asarray(N, dtype=float)
        N_safe = np.maximum(N_array, 1.0)
        gamma_eff = Gamma_coll * (N_safe / max(float(N_ref), 1.0))**beta
        eta = N_array * self.P_single / (const.kB * self.params.temperature * gamma_eff)
        if np.ndim(eta) == 0:
            return float(eta)
        return eta

    def critical_N(self, Gamma_coll: float = 1e6, beta: float = 0.0,
                   N_ref: int = 1):
        if beta >= 1.0:
            return None
        prefactor = const.kB * self.params.temperature * Gamma_coll / (self.P_single * max(N_ref, 1)**beta)
        exponent = 1.0 / (1.0 - beta)
        return int(np.ceil(prefactor**exponent))

class QEDCavityModel:
    def __init__(self, L_MT=25e-6, R_inner=7.5e-9, epsilon=80.0,
                 d_transition=18.7, lambda_probe_nm=280.0):
        self.L = L_MT; self.V = np.pi * R_inner**2 * L_MT
        self.N_dimers = 13 * int(L_MT / 8e-9)
        self.d_trans = d_transition * const.Debye_to_Cm
        self.omega_c_rad = (6e12) * 2 * np.pi; self.epsilon = epsilon
        self.lambda_probe_nm = lambda_probe_nm

    def lambda_0(self) -> float:
        E_vac = np.sqrt(2 * np.pi * const.hbar * self.omega_c_rad / (self.epsilon * const.epsilon_0 * self.V))
        return self.d_trans * E_vac / (const.hbar * 2 * np.pi)

    def collective_rabi_Hz(self) -> float:
        return 2 * self.lambda_0() * np.sqrt(self.N_dimers)

    def splitting_eV(self) -> float:
        return self.collective_rabi_Hz() * const.h / const.eV_to_J

    def splitting_nm_at_280(self) -> float:
        lambda_probe_m = self.lambda_probe_nm * 1e-9
        return lambda_probe_m**2 * self.collective_rabi_Hz() / const.c * 1e9

    def tau_cavity(self) -> float:
        return 5e-7 * (self.epsilon / 80.0)**2

if __name__ == "__main__":
    print("[noneq] Realizando subprocesos de fondo...")
