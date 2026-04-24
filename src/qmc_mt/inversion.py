
"""
Multi-objective Parameter Inversion Engine for MT Coherence.
Enhanced Version (Phase P2.1).

Uses multi-temperature linewidths to break the eta - omega_c degeneracy.
"""

import numpy as np
import sys
from pathlib import Path
from scipy.optimize import least_squares

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raíz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.core import const, TubulinDimer, ExperimentalParameters, DecoherenceModel
import json

class MultiTempInversionEngine:
    def __init__(self, temps=[280, 295, 310]):
        self.temps = temps
        self.dimer_ref = TubulinDimer()
        self.params_lib = [ExperimentalParameters(temperature=T) for T in temps]
        self.default_scales = None
        
    def forward(self, eta, wc, gap_eV):
        """Returns [Gap_310, FWHM_280, FWHM_295, FWHM_310] in meV."""
        results = [gap_eV * 1e3]
        for p in self.params_lib:
            m = DecoherenceModel(self.dimer_ref, p, eta=eta, omega_c=wc)
            r = m.get_all_rates()
            fwhm = const.hbar * r['total_dephasing'] * 1e3 / const.e
            results.append(fwhm)
        return np.array(results)

    def misfit(self, p, target):
        eta, log_wc, gap = p
        pred = self.forward(eta, 10**log_wc, gap)
        scales = self.default_scales
        if scales is None:
            scales = np.maximum(np.abs(target), np.array([1.0, 1e-3, 1e-3, 1e-3]))
        return (pred - target) / scales

    def invert(self, data_vector):
        bounds = ([0.05, 11, 0.01], [1.5, 14, 1.0])
        self.default_scales = np.maximum(np.abs(data_vector), np.array([1.0, 1e-3, 1e-3, 1e-3]))
        initial_guesses = [
            [0.3, 12.6, 0.15],
            [0.15, 12.3, 0.14],
            [0.6, 12.8, 0.18],
            [0.9, 13.0, 0.12],
        ]
        best = None
        for p0 in initial_guesses:
            res = least_squares(self.misfit, p0, bounds=bounds, args=(data_vector,))
            if best is None or res.cost < best.cost:
                best = res
        return best

def run_multi_test():
    engine = MultiTempInversionEngine()
    truth = [0.42, 6.2e12, 0.155]
    print(f"TRUTH (Multi-T): eta={truth[0]}, wc={truth[1]:.2e}")
    
    data = engine.forward(truth[0], truth[1], truth[2])
    print(f"Data Vector (Center, FWHM @ {engine.temps} K): {data}")
    
    res = engine.invert(data)
    eta_fit, log_wc_fit, gap_fit = res.x
    wc_fit = 10**log_wc_fit
    
    print(f"\nRECOVERED (Multi-T):")
    print(f"  eta: {eta_fit:.4f} (Err: {abs(eta_fit - truth[0]):.4e})")
    print(f"  wc: {wc_fit:.2e} (Err: {abs(wc_fit - truth[1]):.4e})")
    print(f"  Success: {res.success}")

if __name__ == "__main__":
    run_multi_test()
