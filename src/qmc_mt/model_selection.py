"""
Statistical model selection for QED cavity detection.
Calculates BIC (Bayesian Information Criterion) for doublet vs singlet models.
"""
import numpy as np
import sys
from pathlib import Path
from scipy.optimize import curve_fit

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raíz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.noneq import QEDCavityModel

def lorentzian(E, E0, gamma):
    return (gamma/2)**2 / ((E - E0)**2 + (gamma/2)**2)

def bic_analysis(n_realizations: int = 25, rng_seed: int = 42, effective_points: int = 20):
    """Doublet (Cavity) vs singlet model selection.

    The synthetic spectrum is sampled on a fine grid for numerical stability,
    but the BIC penalty is evaluated using an ``effective_points`` count to
    reflect the number of statistically independent spectral resolution
    elements rather than the raw interpolation grid.
    """
    # Data from QEDCavityModel (v3.3)
    cav = QEDCavityModel(L_MT=25e-6)
    delta_E = cav.splitting_eV() # ~14.15 meV (0.01415 eV)
    rng = np.random.default_rng(rng_seed)
    
    # Spectroscopy params
    E_center = 4.43 # eV (~280 nm)
    sigma = 0.08    # eV (UV broadening)
    
    # Simulation range (eV)
    x = np.linspace(E_center - 0.5, E_center + 0.5, 1000)
    
    # Models
    def model_singlet(x, A, E0, S):
        return A * lorentzian(x, E0, S)
    
    def model_doublet(x, A, E0, S, dE):
        return A * (lorentzian(x, E0 - dE/2, S) + lorentzian(x, E0 + dE/2, S))

    # True data = Doublet + noise
    y_true = model_doublet(x, 1.0, E_center, sigma, delta_E)
    
    snr_levels = np.logspace(1, 4, 20)
    delta_bic = []
    
    for snr in snr_levels:
        noise = 1.0 / snr
        bic_trials = []
        for _ in range(n_realizations):
            y_obs = y_true + rng.normal(0, noise, len(x))

            p0_s = [1.0, E_center, sigma]
            p0_d = [1.0, E_center, sigma, delta_E]

            popt_s, _ = curve_fit(
                model_singlet, x, y_obs, p0=p0_s,
                bounds=([0.0, E_center - 0.1, 0.01], [5.0, E_center + 0.1, 0.2]),
                maxfev=20000,
            )
            popt_d, _ = curve_fit(
                model_doublet, x, y_obs, p0=p0_d,
                bounds=([0.0, E_center - 0.1, 0.01, 0.0], [5.0, E_center + 0.1, 0.2, 0.05]),
                maxfev=20000,
            )

            rss_s = np.sum((y_obs - model_singlet(x, *popt_s))**2)
            rss_d = np.sum((y_obs - model_doublet(x, *popt_d))**2)

            n_eff = effective_points
            rss_s_eff = rss_s * (n_eff / len(x))
            rss_d_eff = rss_d * (n_eff / len(x))
            bic_s = n_eff * np.log(rss_s_eff / n_eff) + 3 * np.log(n_eff)
            bic_d = n_eff * np.log(rss_d_eff / n_eff) + 4 * np.log(n_eff)
            bic_trials.append(bic_s - bic_d)

        delta_bic.append(np.mean(bic_trials))
        
    return snr_levels, np.array(delta_bic)

if __name__ == "__main__":
    snr, dbic = bic_analysis()
    print(f"Max BIC Difference: {dbic[-1]:.2e}")
