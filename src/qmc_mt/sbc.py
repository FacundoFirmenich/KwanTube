"""
Simulation-Based Calibration (Talts, Betancourt, Simpson, Vehtari, Gelman 2018).

For a correctly-implemented Bayesian sampler:
  θ̃ ~ π(θ);  ỹ ~ p(y|θ̃);  θ*_{1..L} ~ p(θ|ỹ)
  rank(θ̃) = #{l : θ*_l < θ̃}  ~  Uniform{0,1,...,L}

Deviation from discrete uniformity → miscalibration. We report
rank histogram and χ² uniformity test.

Validated on conjugate Normal-Normal (should pass) and on a
deliberately-too-narrow posterior (should be flagged).

Ref: arXiv:1804.06788.
"""
import numpy as np
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Callable
from scipy.stats import chi2 as chi2_dist

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raíz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))


@dataclass
class SBCResult:
    ranks: np.ndarray
    L: int
    n_sim: int
    chi2: float
    dof: int
    p_value: float
    bin_counts: np.ndarray
    bin_edges: np.ndarray


def simulation_based_calibration(
    prior_sampler:     Callable[[np.random.Generator], float],
    data_sampler:      Callable[[float, np.random.Generator], np.ndarray],
    posterior_sampler: Callable[[np.ndarray, int, np.random.Generator], np.ndarray],
    n_sim: int = 1000,
    L: int = 99,
    n_bins: int = 20,
    seed: int = 0,
) -> SBCResult:
    rng = np.random.default_rng(seed)
    ranks = np.empty(n_sim, dtype=int)
    for i in range(n_sim):
        theta_true = float(prior_sampler(rng))
        y = data_sampler(theta_true, rng)
        post = np.asarray(posterior_sampler(y, L, rng))
        ranks[i] = int(np.sum(post < theta_true))

    edges = np.linspace(-0.5, L + 0.5, n_bins + 1)
    counts, _ = np.histogram(ranks, bins=edges)
    expected = n_sim / n_bins
    chi2 = float(((counts - expected) ** 2 / expected).sum())
    dof = n_bins - 1
    return SBCResult(
        ranks=ranks, L=L, n_sim=n_sim,
        chi2=chi2, dof=dof,
        p_value=float(chi2_dist.sf(chi2, dof)),
        bin_counts=counts, bin_edges=edges,
    )


if __name__ == "__main__":
    import json
    tau2, sigma2 = 4.0, 1.0
    def prior(rng): return float(rng.normal(0.0, np.sqrt(tau2)))
    def data(theta, rng): return rng.normal(theta, np.sqrt(sigma2), size=1)
    def posterior_correct(y, L, rng):
        pm = y[0] * tau2 / (sigma2 + tau2)
        pv = sigma2 * tau2 / (sigma2 + tau2)
        return rng.normal(pm, np.sqrt(pv), size=L)
    def posterior_bad(y, L, rng):
        pm = y[0] * tau2 / (sigma2 + tau2)
        pv = 0.25 * sigma2 * tau2 / (sigma2 + tau2)    # 4× too narrow
        return rng.normal(pm, np.sqrt(pv), size=L)

    good = simulation_based_calibration(prior, data, posterior_correct,
                                        n_sim=1000, L=99, n_bins=20, seed=0)
    bad  = simulation_based_calibration(prior, data, posterior_bad,
                                        n_sim=1000, L=99, n_bins=20, seed=1)
    print(json.dumps({
        "conjugate_correct": {
            "chi2": good.chi2, "dof": good.dof, "p_value": good.p_value,
            "calibrated_at_0.05": good.p_value > 0.05,
        },
        "conjugate_too_narrow": {
            "chi2": bad.chi2, "dof": bad.dof, "p_value": bad.p_value,
            "flagged_miscalibrated_at_0.05": bad.p_value < 0.05,
        },
    }, indent=2))
