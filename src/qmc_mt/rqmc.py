"""
Randomized QMC: R independent Owen-scrambled Sobol nets of size N.

  mu_hat = (1/R) sum_r mu_r,   mu_r = (1/N) sum_i f(x_{r,i})
  Var(mu_hat) unbiasedly estimated as s^2({mu_r}) / R.

Gives unbiased estimate + honest SE (not available in deterministic QMC),
while keeping O(N^{-3/2}) RMSE for smooth integrands in low d.

Ref: Owen, SIAM J. Numer. Anal. 34 (1997) 1884-1910.
"""
import numpy as np
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Callable
from scipy.stats import qmc

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))


@dataclass
class RQMCResult:
    estimate: float
    se: float
    n_per_replicate: int
    n_replicates: int
    replicate_means: np.ndarray


def rqmc_integrate(
    f: Callable[[np.ndarray], np.ndarray],
    d: int,
    n: int,
    n_replicates: int = 16,
    seed: int = 0,
) -> RQMCResult:
    seeds = np.random.SeedSequence(seed).spawn(n_replicates)
    m = int(np.ceil(np.log2(max(n, 1))))
    means = np.empty(n_replicates)
    for r in range(n_replicates):
        rng = np.random.default_rng(seeds[r])
        eng = qmc.Sobol(d=d, scramble=True, seed=rng)
        means[r] = float(np.mean(f(eng.random_base2(m=m)[:n])))
    return RQMCResult(
        estimate=float(means.mean()),
        se=float(means.std(ddof=1) / np.sqrt(n_replicates)),
        n_per_replicate=n,
        n_replicates=n_replicates,
        replicate_means=means,
    )


if __name__ == "__main__":
    import json
    # Integral of sin(pi*x)sin(pi*y) dx dy = (2/pi)^2 over [0,1]^2
    f = lambda x: np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1])
    truth = (2 / np.pi) ** 2
    out = []
    for n in (2**8, 2**10, 2**12, 2**14):
        r = rqmc_integrate(f, 2, n, n_replicates=32, seed=0)
        out.append({"n": n, "est": r.estimate, "se": r.se,
                    "truth": truth, "z": (r.estimate - truth) / r.se})
    print(json.dumps(out, indent=2))
