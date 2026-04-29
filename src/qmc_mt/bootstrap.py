"""
Bootstrap confidence intervals.

- Percentile.
- BCa: bias-corrected and accelerated (Efron 1987). Second-order accurate,
  transformation-respecting, handles skewed sampling distributions.
  Acceleration via jackknife.

Validated by empirical coverage of the mean of a Normal sample against
the nominal 95% level.

Ref: Efron & Tibshirani, An Introduction to the Bootstrap (1993), Section14.
"""
from __future__ import annotations
import numpy as np
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Callable
from scipy.stats import norm

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))


@dataclass
class BootstrapResult:
    estimate: float
    se: float
    ci_percentile: tuple[float, float]
    ci_bca: tuple[float, float]
    z0: float
    acceleration: float
    n_boot: int
    level: float


def _jackknife_acceleration(x: np.ndarray, stat: Callable[[np.ndarray], float]) -> float:
    n = x.shape[0]
    jk = np.array([stat(np.delete(x, i, axis=0)) for i in range(n)])
    mean_jk = jk.mean()
    num = float(np.sum((mean_jk - jk) ** 3))
    den = 6.0 * float(np.sum((mean_jk - jk) ** 2)) ** 1.5
    return num / den if den > 0 else 0.0


def bootstrap_ci(
    x: np.ndarray,
    stat: Callable[[np.ndarray], float],
    n_boot: int = 4000,
    level: float = 0.95,
    seed: int = 0,
) -> BootstrapResult:
    x = np.asarray(x)
    n = x.shape[0]
    rng = np.random.default_rng(seed)

    theta_hat = float(stat(x))
    idx = rng.integers(0, n, size=(n_boot, n))
    boots = np.array([stat(x[i]) for i in idx])

    alpha = 1.0 - level
    lo_pct = float(np.quantile(boots, alpha / 2))
    hi_pct = float(np.quantile(boots, 1 - alpha / 2))

    prop_less = float(np.mean(boots < theta_hat))
    prop_less = min(max(prop_less, 1.0 / (n_boot + 1)), 1 - 1.0 / (n_boot + 1))
    z0 = float(norm.ppf(prop_less))
    a = _jackknife_acceleration(x, stat)
    zl, zu = norm.ppf(alpha / 2), norm.ppf(1 - alpha / 2)
    a1 = float(norm.cdf(z0 + (z0 + zl) / (1 - a * (z0 + zl))))
    a2 = float(norm.cdf(z0 + (z0 + zu) / (1 - a * (z0 + zu))))

    return BootstrapResult(
        estimate=theta_hat,
        se=float(boots.std(ddof=1)),
        ci_percentile=(lo_pct, hi_pct),
        ci_bca=(float(np.quantile(boots, a1)), float(np.quantile(boots, a2))),
        z0=z0, acceleration=a,
        n_boot=n_boot, level=level,
    )


if __name__ == "__main__":
    import json
    rng = np.random.default_rng(0)
    n_trials, n_sample = 10000, 40
    cov_pct = cov_bca = 0
    for t in range(n_trials):
        x = rng.normal(0.0, 1.0, size=n_sample)
        r = bootstrap_ci(x, np.mean, n_boot=1000, seed=t)
        cov_pct += int(r.ci_percentile[0] <= 0.0 <= r.ci_percentile[1])
        cov_bca += int(r.ci_bca[0] <= 0.0 <= r.ci_bca[1])
    print(json.dumps({
        "n_trials": n_trials, "n_sample": n_sample, "nominal": 0.95,
        "coverage_percentile": cov_pct / n_trials,
        "coverage_bca": cov_bca / n_trials,
    }, indent=2))
