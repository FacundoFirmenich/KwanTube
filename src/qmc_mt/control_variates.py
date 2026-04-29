"""
Optimal linear control variates for variance reduction.

  mu_cv = E[f] - beta_hat * (E[g] - m_g),    beta_hat = Cov_hat(f,g) Var_hat(g)^{-1}
  Var(mu_cv) / Var(mu) = 1 - rho^2   (multi-dim: 1 - R^2)

Plug-in beta_hat from the same sample; bias O(1/N), variance
reduction factor reported exactly on sample.

Ref: Glasserman, Monte Carlo Methods in Financial Engineering (2004), Section4.1.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass


@dataclass
class CVResult:
    estimate: float
    se: float
    estimate_plain: float
    se_plain: float
    variance_reduction: float
    beta: np.ndarray
    n: int


def control_variates(
    f_samples: np.ndarray,
    g_samples: np.ndarray,          # (N,) or (N,k)
    g_means:   np.ndarray,          # (k,) analytically known
) -> CVResult:
    f = np.asarray(f_samples).ravel()
    G = np.atleast_2d(np.asarray(g_samples))
    if G.shape[0] != f.size:
        G = G.T
    g_means = np.atleast_1d(np.asarray(g_means).ravel())
    N = f.size

    mu_f = float(f.mean())
    var_f = float(f.var(ddof=1))
    Gc = G - G.mean(axis=0, keepdims=True)
    fc = f - mu_f
    Cgg = (Gc.T @ Gc) / (N - 1)
    Cgf = (Gc.T @ fc) / (N - 1)
    beta = np.linalg.solve(Cgg, Cgf)

    mu_cv = mu_f - float(beta @ (G.mean(axis=0) - g_means))
    var_cv = float((fc - Gc @ beta).var(ddof=1))

    return CVResult(
        estimate=mu_cv,
        se=float(np.sqrt(var_cv / N)),
        estimate_plain=mu_f,
        se_plain=float(np.sqrt(var_f / N)),
        variance_reduction=float(1 - var_cv / var_f) if var_f > 0 else 0.0,
        beta=beta,
        n=N,
    )


if __name__ == "__main__":
    import json
    # f(x)=e^x on U[0,1], CV g(x)=x with known mean 1/2.
    rng = np.random.default_rng(0)
    truth = float(np.e - 1)
    out = []
    for N in (200, 2000, 20000):
        x = rng.uniform(size=N)
        r = control_variates(np.exp(x), x.reshape(-1, 1), np.array([0.5]))
        out.append({
            "N": N, "truth": truth,
            "plain": r.estimate_plain, "se_plain": r.se_plain,
            "cv": r.estimate, "se_cv": r.se,
            "var_reduction_pct": 100 * r.variance_reduction,
        })
    print(json.dumps(out, indent=2))
