"""
Quasi-Monte Carlo integration with Sobol and Halton low-discrepancy sequences.

- Sobol via scipy (Joe-Kuo direction numbers, d up to 21201).
- Halton via scipy (van der Corput in successive primes).
- Owen-scrambled by default for unbiasedness.
- Error O(N^{-1}(log N)^d) for smooth f vs MC's O(N^{-1/2}); empirically
  validated with log-log slope fit on closed-form integrals.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Callable, Literal
from scipy.stats import qmc


@dataclass
class QMCResult:
    estimate: float
    n_samples: int
    sequence: str
    d: int


def sobol_points(n: int, d: int, seed: int | None = 0, scramble: bool = True) -> np.ndarray:
    eng = qmc.Sobol(d=d, scramble=scramble, seed=seed)
    m = int(np.ceil(np.log2(max(n, 1))))
    return eng.random_base2(m=m)[:n]


def halton_points(n: int, d: int, seed: int | None = 0, scramble: bool = True) -> np.ndarray:
    eng = qmc.Halton(d=d, scramble=scramble, seed=seed)
    return eng.random(n)


def qmc_integrate(
    f: Callable[[np.ndarray], np.ndarray],
    d: int,
    n: int,
    sequence: Literal["sobol", "halton"] = "sobol",
    seed: int = 0,
) -> QMCResult:
    pts = (sobol_points(n, d, seed=seed) if sequence == "sobol"
           else halton_points(n, d, seed=seed))
    return QMCResult(estimate=float(np.mean(f(pts))),
                     n_samples=n, sequence=sequence, d=d)


def convergence_study(
    f: Callable[[np.ndarray], np.ndarray],
    truth: float,
    d: int,
    ns: list[int],
    sequence: str = "sobol",
    seed: int = 0,
) -> dict:
    errs = [abs(qmc_integrate(f, d, n, sequence=sequence, seed=seed).estimate - truth)
            for n in ns]
    slope = float(np.polyfit(np.log(ns), np.log(np.maximum(errs, 1e-16)), 1)[0])
    return {"ns": list(ns), "errors": errs, "loglog_slope": slope, "sequence": sequence}


if __name__ == "__main__":
    import json
    out = {}
    # ∫_{[0,1]^d} Σ x_i = d/2
    for d in (2, 5, 10):
        f = (lambda x, d=d: x.sum(axis=1))
        out[f"sum_d{d}"] = convergence_study(f, d/2, d, [2**k for k in range(6, 14)])
    # ∫_{[0,1]} e^x = e - 1
    for seq in ("sobol", "halton"):
        out[f"exp_{seq}"] = convergence_study(
            lambda x: np.exp(x[:, 0]), float(np.e - 1), 1,
            [2**k for k in range(6, 14)], sequence=seq)
    # Sanity: slope should be close to -1 for QMC, not -0.5
    print(json.dumps({k: {"slope": v["loglog_slope"],
                          "final_err": v["errors"][-1]}
                      for k, v in out.items()}, indent=2))
