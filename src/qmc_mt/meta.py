"""
Per-study Bayesian evidence for microtubule quantum-coherence claims.

Design decisions (honesty constraints):
  1. No pooled log-ratio meta-analysis: only Babcock2024 reports an
     extractable log-ratio effect in the primary paper. A "meta-analysis"
     of k=1 is not a meta-analysis.
  2. Kalra2024 is reported on its native scale (raw seconds of LORR
     latency delay) with its own BF, not folded into a log-ratio pool.
  3. Bandyopadhyay2013 and Craddock2012 are mechanistic, not meta-analytic:
     we compute no BF for them.
  4. Each BF is computed by nested sampling (qmc_mt.nested_sampling)
     integrating a Gaussian likelihood over an explicit prior,
     against the point null θ=0 with its analytic marginal likelihood.

Output: one BF10 per usable study, plus a combined multiplicative BF
under strict independence assumption (reported with that caveat).
"""
from __future__ import annotations
import math
import numpy as np
from dataclasses import dataclass, asdict
import sys
from pathlib import Path

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raíz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.primary_data import StudyRecord, BABCOCK_2024, KALRA_2024
from qmc_mt.nested_sampling import nested_sample


@dataclass
class StudyEvidence:
    key: str
    scale: str
    effect: float
    se: float
    logZ_H1: float
    logZ_H1_err: float
    logZ_H0: float
    log_BF10: float
    BF10: float
    BF10_analytic: float        # Exact for Gaussian likelihood + Gaussian prior
    prior_sd: float
    n_live: int


def _logZ_H0(effect: float, se: float) -> float:
    """Marginal likelihood at θ=0 (point null), analytic Gaussian."""
    return -0.5 * (effect / se) ** 2 - 0.5 * math.log(2 * math.pi * se * se)


def _logZ_H1_analytic(effect: float, se: float, prior_sd: float) -> float:
    """Exact marginal likelihood for θ ~ N(0, σ_p²), y|θ ~ N(θ, s²)."""
    var_total = se**2 + prior_sd**2
    return -0.5 * (effect**2 / var_total) - 0.5 * math.log(2 * math.pi * var_total)


def _study_bf(
    study: StudyRecord,
    prior_sd: float,
    n_live: int = 600,
    seed: int = 0,
) -> StudyEvidence:
    if study.effect is None or study.se is None:
        raise ValueError(f"{study.key}: missing effect or SE; cannot compute BF.")

    y, s = float(study.effect), float(study.se)

    # H1: θ ~ N(0, prior_sd²); y | θ ~ N(θ, s²).
    def prior_transform(u: np.ndarray) -> np.ndarray:
        # inverse standard-normal CDF (Acklam approximation or simple proxy)
        return np.array([prior_sd * _norm_ppf(float(u[0]))])

    def loglike(theta: np.ndarray) -> float:
        t = float(theta[0])
        return -0.5 * ((y - t) / s) ** 2 - 0.5 * math.log(2 * math.pi * s * s)

    res = nested_sample(
        loglike=loglike,
        prior_transform=prior_transform,
        ndim=1,
        n_live=n_live,
        seed=seed,
    )

    logZ0 = _logZ_H0(y, s)
    log_bf = res.logZ - logZ0
    logZ1_analytic = _logZ_H1_analytic(y, s, prior_sd)
    
    return StudyEvidence(
        key=study.key,
        scale=study.scale,
        effect=y,
        se=s,
        logZ_H1=res.logZ,
        logZ_H1_err=res.logZ_err,
        logZ_H0=logZ0,
        log_BF10=log_bf,
        BF10=math.exp(log_bf),
        BF10_analytic=math.exp(logZ1_analytic - logZ0),
        prior_sd=prior_sd,
        n_live=n_live,
    )


def _norm_ppf(p: float) -> float:
    """Acklam inverse-normal CDF (double precision)."""
    if not (0.0 < p < 1.0):
        p = min(max(p, 1e-16), 1.0 - 1e-16)
    a = [-3.969683028665376e+01,  2.209460984245205e+02,
         -2.759285104469687e+02,  1.383577518672690e+02,
         -3.066479806614716e+01,  2.506628277459239e+00]
    b = [-5.447609879822406e+01,  1.615858368580409e+02,
         -1.556989798598866e+02,  6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
          4.374664141464968e+00,  2.938163982698783e+00]
    d = [ 7.784695709041462e-03,  3.224671290700398e-01,
          2.445134137142996e+00,  3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    if p <= phigh:
        q = p - 0.5
        r = q * q
        return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / \
               (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)
    q = math.sqrt(-2 * math.log(1 - p))
    return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
             ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)


def per_study_evidence(seed: int = 0) -> dict:
    """
    Per-study BF10 with scale-appropriate weakly-informative priors:
      - Babcock (log-ratio): θ ~ N(0, 1), covering ratios from ~0.14 to ~7.4.
      - Kalra (seconds):     θ ~ N(0, 120 s), covering ±2 min of latency shift.
    """
    babcock = _study_bf(BABCOCK_2024, prior_sd=1.0, seed=seed)
    kalra   = _study_bf(KALRA_2024,   prior_sd=120.0, seed=seed)

    log_bf_combined = babcock.log_BF10 + kalra.log_BF10

    return {
        "per_study": {
            babcock.key: asdict(babcock),
            kalra.key:   asdict(kalra),
        },
        "combined_under_independence": {
            "log_BF10": log_bf_combined,
            "BF10": math.exp(log_bf_combined),
            "caveat": ("Multiplicative combination assumes conditional "
                       "independence of the two studies' likelihoods given "
                       "the hypothesis."),
            "not_done": [
                "No pooled log-ratio estimator (only 1 study on log scale).",
                "No random-effects tau2 (undefined for k<2).",
                "No inclusion of Bandyopadhyay2013 or Craddock2012 (no numeric data)."
            ],
        }
    }


if __name__ == "__main__":
    import json
    print(json.dumps(per_study_evidence(seed=0), indent=2, default=float))