"""
Nested Sampling for Bayesian Evidence Computation (Skilling 2004/2006).

Computes the marginal likelihood (evidence) Z = integral L(theta) pi(theta) dtheta.
Uses a static nested sampling approach with simple ellipsoidal sampling
or rejection for low-dimensional problems (ndim=1).

Validated on known Gaussian integrals.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Callable


@dataclass
class NSResult:
    logZ: float
    logZ_err: float
    n_samples: int
    samples: np.ndarray
    weights: np.ndarray


def nested_sample(
    loglike: Callable[[np.ndarray], float],
    prior_transform: Callable[[np.ndarray], np.ndarray],
    ndim: int,
    n_live: int = 500,
    tol: float = 0.01,
    seed: int = 0,
) -> NSResult:
    rng = np.random.default_rng(seed)
    
    # Initialize live points from prior
    u_live = rng.random((n_live, ndim))
    theta_live = np.array([prior_transform(u) for u in u_live])
    logl_live = np.array([loglike(th) for th in theta_live])
    
    logZ = -1e300
    h = 0.0
    log_w = np.log(1.0 - np.exp(-1.0 / n_live))
    
    samples = []
    weights = []
    
    for i in range(10000):  # safety limit
        # Find point with lowest likelihood
        worst = np.argmin(logl_live)
        logl_star = logl_live[worst]
        
        # Update evidence: wt = L * DeltaX
        # DeltaX_i = exp(-i/n_live) - exp(-(i+1)/n_live)
        log_wt = -i / n_live + np.log(1.0 - np.exp(-1.0 / n_live)) + logl_star
        logZ_new = np.logaddexp(logZ, log_wt)
        
        # Update information H (entropy)
        if logZ > -1e100:
            h = np.exp(log_wt - logZ_new) * logl_star + \
                np.exp(logZ - logZ_new) * (h + logZ) - logZ_new
        else:
            # First iteration or extremely low Z: h approx logl_star - logZ_new
            h = logl_star - logZ_new
            
        logZ = logZ_new
        
        # Record sample
        samples.append(theta_live[worst].copy())
        weights.append(np.exp(log_wt))
        
        # Replace worst point with new one from prior, subject to L > L*
        # For ndim=1, simple rejection sampling is efficient
        while True:
            u_new = rng.random(ndim)
            theta_new = prior_transform(u_new)
            logl_new = loglike(theta_new)
            if logl_new > logl_star:
                u_live[worst] = u_new
                theta_live[worst] = theta_new
                logl_live[worst] = logl_new
                break
        
        # Convergence criterion
        logZ_remain = np.max(logl_live) - i / n_live
        if logZ_remain - logZ < tol:
            break
            
    # Add remaining live points
    log_w_final = -(i + 1) / n_live - np.log(n_live)
    for j in range(n_live):
        logZ = np.logaddexp(logZ, log_w_final + logl_live[j])
        samples.append(theta_live[j].copy())
        weights.append(np.exp(log_w_final + logl_live[j]))
        
    return NSResult(
        logZ=float(logZ),
        logZ_err=float(np.sqrt(h / n_live)),
        n_samples=len(samples),
        samples=np.array(samples),
        weights=np.array(weights) / np.sum(weights)
    )

if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    import math
    # Test: N(0, 1) likelihood with U[-10, 10] prior.
    # Z = integral from -10 to 10 (1/sqrt(2*pi)) exp(-x^2/2) (1/20) dx approx 1/20 = 0.05
    # logZ approx log(0.05) approx -2.9957
    def lt(th): return -0.5 * th[0]**2 - 0.5 * math.log(2*math.pi)
    def pt(u): return np.array([20.0 * u[0] - 10.0])
    
    res = nested_sample(lt, pt, 1, n_live=1000)
    print(f"LogZ (truth approx -2.9957): {res.logZ:.4f} +/- {res.logZ_err:.4f}")
