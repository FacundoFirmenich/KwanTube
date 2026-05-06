"""
test_ns_consistency.py -- Statistical stress test for Nested Sampling engine.
Compares NS results against analytic Gaussian-Gaussian solutions across 20 seeds.
"""
import numpy as np
import math
import sys
from pathlib import Path
import unittest

src_path = str(Path(__file__).resolve().parent.parent)
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from qmc_mt.meta import _study_bf
from qmc_mt.primary_data import BABCOCK_2024, KALRA_2024

def run_stress_test(study, prior_sd, n_seeds=21):
    print(f"\n--- Stress Test: {study.key} (prior_sd={prior_sd}) ---")
    deviations = []
    
    for seed in range(n_seeds):
        res = _study_bf(study, prior_sd=prior_sd, n_live=400, seed=seed)
        
        log_bf_ns = res.log_BF10
        log_bf_an = math.log(res.BF10_analytic)
        
        # Standardized deviation: abs(diff) / logZ_err
        # logZ_err is the estimated standard deviation of logZ (and thus of logBF since logZ0 is fixed)
        sigma_dev = abs(log_bf_ns - log_bf_an) / res.logZ_H1_err
        deviations.append(sigma_dev)
        
        status = "OK" if sigma_dev < 3 else "FAIL (>3 sigma)"
        print(f"  Seed {seed:2d}: NS={res.BF10:7.1f} | AN={res.BF10_analytic:7.1f} | Dev={sigma_dev:5.2f} sigma | {status}")

    avg_dev = np.mean(deviations)
    max_dev = np.max(deviations)
    success_rate = np.mean(np.array(deviations) < 3.0) * 100
    
    print(f"\nSummary for {study.key}:")
    print(f"  Success Rate (<3 sigma): {success_rate:3.0f}%")
    print(f"  Mean deviation:     {avg_dev:5.2f} sigma")
    print(f"  Max deviation:      {max_dev:5.2f} sigma")
    
    return success_rate >= 90.0


class TestNestedSamplingConsistency(unittest.TestCase):
    """CI-friendly NS sanity checks against analytic BF on a reduced seed grid."""

    def test_babcock_ns_consistency(self):
        self.assertTrue(run_stress_test(BABCOCK_2024, prior_sd=1.0, n_seeds=5))

    def test_kalra_ns_consistency(self):
        self.assertTrue(run_stress_test(KALRA_2024, prior_sd=120.0, n_seeds=5))

if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    # Optional verbose stress mode when run as script.
    babcock_ok = run_stress_test(BABCOCK_2024, prior_sd=1.0)
    kalra_ok = run_stress_test(KALRA_2024, prior_sd=120.0)
    
    if babcock_ok and kalra_ok:
        print("\n[VERDICT] NS vs Analytic consistency: PASSED.")
    else:
        print("\n[VERDICT] NS vs Analytic consistency: FAILED.")
        exit(1)
