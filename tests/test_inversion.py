
import unittest
import sys
from pathlib import Path
import numpy as np

# Add src to path
src_path = str(Path(__file__).parents[1] / "src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from qmc_mt import MultiTempInversionEngine

class TestInversion(unittest.TestCase):
    def setUp(self):
        self.engine = MultiTempInversionEngine()
        self.truth = [0.42, 6.2e12, 0.155]

    def test_synthetic_recovery(self):
        data = self.engine.forward(self.truth[0], self.truth[1], self.truth[2])
        res = self.engine.invert(data)
        
        eta_fit, log_wc_fit, gap_fit = res.x
        wc_fit = 10**log_wc_fit
        
        self.assertAlmostEqual(eta_fit, self.truth[0], places=4)
        self.assertAlmostEqual(wc_fit / 1e12, self.truth[1] / 1e12, places=4)
        self.assertAlmostEqual(gap_fit, self.truth[2], places=5)
        self.assertTrue(res.success)

if __name__ == "__main__":
    unittest.main()
