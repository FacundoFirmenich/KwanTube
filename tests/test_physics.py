
import unittest
import sys
from pathlib import Path

# Add src to path
src_path = str(Path(__file__).parents[1] / "src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from qmc_mt import TubulinDimer, ExperimentalParameters, DecoherenceModel, EquilibriumBound

class TestPhysics(unittest.TestCase):
    def setUp(self):
        self.dimer = TubulinDimer()
        self.params = ExperimentalParameters(temperature=310.0)
        self.model = DecoherenceModel(self.dimer, self.params, eta=0.3, omega_c=4.5e12)

    def test_positive_rates(self):
        rates = self.model.get_all_rates()
        for name, rate in rates.items():
            if isinstance(rate, (int, float)):
                self.assertGreater(rate, 0, f"Rate {name} should be positive.")

    def test_equilibrium_bound(self):
        bound = EquilibriumBound.tau_max(310.0, 0.3)
        rates = self.model.get_all_rates()
        self.assertLessEqual(1.0/rates['total_dephasing'], bound, 
                             "Coherence time should respect equilibrium bound.")

    def test_high_temp_scaling(self):
        # Scale T to 620 K (twice)
        p2 = ExperimentalParameters(temperature=620.0)
        m2 = DecoherenceModel(self.dimer, p2, eta=0.3, omega_c=4.5e12)
        r1 = self.model.get_all_rates()['vibrational_dephasing']
        r2 = m2.get_all_rates()['vibrational_dephasing']
        # T is doubled, so Gamma should roughly double (scaling eta*T)
        self.assertAlmostEqual(r2/r1, 2.0, delta=0.1)

    def test_vibrational_numeric_matches_analytic_scale(self):
        rates = self.model.get_all_rates()
        ratio = rates['vibrational_dephasing'] / rates['vibrational_analytic']
        self.assertAlmostEqual(ratio, 1.0, delta=0.05)

    def test_ionic_mass_parameter_changes_rate(self):
        na_model = DecoherenceModel(self.dimer, self.params, eta=0.3, omega_c=4.5e12, ion_mass_amu=23.0)
        k_model = DecoherenceModel(self.dimer, self.params, eta=0.3, omega_c=4.5e12, ion_mass_amu=39.0)
        self.assertNotAlmostEqual(na_model.Gamma_ions(), k_model.Gamma_ions(), places=6)

if __name__ == "__main__":
    unittest.main()
