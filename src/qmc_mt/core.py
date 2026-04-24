r"""Core quantum microtubule coherence module.

This module contains the minimal open-system ingredients used throughout the
microtubule framework. The present implementation is still intentionally light,
but now makes explicit several assumptions that were previously only implicit:

- the effective open-system treatment is a phenomenological secular/Markovian
  rate model used for calibrated order-of-magnitude analysis;
- ionic decoherence is parameterized by an effective ion mass, defaulting to
  Na+ but exposed explicitly for sensitivity and physiology-aware studies;
- the vibrational Ohmic integral is regularized analytically in the
  :math:`\omega \to 0` limit instead of silently nulling the integrand below an
  arbitrary cutoff.
"""

import numpy as np
from scipy.integrate import quad


class PhysicalConstants:
    def __init__(self):
        self.hbar = 1.054571817e-34
        self.h = 6.62607015e-34
        self.kB = 1.380649e-23
        self.c = 2.99792458e8
        self.e = 1.602176634e-19
        self.epsilon_0 = 8.8541878128e-12
        self.eV_to_J = 1.602176634e-19
        self.Debye_to_Cm = 3.33564e-30
        self.amu_to_kg = 1.66053906660e-27
        self.cm1_to_Hz = 2.99792458e10
        self.NA = 6.02214076e23


const = PhysicalConstants()


class TubulinDimer:
    def __init__(self, energy_gap=0.15, tunneling=0.01):
        self.energy_gap = energy_gap
        self.tunneling = tunneling
        self.dipole_SI = 4000.0 * const.Debye_to_Cm
        self.omega_0 = energy_gap * const.eV_to_J / const.hbar

    def transition_dipole(self, dipole_type='optical'):
        d = 5.0 if dipole_type == 'optical' else 18.7
        return d * const.Debye_to_Cm


class ExperimentalParameters:
    def __init__(self, temperature=310.0, ionic_strength=0.15, dielectric=80.0):
        self.temperature = temperature
        self.dielectric = dielectric
        self.ionic_strength = ionic_strength
        self.kT = const.kB * temperature
        self.beta = 1.0 / self.kT
        self.debye_length = np.sqrt(
            const.epsilon_0 * dielectric * const.kB * temperature
            / (2 * const.NA * const.e**2 * ionic_strength * 1000.0)
        )


class DecoherenceModel:
    def __init__(
        self,
        dimer: TubulinDimer,
        params: ExperimentalParameters,
        protection_factor=1e-3,
        eta=0.1,
        omega_c=4.5e12,
        ionic_screening_factor=0.005,
        ion_mass_amu=23.0,
        formalism='secular_lindblad',
        low_frequency_cutoff=1e6,
    ):
        self.dimer = dimer
        self.params = params
        self.f_prot = protection_factor
        self.eta = eta
        self.omega_c = omega_c
        self.ionic_screening_factor = ionic_screening_factor
        self.ion_mass_amu = ion_mass_amu
        self.formalism = formalism
        self.low_frequency_cutoff = low_frequency_cutoff

    def Gamma_water(self):
        v_th = np.sqrt(3 * const.kB * self.params.temperature / (18 * const.amu_to_kg))
        return 3.3e28 * 1e-19 * v_th * self.f_prot

    def Gamma_ions(self):
        mu = self.dimer.transition_dipole('conformational')
        ld = self.params.debye_length
        E_rms = const.e / (4 * np.pi * const.epsilon_0 * self.params.dielectric * ld**2)
        E_rms *= self.ionic_screening_factor
        ion_mass_kg = self.ion_mass_amu * const.amu_to_kg
        tau_corr = ld / np.sqrt(3 * const.kB * self.params.temperature / ion_mass_kg)
        return (mu * E_rms / const.hbar) ** 2 * tau_corr

    def Gamma_vibrational_analytic(self):
        """High-temperature Ohmic estimate used for consistency checks."""
        return 2.0 * np.pi * self.eta * self.params.kT / const.hbar

    def Gamma_vibrational(self):
        w_min = self.low_frequency_cutoff

        def integrand(w):
            J = self.eta * (w / self.omega_c) * np.exp(-w / self.omega_c) * const.hbar
            coth = 1.0 / np.tanh(np.clip(const.hbar * w / (2 * self.params.kT), 1e-12, 300))
            return (np.pi / const.hbar) * J * coth

        low_frequency_piece = 2.0 * np.pi * self.eta * self.params.kT / const.hbar
        low_frequency_piece *= (1.0 - np.exp(-w_min / self.omega_c))
        res, _ = quad(integrand, w_min, 50 * self.omega_c, limit=200)
        return low_frequency_piece + res

    def Gamma_radiation(self):
        mu = self.dimer.transition_dipole('optical')
        return self.dimer.omega_0**3 * mu**2 / (3 * np.pi * const.epsilon_0 * const.hbar * const.c**3)

    def get_all_rates(self):
        Gw = self.Gamma_water()
        Gi = self.Gamma_ions()
        Gv = self.Gamma_vibrational()
        Gv_analytic = self.Gamma_vibrational_analytic()
        Gr = self.Gamma_radiation()
        Gtot = Gw + Gi + Gv + Gr
        dominant_source = max(
            {
                'water': Gw,
                'ions': Gi,
                'vibrational': Gv,
                'radiation': Gr,
            }.items(),
            key=lambda item: item[1],
        )[0]
        return {
            'water_collisions': Gw,
            'ionic_fluctuations': Gi,
            'vibrational_dephasing': Gv,
            'vibrational_analytic': Gv_analytic,
            'radiation': Gr,
            'total_dephasing': Gtot,
            'T2_star': 1.0 / Gtot,
            'T2_star_ps': 1e12 / Gtot,
            'dominant_source': dominant_source,
            'formalism': self.formalism,
        }


class EquilibriumBound:
    @staticmethod
    def tau_max(temperature: float, alpha: float = 0.3) -> float:
        """Loose equilibrium upper bound on coherence time."""
        return const.hbar / (alpha * const.kB * temperature)


class CoherenceUtility:
    @staticmethod
    def utility(tau_coh: float, tau_func: float, K: float) -> float:
        return K * tau_coh / tau_func

    @staticmethod
    def K_required(tau_coh: float, tau_func: float) -> float:
        return tau_func / tau_coh
