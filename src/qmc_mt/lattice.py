"""Microtubule lattice collective modes (COMP-6).

Tight-binding exciton Hamiltonian on the canonical 13-protofilament B-lattice
with point-dipole nearest-neighbour couplings:

  J_axial   = (1/4πε₀ε_r) · (−2 μ²/r_axial³)     head-to-tail, attractive (J)
  J_lateral = (1/4πε₀ε_r) · (+1 μ²/r_lateral³)   side-by-side, repulsive (H)

Periodic boundary in the protofilament index, open in the axial index, with
the B-lattice 0.92 nm lateral rise realised as a shift-by-one along the
cylinder seam. Returns the superradiant/subradiant band edges, the spectral
gap, the NN couplings in meV, and the inverse participation ratio (IPR) of
the subradiant eigenstate.
"""
from __future__ import annotations
import numpy as np
import sys
from pathlib import Path

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raíz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.core import const


def summary(n_layers: int = 20, mu_debye: float = 1700.0,
            eps_r: float = 80.0,
            r_axial: float = 8.0e-9, r_lateral: float = 5.2e-9,
            n_protofilaments: int = 13) -> dict:
    mu_SI    = mu_debye * const.Debye_to_Cm
    pref     = 1.0 / (4.0 * np.pi * const.epsilon_0 * eps_r)
    J_ax_J   = pref * (-2.0) * mu_SI**2 / r_axial**3
    J_lat_J  = pref * (+1.0) * mu_SI**2 / r_lateral**3
    J_ax_meV  = J_ax_J  / const.eV_to_J * 1e3
    J_lat_meV = J_lat_J / const.eV_to_J * 1e3

    N_pf = int(n_protofilaments)
    N_L  = int(n_layers)
    N    = N_pf * N_L

    def idx(p: int, l: int) -> int: return p * N_L + l

    H = np.zeros((N, N), dtype=float)
    for p in range(N_pf):
        for l in range(N_L):
            i = idx(p, l)
            # Axial NN (open BC along axis)
            if l + 1 < N_L:
                j = idx(p, l + 1)
                H[i, j] = H[j, i] = J_ax_meV
            # Lateral NN (periodic along cylinder, with seam shift +1)
            pn = (p + 1) % N_pf
            ln = l + 1 if pn == 0 else l        # +1 dimer rise at the seam
            if 0 <= ln < N_L:
                j = idx(pn, ln)
                H[i, j] = H[j, i] = J_lat_meV

    evals, evecs = np.linalg.eigh(H)
    E_sub, E_super = float(evals[0]), float(evals[-1])
    gap_meV        = float(E_super - E_sub)

    v_sub = evecs[:, 0]
    # Inverse participation ratio: 1 ≤ IPR ≤ N (extended when IPR ≫ 1)
    IPR_sub = float(1.0 / np.sum(v_sub**4))

    return {
        "N_dimers":          N,
        "n_layers":          N_L,
        "n_protofilaments":  N_pf,
        "mu_Debye":          mu_debye,
        "eps_r":             eps_r,
        "r_axial_nm":        r_axial * 1e9,
        "r_lateral_nm":      r_lateral * 1e9,
        "nn_axial_meV":      J_ax_meV,
        "nn_lateral_meV":    J_lat_meV,
        "E_super_meV":       E_super,
        "E_sub_meV":         E_sub,
        "gap_meV":           gap_meV,
        "subradiant_IPR":    IPR_sub,
    }


if __name__ == "__main__":
    import pprint; pprint.pprint(summary())