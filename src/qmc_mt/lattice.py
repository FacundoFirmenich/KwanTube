"""Microtubule lattice collective modes (COMP-6).

Tight-binding exciton Hamiltonian on the canonical 13-protofilament B-lattice
with point-dipole nearest-neighbour couplings:

  J_axial   = (1/(4*pi*epsilon_0*epsilon_r)) * (-2 * mu^2 / r_axial^3)     head-to-tail, attractive (J)
  J_lateral = (1/(4*pi*epsilon_0*epsilon_r)) * (+1 * mu^2 / r_lateral^3)   side-by-side, repulsive (H)

Periodic boundary in the protofilament index, open in the axial index, with
the B-lattice 0.92 nm lateral rise realised as a shift-by-one along the
cylinder seam. This module reports excitonic band-edge energies and the IPR of
the lowest-energy eigenmode. Radiative sub/superradiance requires the separate
decay-spectrum analysis and is intentionally not inferred here from the IPR.
"""
from __future__ import annotations
import numpy as np
import sys
from pathlib import Path

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.core import const


DEFAULT_LAYER_COUNTS = (10, 20, 40)


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
    # Inverse participation ratio: 1 <= IPR <= N (extended when IPR >> 1)
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
        "lowest_mode_ipr":   IPR_sub,
        "ipr_over_n":        IPR_sub / max(N, 1),
        "subradiant_IPR":    IPR_sub,
    }


def summary_family(
    layer_counts: tuple[int, ...] = DEFAULT_LAYER_COUNTS,
    mu_debye: float = 1700.0,
    eps_r: float = 80.0,
    r_axial: float = 8.0e-9,
    r_lateral: float = 5.2e-9,
    n_protofilaments: int = 13,
) -> dict:
    family = {}
    for n_layers in layer_counts:
        item = summary(
            n_layers=n_layers,
            mu_debye=mu_debye,
            eps_r=eps_r,
            r_axial=r_axial,
            r_lateral=r_lateral,
            n_protofilaments=n_protofilaments,
        )
        family[f"N{item['N_dimers']}"] = item
    return family


def compare_family(family: dict) -> dict:
    ordered = sorted(
        family.items(),
        key=lambda kv: int(kv[1]["N_dimers"]),
    )
    rows = []
    for label, item in ordered:
        rows.append({
            "label": label,
            "N_dimers": item["N_dimers"],
            "gap_meV": item["gap_meV"],
            "lowest_mode_ipr": item["lowest_mode_ipr"],
            "ipr_over_n": item["ipr_over_n"],
        })

    pairwise = []
    for prev, curr in zip(ordered, ordered[1:]):
        prev_label, prev_item = prev
        curr_label, curr_item = curr
        pairwise.append({
            "from": prev_label,
            "to": curr_label,
            "delta_gap_meV": curr_item["gap_meV"] - prev_item["gap_meV"],
            "delta_gap_pct": 100.0 * (curr_item["gap_meV"] - prev_item["gap_meV"]) / max(abs(prev_item["gap_meV"]), 1e-12),
            "delta_lowest_mode_ipr": curr_item["lowest_mode_ipr"] - prev_item["lowest_mode_ipr"],
            "delta_lowest_mode_ipr_pct": 100.0 * (curr_item["lowest_mode_ipr"] - prev_item["lowest_mode_ipr"]) / max(abs(prev_item["lowest_mode_ipr"]), 1e-12),
            "delta_ipr_over_n": curr_item["ipr_over_n"] - prev_item["ipr_over_n"],
        })

    return {
        "ordered": rows,
        "pairwise": pairwise,
        "gap_range_meV": max((row["gap_meV"] for row in rows), default=0.0) - min((row["gap_meV"] for row in rows), default=0.0),
        "lowest_mode_ipr_range": max((row["lowest_mode_ipr"] for row in rows), default=0.0) - min((row["lowest_mode_ipr"] for row in rows), default=0.0),
    }


if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    import pprint
    family = summary_family()
    for label, item in family.items():
        print(f"--- LATTICE {label} ---")
        pprint.pprint(item)
        print()
    print("--- LATTICE FAMILY COMPARISON ---")
    pprint.pprint(compare_family(family))
