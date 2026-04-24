#!/usr/bin/env python3
"""
jff_nk_stress_probe.py

Minimal probe to compare Nk=1 vs Nk=2 at a fixed NC level for 1JFF.
Designed for rapid verification of the Padé bath truncation.
"""

import argparse
import json
import time
from pathlib import Path
import numpy as np
import qutip as qt
from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver

# Constants
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
LAM_CM, GAM_CM, T_K = 35.0, 53.0, 300.0
T_RADFS = T_K * 0.69503 * CM_TO_RADFS
LAM_RADFS = LAM_CM * CM_TO_RADFS
GAM_RADFS = GAM_CM * CM_TO_RADFS

def run_heom(H_S, NC, Nk, tlist, rho0):
    n_sites = H_S.shape[0]
    baths = [DrudeLorentzPadeBath(qt.basis(n_sites, i) * qt.basis(n_sites, i).dag(), 
                                  lam=LAM_RADFS, gamma=GAM_RADFS, T=T_RADFS, Nk=Nk) 
             for i in range(n_sites)]
    solver = HEOMSolver(H_S, baths, max_depth=NC, options={"store_states": True})
    return solver.run(rho0, tlist).states

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nc", type=int, default=5)
    parser.add_argument("--site", type=int, default=5)
    parser.add_argument("--tmax", type=float, default=80.0)
    parser.add_argument("--sample", type=float, default=10.0)
    args = parser.parse_args()

    # Load 1JFF
    npz_path = Path("H_1JFF.npz")
    if not npz_path.exists():
        npz_path = Path(__file__).parent / "H_1JFF.npz"
    data = np.load(npz_path)
    H_cm = data["H_cm1"]
    H_rad = (H_cm - np.mean(np.diag(H_cm)) * np.eye(H_cm.shape[0])) * CM_TO_RADFS
    H_S = qt.Qobj(H_rad)
    N = H_S.shape[0]
    rho0 = qt.basis(N, args.site) * qt.basis(N, args.site).dag()
    tlist = np.linspace(0, args.tmax, int(args.tmax/args.sample) + 1)

    print(f"Starting Nk Stress Probe (NC={args.nc}, tmax={args.tmax}fs)")
    
    # Run Nk=1
    print("Running Nk=1...")
    states1 = run_heom(H_S, args.nc, 1, tlist, rho0)
    
    # Run Nk=2
    print("Running Nk=2...")
    states2 = run_heom(H_S, args.nc, 2, tlist, rho0)

    # Metrics
    dpop = []
    dfrob = []
    for s1, s2 in zip(states1, states2):
        diff = s1.full() - s2.full()
        dpop.append(np.abs(np.real(np.diag(diff))))
        dfrob.append(np.linalg.norm(diff, ord="fro"))
    
    dpop = np.array(dpop)
    max_pop_idx = np.unravel_index(np.argmax(dpop), dpop.shape)
    max_frob_idx = np.argmax(dfrob)

    results = {
        "dPop_max": float(np.max(dpop)),
        "t_dPop_max": float(tlist[max_pop_idx[0]]),
        "site_dPop_max": int(max_pop_idx[1]),
        "dFrob_max": float(np.max(dfrob)),
        "t_dFrob_max": float(tlist[max_frob_idx]),
        "final_dPop": float(dpop[-1].max()),
        "final_dFrob": float(dfrob[-1])
    }

    # Output
    with open("jff_nk_stress_probe.json", "w") as f:
        json.dump(results, f, indent=2)
    
    with open("jff_nk_stress_probe.txt", "w") as f:
        f.write("=== JFF NK STRESS PROBE RESULTS ===\n")
        for k, v in results.items():
            f.write(f"{k}: {v}\n")

    print("\nProbe Complete.")
    print(f"dPop_max: {results['dPop_max']:.2e} at {results['t_dPop_max']} fs")
    print(f"dFrob_max: {results['dFrob_max']:.2e}")

if __name__ == "__main__":
    main()
