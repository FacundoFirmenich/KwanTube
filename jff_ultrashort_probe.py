
"""
jff_ultrashort_probe.py

Probe ultracompacto para validar (NC, Nk) en 1JFF sin freír el PC.
Diseñado para máquinas saturadas: una corrida baseline y, opcionalmente,
una sola corrida NC+1. No depende de JSON externos si se pasan --nc y --nk.

Uso recomendado:
    python jff_ultrashort_probe.py --nc 5 --nk 1
Opcional:
    python jff_ultrashort_probe.py --nc 5 --nk 1 --with-nc-stress
    python jff_ultrashort_probe.py --nc 5 --nk 1 --tmax 40 --sample 10
"""

from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import qutip as qt
from qutip.solver.heom import DrudeLorentzPadeBath, HEOMSolver

LAM_CM = 35.0
GAM_CM = 53.0
T_K = 300.0
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
T_RADFS = T_K * 0.69503 * CM_TO_RADFS


def locate_project_root() -> Path:
    here = Path(__file__).resolve().parent
    candidates = [
        here,
        here.parent,
        Path.cwd(),
        Path.cwd() / "biofisicaquantiqaCLINE",
        Path.cwd() / "git_repo",
    ]
    for c in candidates:
        if (c / "H_1JFF.npz").exists():
            return c
        if (c / "git_repo" / "H_1JFF.npz").exists():
            return c / "git_repo"
    return here


def load_hamiltonian(project_root: Path) -> qt.Qobj:
    candidates = [
        project_root / "H_1JFF.npz",
        project_root.parent / "H_1JFF.npz",
    ]
    npz_path = next((p for p in candidates if p.exists()), None)
    if npz_path is None:
        raise FileNotFoundError(
            "No encuentro H_1JFF.npz. Busqué en:\n" + "\n".join(str(p) for p in candidates)
        )
    d = np.load(npz_path, allow_pickle=True)
    H_cm = d["H_cm1"]
    H = (H_cm - np.mean(np.diag(H_cm)) * np.eye(H_cm.shape[0])) * CM_TO_RADFS
    return qt.Qobj(H)


def site_projectors(n: int):
    return [qt.basis(n, i) * qt.basis(n, i).dag() for i in range(n)]


def initial_state(n: int, site: int = 0):
    psi = qt.basis(n, site)
    return psi * psi.dag()


def build_solver(H_S, coupling_ops, NC: int, Nk: int) -> HEOMSolver:
    lam_rad = LAM_CM * CM_TO_RADFS
    gam_rad = GAM_CM * CM_TO_RADFS
    baths = [
        DrudeLorentzPadeBath(Q=Q, lam=lam_rad, gamma=gam_rad, T=T_RADFS, Nk=Nk)
        for Q in coupling_ops
    ]
    return HEOMSolver(
        H_S,
        baths,
        max_depth=NC,
        options={
            "nsteps": 100_000,
            "store_states": True,
            "progress_bar": False,
        },
    )


def run_heom(H_S, coupling_ops, NC: int, Nk: int, tmax_fs: float, sample_fs: float, rho0, label: str):
    tlist = np.arange(0.0, tmax_fs + sample_fs, sample_fs)
    solver = build_solver(H_S, coupling_ops, NC, Nk)
    print(f"[run {label}] NC={NC}, Nk={Nk}, tmax={tmax_fs} fs, sample={sample_fs} fs", flush=True)
    t0 = time.time()
    result = solver.run(rho0, tlist)
    wall = time.time() - t0
    print(f"           finished in {wall:.1f} s", flush=True)
    return tlist, result.states, wall


def max_fro_diff(states_A, states_B):
    diffs = [np.linalg.norm((a - b).full(), ord="fro") for a, b in zip(states_A, states_B)]
    idx = int(np.argmax(diffs))
    return float(diffs[idx]), idx


def max_pop_diff(states_A, states_B):
    diffs = []
    for a, b in zip(states_A, states_B):
        pa = np.real(np.diag(a.full()))
        pb = np.real(np.diag(b.full()))
        diffs.append(np.max(np.abs(pa - pb)))
    idx = int(np.argmax(diffs))
    return float(diffs[idx]), idx


def max_coh_diff(states_A, states_B):
    diffs = []
    for a, b in zip(states_A, states_B):
        da = a.full().copy()
        db = b.full().copy()
        np.fill_diagonal(da, 0.0)
        np.fill_diagonal(db, 0.0)
        diffs.append(np.max(np.abs(da - db)))
    idx = int(np.argmax(diffs))
    return float(diffs[idx]), idx


def estimate_30ps_minutes(wall_seconds: float, tmax_fs: float) -> float:
    return (wall_seconds * (30000.0 / tmax_fs)) / 60.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nc", type=int, default=5, help="NC base")
    ap.add_argument("--nk", type=int, default=1, help="Nk base")
    ap.add_argument("--tmax", type=float, default=40.0, help="ventana corta en fs")
    ap.add_argument("--sample", type=float, default=10.0, help="espaciado de salida en fs")
    ap.add_argument("--site", type=int, default=0, help="sitio inicial")
    ap.add_argument("--with-nc-stress", action="store_true", help="correr también NC+1")
    ap.add_argument("--outfile", type=str, default="jff_ultrashort_probe.json")
    args = ap.parse_args()

    project_root = locate_project_root()
    H_S = load_hamiltonian(project_root)
    N = H_S.shape[0]
    S_ops = site_projectors(N)
    rho0 = initial_state(N, site=args.site)

    out = {
        "project_root": str(project_root),
        "nc": args.nc,
        "nk": args.nk,
        "tmax_fs": args.tmax,
        "sample_fs": args.sample,
        "site": args.site,
        "runs": {},
    }

    t_b, st_b, w_b = run_heom(H_S, S_ops, args.nc, args.nk, args.tmax, args.sample, rho0, "baseline")
    out["runs"]["baseline"] = {
        "wall_s": w_b,
        "est_30ps_min": estimate_30ps_minutes(w_b, args.tmax),
        "n_times": len(t_b),
    }

    if args.with_nc_stress:
        t_n, st_n, w_n = run_heom(H_S, S_ops, args.nc + 1, args.nk, args.tmax, args.sample, rho0, "nc_stress")
        dF, iF = max_fro_diff(st_b, st_n)
        dP, iP = max_pop_diff(st_b, st_n)
        dC, iC = max_coh_diff(st_b, st_n)
        out["runs"]["nc_stress"] = {
            "wall_s": w_n,
            "est_30ps_min": estimate_30ps_minutes(w_n, args.tmax),
            "dFrob_max": dF,
            "dPop_max": dP,
            "dCoh_max": dC,
            "time_index_frob": iF,
            "time_index_pop": iP,
            "time_index_coh": iC,
        }

        print("\n" + "=" * 68)
        print("ULTRA-SHORT 1JFF PROBE")
        print(f"baseline       NC={args.nc}, Nk={args.nk}, wall={w_b:.1f} s, est30ps={estimate_30ps_minutes(w_b, args.tmax):.1f} min")
        print(f"NC stress      NC={args.nc+1}, Nk={args.nk}, wall={w_n:.1f} s, est30ps={estimate_30ps_minutes(w_n, args.tmax):.1f} min")
        print(f"dPop max       {dP:.2e}")
        print(f"dCoh max       {dC:.2e}")
        print(f"dFrob max      {dF:.2e}")
    else:
        print("\n" + "=" * 68)
        print("ULTRA-SHORT 1JFF PROBE")
        print(f"baseline only  NC={args.nc}, Nk={args.nk}, wall={w_b:.1f} s, est30ps={estimate_30ps_minutes(w_b, args.tmax):.1f} min")

    with open(args.outfile, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {args.outfile}")


if __name__ == "__main__":
    main()
