"""
heom_benchmark.py -- Converged HEOM benchmark for Redfield validation.

Systems:
  A) 1JFF full         -- 8 sites, 8 independent Drude-Lorentz baths
  B) 6DPU fragment     -- 4 sites along strongest-coupling path from A:346

Bath (synced with redfield_tubulin.py):
  lambda = 35 cm-1, gamma = 53 cm-1, T = 300 K

HEOM numerics (Converged Reference):
  hierarchy depth  NC = 7
  Matsubara terms  Nk = 1  (Pade decomposition)
  window           t_max = 30 ps, dt = 5 fs

Internal units: rad/fs (hbar = 1), time in fs.
"""
from __future__ import annotations
import json, time, pickle
from pathlib import Path
import numpy as np
import sys

# Boilerplate for package-level imports resolution
PROJECT_ROOT = Path(__file__).resolve().parents[2] # back from src/qmc_mt/ to root
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import qutip as qt
from qutip import Qobj, basis

try:
    from qutip.solver.heom import HEOMSolver, DrudeLorentzPadeBath as PadeBath
except ImportError:
    # Fallback for older versions if needed, but we target QuTiP 5
    from qutip.nonmarkov.heom import HEOMSolver, DrudeLorentzBath as PadeBath

# ---- constants ----
CM_TO_RADFS = 2 * np.pi * 2.9979e-5   # rad/fs per cm-1
KB_CM_K     = 0.69503                  # cm-1 / K
T_K         = 300.0
KT_CM       = KB_CM_K * T_K            # ~208.5 cm-1
LAMBDA_CM   = 35.0
GAMMA_B_CM  = 53.0

# HEOM numerics
NC    = 7
NK    = 1
T_MAX = 30000.0   # fs
DT    = 5.0       # fs

FIG_DIR = PROJECT_ROOT / "outputs_data" / "figures_final"
FIG_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
def select_fragment(H_cm, labels, initial_label, n_sites=4):
    """Greedy path: start at initial_label, grow fragment by strongest |V|."""
    N = len(labels)
    i0 = labels.index(initial_label)
    chosen = [i0]
    Vabs = np.abs(H_cm - np.diag(np.diag(H_cm)))
    while len(chosen) < n_sites and len(chosen) < N:
        best_j, best_v = -1, -1.0
        for i in chosen:
            for j in range(N):
                if j in chosen:
                    continue
                if Vabs[i, j] > best_v:
                    best_v, best_j = Vabs[i, j], j
        if best_j < 0:
            break
        chosen.append(best_j)
    idx = sorted(chosen)
    H_sub = H_cm[np.ix_(idx, idx)]
    lab_sub = [labels[i] for i in idx]
    new_i0 = idx.index(labels.index(initial_label))
    return H_sub, lab_sub, new_i0, idx


# --------------------------------------------------------------------------
def run_heom(H_cm, labels, initial_site_idx, tag, nc_override=None):
    nc_run = nc_override if nc_override is not None else NC
    N = len(labels)
    H_rad = H_cm * CM_TO_RADFS
    H = Qobj(H_rad)

    lam_rad = LAMBDA_CM  * CM_TO_RADFS
    gam_rad = GAMMA_B_CM * CM_TO_RADFS
    T_rad   = KT_CM      * CM_TO_RADFS

    baths = []
    L_total = qt.liouvillian(H)
    for i in range(N):
        ket = basis(N, i)
        Q = ket * ket.dag()
        bath = PadeBath(Q, lam=lam_rad, gamma=gam_rad, T=T_rad, Nk=NK)
        baths.append(bath)
        # Add terminator
        _, L_term = bath.terminator()
        L_total += L_term

    rho0 = basis(N, initial_site_idx) * basis(N, initial_site_idx).dag()
    tlist = np.arange(0, T_MAX + DT, DT)
    
    solver = HEOMSolver(L_total, baths, max_depth=nc_run, 
                        options={"nsteps": 100_000, "store_states": True})

    print(f"[{tag}] starting propagation (N={N}, NC={nc_run}, Nk={NK}, ADOs={len(solver.ados.labels)})")
    e_ops = [basis(N, i) * basis(N, i).dag() for i in range(N)]
    t0 = time.time()
    result = solver.run(rho0, tlist, e_ops=e_ops)
    elapsed = time.time() - t0
    print(f"[{tag}] HEOM finished in {elapsed/60:.1f} min")

    P_site = np.array(result.expect).T                          # (nt, N)
    rho_t_q = result.states                                     # List of Qobj
    rho_t  = np.array([s.full() for s in rho_t_q])              # (nt, N, N)
    
    return tlist, P_site, rho_t, elapsed, H, rho_t_q, [], []


# --------------------------------------------------------------------------
def analyze(tlist, P_site, rho_t, H_cm, i0, tag):
    N = P_site.shape[1]

    # tau_transfer (1/e of initial site)
    P0 = P_site[:, i0]
    P_inf = P0[-int(0.1*len(tlist)):].mean()
    target = P_inf + (P0[0] - P_inf) / np.e
    cross = np.where(P0 <= target)[0]
    tau_tr = float(tlist[cross[0]]) if len(cross) else None

    # coherence L1 (off-diagonal)
    off = rho_t.copy()
    for i in range(N):
        off[:, i, i] = 0.0
    coh_L1 = np.sum(np.abs(off), axis=(1, 2))

    peak_idx = int(np.argmax(coh_L1))
    peak_val = float(coh_L1[peak_idx])
    if peak_val < 1e-10:
        tau_coh, status = None, "no_coherence_generated"
    elif peak_idx >= len(tlist) - 5:
        tau_coh, status = None, "peak_at_end_of_window"
    else:
        after = coh_L1[peak_idx:]
        t_after = tlist[peak_idx:]
        cc = np.where(after <= peak_val / np.e)[0]
        if len(cc):
            tau_coh, status = float(t_after[cc[0]] - tlist[peak_idx]), "measured"
        else:
            tau_coh, status = None, "no_decay_in_window"

    # excitonic populations -> KL vs Gibbs
    E, C = np.linalg.eigh(H_cm)
    P_exc_t = np.real(np.einsum('im,tij,jm->tm', C.conj(), rho_t, C))
    P_ss = P_exc_t[-int(0.1*len(tlist)):].mean(axis=0)
    P_bz = np.exp(-(E - E.min()) / KT_CM); P_bz /= P_bz.sum()
    P_ss_c = np.clip(P_ss, 1e-15, 1.0)
    P_bz_c = np.clip(P_bz, 1e-15, 1.0)
    kl = float(np.sum(P_ss_c * np.log(P_ss_c / P_bz_c)))

    tr_err = float(np.max(np.abs(np.einsum('tii->t', rho_t) - 1.0)))
    herm_err = float(np.max(np.abs(rho_t - np.conj(rho_t.transpose(0, 2, 1)))))
    # spot-check min eigenvalue of rho
    sample = range(0, len(tlist), max(1, len(tlist)//50))
    min_eig = min(float(np.linalg.eigvalsh(rho_t[k]).min().real) for k in sample)

    return dict(
        tag=tag, N=N,
        tau_tr_fs=tau_tr,
        tau_coh_fs=tau_coh, tau_coh_status=status,
        coh_peak_time_fs=float(tlist[peak_idx]), coh_peak_value=peak_val,
        kl_tfinal_boltzmann=kl,
        tr_err=tr_err, herm_err=herm_err, min_eig_rho=min_eig,
    )


# --------------------------------------------------------------------------
def _heom_title(tag: str) -> str:
    """Return a publication-safe title for HEOM benchmark panels."""
    label = tag.replace("_", " ")
    return f"HEOM benchmark: {label} (NC={NC}, Nk={NK})"


def save_plot(tlist, P_site, rho_t, labels, i0, tag, outpng):
    N = P_site.shape[1]
    off = rho_t.copy()
    for i in range(N):
        off[:, i, i] = 0.0
    coh = np.sum(np.abs(off), axis=(1, 2))

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), constrained_layout=True)
    ax = axes[0]
    for i in range(N):
        ls = "-" if i == i0 else "--"
        lw = 2.0 if i == i0 else 1.0
        ax.plot(tlist/1000, P_site[:, i].real, ls, lw=lw, label=labels[i])
    ax.set_xlabel("t (ps)"); ax.set_ylabel("P(site)")
    ax.set_title(f"Populations | init={labels[i0]}")
    ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(tlist/1000, coh.real, "k-")
    ax.set_xlabel("t (ps)"); ax.set_ylabel("||rho_off||_1  (L1 coherence)")
    ax.set_title("Total L1 coherence")
    ax.grid(alpha=0.3)

    fig.suptitle(_heom_title(tag), fontsize=12)
    for ext in ("png", "pdf"):
        fname = outpng.with_suffix(f".{ext}")
        fig.savefig(fname, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure -> {outpng.stem}.{{png,pdf}}")


# --------------------------------------------------------------------------
def main():
    import sys
    skip_1jff = "--only-6dpu" in sys.argv
    skip_6dpu = "--only-1jff" in sys.argv
    print(f"Parameters: lam={LAMBDA_CM} cm-1, gamma={GAMMA_B_CM} cm-1, "
          f"T={T_K} K (kT={KT_CM:.1f} cm-1)")
    print(f"HEOM: NC={NC}, Nk={NK} (Pade), t_max={T_MAX:.0f} fs, dt={DT} fs\n")

    results = {}

    # ---- A: 1JFF full ----
    if not skip_1jff:
        d = np.load(PROJECT_ROOT / "outputs_data" / "raw_npz" / "H_1JFF.npz", allow_pickle=True)
        H_cm = d["H_cm1"]; labels = list(d["labels"])
        i0 = labels.index("B:103")
        tlist, P_site, rho_t, elapsed, H_q, rho_t_q, c_ops, lams = run_heom(H_cm, labels, i0, "1JFF_full")
        A = analyze(tlist, P_site, rho_t, H_cm, i0, "1JFF_full")
        A["initial_site"] = "B:103"
        A["labels"] = labels
        A["wallclock_s"] = elapsed
        results["1JFF_full"] = A
        save_plot(tlist, P_site, rho_t, labels, i0, "1JFF_full",
                  FIG_DIR / "heom_1JFF_full.png")
        np.savez(PROJECT_ROOT / "outputs_data" / "raw_npz" / "heom_1JFF_full.npz",
                 t_fs=tlist, P_site=P_site, rho_t=rho_t, labels=np.array(labels))
        with open(PROJECT_ROOT / "outputs_data" / "raw_pkl" / "heom_1JFF_full_trajectory.pkl", "wb") as f:
            pickle.dump({"tlist": tlist, "rho_t": rho_t_q, "H_S": H_q, 
                         "coupling_ops": c_ops, "lam_per_site": lams}, f)

    # ---- B: 6DPU fragment ----
    if not skip_6dpu:
        d = np.load(PROJECT_ROOT / "outputs_data" / "raw_npz" / "H_6DPU.npz", allow_pickle=True)
        H_cm = d["H_cm1"]; labels = list(d["labels"])
        H_sub, lab_sub, i0_sub, idx = select_fragment(H_cm, labels, "A:346", n_sites=4)
        print(f"\n[6DPU_frag] selected sites: {lab_sub}  (indices {idx})")
        tlist, P_site, rho_t, elapsed, H_q, rho_t_q, c_ops, lams = run_heom(H_sub, lab_sub, i0_sub, "6DPU_frag")
        B = analyze(tlist, P_site, rho_t, H_cm[np.ix_(idx, idx)], i0_sub, "6DPU_frag")
        B["initial_site"] = "A:346"
        B["labels"] = lab_sub
        B["wallclock_s"] = elapsed
        B["fragment_selection"] = "greedy strongest-coupling path from A:346"
        results["6DPU_fragment"] = B
        save_plot(tlist, P_site, rho_t, lab_sub, i0_sub, "6DPU_frag",
                  FIG_DIR / "heom_6DPU_frag.png")
        np.savez(PROJECT_ROOT / "outputs_data" / "raw_npz" / "heom_6DPU_frag.npz",
                 t_fs=tlist, P_site=P_site, rho_t=rho_t, labels=np.array(lab_sub))
        with open(PROJECT_ROOT / "outputs_data" / "raw_pkl" / "heom_6DPU_frag_trajectory.pkl", "wb") as f:
            pickle.dump({"tlist": tlist, "rho_t": rho_t_q, "H_S": H_q, 
                         "coupling_ops": c_ops, "lam_per_site": lams}, f)

    results["_params"] = dict(
        lambda_cm=LAMBDA_CM, gamma_bath_cm=GAMMA_B_CM, T_K=T_K, kT_cm=KT_CM,
        NC=NC, Nk=NK, t_max_fs=T_MAX, dt_fs=DT,
        method="HEOM (QuTiP) Drude-Lorentz per site, Pade Matsubara",
        units="times in fs; energies in cm-1 (converted to rad/fs internally)",
    )

    out = PROJECT_ROOT / "outputs_data" / "raw_json" / "heom_summary.json"
    out.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nwrote {out}")

    print("\n=== HEOM Converged Reference Benchmark (NC=7) ===")
    for key in results.keys():
        if key == "_params": continue
        r = results[key]
        print(f"\n[{key}]  init={r['initial_site']}  N={r['N']}  "
              f"wallclock={r.get('wallclock_s',0)/60:.1f} min")
        tr = f"{r['tau_tr_fs']:.1f}" if r['tau_tr_fs'] is not None else "N/A"
        tc = f"{r['tau_coh_fs']:.1f}" if r['tau_coh_fs'] is not None else "N/A"
        print(f"  tau_transfer      = {tr} fs")
        print(f"  tau_coherence     = {tc} fs  ({r['tau_coh_status']})")
        print(f"  KL(P_tfinal || Gibbs) = {r['kl_tfinal_boltzmann']:.4f}")
        print(f"  trace error       = {r['tr_err']:.2e}")
        print(f"  hermiticity err   = {r['herm_err']:.2e}")
        print(f"  min eig(rho)      = {r['min_eig_rho']:+.2e}")


if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    main()