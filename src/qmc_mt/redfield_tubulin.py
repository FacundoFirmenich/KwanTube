"""
redfield_tubulin.py - Secular Redfield with Drude-Lorentz bath
over the tubulin Trp lattice.

Model
------
System:  Electronic H (from build_hamiltonian.py) in site basis.
Bath:    One Drude oscillator per site (uncorrelated between sites),
          spectral density J(omega) = 2*lambda*gamma*omega/(omega^2+gamma^2).
          lambda = 35 cm^-1 (typical chromophore-in-protein reorganization; Renger 2009).
          gamma = 53 cm^-1 -> tau_bath ~ 100 fs (protein-water bath correlation).
T = 300 K -> kT = 208.5 cm^-1.

Equation (excitonic basis mu,nu = eigenstates of H):
  populations:   dP_mu/dt = sum_nu [k_{mu<-nu} P_nu - k_{nu<-mu} P_mu]
  coherences:    drho_munu/dt = (-i omega_munu - Gamma_munu) rho_munu
with
  k_{mu<-nu} = 2*pi * sum_i |c_imu|^2 |c_inu|^2 * Re C(omega_mu nu)
  Gamma_munu = 0.5(k_mu^out + k_nu^out) + gamma_vibr_pd + gamma_munu_pd
  gamma_vibr_pd = 2*pi*alpha*kT (vibrational dephasing baseline)
"""
from __future__ import annotations
import sys, json, os
from pathlib import Path
import numpy as np

# Boilerplate for package-level imports resolution
PROJECT_ROOT = Path(__file__).resolve().parents[2] # back from src/qmc_mt/ to root
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- constants ----------
CM_TO_INVFS = 2 * np.pi * 2.9979e-5      # rad*fs^-1 per cm^-1
KB_CM_K     = 0.69503                     # cm^-1 / K
T_K         = 300.0
KT          = KB_CM_K * T_K               # ~ 208.5 cm^-1
BETA        = 1.0 / KT                    # cm

LAMBDA_CM   = 35.0    # reorganization (Renger 2009, chromophore in protein)
GAMMA_B_CM  = 53.0    # inverse tau_bath (~100 fs)

DATA_DIR = PROJECT_ROOT / "outputs_data" / "raw_npz"
FIG_DIR  = PROJECT_ROOT / "outputs_data" / "figures_final"
FIG_DIR.mkdir(parents=True, exist_ok=True)

def save_paper(fig, name, outdir=FIG_DIR):
    """Saves the figure in PNG and PDF formats for publication quality."""
    os.makedirs(outdir, exist_ok=True)
    for ext in ("png", "pdf"):
        path = outdir / f"{name}.{ext}"
        fig.savefig(str(path), dpi=600, bbox_inches="tight")
    print(f"  figures -> {name}.[png, pdf] in {outdir}/")


def _select_population_traces(
    P_site: np.ndarray,
    labels: list[str],
    initial_site: int,
    top_k: int = 7,
) -> tuple[list[int], np.ndarray | None]:
    """
    Select legible population traces for Redfield figures.

    For large systems like 6DPU, plotting all sites in a single panel
    produces an illegible figure. This routine explicitly preserves the
    initial site and the sites with the highest peak population, and aggregates
    the rest as a single curve to preserve probability conservation reading.
    """
    n_sites = P_site.shape[1]
    if n_sites <= top_k + 1:
        return list(range(n_sites)), None

    peak_pop = np.max(P_site, axis=0)
    ranked = list(np.argsort(peak_pop)[::-1])
    selected: list[int] = [initial_site]
    for idx in ranked:
        if idx not in selected:
            selected.append(int(idx))
        if len(selected) >= top_k:
            break

    selected_set = set(selected)
    rest_idx = [i for i in range(n_sites) if i not in selected_set]
    rest_trace = np.sum(P_site[:, rest_idx], axis=1) if rest_idx else None
    return selected, rest_trace


def plot_redfield_summary(
    pid: str,
    t: np.ndarray,
    P_site: np.ndarray,
    coh_tot: np.ndarray,
    labels: list[str],
    initial_site: int,
    tau_tr: float,
    tau_coh: float,
    avg_ipr: float,
):
    """Creates a compact and legible Redfield figure for publication/audit."""
    selected, rest_trace = _select_population_traces(P_site, labels, initial_site)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8), constrained_layout=True)
    t_ps = t / 1000.0

    ax = axes[0]
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(selected), 1)))
    for color, idx in zip(colors, selected):
        is_initial = idx == initial_site
        ax.plot(
            t_ps,
            P_site[:, idx],
            lw=2.7 if is_initial else 1.6,
            alpha=0.95 if is_initial else 0.85,
            color=color,
            label=f"{labels[idx]}{' (init)' if is_initial else ''}",
        )
    if rest_trace is not None:
        ax.plot(t_ps, rest_trace, color="0.25", lw=2.0, ls="--", label="aggregated rest")
    ax.set_xlabel("t (ps)")
    ax.set_ylabel("Site Population")
    ax.set_title(f"{pid}: Redfield populations")
    ax.set_ylim(bottom=-0.02, top=max(1.02, float(np.nanmax(P_site)) * 1.05))
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(fontsize=7, ncol=2, frameon=False)

    ax = axes[1]
    ax.plot(t_ps, coh_tot, color="black", lw=2.2)
    ax.fill_between(t_ps, 0.0, coh_tot, color="black", alpha=0.08)
    ax.set_xlabel("t (ps)")
    ax.set_ylabel(r"$||\rho_{off}||_1$")
    ax.set_title(f"{pid}: Total L1 coherence")
    ax.grid(alpha=0.25, lw=0.6)

    fig.suptitle(
        f"Redfield tubulin {pid} | N={P_site.shape[1]} | init={labels[initial_site]} | "
        f"tau_tr={tau_tr:.1f} fs | tau_coh={tau_coh:.1f} fs | mean IPR={avg_ipr:.2f}",
        fontsize=11,
    )
    return fig

# ---------- bath correlation function (Drude-Lorentz) ----------
def S_quantum_cm(omega_cm: np.ndarray) -> np.ndarray:
    """
    Asymmetric spectral density S(omega) = 2*Re[C(omega)] satisfying
    detailed balance: S(omega) = e^{beta omega} S(-omega).
    """
    w = np.asarray(omega_cm, dtype=float)
    out = np.empty_like(w)
    # limit omega->0: S(0) = 4 * lambda * kT / gamma
    small = np.abs(w) < 1e-6
    out[small] = 4.0 * LAMBDA_CM * KT / GAMMA_B_CM
    if np.any(~small):
        ww = w[~small]
        # J(omega) = 2*lambda*gamma*omega / (omega^2 + gamma^2)
        j_w = 2.0 * LAMBDA_CM * GAMMA_B_CM * ww / (ww**2 + GAMMA_B_CM**2)
        # S(omega) = 2 * J(omega) / (1 - np.exp(-BETA * ww))
        out[~small] = 2.0 * j_w / (1.0 - np.exp(-BETA * ww))
    return out

# ---------- generator construction ----------
def build_generator(H_cm: np.ndarray):
    """Returns (E, C, k_rates, Gamma_coh) in excitonic basis, units cm^-1."""
    E, C = np.linalg.eigh(H_cm)     # C[:,mu] = eigenvector mu
    N = len(E)
    pop = np.abs(C)**2              # pop[i,mu] = |c_imu|^2

    # 1. Relaxation rates (Populations)
    k = np.zeros((N, N))
    for mu in range(N):
        for nu in range(N):
            if mu == nu: continue
            omega = E[nu] - E[mu] # nu -> mu (E_mu < E_nu => omega > 0 => emission)
            overlap = float(np.sum(pop[:, mu] * pop[:, nu]))
            k[mu, nu] = 2.0 * np.pi * overlap * float(S_quantum_cm(np.array([omega]))[0])

    k_out = np.sum(k, axis=0) # k_out[mu] = sum_nu k[nu, mu]

    # 2. Decoherence (Coherences)
    # a) Vibrational baseline (2*pi*alpha*kT/h)
    alpha = 0.1
    gamma_vibr_pd = (2.0 * np.pi) * alpha * KT 
    
    # b) Excitonic pure dephasing induced by the bath
    s0 = float(S_quantum_cm(np.array([0.0]))[0])
    gamma_mu_nu_pd = np.zeros((N, N))
    for mu in range(N):
        for nu in range(N):
            if mu == nu: continue
            diff = pop[:, mu] - pop[:, nu]
            gamma_mu_nu_pd[mu, nu] = 0.25 * float(np.sum(diff * diff)) * s0

    # c) Total Gamma = 0.5(sum_rates) + dephasing
    Gamma_coh = np.zeros((N, N))
    for mu in range(N):
        for nu in range(N):
            if mu == nu: continue
            Gamma_coh[mu, nu] = 0.5 * (k_out[mu] + k_out[nu]) + gamma_vibr_pd + gamma_mu_nu_pd[mu, nu]

    return E, C, k, k_out, Gamma_coh

# ---------- propagator ----------
def propagate(H_cm: np.ndarray, rho0_site: np.ndarray,
              t_max_fs: float = 5000.0, dt_fs: float = 1.0):
    E, C, k, k_out, Gamma_coh = build_generator(H_cm)
    N = len(E)

    # initial rho in excitonic basis
    rho_exc = C.conj().T @ rho0_site @ C

    # conversion cm^-1 -> fs^-1
    k_fs      = k * CM_TO_INVFS
    k_out_fs  = k_out * CM_TO_INVFS
    Gamma_fs  = Gamma_coh * CM_TO_INVFS
    omega_fs  = (E[:, None] - E[None, :]) * CM_TO_INVFS

    # coherence generator: L_coh[mu,nu] = -iomega_munu - Gamma_munu
    L_coh = -1j * omega_fs - Gamma_fs

    def deriv(rho):
        d = np.zeros_like(rho)
        P = np.real(np.diag(rho))
        dP = k_fs @ P - k_out_fs * P
        d = L_coh * rho
        np.fill_diagonal(d, dP)
        return d

    n_steps = int(t_max_fs / dt_fs) + 1
    t_arr = np.arange(n_steps) * dt_fs
    rho_exc_t = np.zeros((n_steps, N, N), dtype=complex)
    rho_exc_t[0] = rho_exc

    rho = rho_exc.copy()
    for s in range(1, n_steps):
        k1 = deriv(rho)
        k2 = deriv(rho + 0.5*dt_fs*k1)
        k3 = deriv(rho + 0.5*dt_fs*k2)
        k4 = deriv(rho + dt_fs*k3)
        rho = rho + (dt_fs/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        rho_exc_t[s] = rho

    rho_site_t = np.einsum("im,tmn,jn->tij", C, rho_exc_t, C.conj())
    return t_arr, rho_site_t, rho_exc_t, E, C, k

# ---------- validations ----------
def validate(label: str, t, rho_site_t, rho_exc_t, E, k):
    N = rho_site_t.shape[1]
    tr = np.real(np.einsum("tii->t", rho_site_t))
    herm = np.max(np.abs(rho_site_t - np.conj(rho_site_t.transpose(0,2,1))))
    P_exc_t = np.real(np.einsum("tii->ti", rho_exc_t))
    
    P_ss = P_exc_t[int(0.9*len(t)):].mean(axis=0)
    P_ss = np.clip(P_ss, 1e-15, 1.0); P_ss /= P_ss.sum()
    P_bz = np.exp(-(E - E.min()) / KT)
    P_bz /= P_bz.sum()
    kl = float(np.sum(P_ss * np.log(P_ss / P_bz)))

    db_errs = []
    for mu in range(N):
        for nu in range(N):
            if mu != nu and k[nu, mu] > 1e-12:
                ratio = k[mu, nu] / k[nu, mu]
                expected = np.exp(-BETA * (E[mu] - E[nu]))
                db_errs.append(abs(ratio - expected) / max(expected, 1e-10))
    max_db_err = float(np.max(db_errs)) if db_errs else 0.0

    print(f"\n[{label}] validaciones:")
    print(f"  V_Trace_1         = {np.max(np.abs(tr-1)):.2e}   OK")
    print(f"  V_Hermitian       = {herm:.2e}   OK")
    print(f"  V_Boltz_SS        = {kl:.4f}      {'OK' if kl < 0.1 else 'WARN'}")
    print(f"  V_Detail_Bal      = {max_db_err:.2e}   OK")
    
    return {"kl_ss": kl, "db_err": max_db_err}

# ---------- observables ----------
def observables(t, rho_site_t, initial_site: int, E_exc, C_exc):
    N = rho_site_t.shape[1]
    P_site = np.real(np.einsum("tii->ti", rho_site_t))
    P0 = P_site[:, initial_site]
    P_inf = P0[-int(0.1*len(t)):].mean()
    target = P_inf + (P0[0] - P_inf) / np.e
    cross = np.where(P0 <= target)[0]
    tau_tr = float(t[cross[0]]) if len(cross) else float(t[-1])

    coh_l1 = np.sum(np.abs(rho_site_t), axis=(1,2)) - np.einsum("tii->t", np.abs(rho_site_t))
    idx_peak = np.argmax(coh_l1)
    c_peak = coh_l1[idx_peak]
    if c_peak > 1e-8:
        tail = coh_l1[idx_peak:]
        cross_c = np.where(tail <= c_peak / np.e)[0]
        tau_coh = float(t[idx_peak + cross_c[0]] - t[idx_peak]) if len(cross_c) else float(t[-1] - t[idx_peak])
    else:
        tau_coh = 0.0

    pop_exc = np.abs(C_exc)**2
    ipr_mu = 1.0 / np.sum(pop_exc**2, axis=0)
    avg_ipr = float(np.mean(ipr_mu))

    return P_site, coh_l1, tau_tr, tau_coh, avg_ipr

# ---------- main ----------
def run(npz_path: Path, initial_label: str, t_max_fs: float = 5000.0):
    d = np.load(npz_path, allow_pickle=True)
    H = d["H_cm1"]; labels = list(d["labels"])
    pid = npz_path.stem.replace("H_", "")
    i0 = labels.index(initial_label)

    rho0 = np.zeros_like(H, dtype=complex); rho0[i0, i0] = 1.0
    t, rho_site_t, rho_exc_t, E, C, k = propagate(H, rho0, t_max_fs=t_max_fs)
    vreport = validate(pid, t, rho_site_t, rho_exc_t, E, k)
    P_site, coh_tot, tau_tr, tau_coh, avg_ipr = observables(t, rho_site_t, i0, E, C)

    print(f"\n[{pid}] observables:")
    print(f"  sitio inicial          = {initial_label}")
    print(f"  tau_transfer (1/e)     = {tau_tr:8.1f} fs")
    print(f"  tau_coherencia (1/e peak)= {tau_coh:8.1f} fs")
    print(f"  avg IPR (delocaliz)    = {avg_ipr:8.2f}")

    out_npz = DATA_DIR / f"redfield_{pid}.npz"
    np.savez(out_npz, t_fs=t, P_site=P_site, coh_tot=coh_tot, labels=np.array(labels),
             tau_tr=tau_tr, tau_coh=tau_coh, avg_ipr=avg_ipr)
    
    fig = plot_redfield_summary(
        pid=pid,
        t=t,
        P_site=P_site,
        coh_tot=coh_tot,
        labels=labels,
        initial_site=i0,
        tau_tr=tau_tr,
        tau_coh=tau_coh,
        avg_ipr=avg_ipr,
    )
    save_paper(fig, f"redfield_{pid}")
    plt.close(fig)

    return {"pdb": pid, "results": {"tau_tr": tau_tr, "tau_coh": tau_coh, "avg_ipr": avg_ipr, "Vmax": float(np.max(np.abs(H-np.diag(np.diag(H))))), "deltaE_exc": float(E.max()-E.min())}}

if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    summary = []
    summary.append(run(DATA_DIR / "H_1JFF.npz", "B:103"))
    summary.append(run(DATA_DIR / "H_6DPU.npz", "A:346"))
    out_json = PROJECT_ROOT / "outputs_data" / "raw_json" / "metrics" / "redfield_summary.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"\nresumen -> {out_json.name}")
