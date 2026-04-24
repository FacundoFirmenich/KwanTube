"""
redfield_tubulin.py — Redfield secular con baño Drude-Lorentz
sobre la red de Trps de tubulina.

Modelo
------
Sistema:  H electrónico (desde build_hamiltonian.py)  en base de sitios.
Baño:     un oscilador de Drude por sitio (no correlacionado entre sitios),
          densidad espectral J(ω) = 2λγω/(ω²+γ²).
          λ = 35 cm⁻¹ (reorganización típica cromóforo en proteína; Renger 2009).
          γ = 53 cm⁻¹ ↔ τ_bath ≈ 100 fs (correlación protein-water bath).
T = 300 K  ⇒  kT = 208.5 cm⁻¹.

Ecuación (base excitónica μ,ν = autoestados de H):
  poblaciones:   dP_μ/dt = Σ_ν [k_{μ←ν} P_ν  −  k_{ν←μ} P_μ]
  coherencias:   dρ_μν/dt = (−i ω_μν − Γ_μν) ρ_μν
con
  k_{μ←ν} = 2π Σ_i |c_iμ|² |c_iν|² · Re C(ω_μν)
  Γ_μν    = ½(k_μ^out + k_ν^out) + γ_vibr_pd + γ_μν^pd
  γ_vibr_pd = 2π·α·kT (vibrational dephasing baseline)
"""
from __future__ import annotations
import sys, json, os
from pathlib import Path
import numpy as np

# Boilerplate para resolver importaciones y rutas desde la raíz del proyecto
PROJECT_ROOT = Path(__file__).resolve().parents[3] # retrocede desde src/qmc_mt/ a la raíz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- constantes ----------
CM_TO_INVFS = 2 * np.pi * 2.9979e-5      # rad·fs⁻¹ por cm⁻¹
KB_CM_K     = 0.69503                     # cm⁻¹ / K
T_K         = 300.0
KT          = KB_CM_K * T_K               # ≈ 208.5 cm⁻¹
BETA        = 1.0 / KT                    # cm

LAMBDA_CM   = 35.0    # reorganización (Renger 2009, cromóforo en proteína)
GAMMA_B_CM  = 53.0    # inverso tau_bath (~100 fs)

DATA_DIR = PROJECT_ROOT
FIG_DIR  = PROJECT_ROOT / "figures_final"
FIG_DIR.mkdir(exist_ok=True)

def save_paper(fig, name, outdir=FIG_DIR):
    """Guarda la figura en formatos PNG y PDF para calidad de publicación."""
    os.makedirs(outdir, exist_ok=True)
    for ext in ("png", "pdf"):
        path = outdir / f"{name}.{ext}"
        fig.savefig(str(path), dpi=600, bbox_inches="tight")
    print(f"  figuras -> {name}.[png, pdf] en {outdir.name}/")

# ---------- función correlación baño (Drude-Lorentz) ----------
def S_quantum_cm(omega_cm: np.ndarray) -> np.ndarray:
    """
    Densidad espectral asimétrica S(ω) = 2*Re[C(ω)] que satisface
    el balance detallado: S(ω) = e^{βω} S(-ω).
    """
    w = np.asarray(omega_cm, dtype=float)
    out = np.empty_like(w)
    # límite ω->0: S(0) = 4 * λ * kT / γ
    small = np.abs(w) < 1e-6
    out[small] = 4.0 * LAMBDA_CM * KT / GAMMA_B_CM
    if np.any(~small):
        ww = w[~small]
        # J(ω) = 2λγω / (ω² + γ²)
        j_w = 2.0 * LAMBDA_CM * GAMMA_B_CM * ww / (ww**2 + GAMMA_B_CM**2)
        # S(ω) = 2 * J(ω) / (1 - np.exp(-BETA * ww))
        out[~small] = 2.0 * j_w / (1.0 - np.exp(-BETA * ww))
    return out

# ---------- construcción del generador ----------
def build_generator(H_cm: np.ndarray):
    """Devuelve (E, C, k_rates, Gamma_coh) en base excitónica, unidades cm⁻¹."""
    E, C = np.linalg.eigh(H_cm)     # C[:,μ] = autovector μ
    N = len(E)
    pop = np.abs(C)**2              # pop[i,μ] = |c_iμ|²

    # 1. Tasas de relajación (Poblaciones)
    k = np.zeros((N, N))
    for mu in range(N):
        for nu in range(N):
            if mu == nu: continue
            omega = E[nu] - E[mu] # nu -> mu (E_mu < E_nu => omega > 0 => emisión)
            overlap = float(np.sum(pop[:, mu] * pop[:, nu]))
            k[mu, nu] = 2.0 * np.pi * overlap * float(S_quantum_cm(np.array([omega]))[0])

    k_out = np.sum(k, axis=0) # k_out[mu] = sum_nu k[nu, mu]

    # 2. Decoherencia (Coherencias)
    # a) Baseline vibracional (2*pi*alpha*kT/h)
    alpha = 0.1
    gamma_vibr_pd = (2.0 * np.pi) * alpha * KT 
    
    # b) Pure dephasing excitónico inducido por el baño
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

# ---------- propagador ----------
def propagate(H_cm: np.ndarray, rho0_site: np.ndarray,
              t_max_fs: float = 5000.0, dt_fs: float = 1.0):
    E, C, k, k_out, Gamma_coh = build_generator(H_cm)
    N = len(E)

    # rho inicial en base excitónica
    rho_exc = C.conj().T @ rho0_site @ C

    # conversión cm⁻¹ -> fs⁻¹
    k_fs      = k * CM_TO_INVFS
    k_out_fs  = k_out * CM_TO_INVFS
    Gamma_fs  = Gamma_coh * CM_TO_INVFS
    omega_fs  = (E[:, None] - E[None, :]) * CM_TO_INVFS

    # generador coherencias: L_coh[μ,ν] = -iω_μν - Γ_μν
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

# ---------- validaciones ----------
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

    out_npz = FIG_DIR / f"redfield_{pid}.npz"
    np.savez(out_npz, t_fs=t, P_site=P_site, coh_tot=coh_tot, labels=np.array(labels),
             tau_tr=tau_tr, tau_coh=tau_coh, avg_ipr=avg_ipr)
    
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for i in range(P_site.shape[1]):
        axes[0].plot(t/1000, P_site[:, i], label=labels[i], alpha=0.7)
    axes[0].set_title(f"{pid} - Poblaciones"); axes[0].legend(fontsize=6)
    axes[1].plot(t/1000, coh_tot, "k-")
    axes[1].set_title(f"{pid} - Coherencia L1"); fig.tight_layout()
    save_paper(fig, f"redfield_{pid}")
    plt.close(fig)

    return {"pdb": pid, "results": {"tau_tr": tau_tr, "tau_coh": tau_coh, "avg_ipr": avg_ipr, "Vmax": float(np.max(np.abs(H-np.diag(np.diag(H))))), "deltaE_exc": float(E.max()-E.min())}}

if __name__ == "__main__":
    summary = []
    summary.append(run(DATA_DIR / "H_1JFF.npz", "B:103"))
    summary.append(run(DATA_DIR / "H_6DPU.npz", "A:346"))
    out_json = FIG_DIR / "redfield_summary.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"\nresumen -> {out_json.name}")