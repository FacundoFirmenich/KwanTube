"""
diagnose_ss_and_meanforce.py
Diagnóstico del "fallo" KL(ρ_∞ ‖ ρ_Gibbs)=0.93 en 1JFF.

Hipótesis a discriminar:
  H1: El HEOM no llegó al steady state en t_max=10 ps.
  H2: El steady state físico es ρ_mean-force ≠ ρ_Gibbs_bare.
  H3: Ambas.

Uso: requiere tener guardados rho_t (lista de Qobj) y tlist en archivos .pkl.
"""

import numpy as np
import qutip as qt
import sys
import pickle, json
from pathlib import Path

# Boilerplate para resolver importaciones desde la raíz del paquete
PROJECT_ROOT = Path(__file__).resolve().parent.parent # La raíz es biofisicaquantiqaCLINE
if str(PROJECT_ROOT / "git_repo" / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "git_repo" / "src"))

# ---------- INPUTS (ajustar rutas al repo) ----------
LAM_CM    = 35.0       # cm^-1, reorganization energy per site
GAM_CM    = 53.0       # cm^-1, Drude cutoff
T_K       = 300.0
KT_CM     = 208.5      # cm^-1
CM_TO_FS  = 5308.8     # conversion factor omega[cm^-1] * t[fs] / CM_TO_FS = phase
# beta in units where H is in cm^-1:
BETA = 1.0 / KT_CM

# Paths — sincronizados con la raíz del proyecto
HEOM_DUMP_1JFF = PROJECT_ROOT / "heom_1JFF_full_trajectory.pkl"
HEOM_DUMP_6DPU = PROJECT_ROOT / "heom_6DPU_frag_trajectory.pkl"


def load_dump(path):
    with open(path, "rb") as f:
        d = pickle.load(f)
    return d


# ---------- (A) Diagnóstico de convergencia al steady state ----------
def steady_state_drift(tlist, rho_t, window_fs=1000.0):
    """
    Mide ||dρ/dt|| promediado sobre los últimos `window_fs` femtosegundos.
    Si el drift es >> 1/t_max, NO estás en steady state.
    """
    tlist = np.asarray(tlist)
    dt = tlist[1] - tlist[0]
    mask = tlist >= (tlist[-1] - window_fs)
    idx = np.where(mask)[0]
    drifts = []
    for k in idx[:-1]:
        drho = rho_t[k+1] - rho_t[k]
        # Frobenius norm
        drifts.append(np.linalg.norm(drho.full(), ord='fro') / dt)
    drift_rate = np.mean(drifts)           # norm/fs
    
    # Extrapolación lineal del steady state asumiendo relajación exponencial
    k_mid = len(rho_t) - len(idx)//2
    k_end = len(rho_t) - 1
    rho_mid = rho_t[k_mid]
    rho_end = rho_t[k_end]
    delta = rho_end - rho_mid
    delta_norm = np.linalg.norm(delta.full(), ord='fro')
    
    return {
        "drift_rate_per_fs": float(drift_rate),
        "delta_last_window_frobenius": float(delta_norm),
        "t_end_fs": float(tlist[-1]),
        "verdict": (
            "NOT_STEADY" if delta_norm > 1e-3 else
            "MARGINAL"   if delta_norm > 1e-4 else
            "STEADY"
        ),
    }


# ---------- (B) Referencias termodinámicas ----------
def bare_gibbs(H_S, beta=BETA):
    """ρ_0 = exp(-βH_S)/Z."""
    E = H_S.eigenenergies()
    E0 = E.min()
    op = (-beta * (H_S - E0 * qt.qeye(H_S.dims[0][0]))).expm()
    return op / op.tr()


def polaron_shifted_gibbs(H_S, coupling_ops, lam_list, beta=BETA):
    """Reference #2: 'polaron' reference."""
    H_pol = H_S.copy()
    for Sn, lam in zip(coupling_ops, lam_list):
        H_pol = H_pol - lam * (Sn * Sn)
    return bare_gibbs(H_pol, beta=beta)


def mean_force_gibbs_2nd_order(H_S, coupling_ops, lam_list, gam_cm=GAM_CM, beta=BETA):
    """Mean-force Gibbs a 2º orden en acoplamiento (alta T, Drude-Lorentz)."""
    H_MF = H_S.copy()
    # Primer término: shift de reorganización
    for Sn, lam in zip(coupling_ops, lam_list):
        H_MF = H_MF - lam * (Sn * Sn)
    # Segundo término: conmutador anidado
    f_mod = beta / (beta * gam_cm + 2.0) / gam_cm
    for Sn, lam in zip(coupling_ops, lam_list):
        nested = qt.commutator(Sn, qt.commutator(Sn, H_S))
        H_MF = H_MF - 0.5 * lam * nested * f_mod
    return bare_gibbs(H_MF, beta=beta)


# ---------- (C) KL entre estados mixtos ----------
def kl_quantum(rho, sigma, eps=1e-12):
    """S(ρ ‖ σ) = Tr[ρ log ρ] - Tr[ρ log σ] (en nats)."""
    w_r, V_r = np.linalg.eigh(rho.full())
    w_r = np.clip(w_r, eps, None)
    log_rho = V_r @ np.diag(np.log(w_r)) @ V_r.conj().T
    
    w_s, V_s = np.linalg.eigh(sigma.full())
    w_s = np.clip(w_s, eps, None)
    log_sigma = V_s @ np.diag(np.log(w_s)) @ V_s.conj().T
    
    rho_np = rho.full()
    kl = np.real(np.trace(rho_np @ log_rho) - np.trace(rho_np @ log_sigma))
    return float(kl)


# ---------- (D) Pipeline ----------
def diagnose(dump_path, label):
    print(f"\n{'='*60}\n  {label}\n{'='*60}")
    d = load_dump(dump_path)
    tlist   = d["tlist"]
    rho_t   = d["rho_t"]
    H_S     = d["H_S"]
    S_ops   = d["coupling_ops"]
    lams    = d["lam_per_site"]
    
    ss = steady_state_drift(tlist, rho_t, window_fs=1000.0)
    print(f"[SS convergence]")
    for k, v in ss.items():
        print(f"  {k:35s} = {v}")
    
    rho_end = rho_t[-1]
    
    rho_bare   = bare_gibbs(H_S)
    rho_pol    = polaron_shifted_gibbs(H_S, S_ops, lams)
    rho_mf     = mean_force_gibbs_2nd_order(H_S, S_ops, lams)
    
    kl_bare = kl_quantum(rho_end, rho_bare)
    kl_pol  = kl_quantum(rho_end, rho_pol)
    kl_mf   = kl_quantum(rho_end, rho_mf)
    kl_refs = kl_quantum(rho_bare, rho_mf)
    
    print(f"\n[KL divergences at t={tlist[-1]:.0f} fs]")
    print(f"  KL(rho_HEOM || rho_bare_Gibbs)     = {kl_bare:.4f} nats")
    print(f"  KL(rho_HEOM || rho_polaron_Gibbs)  = {kl_pol:.4f} nats")
    print(f"  KL(rho_HEOM || rho_meanforce_2nd)  = {kl_mf:.4f} nats")
    print(f"  KL(rho_bare || rho_meanforce_2nd)  = {kl_refs:.4f} nats")
    
    threshold = 0.05
    best_kl = min(kl_bare, kl_pol, kl_mf)
    best_ref = ["bare", "polaron", "meanforce"][np.argmin([kl_bare, kl_pol, kl_mf])]
    print(f"\n[Verdict section 5]")
    print(f"  Best reference: rho_{best_ref}  (KL={best_kl:.4f})")
    print(f"  Threshold:      KL < {threshold}")
    print(f"  Status:         {'PASS' if best_kl < threshold else 'FAIL'}")
    
    return {
        "label": label,
        "ss_diagnosis": ss,
        "kl_bare": kl_bare,
        "kl_polaron": kl_pol,
        "kl_meanforce": kl_mf,
        "kl_ref_gap": kl_refs,
        "best_reference": best_ref,
        "best_kl": best_kl,
        "passes_threshold": best_kl < threshold,
    }


if __name__ == "__main__":
    results = {}
    for path, label in [(HEOM_DUMP_1JFF, "1JFF_full"),
                        (HEOM_DUMP_6DPU, "6DPU_frag")]:
        if not path.exists():
            print(f"[skip] {path} no existe")
            continue
        results[label] = diagnose(path, label)
    
    with open("meanforce_diagnosis.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nwrote meanforce_diagnosis.json")
