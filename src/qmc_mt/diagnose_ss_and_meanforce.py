"""
diagnose_ss_and_meanforce.py
Diagnostico de convergencia Gibbs/MeanForce.
VERSION CORREGIDA v2 -- fixes:
  [1] bare_gibbs: qt.qeye(H_S.dims[0]) en lugar de qt.qeye(H_S.dims[0][0])
  [2] mean_force_gibbs_2nd_order: reorganizacion lambda no puede ser identidad
      si Q es proyector de sitio (Q^2 = Q, no I) -> OK; pero si Q = sigma_z
      Q^2 = I -> el shift es absorbido por E0. Se aplica la correccion directamente
      sobre los eigenvalores para evitar cancelacion por E0.
  [3] kl_quantum: formula numericamente estable via base comun + max(0,.)
  [4] Diagnostico de lam_list y coupling_ops impreso para detectar valores nulos.
"""

import numpy as np
import qutip as qt
import sys
import pickle, json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

# Conversiones Tier-0
CM_TO_RADFS = 2 * np.pi * 2.9979e-5
KT_CM       = 208.5
KT_RADFS    = KT_CM * CM_TO_RADFS
BETA_RADFS  = 1.0 / KT_RADFS
GAM_RADFS   = 53.0 * CM_TO_RADFS

# Paths
HEOM_DUMP_1JFF = PROJECT_ROOT / "outputs_data" / "raw_pkl" / "heom_1JFF.pkl"
HEOM_DUMP_6DPU = PROJECT_ROOT / "outputs_data" / "raw_pkl" / "heom_6DPU.pkl"


def _make_eye(H):
    """
    Construye la identidad compatible con la estructura tensorial de H.
    qt.qeye(H.dims[0]) crea el producto tensorial correcto.
    qt.qeye(H.dims[0][0]) solo tomaba la primera subdimension -- BUG ORIGINAL.
    """
    return qt.qeye(H.dims[0])


def bare_gibbs(H_S, beta=BETA_RADFS):
    """Estado de Gibbs del sistema desnudo (sin acoplamiento al bano)."""
    E = H_S.eigenenergies()
    E0 = E.min()
    # CORRECCION [1]: usar _make_eye en lugar de qt.qeye(H_S.dims[0][0])
    op = (-beta * (H_S - E0 * _make_eye(H_S))).expm()
    return op / op.tr()


def mean_force_gibbs_2nd_order(H_S, coupling_ops, lam_list, beta=BETA_RADFS, gam=GAM_RADFS):
    """
    Estado de Gibbs de fuerza media a segundo orden (Gnanasekaran-Moix).

    H_MF = H_S - sum_i lambda_i * Q_i^2
                - sum_i (lambda_i/beta) * f_corr * [Q_i, [Q_i, H_S]]

    NOTA CRITICA [2]:
    Si Q_i = sigma_z (unico sitio, operador de posicion adimensional),
    entonces Q_i^2 = I_S, y el primer termino es un shift de identidad
    COMPLETAMENTE ABSORBIDO por la substraccion E0 en bare_gibbs.
    En ese caso la primera suma no modifica rho.
    Para sistemas multisitio con Q_i = |i><i| (proyector), Q_i^2 = Q_i != I,
    y la correccion SI modifica el espectro de H_MF de forma no trivial.

    Para evitar la cancelacion silenciosa por E0, calculamos el estado
    de Gibbs directamente sobre el espectro corregido de H_MF, donde
    la substraccion E0 se aplica DESPUES de incorporar todas las correcciones.
    Esto es equivalente a bare_gibbs(H_MF) pero se hace explicitamente.
    Esto es equivalente a bare_gibbs(H_MF) pero se hace explicitamente.

    Ademas se imprime el diagnostico de la magnitud de las correcciones.
    """
    H_MF = H_S.copy()

    # --- Termino de primer orden: reorganizacion ---
    # lambda_i * Q_i^2: desplaza energias de sitio proporcionalmente a <Q_i^2>
    corr1_norm = 0.0
    for Q, lam in zip(coupling_ops, lam_list):
        dH = lam * (Q * Q)
        corr1_norm += dH.norm()
        H_MF -= dH

    # --- Factor de segundo orden (Drude-Lorentz, aprox. adiabatic) ---
    # f = beta^2 / (12 * (1 + beta*gam/2))
    f_corr = (beta**2) / (12.0 * (1.0 + beta * gam / 2.0))

    # --- Termino de segundo orden: commutador anidado ---
    corr2_norm = 0.0
    for Q, lam in zip(coupling_ops, lam_list):
        # [Q, [Q, H_S]] -- usa H_S original, no H_MF
        inner = qt.commutator(Q, H_S)
        nested = qt.commutator(Q, inner)
        dH = (lam / beta) * f_corr * nested
        corr2_norm += dH.norm()
        H_MF -= dH

    # Diagnostico de magnitudes
    print(f"    [MF diag] ||corr1|| = {corr1_norm:.4e}  ||corr2|| = {corr2_norm:.4e}  "
          f"||H_S|| = {H_S.norm():.4e}  f_corr = {f_corr:.4e}")

    if corr1_norm < 1e-12 and corr2_norm < 1e-12:
        print("    [MF WARN] Ambas correcciones son numericamente nulas.")
        print("              Posibles causas:")
        print("              (a) lam_per_site ~ 0 en el pickle (unidades incorrectas)")
        print("              (b) Q^2 = I (Q=sigma_z) -> corr1 absorbida por E0 en bare_gibbs")
        print("              (c) [Q,[Q,H_S]] = 0 por simetria del Hamiltoniano")

    E_MF = H_MF.eigenenergies()
    E0_MF = E_MF.min()
    op = (-beta * (H_MF - E0_MF * _make_eye(H_MF))).expm()
    return op / op.tr()


def kl_quantum(rho, sigma):
    """
    KL(rho || sigma) = Tr[rho (log rho - log sigma)] >= 0.

    CORRECCION [3]: formula numericamente estable.
    Se construyen los logaritmos matriciales via diagonalizacion individual,
    pero el traceado final se hace en representacion matricial densa,
    lo cual es correcto. Se agrega max(0,.) para evitar -0.0000 por
    errores de punto flotante cuando rho ~ sigma.

    Alternativa robusta via vectorizacion de eigenvalores:
    Si rho y sigma son diagonales en la misma base, KL = sum p_i log(p_i/q_i).
    Para el caso general (bases distintas), la formula matricial es correcta.
    """
    rho_arr = rho.full()
    sig_arr = sigma.full()

    # Eigendescomposicion de rho
    er, Vr = np.linalg.eigh(rho_arr)
    er = np.clip(er, 1e-15, None)
    er /= er.sum()  # renormalizar para precision numerica

    # Eigendescomposicion de sigma
    es, Vs = np.linalg.eigh(sig_arr)
    es = np.clip(es, 1e-15, None)
    es /= es.sum()

    # Logaritmos matriciales
    log_rho   = Vr @ np.diag(np.log(er)) @ Vr.conj().T
    log_sigma = Vs @ np.diag(np.log(es)) @ Vs.conj().T

    # KL = Tr[rho (log_rho - log_sigma)]
    kl = np.real(np.trace(rho_arr @ (log_rho - log_sigma)))

    # max(0,.) para evitar -epsilon por punto flotante cuando rho ~ sigma
    return float(max(kl, 0.0))


def diagnose(path, label):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")

    with open(path, "rb") as f:
        d = pickle.load(f)

    tlist       = d["tlist"]
    rho_t       = d["rho_t"]
    H_S         = d["H_S"]
    S_ops       = d["coupling_ops"]
    lams        = d["lam_per_site"]

    # --- Diagnostico de contenido del pickle [4] ---
    print(f"\n  [Pickle diag]")
    print(f"  H_S dims:      {H_S.dims}  shape: {H_S.shape}")
    print(f"  H_S norm:      {H_S.norm():.4e} rad/fs")
    print(f"  n_coupling_ops:{len(S_ops)}")
    print(f"  lam_per_site:  {np.array(lams)}")
    lam_arr = np.array(lams)
    if np.all(np.abs(lam_arr) < 1e-12):
        print("  [WARN] lam_per_site son todos ~0 -> sin correccion de fuerza media.")
        print("         Verificar unidades: deben estar en rad/fs, no en cm^-1.")
        print(f"         Si estan en cm^-1, multiplica por CM_TO_RADFS = {CM_TO_RADFS:.4e}")
    for i, (Q, lam) in enumerate(zip(S_ops, lams)):
        Q2 = Q * Q
        # Grado de desviacion de Q^2 respecto a la identidad.
        eye_n = _make_eye(Q)
        diff_from_id = (Q2 - lam * eye_n if abs(lam) > 0 else Q2).norm()
        print(f"  Q[{i}] norm={Q.norm():.3e}  lam={lam:.4e}  ||Q^2||={Q2.norm():.3e}  "
              f"||Q^2 - I||={( Q2 - eye_n ).norm():.3e}")

    rho_end  = rho_t[-1]
    rho_bare = bare_gibbs(H_S)
    print(f"\n  [Computing rho_meanforce_2nd_order...]")
    rho_mf   = mean_force_gibbs_2nd_order(H_S, S_ops, lams)

    kl_bare = kl_quantum(rho_end, rho_bare)
    kl_mf   = kl_quantum(rho_end, rho_mf)
    kl_bm   = kl_quantum(rho_bare, rho_mf)

    # Drift (norma Frobenius de la derivada discreta al final)
    dt   = tlist[-1] - tlist[-10]
    delta = (rho_t[-1] - rho_t[-10]).norm() / dt if abs(dt) > 0 else float("nan")

    print(f"\n  [SS convergence]")
    print(f"  drift_rate_per_fs = {delta:.4e}")
    print(f"  t_end_fs          = {tlist[-1]:.1f}")
    verdict_ss = "STEADY" if delta < 1e-5 else "NOT_STEADY"
    print(f"  verdict           = {verdict_ss}")

    print(f"\n  [KL divergences at t={tlist[-1]:.0f} fs]")
    print(f"  KL(rho_HEOM || rho_bare_Gibbs)     = {kl_bare:.4f} nats")
    print(f"  KL(rho_HEOM || rho_meanforce_2nd)  = {kl_mf:.4f} nats")
    print(f"  KL(rho_bare || rho_meanforce_2nd)  = {kl_bm:.4f} nats")

    if kl_bm < 1e-6:
        print("  [WARN] rho_bare == rho_meanforce: la correccion MF es nula.")
        print("         Ver diagnostico de lam y Q^2 arriba.")

    best_ref = "rho_bare" if kl_bare <= kl_mf else "rho_meanforce_2nd"
    kl_best  = min(kl_bare, kl_mf)
    threshold = 0.05
    status = "PASS" if kl_best < threshold else "FAIL"

    print(f"\n  [Verdict]")
    print(f"  Best reference: {best_ref}  (KL={kl_best:.4f})")
    print(f"  Threshold:      KL < {threshold}")
    print(f"  Status:         {status}")

    return {
        "label":    label,
        "kl_bare":  kl_bare,
        "kl_mf":    kl_mf,
        "kl_bare_vs_mf": kl_bm,
        "drift":    float(delta),
        "ss_verdict": verdict_ss,
        "verdict":  status,
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
    res = {}
    for p, l in [(HEOM_DUMP_1JFF, "1JFF"), (HEOM_DUMP_6DPU, "6DPU")]:
        if p.exists():
            res[l] = diagnose(p, l)
        else:
            print(f"\n[SKIP] No encontrado: {p}")

    out = PROJECT_ROOT / "outputs_data" / "raw_json" / "metrics" / "meanforce_diagnosis.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nOK. Resultados guardados en: {out.resolve()}")