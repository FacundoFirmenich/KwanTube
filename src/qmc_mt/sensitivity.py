"""
Sobol sensitivity analysis for the DecohereceModel parameters.
Uses Saltelli sampling and Jansen estimator for variance decomposition.
"""

import numpy as np
import sys
import json
from pathlib import Path

# Boilerplate for package-level imports resolution
PROJECT_ROOT = Path(__file__).resolve().parents[2] # back from src/qmc_mt/ to root
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.core import (TubulinDimer, ExperimentalParameters, DecoherenceModel, const)


def sobol_indices(n_samples: int = 50000, seed: int = 42):
    """
    Compatibility adapter expected by `reproduce_paper_results.py`.

    Returns a compact dictionary with first/total Sobol indices and basic
    T2* summary statistics, derived from the bootstrap estimator implemented
    in this module.
    """
    # Use 200 bootstrap iterations for Tier-0 reporting
    n_boot = 200
    rep = sobol_indices_bootstrap(n_samples=int(n_samples), n_boot=int(n_boot))

    params = [r["parameter"] for r in rep["results"]]
    s1 = [float(r["S1"]["mean"]) for r in rep["results"]]
    st = [float(r["ST"]["mean"]) for r in rep["results"]]
    s1_ci95 = [list(map(float, r["S1"]["ci95"])) for r in rep["results"]]
    st_ci95 = [list(map(float, r["ST"]["ci95"])) for r in rep["results"]]

    # Internal MC to provide T2* summary used in SI rendering
    rng = np.random.default_rng(seed)
    X = rng.uniform(0, 1, (max(512, n_samples // 2), 4))
    X[:, 0] = 0.1 + 0.9 * X[:, 0]               # eta: 0.1 - 1.0
    X[:, 1] = 1e12 + 9e12 * X[:, 1]             # wc: 1e12 - 1e13
    X[:, 2] = 280.0 + 40.0 * X[:, 2]            # Temp: 280 - 320
    X[:, 3] = 1e-4 * (1e-2 / 1e-4) ** X[:, 3]   # f_prot: log-unif 1e-4 - 1e-2

    dimer = TubulinDimer()
    t2_ps = []
    for eta, wc, temp, f_prot in X:
        p = ExperimentalParameters(temperature=float(temp))
        m = DecoherenceModel(dimer, p, protection_factor=float(f_prot), eta=float(eta), omega_c=float(wc))
        t2_ps.append(1e12 / m.get_all_rates()["total_dephasing"])
    t2_ps = np.asarray(t2_ps, dtype=float)

    return {
        "parameters": params,
        "n_samples": int(rep["n_samples"]),
        "n_boot": int(rep["n_boot"]),
        "confidence_level": 0.95,
        "first_order": s1,
        "first_order_ci95": s1_ci95,
        "total_order": st,
        "total_order_ci95": st_ci95,
        "T2_ps_mean": float(np.mean(t2_ps)),
        "T2_ps_range": [float(np.min(t2_ps)), float(np.max(t2_ps))],
    }

def sobol_indices_bootstrap(n_samples: int = 50000, n_boot: int = 200):
    """
    Executes Sobol sensitivity analysis (Saltelli) with bootstrapping for 95% CI.
    """
    dimer = TubulinDimer()
    rng = np.random.default_rng(42)
    dim = 4
    params_names = ['eta', 'omega_c', 'Temperature', 'f_prot']

    def model(x):
        # Vectorized or batch evaluation of the decoherence model
        y = np.zeros(len(x))
        for i in range(len(x)):
            eta, wc, temp, f_prot = x[i]
            p = ExperimentalParameters(temperature=temp)
            m = DecoherenceModel(dimer, p, protection_factor=f_prot, eta=eta, omega_c=wc)
            y[i] = 1e12 / m.get_all_rates()['total_dephasing'] # T2 in ps
        return y

    # 1. Generate A and B matrices
    A = rng.uniform(0, 1, (n_samples, dim))
    B = rng.uniform(0, 1, (n_samples, dim))
    
    def scale(mat):
        res = mat.copy()
        res[:, 0] = 0.1 + 0.9 * mat[:, 0]               # eta: 0.1 - 1.0
        res[:, 1] = 1e12 + 9e12 * mat[:, 1]            # wc: 1e12 - 1e13
        res[:, 2] = 280.0 + 40.0 * mat[:, 2]           # Temp: 280 - 320
        res[:, 3] = 1e-4 * (1e-2/1e-4)**mat[:, 3]      # f_prot: log-unif 1e-4 - 1e-2
        return res

    A_s, B_s = scale(A), scale(B)
    yA, yB = model(A_s), model(B_s)
    
    # Ci matrices for each parameter
    yCi = []
    for i in range(dim):
        Ci = A_s.copy()
        Ci[:, i] = B_s[:, i]
        yCi.append(model(Ci))

    def compute_indices(ya, yb, yci_list):
        var_y = np.var(np.concatenate([ya, yb]))
        s1, st = [], []
        for i in range(dim):
            # Saltelli/Jansen estimator
            s1_i = (np.mean(yb * (yci_list[i] - ya))) / var_y
            st_i = (np.mean((ya - yci_list[i])**2)) / (2 * var_y)
            s1.append(s1_i)
            st.append(st_i)
        return np.array(s1), np.array(st)

    # 2. Bootstrap for confidence intervals (95%)
    boot_s1 = np.zeros((n_boot, dim))
    boot_st = np.zeros((n_boot, dim))
    
    all_idx = np.arange(n_samples)
    for b in range(n_boot):
        idx = rng.choice(all_idx, size=n_samples, replace=True)
        s1, st = compute_indices(yA[idx], yB[idx], [y[idx] for y in yCi])
        boot_s1[b], boot_st[b] = s1, st

    # 3. Consolidate results
    s1_mu = np.mean(boot_s1, axis=0)
    s1_lo = np.percentile(boot_s1, 2.5, axis=0)
    s1_hi = np.percentile(boot_s1, 97.5, axis=0)
    
    st_mu = np.mean(boot_st, axis=0)
    st_lo = np.percentile(boot_st, 2.5, axis=0)
    st_hi = np.percentile(boot_st, 97.5, axis=0)

    report = []
    for i in range(dim):
        report.append({
            "parameter": params_names[i],
            "S1": {"mean": float(s1_mu[i]), "ci95": [float(s1_lo[i]), float(s1_hi[i])]},
            "ST": {"mean": float(st_mu[i]), "ci95": [float(st_lo[i]), float(st_hi[i])]}
        })

    return {
        "n_samples": n_samples,
        "n_boot": n_boot,
        "results": report,
        "is_one_parameter_dominated": s1_mu[0] > 0.9 and all(s < 0.05 for s in s1_mu[1:])
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
    print(f"Starting high-precision Sobol analysis (Saltelli + Bootstrap n=200)...")
    summary = sobol_indices_bootstrap(n_samples=50000, n_boot=200) 
    
    print("\n=== SENSITIVITY REPORT (95% CI) ===")
    for p in summary["results"]:
        print(f"{p['parameter']:12}: S1={p['S1']['mean']:.4f} [{p['S1']['ci95'][0]:.4f}, {p['S1']['ci95'][1]:.4f}] | "
              f"ST={p['ST']['mean']:.4f} [{p['ST']['ci95'][0]:.4f}, {p['ST']['ci95'][1]:.4f}]")
    
    if summary["is_one_parameter_dominated"]:
        print("\nVERDICT: The model is effectively one-parameter-dominated (eta).")
    
    out_path = PROJECT_ROOT / "outputs_data" / "raw_json" / "metrics" / "sensitivity_sobol_final.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults saved to: {out_path}")
