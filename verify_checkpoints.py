import pickle
import numpy as np
import qutip as qt
from pathlib import Path

def verify_checkpoint(path):
    p = Path(path)
    if not p.exists():
        print(f"File not found: {path}")
        return
    
    with open(p, "rb") as f:
        data = pickle.load(f)
    
    nc = data["NC"]
    nk = data["Nk"]
    rho_t = data["rho_t"]
    stats = data["stats"]
    H_S = data["H_S"]
    
    n_sites = H_S.shape[0]
    diag_init = np.real(np.diag(rho_t[0].full()))
    diag_final = np.real(np.diag(rho_t[-1].full()))
    trace_final = rho_t[-1].tr().real
    herm_err = (rho_t[-1] - rho_t[-1].dag()).norm()
    
    print(f"\n--- Checkpoint: {p.name} ---")
    print(f"NC: {nc}, Nk: {nk}")
    print(f"System size: {n_sites}")
    print(f"Initial diagonal: {diag_init}")
    print(f"Final diagonal:   {diag_final}")
    print(f"Final trace:      {trace_final:.8f}")
    print(f"Hermiticity error: {herm_err:.2e}")

if __name__ == "__main__":
    verify_checkpoint("pade_refined_ckpt_NC5_Nk1.pkl")
    verify_checkpoint("pade_refined_ckpt_NC5_Nk2.pkl")
