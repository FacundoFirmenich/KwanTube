"""
recompute_tau_coh.py -- Precise tau_coherence calculation from Redfield outputs.

Method: 
1/e decay measured from transient peak, reporting non-decaying cases explicitly.
This avoids artefacts from initial zero-coherence states.
"""
import json
import numpy as np
from pathlib import Path

ROOT = Path(r"C:\Users\User\3D Objects\biofisicaquantiqaCLINE")


def tau_coh_from_peak(t, coh):
    """
    Calculates coherence time by finding the transient peak and measuring
    the 1/e decay time from that point onwards.
    """
    coh = np.asarray(coh)
    peak_idx = int(np.argmax(coh))
    peak_val = float(coh[peak_idx])
    if peak_val < 1e-10:
        return None, "no_coherence_generated", 0.0, 0.0
    if peak_idx >= len(t) - 5:
        return None, "peak_at_end_of_window", float(t[peak_idx]), peak_val
    after = coh[peak_idx:]
    t_after = np.asarray(t)[peak_idx:]
    cross = np.where(after <= peak_val / np.e)[0]
    if len(cross):
        tau = float(t_after[cross[0]] - t[peak_idx])
        return tau, "measured", float(t[peak_idx]), peak_val
    return None, "no_decay_in_window", float(t[peak_idx]), peak_val


def main():
    updates = {}
    for pid in ("1JFF", "6DPU"):
        f = ROOT / f"redfield_{pid}.npz"
        if not f.exists():
            print(f"[{pid}] skip -- {f.name} not found")
            continue
        d = np.load(f, allow_pickle=True)
        t = d["t_fs"]
        coh = d["coh_tot"]
        tau, status, peak_t, peak_v = tau_coh_from_peak(t, coh)
        updates[pid] = dict(
            tau_coh_fs=tau,
            tau_coh_status=status,
            coh_peak_time_fs=peak_t,
            coh_peak_value=peak_v,
            window_fs=float(t[-1]),
        )
        tau_str = f"{tau:8.1f} fs" if tau is not None else "   N/A   "
        print(f"[{pid}] tau_coh = {tau_str}  ({status})   "
              f"peak at {peak_t:7.1f} fs, value={peak_v:.3e}")

    out = ROOT / "redfield_tau_coh.json"
    out.write_text(json.dumps(updates, indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
