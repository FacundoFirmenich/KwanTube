#!/usr/bin/env python3
"""
heom_production_driver.py  — vFinal
Hierarchical Equations of Motion (HEOM) production driver optimized for 
memory-constrained architectures (e.g., consumer-grade hardware). 

Implements strict ODE solver tolerances, Write-Ahead Logging (WAL) fault 
tolerance, and dynamic memory monitoring to ensure data integrity during 
extended-time propagation of the 1JFF tubulin dimer system (Nc=7, Nk=1).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import tempfile
import threading
from pathlib import Path

import numpy as np
import qutip as qt
from qutip.solver.heom import (
    DrudeLorentzPadeBath, 
    HEOMSolver,
)

try:
    from qutip.solver.heom import HierarchyADOsState as _HierarchyState
except ImportError:  # QuTiP < 5 compatibility
    try:
        from qutip.solver.heom import HierarchyADOState as _HierarchyState
    except ImportError:
        _HierarchyState = object

# Optional dependency for system memory monitoring
try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False

# ---------------------------------------------------------------------------
# Physical constants  (cm⁻¹ → rad fs⁻¹)
# ---------------------------------------------------------------------------
CM_TO_RADFS: float = 2.0 * np.pi * 2.9979e-5
LAM_RADFS:   float = 35.0  * CM_TO_RADFS
GAM_RADFS:   float = 53.0  * CM_TO_RADFS
T_RADFS:     float = 300.0 * 0.69503 * CM_TO_RADFS

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
def resolve_project_root(explicit: str | None = None) -> Path:
    """Resolves the base directory of the repository."""
    if explicit:
        return Path(explicit).expanduser().resolve()
    here = Path(__file__).resolve()
    root = here.parents[3]
    if (root / "outputs_data").exists(): return root
    cwd = Path.cwd().resolve()
    if (cwd / "outputs_data").exists(): return cwd
    return root

# ---------------------------------------------------------------------------
# System loader
# ---------------------------------------------------------------------------
def load_1jff(project_root: Path) -> tuple[qt.Qobj, list[str]]:
    """Loads the 1JFF Hamiltonian and site labels from archived NPZ data."""
    npz_path = project_root / "outputs_data" / "raw_npz" / "H_1JFF.npz"
    if not npz_path.exists():
        raise FileNotFoundError(f"H_1JFF.npz not found in: {npz_path}")
    data   = np.load(npz_path, allow_pickle=True)
    H_cm   = data["H_cm1"]
    labels = list(data["labels"])
    # Subtract mean trace to remove global energy offset (gauge freedom)
    H_rad  = (H_cm - np.mean(np.diag(H_cm)) * np.eye(H_cm.shape[0])) * CM_TO_RADFS
    return qt.Qobj(H_rad), labels

# ---------------------------------------------------------------------------
# Solver builder
# ---------------------------------------------------------------------------
def build_solver(
    H_S: qt.Qobj, Q_ops: list[qt.Qobj], nk: int, nc: int, 
    store_ados: bool, atol: float, rtol: float
) -> HEOMSolver:
    """Constructs the HEOM solver with strictly enforced numerical tolerances."""
    baths = [
        DrudeLorentzPadeBath(Q, lam=LAM_RADFS, gamma=GAM_RADFS, T=T_RADFS, Nk=nk)
        for Q in Q_ops
    ]
    return HEOMSolver(
        H_S, baths, max_depth=nc,
        options={
            "nsteps": 200_000, 
            "store_ados": True,
            "store_states": store_ados,
            "store_final_state": True,
            "progress_bar": False,
            "atol": atol,     # Required for convergence verification
            "rtol": rtol,     # Required for stability in intermediate coupling
        },
    )

# ---------------------------------------------------------------------------
# Crash-safe I/O (WAL) and Serialization
# ---------------------------------------------------------------------------
def _atomic_replace(path: Path, writer_fn) -> None:
    """Executes writer function to a temporary file, then atomically renames."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=f".{path.stem}.", suffix=path.suffix)
    try:
        os.close(fd)
        writer_fn(Path(tmp))
        os.replace(tmp, path)
    except Exception:
        if os.path.exists(tmp): os.unlink(tmp)
        raise

def save_checkpoint(path: Path, state, meta: dict) -> None:
    """Serializes HEOM state and metadata into separate binary and JSON files."""
    state_path = path.with_suffix(".ados.npz")
    meta_path = path.with_suffix(".meta.json")
    
    # Compatibility layer for QuTiP ADO serialization
    if isinstance(state, qt.Qobj):
        save_npz_atomic(state_path, rho=np.asarray(state.full()))
    else:
        ado_payload = getattr(state, "_ado_state", None)
        rho_payload = state.rho.full() if hasattr(state, "rho") else None
        if ado_payload is None or rho_payload is None:
            raise TypeError(f"Unsupported HEOM checkpoint state type: {type(state)!r}")
        save_npz_atomic(state_path, ado_state=np.asarray(ado_payload), rho=np.asarray(rho_payload))
        
    _atomic_replace(meta_path, lambda p: p.write_text(json.dumps(meta, indent=2), encoding="utf-8"))
    write_sha256_sidecar(meta_path)

def load_checkpoint(path: Path, solver: HEOMSolver) -> tuple[qt.Qobj, dict]:
    """Restores system state from checkpoint files."""
    state_path = path.with_suffix(".ados.npz")
    meta_path = path.with_suffix(".meta.json")
    if not state_path.exists() or not meta_path.exists():
        raise FileNotFoundError("Checkpoint data or metadata file not found.")
    with open(meta_path, "r") as f: meta = json.load(f)
    data = np.load(state_path, allow_pickle=False)
    rho = qt.Qobj(data["rho"])
    return rho, meta

def save_npz_atomic(path: Path, **arrays) -> None:
    """Writes NumPy arrays to disk using atomic file replacement."""
    def _write_npz(p: Path) -> None:
        with p.open("wb") as handle:
            np.savez_compressed(handle, **arrays)
    _atomic_replace(path, _write_npz)
    write_sha256_sidecar(path)

def compute_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Generates a SHA-256 cryptographic hash for file integrity verification."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()

def write_sha256_sidecar(path: Path) -> None:
    """Writes a checksum sidecar file for validation pipelines."""
    sidecar = path.with_suffix(path.suffix + ".sha256")
    sidecar.write_text(f"{compute_sha256(path)}  {path.name}\n", encoding="utf-8")

# ---------------------------------------------------------------------------
# Hardware Monitoring
# ---------------------------------------------------------------------------
def check_ram_safety(limit_gb: float) -> None:
    """Terminates the process safely if available system memory is below threshold."""
    if not HAS_PSUTIL:
        return
    available_gb = psutil.virtual_memory().available / (1024 ** 3)
    if available_gb < limit_gb:
        raise MemoryError(
            f"MEMORY THRESHOLD BREACHED: Available system RAM ({available_gb:.2f} GB) "
            f"is below the defined safety limit ({limit_gb:.2f} GB). "
            f"Process terminated gracefully to prevent kernel-level out-of-memory termination."
        )

# ---------------------------------------------------------------------------
# Quantum Metrics
# ---------------------------------------------------------------------------
def make_tlist(t_start: float, t_end: float, dt: float) -> np.ndarray:
    """Generates a uniformly spaced time array to prevent floating-point drift."""
    n = max(1, round((t_end - t_start) / dt))
    return np.linspace(t_start, t_end, n + 1)

def purity_lower_bound(pop: np.ndarray) -> np.ndarray:
    """Calculates Tr(ρ²) strictly from diagonal populations."""
    return np.sum(pop ** 2, axis=0)

def participation_ratio(pop: np.ndarray) -> np.ndarray:
    """Calculates the effective number of populated states."""
    return 1.0 / np.maximum(np.sum(pop ** 2, axis=0), 1e-15)

def vn_entropy_diag(pop: np.ndarray) -> np.ndarray:
    """Calculates von Neumann entropy lower bound from populations."""
    p = np.clip(pop, 1e-15, 1.0)
    return -np.sum(p * np.log(p), axis=0)

# ---------------------------------------------------------------------------
# Argument Parser
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="HEOM Production Driver for Memory-Constrained Environments", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--project-root", type=str, default=None)
    p.add_argument("--nc", type=int, default=7)
    p.add_argument("--nk", type=int, default=1)
    
    p.add_argument("--atol", type=float, default=1e-10, help="Absolute ODE solver tolerance")
    p.add_argument("--rtol", type=float, default=1e-8,  help="Relative ODE solver tolerance")
    
    p.add_argument("--target-ps", type=float, default=30.0)
    p.add_argument("--window-ps", type=float, default=1.0, help="Integration chunk size (ps). Reduced to minimize peak RAM footprint.") 
    p.add_argument("--dt-fs", type=float, default=10.0)
    p.add_argument("--start-site-label", type=str, default="B:103")
    p.add_argument("--run-tag", type=str, default="nc7_nk1_i5_survival")
    
    p.add_argument("--store-ados", action="store_true", 
                   help="Retain complete ADO matrices. Substantially increases RAM requirements.")
    p.add_argument("--resume", action="store_true", help="Resume execution from existing checkpoint.")
    p.add_argument("--non-interactive", action="store_true", help="Suppress interactive prompts between windows.")
    p.add_argument("--prompt-timeout-sec", type=int, default=60, help="Timeout for interactive hold.")
    p.add_argument("--max-windows-per-run", type=int, default=0, help="Limit discrete execution length (0=unlimited).")
    p.add_argument("--telemetry-sec", type=int, default=30, help="Interval for periodic system state logging.")
    
    p.add_argument("--min-ram-gb", type=float, default=3.0, help="Minimum free RAM required to proceed (GB).")
    
    p.add_argument("--save-coherences", action="store_true", help="Evaluate off-diagonal density matrix elements.")
    p.add_argument("--rho-snapshot-every", type=int, default=0, help="Frequency of full density matrix archival.")
    return p.parse_args()

def get_final_ado_state(result):
    """Extracts the terminal ADO state across QuTiP API versions."""
    if hasattr(result, "final_ado_state"):
        state = result.final_ado_state
        if state is not None: return state
    if hasattr(result, "final_ados_state"):
        state = result.final_ados_state
        if state is not None: return state
    if hasattr(result, "ado_states") and result.ado_states:
        return result.ado_states[-1]
    available = ", ".join(name for name in dir(result) if "ado" in name.lower())
    raise AttributeError(f"HEOMResult lacks a valid terminal ADO state. Detected attributes: {available}")

def get_system_telemetry() -> str:
    """Formats current system memory status for logging."""
    if not HAS_PSUTIL:
        return "RAM_TELEMETRY=UNAVAILABLE"
    vm = psutil.virtual_memory()
    return f"RAM_AVAILABLE={vm.available / (1024 ** 3):.2f}GB USED={vm.percent:.1f}%"

def start_telemetry_thread(log, stop_event: threading.Event, interval_sec: int, context_fn) -> threading.Thread | None:
    """Initiates a background thread for periodic hardware state logging."""
    if interval_sec <= 0:
        return None
    def _loop() -> None:
        while not stop_event.wait(interval_sec):
            log(f"PERIODIC TELEMETRY | {context_fn()} | {get_system_telemetry()}")
    thread = threading.Thread(target=_loop, daemon=True)
    thread.start()
    return thread

def prompt_continue_with_timeout(timeout_sec: int, log) -> bool:
    """Handles interactive continuation prompts with a hard timeout."""
    if timeout_sec <= 0:
        log("Interactive hold skipped (timeout <= 0). Defaulting to continuation.")
        return True
    prompt = f"Proceed with next integration window? [Y/n] (timeout {timeout_sec}s, default Y): "
    print(prompt, end="", flush=True)
    
    # OS-specific non-blocking input handling
    if os.name == "nt":
        import msvcrt
        chars: list[str] = []
        deadline = time.time() + timeout_sec
        while time.time() < deadline:
            if msvcrt.kbhit():
                ch = msvcrt.getwch()
                if ch in ("\r", "\n"): break
                chars.append(ch)
                print(ch, end="", flush=True)
            time.sleep(0.05)
        print()
        answer = "".join(chars).strip().lower()
    else:
        import select
        readable, _, _ = select.select([sys.stdin], [], [], timeout_sec)
        answer = sys.stdin.readline().strip().lower() if readable else ""
        
    if answer in {"n", "no"}:
        log("Interactive continuation halted by user input at checkpoint. Process terminating cleanly.")
        return False
    if answer == "":
        log(f"No input received within {timeout_sec}s timeout. Defaulting to continuation.")
    else:
        log(f"User input confirmed: '{answer}'. Proceeding.")
    return True

def write_json_atomic(path: Path, payload: dict) -> None:
    _atomic_replace(path, lambda p: p.write_text(json.dumps(payload, indent=2), encoding="utf-8"))
    write_sha256_sidecar(path)

def write_resume_instructions(path: Path, project_root: Path, args: argparse.Namespace) -> None:
    """Generates a file containing the exact CLI command required to resume execution."""
    cmd = (
        f'python "{project_root / "src" / "scripts" / "heom" / "heom_production_driver.py"}" '
        f'--resume --run-tag {args.run_tag} --nc {args.nc} --nk {args.nk} '
        f'--atol {args.atol} --rtol {args.rtol} --target-ps {args.target_ps} '
        f'--window-ps {args.window_ps} --dt-fs {args.dt_fs}'
    )
    path.write_text("Resume command:\n" + cmd + "\n", encoding="utf-8")

# ---------------------------------------------------------------------------
# Main Execution Loop
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    
    if args.save_coherences or args.rho_snapshot_every > 0:
        print("[CONFIGURATION WARNING] Off-diagonal observables and/or density matrix snapshots enabled. Estimated memory footprint may exceed safe operational margins for constrained hardware.", flush=True)

    project_root = resolve_project_root(args.project_root)
    out_dir     = project_root / "outputs_data" / "production" / f"heom_prod_{args.run_tag}"
    windows_dir = out_dir / "windows"
    metrics_dir = out_dir / "metrics"
    for d in (windows_dir, metrics_dir): d.mkdir(parents=True, exist_ok=True)

    ckpt_base    = out_dir / "checkpoint"
    log_file     = out_dir / "production_progress.log"
    manifest_file = out_dir / "progress_manifest.json"
    resume_file = out_dir / "resume_instructions.txt"

    def log(msg: str):
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
        print(line, flush=True)
        with open(log_file, "a", encoding="utf-8") as f: f.write(line + "\n")

    log(f"=== INITIATING HEOM PROPAGATION SEQUENCE ===")
    write_resume_instructions(resume_file, project_root, args)

    # 1. SYSTEM INITIALIZATION
    H_S, labels = load_1jff(project_root)
    N = len(labels)
    Q_ops = [qt.basis(N, i) * qt.basis(N, i).dag() for i in range(N)]
    
    # Restrict expectation operators to populations to minimize RAM allocation
    e_ops = Q_ops if not args.save_coherences else (
        Q_ops + [qt.basis(N, i) * qt.basis(N, j).dag() for i in range(N) for j in range(i + 1, N)]
    )
    n_pop = N

    solver = build_solver(H_S, Q_ops, nk=args.nk, nc=args.nc, 
                          store_ados=args.store_ados, atol=args.atol, rtol=args.rtol)

    # 2. STATE INITIALIZATION / RESUMPTION
    if ckpt_base.with_suffix(".ados.npz").exists():
        log("Restoring system state from WAL checkpoint...")
        current_state, ckpt_meta = load_checkpoint(ckpt_base, solver)
        t_start_fs = float(ckpt_meta["t_fs"])
        window_index = int(ckpt_meta["window_index"])
        wall_history = ckpt_meta.get("wall_history", [])
        log(f"Checkpoint restored successfully. Resuming propagation at t={t_start_fs:.1f} fs")
    else:
        try: i0 = labels.index(args.start_site_label)
        except ValueError: i0 = 0
        rho0 = qt.basis(N, i0) * qt.basis(N, i0).dag()
        current_state = rho0
        t_start_fs = 0.0
        window_index = 0
        wall_history = []
        log(f"Initializing fresh state: Excited site {labels[i0]}")

    target_fs = args.target_ps * 1000.0
    window_fs = args.window_ps * 1000.0
    total_windows_est = int(np.ceil(max(target_fs - t_start_fs, 0.0) / max(window_fs, 1e-12)))
    windows_this_run = 0

    # 3. INTEGRATION LOOP
    while t_start_fs < target_fs - 1e-12:
        
        check_ram_safety(args.min_ram_gb)

        w_start = t_start_fs
        w_end   = min(target_fs, w_start + window_fs)
        tlist   = make_tlist(w_start, w_end, args.dt_fs)

        log(f"INTEGRATION WINDOW {window_index + 1:03d} | t_span=[{w_start:.1f}, {w_end:.1f}] fs")
        
        t0 = time.perf_counter()
        telemetry_stop = threading.Event()
        context = lambda: f"window={window_index + 1:03d} span=[{w_start:.1f},{w_end:.1f}]fs elapsed={(time.perf_counter() - t0)/60:.2f}min"
        telemetry_thread = start_telemetry_thread(log, telemetry_stop, args.telemetry_sec, context)
        
        try:
            result = solver.run(current_state, tlist, e_ops=e_ops)
        finally:
            telemetry_stop.set()
            if telemetry_thread is not None:
                telemetry_thread.join(timeout=1.0)
                
        wall = time.perf_counter() - t0
        wall_history.append(wall)

        expect_raw = np.array(result.expect)
        pop = np.real(expect_raw[:n_pop])

        # Persist window data
        save_npz_atomic(windows_dir / f"window_{window_index:03d}.npz", tlist=tlist, populations=pop)

        # Calculate and persist physical metrics
        pur = purity_lower_bound(pop)
        pr  = participation_ratio(pop)
        ent = vn_entropy_diag(pop)
        save_npz_atomic(metrics_dir / f"metrics_{window_index:03d}.npz", tlist=tlist, purity=pur, pr=pr, entropy=ent)

        # Update propagation state
        current_state = get_final_ado_state(result)
        t_start_fs = w_end
        window_index += 1
        windows_this_run += 1

        ms_per_fs = wall / max(w_end - w_start, 1e-12) * 1000.0
        eta_min = (target_fs - t_start_fs) / max(w_end - w_start, 1e-12) * float(np.mean(wall_history)) / 60.0
        
        log(f"  Wall time: {wall/60:.2f} min | Speed: {ms_per_fs:.1f} ms/fs | "
            f"Purity: {pur[-1]:.4f} | ETA: {eta_min:.1f} min")

        # Execute WAL Checkpoint
        meta = {
            "t_fs": t_start_fs, "window_index": window_index,
            "nc": args.nc, "nk": args.nk, "atol": args.atol, "rtol": args.rtol,
            "wall_history": wall_history
        }
        save_checkpoint(ckpt_base, current_state, meta)
        
        manifest = {
            **meta,
            "last_completed_window": window_index - 1,
            "target_fs": target_fs,
            "window_fs": window_fs,
            "completed_fraction": float(t_start_fs / max(target_fs, 1e-12)),
            "last_window": {
                "window_npz": str((windows_dir / f"window_{window_index - 1:03d}.npz").resolve()),
                "metrics_npz": str((metrics_dir / f"metrics_{window_index - 1:03d}.npz").resolve()),
                "purity_final": float(pur[-1]),
                "participation_ratio_final": float(pr[-1]),
                "entropy_final": float(ent[-1]),
                "wall_s": float(wall),
            },
            "resume_command_file": str(resume_file.resolve()),
            "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
        write_json_atomic(manifest_file, manifest)
        
        log(
            f"CHECKPOINT SECURED | Completed: {window_index - 1:03d}/{max(window_index + total_windows_est - 1, window_index):03d} | "
            f"Simulated Time: {t_start_fs:.1f}/{target_fs:.1f} fs | "
            f"Purity: {pur[-1]:.6f} | IPR: {pr[-1]:.3f} | S_vN: {ent[-1]:.3f} | "
            f"Mean Wall Time: {np.mean(wall_history)/60:.2f} min | ETA: {eta_min:.1f} min"
        )
        log(f"Resume configuration archived at: {resume_file.resolve()}")

        if args.max_windows_per_run > 0 and windows_this_run >= args.max_windows_per_run:
            log(f"Reached defined execution limit (--max-windows-per-run={args.max_windows_per_run}). Halting cleanly post-checkpoint.")
            break
        if not args.non_interactive and t_start_fs < target_fs - 1e-12:
            if not prompt_continue_with_timeout(args.prompt_timeout_sec, log):
                break

    if t_start_fs >= target_fs - 1e-12:
        log("=== HEOM PROPAGATION SEQUENCE COMPLETE ===")
        log("Proceed to execute the temporal concatenation script to assemble the master results array.")
    else:
        log("=== PROPAGATION SEQUENCE PAUSED ===")
        log(f"To resume execution, refer to instructions at: {resume_file.resolve()}")

if __name__ == "__main__":
    main()