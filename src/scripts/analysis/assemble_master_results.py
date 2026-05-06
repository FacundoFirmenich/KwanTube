#!/usr/bin/env python3
"""
assemble_master_results.py
Post-processing driver for HEOM production sequences.

This script performs the temporal concatenation of discretized HEOM propagation
windows into a single continuous master array. It enforces data integrity 
verification via SHA-256 sidecar validation prior to assembly, ensuring 
the final dataset is cryptographically sound and suitable for publication-grade 
figure generation (e.g., time-dependent coherence landscapes).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Cryptographic Integrity Verification
# ---------------------------------------------------------------------------
def compute_sha256(file_path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Generates a SHA-256 cryptographic hash for robust file integrity checks."""
    digest = hashlib.sha256()
    with file_path.open('rb') as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()

def verify_file_integrity(file_path: Path) -> bool:
    """
    Validates a data file against its accompanying SHA-256 sidecar.
    
    Returns:
        True if the file matches its sidecar, False if the sidecar is 
        missing or if a hash mismatch is detected.
    """
    sidecar_path = file_path.with_suffix(file_path.suffix + ".sha256")
    if not sidecar_path.exists():
        return False
        
    with open(sidecar_path, 'r', encoding='utf-8') as f:
        expected_hash = f.read().split()[0]
        
    actual_hash = compute_sha256(file_path)
    
    if actual_hash != expected_hash:
        print(f"[INTEGRITY FAILURE] Hash mismatch detected for {file_path.name}.", file=sys.stderr)
        print(f"  Expected: {expected_hash}", file=sys.stderr)
        print(f"  Actual:   {actual_hash}", file=sys.stderr)
        return False
        
    return True

# ---------------------------------------------------------------------------
# Temporal Concatenation Engine
# ---------------------------------------------------------------------------
def assemble_windows(windows_dir: Path, expected_window_count: int, tolerance_fs: float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    """
    Concatenates sequential propagation windows into a single time series.
    
    Overlapping time-points at window boundaries are identified via numerical 
    tolerance and removed to prevent discontinuities in the master array.
    
    Args:
        windows_dir: Directory containing window_XXX.npz files.
        expected_window_count: Exact number of windows expected (from manifest).
        tolerance_fs: Numerical tolerance for identifying overlapping time points.
        
    Returns:
        Tuple of (tlist_master, populations_master) numpy arrays.
    """
    t_master_list = []
    pop_master_list = []
    
    processed_count = 0
    
    for i in range(expected_window_count):
        window_path = windows_dir / f"window_{i:03d}.npz"
        
        if not window_path.exists():
            raise FileNotFoundError(
                f"Sequence integrity broken: Expected {expected_window_count} windows, "
                f"but {window_path.name} is missing."
            )
            
        # 1. Verify cryptographic integrity before loading
        if not verify_file_integrity(window_path):
            raise ValueError(
                f"Aborting assembly: Data corruption detected in {window_path.name}. "
                f"Verify system stability or re-run the production driver."
            )
        
        # 2. Load data (restrict deserialization to raw numerical types)
        data = np.load(window_path, allow_pickle=False)
        t_window = data['tlist']
        pop_window = data['populations']
        
        # 3. Handle temporal overlap between consecutive windows
        if processed_count > 0:
            time_gap = abs(t_window[0] - t_master_list[-1][-1])
            if time_gap < tolerance_fs:
                # Overlap detected; discard the first point of the new window
                t_master_list.append(t_window[1:])
                pop_master_list.append(pop_window[:, 1:])
            else:
                # Gap detected; log a warning but append raw to prevent data loss
                print(f"[WARNING] Temporal gap detected at window {i:03d} "
                      f"({time_gap:.2f} fs). Appending raw data.", file=sys.stderr)
                t_master_list.append(t_window)
                pop_master_list.append(pop_window)
        else:
            t_master_list.append(t_window)
            pop_master_list.append(pop_window)
            
        processed_count += 1
        # Print progress indicator for long runs
        if processed_count % 10 == 0 or processed_count == expected_window_count:
            print(f"[PROGRESS] Assembled {processed_count}/{expected_window_count} windows...")

    t_master = np.concatenate(t_master_list)
    pop_master = np.concatenate(pop_master_list, axis=1)
    
    return t_master, pop_master

# ---------------------------------------------------------------------------
# Main Execution Block
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Assembles HEOM window segments into a single master dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--run-tag", 
        type=str, 
        default="nc7_nk1_i5_survival",
        help="Identifier of the production run (default: nc7_nk1_i5_survival)."
    )
    parser.add_argument(
        "--project-root", 
        type=str, 
        default=None,
        help="Override base directory of the repository."
    )
    args = parser.parse_args()

    start_time = time.time()
    
    # Resolve paths
    project_root = resolve_project_root(args.project_root)
    run_dir = project_root / "outputs_data" / "production" / f"heom_prod_{args.run_tag}"
    windows_dir = run_dir / "windows"
    manifest_path = run_dir / "progress_manifest.json"
    
    # Canonical output location for the master dataset (inside KwanTube)
    output_dir = project_root / "outputs_data" / "raw_npz"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "master_results.npz"
    
    # Validate environment
    if not manifest_path.exists():
        print(f"[ERROR] Manifest file not found at {manifest_path}. "
              f"Cannot verify expected window count.", file=sys.stderr)
        sys.exit(1)
        
    if not windows_dir.exists():
        print(f"[ERROR] Windows directory not found at {windows_dir}.", file=sys.stderr)
        sys.exit(1)

    print("="*70)
    print("HEOM MASTER DATASET ASSEMBLY PROTOCOL")
    print("="*70)
    
    # Load metadata to enforce strict sequence matching
    with open(manifest_path, 'r', encoding='utf-8') as f:
        manifest = json.load(f)
        
    expected_windows = int(manifest.get('window_index', 0))
    target_fs = manifest.get('target_fs', 0)
    
    print(f"Run Tag:           {args.run_tag}")
    print(f"Expected Windows:  {expected_windows}")
    print(f"Target Time (fs):  {target_fs}")
    print("-" * 70)

    try:
        # Execute core assembly
        t_final, pop_final = assemble_windows(windows_dir, expected_windows)
        
        # Final persistence
        np.savez_compressed(output_path, tlist=t_final, populations=pop_final)
        
        # Cryptographically seal the master dataset
        master_hash = compute_sha256(output_path)
        sidecar_path = output_path.with_suffix(".npz.sha256")
        sidecar_path.write_text(f"{master_hash}  {output_path.name}\n", encoding="utf-8")
        
        elapsed_time = time.time() - start_time
        
        print("-" * 70)
        print("[SUCCESS] Assembly sequence completed without errors.")
        print(f"Output File:      {output_path.name}")
        print(f"SHA-256 Hash:      {master_hash}")
        print(f"Time Axis Shape:  {t_final.shape}")
        print(f"Population Shape: {pop_final.shape}")
        print(f"Final Time (fs):  {t_final[-1]:.1f}")
        print(f"Execution Time:   {elapsed_time:.2f} seconds")
        print("="*70)
        print("Master dataset is now ready for downstream figure generation.\n")
        
    except Exception as e:
        print(f"\n[FATAL ERROR] Assembly protocol aborted.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        sys.exit(1)

# Helper to resolve project root if run standalone outside repo context
def resolve_project_root(explicit: str | None) -> Path:
    if explicit:
        return Path(explicit).expanduser().resolve()
    here = Path(__file__).resolve()
    # Adjusted to parents[3] to reach KwanTube/ from src/scripts/analysis/
    root = here.parents[3]
    if (root / "outputs_data").exists(): return root
    cwd = Path.cwd().resolve()
    if (cwd / "outputs_data").exists(): return cwd
    return root

if __name__ == "__main__":
    main()