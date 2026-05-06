#!/usr/bin/env python3
"""Compute collective radiative decay rates for microtubule lattice modes.

The existing lattice diagnostic reports an inverse participation ratio (IPR),
which measures energetic delocalization. Reviewers correctly noted that IPR is
not equivalent to optical subradiance. This script closes that gap by computing
the collective radiative decay matrix using a Lehmberg/Agarwal free-space kernel
and then reporting:

- the radiative-rate expectation value of each excitonic Hamiltonian eigenmode;
- the eigenvalue spectrum of the radiative decay matrix itself;
- IPR-vs-decay-rate scatter data and publication-ready figures.

The model uses the free-space dyadic Green's function.  In the actual MT
geometry (layered dielectric cylinder: inner lumen eps~80, protein wall
eps~2-4, outer water eps~80), local field effects and leaky-wave contributions
would redistribute these rates.  Modes identified as subradiant here therefore
represent an upper bound on protection; true subradiance requires the formalism
of dipole emission in microcylinders.  This caveat is printed in every output.

Usage
-----
    # From the project root (KwanTube/):
    python paper/compute_subradiant_decay_spectrum.py

    # With explicit options:
    python paper/compute_subradiant_decay_spectrum.py \\
        --project-root /path/to/KwanTube \\
        --n-layers 20 --n-protofilaments 13 \\
        --dipole-orientation tangential \\
        --wavelength-nm 280

Notes
-----
- J_axial and J_lateral default reproduce the paper values
  (J_axial=-88.08 meV, J_lateral=160.36 meV).
- The script runs fast enough for N<=520 sites on a laptop; for larger N
  the O(N^2) radiative matrix build dominates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def resolve_project_root(explicit: Optional[str] = None) -> Path:
    """Resolve the KwanTube project root."""
    if explicit:
        return Path(explicit).expanduser().resolve()
    here = Path(__file__).resolve()
    candidates = [
        here.parents[3],
        Path.cwd().resolve(),
        Path.cwd().resolve() / "KwanTube",
    ]
    for candidate in candidates:
        if (candidate / "outputs_data").exists() and (candidate / "src").exists():
            return candidate.resolve()
    return here.parents[3].resolve()


def audit_log(project_root: Path, message: str) -> None:
    """Append a timestamped execution line to the canonical memory log."""
    log_path = (
        project_root
        / "outputs_data"
        / "raw_txt+md"
        / "logs"
        / "execution_memory.log.txt"
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(
            f"[{timestamp}] [compute_subradiant_decay_spectrum.py] {message}\n"
        )


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Compute a SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# Geometry builders
# ---------------------------------------------------------------------------

def infer_radius_nm(
    r_lateral_nm: float,
    seam_rise_nm: float,
    n_protofilaments: int,
) -> float:
    """Infer a cylinder radius matching the requested nearest-neighbour distance."""
    chord = math.sqrt(max(r_lateral_nm**2 - seam_rise_nm**2, 1e-12))
    return chord / (2.0 * math.sin(math.pi / n_protofilaments))


def build_b_lattice_geometry(
    n_layers: int,
    n_protofilaments: int,
    axial_spacing_nm: float,
    seam_rise_nm: float,
    radius_nm: float,
    dipole_orientation: str,
    seed: int,
) -> tuple:
    """Build cylindrical B-lattice positions and unit dipole orientations.

    Returns:
        (positions_m, dipoles, labels) where labels are (protofilament, layer)
        tuples for each site.
    """
    rng = np.random.default_rng(seed)
    positions_nm: list = []
    dipoles: list = []
    labels: list = []

    for p in range(n_protofilaments):
        phi = 2.0 * math.pi * p / n_protofilaments
        radial = np.array([math.cos(phi), math.sin(phi), 0.0], dtype=float)
        tangential = np.array([-math.sin(phi), math.cos(phi), 0.0], dtype=float)
        axial_vec = np.array([0.0, 0.0, 1.0], dtype=float)
        for layer in range(n_layers):
            z_nm = layer * axial_spacing_nm + p * seam_rise_nm / n_protofilaments
            positions_nm.append([
                radius_nm * radial[0],
                radius_nm * radial[1],
                z_nm,
            ])
            if dipole_orientation == "radial":
                mu_hat = radial.copy()
            elif dipole_orientation == "tangential":
                mu_hat = tangential.copy()
            elif dipole_orientation == "axial":
                mu_hat = axial_vec.copy()
            elif dipole_orientation == "helical45":
                mu_hat = (
                    tangential * math.cos(math.pi / 4.0)
                    + axial_vec * math.sin(math.pi / 4.0)
                )
            elif dipole_orientation == "random":
                raw = rng.normal(size=3)
                mu_hat = raw / np.linalg.norm(raw)
            else:
                raise ValueError(
                    f"Unsupported dipole orientation: {dipole_orientation!r}. "
                    "Choose from: radial, tangential, axial, helical45, random."
                )
            norm = np.linalg.norm(mu_hat)
            dipoles.append(mu_hat / max(norm, 1e-15))
            labels.append((p, layer))

    return (
        np.asarray(positions_nm, dtype=float) * 1e-9,
        np.asarray(dipoles, dtype=float),
        labels,
    )


# ---------------------------------------------------------------------------
# Excitonic Hamiltonian
# ---------------------------------------------------------------------------

def build_tight_binding_hamiltonian(
    n_layers: int,
    n_protofilaments: int,
    j_axial_mev: float,
    j_lateral_mev: float,
) -> np.ndarray:
    """Build the nearest-neighbour B-lattice excitonic Hamiltonian.

    Uses the same connectivity as the existing lattice diagnostic.

    Returns:
        (n_sites, n_sites) float array in units of meV.
    """
    n_sites = n_layers * n_protofilaments

    def idx(p: int, layer: int) -> int:
        return p * n_layers + layer

    hamiltonian = np.zeros((n_sites, n_sites), dtype=float)
    for p in range(n_protofilaments):
        for layer in range(n_layers):
            i = idx(p, layer)
            # Axial neighbour
            if layer + 1 < n_layers:
                j = idx(p, layer + 1)
                hamiltonian[i, j] = hamiltonian[j, i] = j_axial_mev
            # Lateral neighbour (B-lattice seam shift)
            p_next = (p + 1) % n_protofilaments
            layer_next = layer + 1 if p_next == 0 else layer
            if 0 <= layer_next < n_layers:
                j = idx(p_next, layer_next)
                hamiltonian[i, j] = hamiltonian[j, i] = j_lateral_mev
    return hamiltonian


# ---------------------------------------------------------------------------
# Radiative decay matrix (free-space Lehmberg/Agarwal kernel)
# ---------------------------------------------------------------------------

def lehmberg_rate_ratio(
    kr: float,
    mu_i: np.ndarray,
    mu_j: np.ndarray,
    r_hat: np.ndarray,
) -> float:
    """Return Gamma_ij / Gamma0 for two dipoles using the real Lehmberg kernel.

    This is the free-space result. Dielectric geometry corrections are NOT
    included (see module docstring).
    """
    mu_dot = float(np.dot(mu_i, mu_j))
    mi_r = float(np.dot(mu_i, r_hat))
    mj_r = float(np.dot(mu_j, r_hat))

    if abs(kr) < 1e-6:
        # Near-field / same-site limit -> mu_dot
        return mu_dot

    term_far = (mu_dot - mi_r * mj_r) * math.sin(kr) / kr
    term_near = (mu_dot - 3.0 * mi_r * mj_r) * (
        math.cos(kr) / (kr**2) - math.sin(kr) / (kr**3)
    )
    return 1.5 * (term_far + term_near)


def build_radiative_matrix(
    positions_m: np.ndarray,
    dipoles: np.ndarray,
    wavelength_nm: float,
    refractive_index: float,
) -> np.ndarray:
    """Build the collective radiative decay matrix Gamma_ij / Gamma0.

    Diagonal entries are 1 (self-decay rate in units of Gamma0).
    Off-diagonal entries encode collective emission coupling.

    Args:
        positions_m: (N, 3) array of site positions in metres.
        dipoles:     (N, 3) array of unit dipole vectors.
        wavelength_nm: Transition wavelength in nm (default 280 for Trp).
        refractive_index: Effective refractive index of the medium.

    Returns:
        (N, N) symmetric float array.
    """
    n_sites = positions_m.shape[0]
    gamma = np.eye(n_sites, dtype=float)
    k_medium = 2.0 * math.pi * refractive_index / (wavelength_nm * 1e-9)

    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            delta = positions_m[i] - positions_m[j]
            distance = float(np.linalg.norm(delta))
            if distance <= 0.0:
                value = float(np.dot(dipoles[i], dipoles[j]))
            else:
                r_hat = delta / distance
                value = lehmberg_rate_ratio(
                    k_medium * distance, dipoles[i], dipoles[j], r_hat
                )
            gamma[i, j] = gamma[j, i] = value

    return gamma


# ---------------------------------------------------------------------------
# Mode analysis helpers
# ---------------------------------------------------------------------------

def inverse_participation_ratio(vectors: np.ndarray) -> np.ndarray:
    """Compute IPR for column eigenvectors in the site basis.

    IPR_k = 1 / sum_i |psi_k(i)|^4.

    High IPR -> delocalized; low IPR -> localized.
    """
    return 1.0 / np.maximum(np.sum(np.abs(vectors) ** 4, axis=0), 1e-15)


def mode_radiative_rates(
    h_evecs: np.ndarray, gamma_matrix: np.ndarray
) -> np.ndarray:
    """Compute the radiative rate <psi_k | Gamma | psi_k> for each eigenmode.

    Args:
        h_evecs:      (N, N) column eigenvectors of the Hamiltonian.
        gamma_matrix: (N, N) radiative decay matrix (Gamma_ij / Gamma0).

    Returns:
        (N,) array of Gamma_k / Gamma0 values.
    """
    # <psi_k|Gamma|psi_k> = evec_k^T . gamma . evec_k
    # Vectorised: diag( evecs^T . gamma . evecs )
    projected = h_evecs.T @ gamma_matrix @ h_evecs
    return np.real(np.diag(projected))


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_modes_csv(
    path: Path,
    h_evals: np.ndarray,
    h_ipr: np.ndarray,
    h_gamma: np.ndarray,
    gamma_evals: np.ndarray,
    gamma_ipr: np.ndarray,
    gamma_energy: np.ndarray,
) -> None:
    """Write excitonic and radiative mode diagnostics to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "basis", "mode_index", "energy_meV",
            "ipr", "gamma_over_gamma0",
        ])
        for idx, (energy, ipr, gval) in enumerate(
            zip(h_evals, h_ipr, h_gamma)
        ):
            writer.writerow([
                "hamiltonian_eigenmode",
                idx,
                f"{float(energy):.12g}",
                f"{float(ipr):.12g}",
                f"{float(gval):.12g}",
            ])
        for idx, (gval, ipr, energy) in enumerate(
            zip(gamma_evals, gamma_ipr, gamma_energy)
        ):
            writer.writerow([
                "radiative_eigenmode",
                idx,
                f"{float(energy):.12g}",
                f"{float(ipr):.12g}",
                f"{float(gval):.12g}",
            ])


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_spectrum(
    fig_prefix: Path,
    h_ipr: np.ndarray,
    h_gamma: np.ndarray,
    gamma_evals: np.ndarray,
    h_evals: np.ndarray,
) -> None:
    """Generate publication-ready histogram and IPR-vs-rate figures.

    Three panels:
        [0] Histogram of log10(Gamma_i/Gamma0) for radiative eigenmodes.
        [1] Scatter: Hamiltonian-mode IPR vs radiative rate, coloured by energy.
        [2] Energy spectrum coloured by radiative rate on log scale.
    """
    import matplotlib.pyplot as plt  # lazy import

    fig_prefix.parent.mkdir(parents=True, exist_ok=True)

    positive_h_gamma = np.maximum(h_gamma, 1e-15)
    positive_gamma_evals = np.maximum(gamma_evals, 1e-15)

    plt.rcParams.update({
        "font.size": 9,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "figure.dpi": 300,
        "axes.linewidth": 0.8,
    })
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.9))

    # Panel 0: histogram of radiative eigenvalue spectrum
    axes[0].hist(
        np.log10(positive_gamma_evals),
        bins=40,
        color="#3b7ea1",
        alpha=0.85,
        edgecolor="white",
        linewidth=0.3,
    )
    axes[0].axvline(0.0, color="black", lw=1.0, ls="--", label=r"$\Gamma_0$")
    axes[0].set_xlabel(r"$\log_{10}(\Gamma_i/\Gamma_0)$ — radiative eigenmodes")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Radiative decay spectrum")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(True, alpha=0.3, linestyle=":")

    # Panel 1: IPR vs <Gamma/Gamma0> scatter, coloured by Hamiltonian energy
    sc = axes[1].scatter(
        h_ipr,
        positive_h_gamma,
        c=h_evals,
        s=18,
        cmap="viridis",
        alpha=0.85,
    )
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Hamiltonian-mode IPR")
    axes[1].set_ylabel(r"$\langle\Gamma\rangle_k / \Gamma_0$")
    axes[1].set_title("IPR vs radiative rate\n(colour = energy / meV)")
    axes[1].axhline(1.0, color="black", lw=0.8, ls="--",
                    label=r"$\Gamma_0$ reference")
    axes[1].legend(frameon=False, fontsize=7)
    axes[1].grid(True, alpha=0.3, linestyle=":")
    cb1 = fig.colorbar(sc, ax=axes[1], fraction=0.046, pad=0.04)
    cb1.set_label("Energy (meV)", fontsize=8)

    # Subradiance annotation: modes with Gamma < 0.1 Gamma0
    n_subradiant = int(np.sum(positive_h_gamma < 0.1))
    n_superradiant = int(np.sum(positive_h_gamma > 10.0))
    axes[1].text(
        0.03, 0.03,
        f"Sub: {n_subradiant}  Super: {n_superradiant}",
        transform=axes[1].transAxes,
        fontsize=7,
        color="#333333",
        verticalalignment="bottom",
    )

    # Panel 2: energy spectrum, coloured by log10(Gamma/Gamma0)
    mode_idx = np.arange(h_evals.size)
    sc2 = axes[2].scatter(
        mode_idx,
        h_evals,
        c=np.log10(positive_h_gamma),
        s=14,
        cmap="RdYlBu_r",
        alpha=0.85,
    )
    axes[2].set_xlabel("Mode index")
    axes[2].set_ylabel("Energy (meV)")
    axes[2].set_title("Excitonic spectrum\n" r"(colour = $\log_{10}\Gamma/\Gamma_0$)")
    axes[2].grid(True, alpha=0.3, linestyle=":")
    cb2 = fig.colorbar(sc2, ax=axes[2], fraction=0.046, pad=0.04)
    cb2.set_label(r"$\log_{10}(\Gamma/\Gamma_0)$", fontsize=8)

    fig.tight_layout()
    fig.savefig(str(fig_prefix) + ".png", bbox_inches="tight")
    fig.savefig(str(fig_prefix) + ".pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description=(
            "Compute collective radiative decay rates for MT tryptophan lattice "
            "modes and generate IPR-vs-rate figures for the paper SI."
        )
    )
    parser.add_argument(
        "--project-root", default=None,
        help="Optional KwanTube project root override.",
    )
    # Geometry parameters
    parser.add_argument(
        "--n-layers", type=int, default=20,
        help="Number of axial layers (dimer rows per protofilament).",
    )
    parser.add_argument(
        "--n-protofilaments", type=int, default=13,
        help="Number of protofilaments (default: 13, B-lattice MT).",
    )
    parser.add_argument(
        "--axial-spacing-nm", type=float, default=4.0,
        help="Axial dimer spacing in nm (default: 4.0 nm).",
    )
    parser.add_argument(
        "--seam-rise-nm", type=float, default=0.92,
        help="B-lattice seam rise per protofilament in nm.",
    )
    parser.add_argument(
        "--radius-nm", type=float, default=None,
        help=(
            "Cylinder radius in nm. If omitted, inferred from "
            "--r-lateral-nm and --seam-rise-nm."
        ),
    )
    parser.add_argument(
        "--r-lateral-nm", type=float, default=5.2,
        help="Lateral nearest-neighbour distance in nm (used to infer radius).",
    )
    parser.add_argument(
        "--dipole-orientation", default="tangential",
        choices=["radial", "tangential", "axial", "helical45", "random"],
        help="Dipole orientation for all sites (default: tangential).",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for 'random' dipole orientation.",
    )
    # Hamiltonian parameters (meV)
    parser.add_argument(
        "--j-axial-mev", type=float, default=-88.08,
        help="Axial nearest-neighbour coupling in meV (default: -88.08).",
    )
    parser.add_argument(
        "--j-lateral-mev", type=float, default=160.36,
        help="Lateral nearest-neighbour coupling in meV (default: 160.36).",
    )
    # Radiative matrix parameters
    parser.add_argument(
        "--wavelength-nm", type=float, default=280.0,
        help="Transition wavelength in nm (default: 280 for Trp UV band).",
    )
    parser.add_argument(
        "--refractive-index", type=float, default=1.33,
        help="Effective refractive index of medium (default: 1.33, water).",
    )
    # Output paths
    parser.add_argument(
        "--output-json",
        default=(
            "outputs_data/raw_json/metrics/"
            "subradiant_decay_spectrum.json"
        ),
        help="Output JSON path relative to project root.",
    )
    parser.add_argument(
        "--output-csv",
        default="outputs_data/raw_csv/subradiant_modes.csv",
        help="Output CSV path relative to project root.",
    )
    parser.add_argument(
        "--fig-prefix",
        default="outputs_data/figures_final/subradiant_spectrum",
        help="Figure prefix relative to project root (no extension).",
    )
    args = parser.parse_args()

    project_root = resolve_project_root(args.project_root)
    audit_log(
        project_root,
        "[RUN_AUDIT] START script=compute_subradiant_decay_spectrum.py",
    )
    start = time.time()

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------
    if args.radius_nm is not None:
        radius_nm = args.radius_nm
    else:
        radius_nm = infer_radius_nm(
            args.r_lateral_nm, args.seam_rise_nm, args.n_protofilaments
        )

    n_sites = args.n_layers * args.n_protofilaments
    print(
        f"[INFO] Building B-lattice: {args.n_protofilaments} pf x "
        f"{args.n_layers} layers = {n_sites} sites, "
        f"radius={radius_nm:.2f} nm, orientation={args.dipole_orientation}"
    )

    positions_m, dipoles, labels = build_b_lattice_geometry(
        n_layers=args.n_layers,
        n_protofilaments=args.n_protofilaments,
        axial_spacing_nm=args.axial_spacing_nm,
        seam_rise_nm=args.seam_rise_nm,
        radius_nm=radius_nm,
        dipole_orientation=args.dipole_orientation,
        seed=args.seed,
    )

    # ------------------------------------------------------------------
    # Excitonic Hamiltonian
    # ------------------------------------------------------------------
    print(
        f"[INFO] Building tight-binding Hamiltonian: "
        f"J_axial={args.j_axial_mev:.2f} meV, "
        f"J_lateral={args.j_lateral_mev:.2f} meV"
    )
    H = build_tight_binding_hamiltonian(
        n_layers=args.n_layers,
        n_protofilaments=args.n_protofilaments,
        j_axial_mev=args.j_axial_mev,
        j_lateral_mev=args.j_lateral_mev,
    )
    h_evals, h_evecs = np.linalg.eigh(H)
    h_ipr = inverse_participation_ratio(h_evecs)

    spectral_gap = float(h_evals[-1] - h_evals[0])
    ipr_max = float(np.max(h_ipr))
    ipr_min = float(np.min(h_ipr))
    print(
        f"[INFO] Hamiltonian: spectral gap={spectral_gap:.2f} meV, "
        f"IPR max={ipr_max:.1f}, IPR min={ipr_min:.1f}"
    )

    # ------------------------------------------------------------------
    # Radiative decay matrix
    # ------------------------------------------------------------------
    print(
        f"[INFO] Building radiative matrix ({n_sites}x{n_sites}): "
        f"lambda={args.wavelength_nm:.0f} nm, n={args.refractive_index:.2f}"
    )
    gamma_matrix = build_radiative_matrix(
        positions_m=positions_m,
        dipoles=dipoles,
        wavelength_nm=args.wavelength_nm,
        refractive_index=args.refractive_index,
    )

    # Radiative rates for Hamiltonian eigenmodes
    h_gamma = mode_radiative_rates(h_evecs, gamma_matrix)

    # Eigendecomposition of radiative matrix itself
    gamma_evals, gamma_evecs = np.linalg.eigh(gamma_matrix)
    gamma_ipr = inverse_participation_ratio(gamma_evecs)

    # Associate radiative eigenmodes with Hamiltonian energy expectation
    # <H>_k = gamma_evec_k^T . H . gamma_evec_k
    gamma_energy = np.real(
        np.diag(gamma_evecs.T @ H @ gamma_evecs)
    )

    # Summary statistics
    n_subradiant = int(np.sum(h_gamma < 0.1))
    n_superradiant = int(np.sum(h_gamma > 10.0))
    gamma_min = float(np.min(gamma_evals))
    gamma_max = float(np.max(gamma_evals))

    print(
        f"[INFO] Radiative rates: "
        f"sub (Gamma<0.1): {n_subradiant}/{n_sites}, "
        f"super (Gamma>10): {n_superradiant}/{n_sites}"
    )
    print(
        f"[INFO] Gamma eigenvalue range: "
        f"[{gamma_min:.4g}, {gamma_max:.4g}] x Gamma0"
    )

    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------
    output_csv = project_root / args.output_csv
    write_modes_csv(
        path=output_csv,
        h_evals=h_evals,
        h_ipr=h_ipr,
        h_gamma=h_gamma,
        gamma_evals=gamma_evals,
        gamma_ipr=gamma_ipr,
        gamma_energy=gamma_energy,
    )

    fig_prefix = project_root / args.fig_prefix
    try:
        plot_spectrum(
            fig_prefix=fig_prefix,
            h_ipr=h_ipr,
            h_gamma=h_gamma,
            gamma_evals=gamma_evals,
            h_evals=h_evals,
        )
        print(f"[INFO] Figures saved: {fig_prefix}.png / .pdf")
    except ImportError:
        print("[WARN] matplotlib not available; figures skipped.")

    # Build JSON report
    free_space_caveat = (
        "Collective decay rates are computed using the free-space dyadic "
        "Green's function. In the actual MT geometry (layered dielectric "
        "cylinder: inner lumen eps~80, protein wall eps~2-4, outer water "
        "eps~80), local field effects and leaky-wave contributions would "
        "redistribute these rates. Modes identified as subradiant here "
        "represent an UPPER BOUND on protection; true subradiance requires "
        "the formalism of dipole emission in microcylinders."
    )

    report = {
        "script": "compute_subradiant_decay_spectrum.py",
        "version": "1.0.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(project_root),
        "free_space_caveat": free_space_caveat,
        "geometry": {
            "n_layers": args.n_layers,
            "n_protofilaments": args.n_protofilaments,
            "n_sites": n_sites,
            "axial_spacing_nm": args.axial_spacing_nm,
            "seam_rise_nm": args.seam_rise_nm,
            "radius_nm": radius_nm,
            "dipole_orientation": args.dipole_orientation,
            "seed": args.seed,
        },
        "hamiltonian": {
            "j_axial_mev": args.j_axial_mev,
            "j_lateral_mev": args.j_lateral_mev,
            "spectral_gap_mev": spectral_gap,
            "energy_min_mev": float(h_evals[0]),
            "energy_max_mev": float(h_evals[-1]),
            "ipr_max": ipr_max,
            "ipr_min": ipr_min,
            "ipr_mean": float(np.mean(h_ipr)),
        },
        "radiative_matrix": {
            "wavelength_nm": args.wavelength_nm,
            "refractive_index": args.refractive_index,
            "gamma_eigenvalue_min": gamma_min,
            "gamma_eigenvalue_max": gamma_max,
        },
        "mode_classification": {
            "n_subradiant_gamma_lt_0p1": n_subradiant,
            "n_superradiant_gamma_gt_10": n_superradiant,
            "n_intermediate": n_sites - n_subradiant - n_superradiant,
            "fraction_subradiant": n_subradiant / max(n_sites, 1),
            "fraction_superradiant": n_superradiant / max(n_sites, 1),
            "mean_h_gamma": float(np.mean(h_gamma)),
            "median_h_gamma": float(np.median(h_gamma)),
        },
        "outputs": {
            "csv": str(output_csv),
            "figure_png": str(fig_prefix) + ".png",
            "figure_pdf": str(fig_prefix) + ".pdf",
        },
        "elapsed_seconds": time.time() - start,
    }

    output_json = project_root / args.output_json
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------
    elapsed = time.time() - start
    print("=" * 78)
    print("SUBRADIANT DECAY SPECTRUM COMPLETE")
    print("=" * 78)
    print(f"  Sites           : {n_sites}")
    print(f"  Spectral gap    : {spectral_gap:.2f} meV")
    print(f"  IPR (max/mean)  : {ipr_max:.1f} / {np.mean(h_ipr):.1f}")
    print(
        f"  Subradiant modes (Gamma<0.1*Gamma0) : "
        f"{n_subradiant}/{n_sites} "
        f"({100.0*n_subradiant/max(n_sites,1):.1f}%)"
    )
    print(
        f"  Superradiant modes (Gamma>10*Gamma0): "
        f"{n_superradiant}/{n_sites} "
        f"({100.0*n_superradiant/max(n_sites,1):.1f}%)"
    )
    print(f"  CAVEAT: {free_space_caveat[:90]}...")
    print(f"  JSON  : {output_json}")
    print(f"[compute_subradiant_decay_spectrum] spectral_gap_meV={spectral_gap:.4f}")
    print(f"[compute_subradiant_decay_spectrum] subradiant_fraction={n_subradiant/max(n_sites,1):.4f}")
    print(f"[compute_subradiant_decay_spectrum] JSON -> {output_json}")
    print(f"  CSV   : {output_csv}")
    print(f"  FIG   : {fig_prefix}.png")
    print(f"  Time  : {elapsed:.1f} s")

    audit_log(
        project_root,
        f"[RUN_AUDIT] END status=ok n_sites={n_sites} "
        f"n_sub={n_subradiant} n_super={n_superradiant} "
        f"elapsed={elapsed:.1f}s output={output_json.name}",
    )
    return 0


if __name__ == "__main__":
    from pathlib import Path as _P
    import sys as _sys
    for _parent in _P(__file__).resolve().parents:
        if (_parent / "qmc_mt" / "run_audit.py").exists():
            _sys.path.insert(0, str(_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    _sys.exit(main())
