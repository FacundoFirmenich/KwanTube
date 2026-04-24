"""
Generate all figures for the paper (Tier 1 Production Version).
Aesthetic Override: Blue/Purple/Silver Palette (No Red/Green/Yellow/Black).
Version: 3.5.0 (Tier-A Production)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import sys
import os
import json

# Robust path resolution relative to repository root
REPO_ROOT = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

# Forzar backend sin interfaz para servidores/terminal
import matplotlib
matplotlib.use("Agg")

from qmc_mt.core import (
    TubulinDimer,
    ExperimentalParameters,
    DecoherenceModel,
    const,
    CoherenceUtility
)
from qmc_mt.noneq import FrohlichCondensation, QEDCavityModel
from qmc_mt.model_selection import bic_analysis

# IEEE / Nature-style global parameters
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'figure.dpi': 200,
    'savefig.dpi': 600,
    'lines.linewidth': 2.5,
    'axes.linewidth': 1.5,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'legend.frameon': True,
    'legend.fontsize': 10,
    'axes.edgecolor': '#003366',
    'xtick.color': '#003366',
    'ytick.color': '#003366',
    'text.color': '#003366',
})

# PREMIUM COLOR PALETTE
COLORS = {
    'navy': '#003366', 'royal': '#0074D9', 'sky': '#00BFFF',
    'murex': '#660066', 'magenta': '#DF00FF', 'violet': '#8A2BE2',
    'silver': '#C0C0C0', 'platinum': '#E5E4E2', 'pearl': '#F8F9FA'
}

MAP = {
    'baseline': COLORS['navy'], 'equilibrium': COLORS['murex'],
    'frohlich': COLORS['violet'], 'cavity': COLORS['royal'],
    'biology': COLORS['magenta'], 'reference': COLORS['silver']
}

OUT = REPO_ROOT / "figures_final"
OUT.mkdir(exist_ok=True)

def figure1_landscape():
    print("Fig 1: Landscape...")
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xscale('log')
    # Load real Redfield results if available for Fig 1
    try:
        # Search for redfield results in multiple locations
        rf_path = REPO_ROOT / "redfield_summary.json"
        if not rf_path.exists():
            rf_path = REPO_ROOT / "outputs" / "redfield_summary.json"
            
        with open(rf_path, 'r') as f:
            rf_data = json.load(f)
        t2_eq = next(r['results']['tau_coh'] for r in rf_data if r['pdb'] == '1JFF') * 1e-15
        lab_eq = 'Tubulin equilibrium (This work)'
    except Exception:
        t2_eq = 3.91e-14  # Redfield secular fallback (approx 39.1 fs)
        lab_eq = 'Tubulin equilibrium (Redfield)'

    data = [
        (60e-15,  'FMO electronic (Duan 2017)', MAP['reference'], 'o'),
        (t2_eq,   lab_eq, MAP['equilibrium'], 'D'),
        (5e-12,   'FMO vibrational (Thyrhaug 2018)', MAP['reference'], 's'),
        (50e-12,  'Tubulin Fröhlich (N=1000)', MAP['frohlich'], '^'),
        (5e-7,    'MT QED cavity (Mavromatos 2025)', MAP['cavity'], 'p'),
        (5e-6,    'Cryptochrome spin lifetime', MAP['reference'], 'v'),
        (25e-3,   'Neural Gamma (40 Hz)', MAP['biology'], '*'),
    ]
    for i, (tau, lab, col, mk) in enumerate(data):
        y_pos = i + 1
        ax.scatter([tau], [y_pos], s=300, c=col, marker=mk, edgecolor=COLORS['pearl'], lw=2, zorder=5)
        ax.axvline(tau, c=col, ls=':', alpha=0.4, zorder=1)
        ax.text(tau * 1.3, y_pos, lab, va='center', fontsize=11, fontweight='bold' if 'Tubulin' in lab else 'normal', color=col)

    ax.axvspan(1e-15, 1e-12, color=COLORS['silver'], alpha=0.15)
    ax.axvspan(1e-12, 1e-6, color=COLORS['royal'], alpha=0.08)
    ax.axvspan(1e-6, 1e-1, color=COLORS['magenta'], alpha=0.08)
    ax.annotate('', xy=(4.08e-14, 7.5), xytext=(25e-3, 7.5), arrowprops=dict(arrowstyle='<->', color=MAP['equilibrium'], lw=3))
    ax.text(3e-8, 7.8, 'AMPLIFICATION GAP Δ ≈ 11.8', color=MAP['equilibrium'], fontweight='bold', ha='center', fontsize=13,
            bbox=dict(boxstyle='round,pad=0.4', fc=COLORS['pearl'], ec=MAP['equilibrium'], lw=2, alpha=0.9))
    ax.set_xlabel('Timescale (seconds)')
    ax.set_xlim(1e-16, 1e-1); ax.set_ylim(0, 8.5); ax.set_yticks([])
    ax2 = ax.twiny(); ax2.set_xscale('log'); ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks([1e-15, 1e-12, 1e-9, 1e-6, 1e-3])
    ax2.set_xticklabels(['fs', 'ps', 'ns', 'μs', 'ms'], fontweight='bold', color=COLORS['navy'])
    plt.savefig(OUT / 'fig1_landscape.pdf', bbox_inches='tight')
    plt.savefig(OUT / 'fig1_landscape.png', bbox_inches='tight')
    plt.close()

def figure2_signatures():
    print("Fig 2: Signatures...")
    fig = plt.figure(figsize=(15, 10), facecolor=COLORS['pearl'])
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.25)
    cav = QEDCavityModel(L_MT=25e-6)
    dE = cav.splitting_eV()
    
    ax_a = fig.add_subplot(gs[0, 0])
    omega_0 = 6.0 # THz
    det = np.linspace(-3, 3, 1000)
    l0_THz = cav.lambda_0() / 1e12
    N = cav.N_dimers
    br_p = omega_0 - det/2 + 0.5 * np.sqrt(det**2 + 4 * N * l0_THz**2)
    br_m = omega_0 - det/2 - 0.5 * np.sqrt(det**2 + 4 * N * l0_THz**2)
    ax_a.plot(det, br_p, color=MAP['cavity'], lw=4, label='Upper Polariton')
    ax_a.plot(det, br_m, color=MAP['equilibrium'], lw=4, label='Lower Polariton')
    ax_a.set_xlabel('Detuning Δ (THz)'); ax_a.set_ylabel('Freq (THz)'); ax_a.legend()

    ax_b = fig.add_subplot(gs[0, 1])
    lam = np.linspace(270, 290, 1000)
    E_eV = const.h * const.c / (lam * 1e-9) / const.eV_to_J
    E0 = 4.43 # 280 nm
    S = 0.08
    def L(E, e0): return (S/2)**2 / ((E-e0)**2 + (S/2)**2)
    A = L(E_eV, E0 - dE/2) + L(E_eV, E0 + dE/2)
    d2A = -np.gradient(np.gradient(A, lam), lam)
    ax_b.plot(lam, d2A/d2A.max(), color=MAP['cavity'], lw=4)
    ax_b.set_xlabel('Wavelength (nm)'); ax_b.set_ylabel('d²A/dλ²'); ax_b.set_title(f'Splitting ≈ {dE*1e3:.1f} meV')

    ax_c = fig.add_subplot(gs[1, :])
    snr, dbic = bic_analysis(n_realizations=10, rng_seed=42, effective_points=20)
    ax_c.semilogx(snr, dbic, color=MAP['frohlich'], lw=4)
    ax_c.axhline(10, ls='--', color=MAP['biology'])
    crossing = next((float(s) for s, v in zip(snr, dbic) if v > 10), None)
    if crossing is not None:
        ax_c.axvline(crossing, ls=':', color=MAP['equilibrium'])
        ax_c.text(crossing * 1.05, 12, f'SNR* ≈ {crossing:.0f}', color=MAP['equilibrium'])
    ax_c.set_xlabel('SNR'); ax_c.set_ylabel('ΔBIC (Decisive evidence > 10)')
    plt.savefig(OUT / 'fig2_signatures.pdf', bbox_inches='tight')
    plt.savefig(OUT / 'fig2_signatures.png', bbox_inches='tight')
    plt.close()

def figure3_frohlich():
    print("Fig 3: Frohlich...")
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6), facecolor=COLORS['pearl'])
    fr = FrohlichCondensation(TubulinDimer(), ExperimentalParameters())
    N = np.logspace(0, 6, 200)
    eta = fr.pumping_parameter(N)
    ax_a.loglog(N, eta, color=MAP['frohlich'], lw=4)
    ax_a.axhline(1, ls='--', color=MAP['equilibrium'])
    ax_a.set_xlabel('N dimers'); ax_a.set_ylabel('η pumping'); ax_a.grid(alpha=0.1)
    
    P_th = const.kB * 310 * 1e6
    P_met = N * 1.0 * 0.42 * const.eV_to_J
    ax_b.loglog(N, P_met, label='Metabolic Power')
    ax_b.axhline(P_th, ls='--', color=MAP['equilibrium'], label='Thermal Floor (1 MHz)')
    ax_b.set_xlabel('N'); ax_b.set_ylabel('Power (W)'); ax_b.legend()
    plt.savefig(OUT / 'fig3_frohlich.pdf', bbox_inches='tight')
    plt.savefig(OUT / 'fig3_frohlich.png', bbox_inches='tight')
    plt.close()

def figure4_scaling():
    print("Fig 4: Scaling...")
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6), facecolor=COLORS['pearl'])
    L = np.linspace(1, 50, 50)
    sp = [QEDCavityModel(L_MT=l*1e-6).splitting_eV()*1e3 for l in L]
    ax_a.plot(L, sp, color=MAP['cavity'], lw=4)
    ax_a.set_xlabel('L_MT (μm)'); ax_a.set_ylabel('Splitting (meV)'); ax_a.set_ylim(13.5, 14.5)

    eps = np.linspace(30, 120, 50)
    tau = [5e-7 * (e/80)**2 * 1e6 for e in eps]
    ax_b.plot(eps, tau, color=MAP['cavity'], lw=4)
    ax_b.set_xlabel('Dielectric ε'); ax_b.set_ylabel('τ_cavity (μs)'); ax_b.set_title('τ ∝ ε² dependence')
    plt.savefig(OUT / 'fig4_scaling.pdf', bbox_inches='tight')
    plt.savefig(OUT / 'fig4_scaling.png', bbox_inches='tight')
    plt.close()

def main():
    figure1_landscape(); figure2_signatures(); figure3_frohlich(); figure4_scaling()
    print("\nVisual excellence v3.5.0: OK (600 DPI Dual Export)")

if __name__ == "__main__": main()
