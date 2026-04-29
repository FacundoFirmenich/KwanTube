import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# IEEE / Nature-style
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 600,
    'lines.linewidth': 1.5,
    'axes.linewidth': 1.2,
    'axes.edgecolor': '#003366',
    'xtick.color': '#003366',
    'ytick.color': '#003366',
    'text.color': '#003366',
})

COLORS = ['#003366', '#0074D9', '#00BFFF', '#8A2BE2', '#DF00FF', '#660066', '#FF4136', '#FF851B']

def main():
    root = Path(r"c:\Users\User\3D Objects\biofisicaquantiqaCLINE\KwanTube")
    data_path = root / "outputs_data" / "raw_npz" / "heom_1JFF_full.npz"
    out_path = root / "outputs_data" / "figures_final" / "heom_benchmark_1JFF.pdf"
    
    if not data_path.exists():
        print(f"ERROR: Data not found at {data_path}")
        return

    print(f"Loading data from {data_path}...")
    d = np.load(data_path, allow_pickle=True)
    t = d['t_fs']
    P = d['P_site'] # (nt, N)
    labels = d['labels']
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    for i in range(P.shape[1]):
        ax.plot(t, P[:, i], label=labels[i], color=COLORS[i % len(COLORS)])
        
    ax.set_xlabel('Time (fs)', fontweight='bold')
    ax.set_ylabel('Site Population', fontweight='bold')
    ax.set_title('HEOM Benchmark (1JFF, NC=7, Nk=1)', fontsize=14, pad=15)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend(frameon=True, facecolor='white', framealpha=0.8, loc='upper right', fontsize=9)
    
    # Inset for early dynamics
    ax_ins = ax.inset_axes([0.15, 0.45, 0.35, 0.4])
    for i in range(P.shape[1]):
        ax_ins.plot(t[:200], P[:200, i], color=COLORS[i % len(COLORS)])
    ax_ins.set_xlim(0, 200)
    ax_ins.set_title('Short-time', fontsize=9)
    ax_ins.tick_params(labelsize=8)
    
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    plt.savefig(str(out_path).replace('.pdf', '.png'), bbox_inches='tight') # Also update PNG
    print(f"SUCCESS: Figure saved as {out_path} (PDF) and PNG.")

if __name__ == "__main__":
    main()
