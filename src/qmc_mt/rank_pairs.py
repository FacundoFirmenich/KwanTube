"""
rank_pairs.py — lee pdb_tubulin_analysis.json y tabula:
  (1) top-10 pares Trp-Trp por G = κ²/r⁶  para cada estructura
  (2) clasificación intra-cadena vs inter-cadena
  (3) comparación 1JFF vs 1TUB (pares compartidos)
"""
from __future__ import annotations
import json
from pathlib import Path

JSON_PATH = Path(r"C:\Users\User\3D Objects\biofisicaquantiqaCLINE\pdb_tubulin_analysis.json")

def load():
    return json.loads(JSON_PATH.read_text())

def classify(pair_str: str) -> str:
    # pair_str = "A:388—B:101"
    left, right = pair_str.split("—")
    cL = left.split(":")[0]
    cR = right.split(":")[0]
    return "intra" if cL == cR else "inter"

def top_table(pid: str, data: dict, n: int = 10):
    r = data[pid]
    if "summary" not in r:
        print(f"[{pid}] sin datos")
        return []
    pairs = r["pairs_within_cutoff"]
    pairs_sorted = sorted(pairs, key=lambda p: -p["G_A-6"])[:n]
    print(f"\n=== {pid}  (nTrp={r['n_trp']}, ⟨κ²⟩={r['summary']['kappa2']['mean']:.3f}) ===")
    print(f"{'rank':>4}  {'pair':<20}  {'r(Å)':>7}  {'κ²':>6}  {'G(Å⁻⁶) ×1e9':>14}  {'tipo':>5}")
    out = []
    for k, p in enumerate(pairs_sorted, 1):
        pair = f"{p['donor']}—{p['acceptor']}"
        tipo = classify(pair)
        print(f"{k:>4}  {pair:<20}  {p['r_A']:>7.2f}  {p['kappa2']:>6.3f}  "
              f"{p['G_A-6']*1e9:>14.3f}  {tipo:>5}")
        out.append((pair, tipo, p["r_A"], p["kappa2"], p["G_A-6"]))
    return out

def intra_inter_balance(pid: str, data: dict):
    pairs = data[pid]["pairs_within_cutoff"]
    G_intra = sum(p["G_A-6"] for p in pairs
                  if classify(f"{p['donor']}—{p['acceptor']}") == "intra")
    G_inter = sum(p["G_A-6"] for p in pairs
                  if classify(f"{p['donor']}—{p['acceptor']}") == "inter")
    total = G_intra + G_inter
    if total == 0:
        print(f"[{pid}] sin pares")
        return
    print(f"[{pid}] peso FRET intra-cadena = {100*G_intra/total:5.1f}%   "
          f"inter-cadena = {100*G_inter/total:5.1f}%")

def pair_key(p):
    # clave robusta: resseq del donor y acceptor ordenados (ignora cadena)
    a = p["donor"].split(":")[1]
    b = p["acceptor"].split(":")[1]
    return tuple(sorted((a, b)))

def compare_structures(pid1: str, pid2: str, data: dict, n: int = 10):
    s1 = {pair_key(p): p for p in data[pid1]["pairs_within_cutoff"]}
    s2 = {pair_key(p): p for p in data[pid2]["pairs_within_cutoff"]}
    top1 = sorted(data[pid1]["pairs_within_cutoff"],
                  key=lambda p: -p["G_A-6"])[:n]
    print(f"\n=== top-{n} de {pid1} — ¿aparecen en {pid2}? ===")
    print(f"{'pair':<20}  {'G1(×1e9)':>10}  {'G2(×1e9)':>10}  {'ratio':>7}")
    for p in top1:
        k = pair_key(p)
        g1 = p["G_A-6"] * 1e9
        if k in s2:
            g2 = s2[k]["G_A-6"] * 1e9
            ratio = g2 / g1 if g1 else float("nan")
            marca = f"{ratio:>7.2f}"
        else:
            g2 = 0.0
            marca = "   ---"
        pair = f"{p['donor']}—{p['acceptor']}"
        print(f"{pair:<20}  {g1:>10.3f}  {g2:>10.3f}  {marca}")

if __name__ == "__main__":
    data = load()

    # (1) top-10 por estructura
    for pid in ("1JFF", "1TUB", "6DPU"):
        if pid in data and "pairs_within_cutoff" in data[pid]:
            top_table(pid, data)

    # (2) balance intra/inter
    print("\n--- balance intra-cadena vs inter-cadena (suma de G) ---")
    for pid in ("1JFF", "1TUB", "6DPU"):
        if pid in data and "pairs_within_cutoff" in data[pid]:
            intra_inter_balance(pid, data)

    # (3) robustez: 1JFF vs 1TUB
    if "1JFF" in data and "1TUB" in data:
        compare_structures("1JFF", "1TUB", data)