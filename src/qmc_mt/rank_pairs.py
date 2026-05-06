"""
rank_pairs.py - lee pdb_tubulin_analysis.json y tabula:
  (1) top-10 pares Trp-Trp por G = kappa^2/r^6 para cada estructura
  (2) clasificacion intra-cadena vs inter-cadena
  (3) comparacion 1JFF vs 1TUB (pares compartidos)
Corregido: Encoding ASCII-safe para terminales Windows.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

JSON_PATH = PROJECT_ROOT / "outputs_data" / "raw_json" / "metrics" / "pdb_tubulin_analysis.json"

def load():
    return json.loads(JSON_PATH.read_text())

def classify(pair_str: str) -> str:
    left, right = pair_str.split("-")
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
    # ASCII-safe: kappa^2 instead of greek
    print(f"\n=== {pid}  (nTrp={r['n_trp']}, <kappa2>={r['summary']['kappa2']['mean']:.3f}) ===")
    print(f"{'rank':>4}  {'pair':<20}  {'r(A)':>7}  {'kappa2':>6}  {'G(A^-6) *1e9':>14}  {'tipo':>5}")
    out = []
    for k, p in enumerate(pairs_sorted, 1):
        pair = f"{p['donor']}-{p['acceptor']}"
        tipo = classify(pair)
        print(f"{k:>4}  {pair:<20}  {p['r_A']:>7.2f}  {p['kappa2']:>6.3f}  "
              f"{p['G_A-6']*1e9:>14.3f}  {tipo:>5}")
        out.append((pair, tipo, p["r_A"], p["kappa2"], p["G_A-6"]))
    return out

def intra_inter_balance(pid: str, data: dict):
    pairs = data[pid]["pairs_within_cutoff"]
    G_intra = sum(p["G_A-6"] for p in pairs if classify(f"{p['donor']}-{p['acceptor']}") == "intra")
    G_inter = sum(p["G_A-6"] for p in pairs if classify(f"{p['donor']}-{p['acceptor']}") == "inter")
    total = G_intra + G_inter
    if total == 0:
        print(f"[{pid}] sin pares")
        return
    print(f"[{pid}] peso FRET intra-cadena = {100*G_intra/total:5.1f}%   "
          f"inter-cadena = {100*G_inter/total:5.1f}%")

def pair_key(p):
    a = p["donor"].split(":")[1]
    b = p["acceptor"].split(":")[1]
    return tuple(sorted((a, b)))

def compare_structures(pid1: str, pid2: str, data: dict, n: int = 10):
    s1 = {pair_key(p): p for p in data[pid1]["pairs_within_cutoff"]}
    s2 = {pair_key(p): p for p in data[pid2]["pairs_within_cutoff"]}
    top1 = sorted(data[pid1]["pairs_within_cutoff"], key=lambda p: -p["G_A-6"])[:n]
    print(f"\n=== top-{n} de {pid1} - aparecen en {pid2} ===")
    print(f"{'pair':<20}  {'G1(*1e9)':>10}  {'G2(*1e9)':>10}  {'ratio':>7}")
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
        pair = f"{p['donor']}-{p['acceptor']}"
        print(f"{pair:<20}  {g1:>10.3f}  {g2:>10.3f}  {marca}")

if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    try:
        data = load()
        for pid in ("1JFF", "1TUB", "6DPU"):
            if pid in data and "pairs_within_cutoff" in data[pid]:
                top_table(pid, data)
        print("\n--- balance intra-cadena vs inter-cadena (suma de G) ---")
        for pid in ("1JFF", "1TUB", "6DPU"):
            if pid in data and "pairs_within_cutoff" in data[pid]:
                intra_inter_balance(pid, data)
        if "1JFF" in data and "1TUB" in data:
            compare_structures("1JFF", "1TUB", data)
    except FileNotFoundError:
        print(f"Error: {JSON_PATH} no encontrado.")