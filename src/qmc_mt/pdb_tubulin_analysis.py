"""
pdb_tubulin_analysis.py - real Trp geometry & Foerster couplings from RCSB.

Pipeline
--------
1. Download mmCIF files from files.rcsb.org (GET, cached on disk).
2. Pull entry metadata (title, resolution) from data.rcsb.org REST v1.
3. Parse with Biopython's MMCIFParser; for every TRP residue collect the
   nine indole heavy atoms (CG, CD1, NE1, CE2, CD2, CE3, CZ2, CZ3, CH2).
4. SVD-fit the indole plane; build in-plane orthonormal frame (x_hat,
   y_hat) with x_hat approx long axis (NE1 -> midpoint(CZ3, CH2)) and y_hat
   oriented such that CE3 lies on the +y side (fixes rotation sign).
5. Place the 1La and 1Lb transition dipoles per Callis (1997):
      mu_1La = cos(+38 deg)*x_hat + sin(+38 deg)*y_hat       (in-plane)
      mu_1Lb = mu_1La rotated +90 deg in-plane               (orthogonal proxy)
   The +38 deg convention and long-axis definition follow Callis, Methods
   Enzymol. 278 (1997) 113-150, Fig. 2.
6. For every ordered pair of Trps compute
      r_vec = R_A - R_D,    r = |r_vec|
      kappa^2   = (mu_D dot mu_A - 3 * (mu_D dot r_hat) * (mu_A dot r_hat))^2
      G    = kappa^2 / r^6        [A^-6, the structure-only FRET factor]
   A FRET rate requires also the donor lifetime and the spectral overlap
   integral; those are residue-environment-dependent and are supplied
   downstream. Here we expose only the geometric factor, which is the
   quantity that differs between tubulin isoforms at fixed photophysics.
7. Self-validate: Monte-Carlo isotropic average of kappa^2 must recover 2/3.

References
----------
Callis, P. R. (1997) "1La and 1Lb transitions of tryptophan: applications
of theory and experimental observations to fluorescence of proteins."
Methods Enzymol. 278, 113-150. doi:10.1016/S0076-6879(97)78009-1
RCSB PDB REST API: https://data.rcsb.org/redoc/index.html
"""
from __future__ import annotations
import json
import sys
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import requests
from Bio.PDB import MMCIFParser

# Boilerplate para resolver importaciones desde la raiz del paquete
PROJECT_ROOT = Path(__file__).resolve().parents[2] # retrocede desde src/qmc_mt/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

# ------------------------------------------------------------------ constants
INDOLE_ATOMS    = ("CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2")
ANGLE_1LA_DEG   = 38.0              # Callis 1997: 1La from long axis
ANGLE_1LB_DEG   = ANGLE_1LA_DEG + 90.0
RCSB_FILE_URL   = "https://files.rcsb.org/download/{pdb}.cif"
RCSB_META_URL   = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"
DEFAULT_CUTOFF_A = 40.0             # report pair detail within this range

# ------------------------------------------------------------------ HTTP
def fetch_cif(pdb_id: str, cache: Path, timeout: float = 60.0) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    out = cache / f"{pdb_id.upper()}.cif"
    if out.exists() and out.stat().st_size > 0:
        return out
    r = requests.get(RCSB_FILE_URL.format(pdb=pdb_id.upper()), timeout=timeout)
    r.raise_for_status()
    out.write_bytes(r.content)
    return out

def fetch_metadata(pdb_id: str, timeout: float = 30.0) -> dict:
    r = requests.get(RCSB_META_URL.format(pdb=pdb_id.upper()), timeout=timeout)
    r.raise_for_status()
    return r.json()

# ------------------------------------------------------------------ geometry
@dataclass
class TrpDipole:
    pdb_id: str
    chain: str
    resseq: int
    center: np.ndarray
    long_axis: np.ndarray
    normal: np.ndarray
    mu_1La: np.ndarray
    mu_1Lb: np.ndarray
    planarity_rmsd_A: float

    def to_dict(self):
        d = asdict(self)
        for k in ("center", "long_axis", "normal", "mu_1La", "mu_1Lb"):
            d[k] = [float(x) for x in getattr(self, k)]
        d["planarity_rmsd_A"] = float(self.planarity_rmsd_A)
        return d

def _best_fit_plane(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    c = xyz.mean(axis=0)
    _, s, vt = np.linalg.svd(xyz - c, full_matrices=False)
    n = vt[-1] / np.linalg.norm(vt[-1])
    rmsd = float(np.sqrt((s[-1] ** 2) / xyz.shape[0]))
    return c, n, rmsd

def _trp_dipole(coords: dict, pdb_id: str, chain: str, resseq: int):
    if not all(a in coords for a in INDOLE_ATOMS):
        return None
    xyz = np.array([coords[a] for a in INDOLE_ATOMS])
    center, normal, planarity = _best_fit_plane(xyz)

    # In-plane long axis: NE1 -> midpoint(CZ3, CH2), projected onto plane
    raw = 0.5 * (coords["CZ3"] + coords["CH2"]) - coords["NE1"]
    raw -= np.dot(raw, normal) * normal
    if np.linalg.norm(raw) < 1e-6:
        return None
    x_hat = raw / np.linalg.norm(raw)

    # In-plane y-axis: orient so CE3 lies on +y side (Callis sign convention)
    y_hat = np.cross(normal, x_hat)
    ce3 = coords["CE3"] - center
    ce3 -= np.dot(ce3, normal) * normal
    if np.dot(ce3, y_hat) < 0:
        y_hat = -y_hat
        normal = -normal

    a, b = np.deg2rad(ANGLE_1LA_DEG), np.deg2rad(ANGLE_1LB_DEG)
    mu_1La = np.cos(a) * x_hat + np.sin(a) * y_hat
    mu_1Lb = np.cos(b) * x_hat + np.sin(b) * y_hat

    return TrpDipole(pdb_id=pdb_id, chain=chain, resseq=resseq,
                     center=center, long_axis=x_hat, normal=normal,
                     mu_1La=mu_1La, mu_1Lb=mu_1Lb,
                     planarity_rmsd_A=planarity)

def parse_trps(cif_path: Path, pdb_id: str) -> list[TrpDipole]:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(cif_path))
    model = next(structure.get_models())
    dipoles: list[TrpDipole] = []
    for chain in model:
        for res in chain:
            if res.get_resname().strip() != "TRP":
                continue
            # Hetero-flag filter: keep standard residues only
            if res.get_id()[0].strip() not in ("", "H_TRP"):
                # keep blank hetflag (standard aa); skip waters/ligands
                if res.get_id()[0].strip() != "":
                    continue
            coords: dict[str, np.ndarray] = {}
            for atom in res:
                alt = atom.get_altloc()
                # accept blank/space (no altloc) or primary conformer 'A'
                if alt not in ("", " ", "A"):
                    continue
                name = atom.get_name()
                # first occurrence wins (avoids overwriting 'A' with 'B')
                if name not in coords:
                    coords[name] = np.array(atom.get_coord(), dtype=float)
            d = _trp_dipole(coords, pdb_id, chain.id, res.get_id()[1])
            if d is not None:
                dipoles.append(d)
    return dipoles

# ------------------------------------------------------------------ FRET
def kappa_squared(mu_d: np.ndarray, mu_a: np.ndarray, r_vec: np.ndarray) -> float:
    r_hat = r_vec / np.linalg.norm(r_vec)
    return float((np.dot(mu_d, mu_a) - 3 * np.dot(mu_d, r_hat) * np.dot(mu_a, r_hat)) ** 2)

def pairwise(dipoles: list[TrpDipole], cutoff_A: float) -> tuple[list[dict], dict]:
    n = len(dipoles)
    all_k2, all_r, all_G = [], [], []
    pairs_within = []
    for i in range(n):
        for j in range(i + 1, n):
            rv = dipoles[j].center - dipoles[i].center
            r = float(np.linalg.norm(rv))
            if r < 1e-6:
                continue
            k2 = kappa_squared(dipoles[i].mu_1La, dipoles[j].mu_1La, rv)
            G = k2 / r**6
            all_k2.append(k2); all_r.append(r); all_G.append(G)
            if r <= cutoff_A:
                pairs_within.append({
                    "donor": f"{dipoles[i].chain}:{dipoles[i].resseq}",
                    "acceptor": f"{dipoles[j].chain}:{dipoles[j].resseq}",
                    "r_A": r, "kappa2": k2, "G_A-6": G,
                })
    pairs_within.sort(key=lambda p: -p["G_A-6"])

    # Always-present skeleton so downstream code never KeyErrors
    summary = {
        "n_trp": n,
        "n_pairs": len(all_k2),
        "n_pairs_within_cutoff": len(pairs_within),
        "cutoff_A": cutoff_A,
        "r_A": {"min": None, "median": None, "max": None},
        "kappa2": {"mean": None, "std": None, "isotropic_reference": 2/3},
        "top5_by_G_A-6": [],
    }
    if all_k2:
        k2 = np.asarray(all_k2); rs = np.asarray(all_r)
        summary["r_A"] = {"min": float(rs.min()),
                          "median": float(np.median(rs)),
                          "max": float(rs.max())}
        summary["kappa2"] = {
            "mean": float(k2.mean()),
            "std":  float(k2.std(ddof=1)) if k2.size > 1 else 0.0,
            "isotropic_reference": 2/3,
        }
        summary["top5_by_G_A-6"] = [
            {"pair": f"{p['donor']}-{p['acceptor']}", "r_A": p["r_A"],
             "kappa2": p["kappa2"], "G_A-6": p["G_A-6"]}
            for p in pairs_within[:5]
        ]
    return pairs_within, summary

# ------------------------------------------------------------------ driver
def analyze(pdb_id: str, cache: Path, cutoff_A: float = DEFAULT_CUTOFF_A) -> dict:
    cif = fetch_cif(pdb_id, cache)
    try:
        meta = fetch_metadata(pdb_id)
        title = meta.get("struct", {}).get("title", "")
        res   = (meta.get("rcsb_entry_info", {})
                     .get("resolution_combined", [None]) or [None])[0]
        method = (meta.get("exptl", [{}])[0] or {}).get("method", "")
    except Exception as e:
        title, res, method = f"(metadata failed: {e})", None, ""

    dipoles = parse_trps(cif, pdb_id)
    pairs, summary = pairwise(dipoles, cutoff_A=cutoff_A)
    return {
        "pdb_id": pdb_id.upper(),
        "title": title,
        "method": method,
        "resolution_A": res,
        "n_trp": len(dipoles),
        "trp": [d.to_dict() for d in dipoles],
        "pairs_within_cutoff": pairs,
        "summary": summary,
    }

# ------------------------------------------------------------------ self-check
def _validate_kappa2_isotropic(n: int = 200_000, seed: int = 0) -> float:
    rng = np.random.default_rng(seed)
    def _u(k):
        v = rng.normal(size=(k, 3))
        return v / np.linalg.norm(v, axis=1, keepdims=True)
    d, a, r = _u(n), _u(n), _u(n)
    return float((((d * a).sum(1) - 3 * (d * r).sum(1) * (a * r).sum(1)) ** 2).mean())

# ------------------------------------------------------------------ main
if __name__ == "__main__":
    cache = Path(tempfile.gettempdir()) / "qmc_mt_pdb_cache"
    targets = ["1JFF", "1TUB", "6DPU"]
    if "--with-3j6f" in sys.argv:
        targets.append("3J6F")

    # 1. kappa^2 isotropic sanity check (must equal 2/3 within MC noise)
    k2_iso = _validate_kappa2_isotropic()
    print(f"self-check  <kappa^2>_iso = {k2_iso:.4f}   expected = 0.6667   "
          f"Delta = {abs(k2_iso - 2/3):.4f}")
    assert abs(k2_iso - 2/3) < 5e-3, "kappa^2 isotropic average failed"

    # 2. structures
    results = {}
    for pid in targets:
        try:
            r = analyze(pid, cache)
            results[pid] = r
            s = r["summary"]
            k2_mean = s["kappa2"]["mean"]
            rmin    = s["r_A"]["min"]
            print(f"[{pid}] {r['method']:20s}  res={r['resolution_A']}  "
                  f"nTrp={r['n_trp']:3d}  "
                  f"<kappa^2>={k2_mean if k2_mean is None else f'{k2_mean:.3f}'}  "
                  f"r_min={rmin if rmin is None else f'{rmin:.2f} A'}  "
                  f"pairs (r < {DEFAULT_CUTOFF_A:.0f} A) = {s['n_pairs_within_cutoff']}")
        except Exception as e:
            results[pid] = {"error": repr(e)}
            print(f"[{pid}] ERROR: {e}", file=sys.stderr)

    out_path = PROJECT_ROOT / "outputs_data" / "raw_json" / "pdb_tubulin_analysis.json"
    out_path.write_text(json.dumps(results, indent=2, default=float))
    print(f"\nwrote {out_path.resolve()}")