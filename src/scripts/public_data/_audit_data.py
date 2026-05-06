"""Audit script: what quantitative data do we actually have.

Resolves PROJECT_ROOT dynamically and auto-detects the latest
timestamped snapshot folder inside data_downloaded_public_repos/raw/public_api/.
"""
import json
import re
import statistics
import sys
from pathlib import Path

# Force UTF-8 output — prevents UnicodeEncodeError on Windows cp1252 consoles
# when printing titles with non-ASCII characters (e.g., '\u2010', '–', '·').
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ---------------------------------------------------------------------------
# Path resolution — robust regardless of CWD
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve()

# Walk up from src/scripts/data/ to find the workspace root (contains KwanTube/)
def _find_project_root(start: Path) -> Path:
    for parent in start.parents:
        if (parent / "data_downloaded_public_repos").exists():
            return parent
        if (parent / "KwanTube").exists() and (parent / "KwanTube" / "pyproject.toml").exists():
            return parent
    raise FileNotFoundError(
        f"Cannot locate project root (data_downloaded_public_repos) from {start}"
    )

WORKSPACE_ROOT = _find_project_root(_HERE)
API_SNAPSHOTS = WORKSPACE_ROOT / "data_downloaded_public_repos" / "raw" / "public_api"

# Auto-detect latest timestamped folder (lexicographic sort = chronological for ISO8601)
def _latest_snapshot(base: Path) -> Path:
    if not base.exists():
        raise FileNotFoundError(f"Snapshots directory not found: {base}")
    candidates = sorted(
        [d for d in base.iterdir() if d.is_dir()],
        key=lambda d: d.name,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(f"No snapshot folders found in {base}")
    return candidates[0]

p = _latest_snapshot(API_SNAPSHOTS)
print(f"[_audit_data] Using snapshot: {p.name}  ({p})")

# ---------------------------------------------------------------------------
# RCSB: Wilson B, resolution, method
# ---------------------------------------------------------------------------
core = json.loads((p / "rcsb_tubulin_core_entries.json").read_text(encoding="utf-8"))
records = [r for r in core.get("records", []) if r.get("status") == "ok"]
print(f"RCSB records OK: {len(records)}")

wilson_b: list[float] = []
resolutions: list[float] = []
methods: dict[str, int] = {}
for r in records:
    payload = r.get("payload", {})
    vrpt = payload.get("pdbx_vrpt_summary_diffraction", [])
    if isinstance(vrpt, list) and vrpt:
        wb = vrpt[0].get("Wilson_B_estimate")
        if wb is not None:
            try:
                wilson_b.append(float(wb))
            except Exception:
                pass
    rcsb_info = payload.get("rcsb_entry_info", {})
    res = rcsb_info.get("resolution_combined", [])
    if res:
        try:
            resolutions.append(float(res[0]))
        except Exception:
            pass
    exptl = payload.get("exptl", [])
    if exptl:
        m = exptl[0].get("method", "UNKNOWN")
        methods[m] = methods.get(m, 0) + 1

if len(wilson_b) >= 2:
    print(
        f"  Wilson B (n={len(wilson_b)}): mean={statistics.mean(wilson_b):.1f}, "
        f"median={statistics.median(wilson_b):.1f}, stdev={statistics.stdev(wilson_b):.1f}"
    )
else:
    print(f"  Wilson B: insufficient data (n={len(wilson_b)})")

if len(resolutions) >= 2:
    print(
        f"  Resolution (n={len(resolutions)}): mean={statistics.mean(resolutions):.2f} Ang, "
        f"median={statistics.median(resolutions):.2f}"
    )
else:
    print(f"  Resolution: insufficient data (n={len(resolutions)})")

print(f"  Methods: {methods}")

# ---------------------------------------------------------------------------
# OpenAlex
# ---------------------------------------------------------------------------
oa = json.loads(
    (p / "openalex_microtubule_quantum.json").read_text(encoding="utf-8", errors="replace")
)
pages = oa.get("pages", [oa])
works_oa: list = []
for pg in pages:
    works_oa.extend(pg.get("results", []))
print(f"\nOpenAlex works: {len(works_oa)}")
for w in works_oa[:8]:
    print(f"  - {str(w.get('title', ''))[:90]}")

# ---------------------------------------------------------------------------
# EuropePMC
# ---------------------------------------------------------------------------
epm = json.loads(
    (p / "europepmc_microtubule_spectroscopy.json").read_text(encoding="utf-8", errors="replace")
)
pages2 = epm.get("pages", [epm])
results_epm: list = []
for pg in pages2:
    results_epm.extend(pg.get("resultList", {}).get("result", []))
print(f"\nEuropePMC results: {len(results_epm)}")
for r in results_epm[:8]:
    print(f"  - {str(r.get('title', ''))[:90]}")

# ---------------------------------------------------------------------------
# Crossref
# ---------------------------------------------------------------------------
cr = json.loads(
    (p / "crossref_microtubule_spectroscopy.json").read_text(encoding="utf-8", errors="replace")
)
pages3 = cr.get("pages", [cr])
items_cr: list = []
for pg in pages3:
    items_cr.extend(pg.get("message", {}).get("items", []))
print(f"\nCrossref items: {len(items_cr)}")
for it in items_cr[:8]:
    title = it.get("title", [""])[0] if it.get("title") else ""
    print(f"  - {title[:90]}")

# ---------------------------------------------------------------------------
# PubChem: TOC structure of first compound file
# ---------------------------------------------------------------------------
pubchem_files = sorted(p.glob("pubchem_compound_*.json"))
print(f"\nPubChem files: {len(pubchem_files)}")
if pubchem_files:
    sample = json.loads(pubchem_files[0].read_text(encoding="utf-8", errors="ignore"))
    toc_headings: list[str] = []

    def get_toc(obj: object, depth: int = 0) -> None:
        if depth > 4:
            return
        if isinstance(obj, dict):
            h = obj.get("TOCHeading")
            if h and depth <= 2:
                toc_headings.append(h)
            for v in obj.values():
                get_toc(v, depth + 1)
        elif isinstance(obj, list):
            for item in obj:
                get_toc(item, depth + 1)

    get_toc(sample)
    unique_headings = list(dict.fromkeys(toc_headings))
    print(f"  TOC sections in sample PubChem file: {unique_headings[:30]}")
else:
    print("  No PubChem files found in snapshot.")

print("\n[_audit_data] DONE")
