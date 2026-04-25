"""Audit script: what quantitative data do we actually have."""
import json, re, statistics
from pathlib import Path

p = Path("data/raw/public_api/20260425T115542Z")

# --- RCSB: Wilson B, resolution, method ---
core = json.loads((p / "rcsb_tubulin_core_entries.json").read_text())
records = [r for r in core.get("records", []) if r.get("status") == "ok"]
print(f"RCSB records OK: {len(records)}")

wilson_b = []
resolutions = []
methods = {}
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

print(f"  Wilson B (n={len(wilson_b)}): mean={statistics.mean(wilson_b):.1f}, "
      f"median={statistics.median(wilson_b):.1f}, stdev={statistics.stdev(wilson_b):.1f}")
print(f"  Resolution (n={len(resolutions)}): mean={statistics.mean(resolutions):.2f} Ang, "
      f"median={statistics.median(resolutions):.2f}")
print(f"  Methods: {methods}")

# --- OpenAlex ---
oa = json.loads((p / "openalex_microtubule_quantum.json").read_text(encoding="utf-8", errors="replace"))
pages = oa.get("pages", [oa])
works_oa = []
for pg in pages:
    works_oa.extend(pg.get("results", []))
print(f"\nOpenAlex works: {len(works_oa)}")
for w in works_oa[:8]:
    print(f"  - {str(w.get('title',''))[:90]}")

# --- EuropePMC ---
epm = json.loads((p / "europepmc_microtubule_spectroscopy.json").read_text(encoding="utf-8", errors="replace"))
pages2 = epm.get("pages", [epm])
results_epm = []
for pg in pages2:
    results_epm.extend(pg.get("resultList", {}).get("result", []))
print(f"\nEuropePMC results: {len(results_epm)}")
for r in results_epm[:8]:
    print(f"  - {str(r.get('title',''))[:90]}")

# --- Crossref ---
cr = json.loads((p / "crossref_microtubule_spectroscopy.json").read_text(encoding="utf-8", errors="replace"))
pages3 = cr.get("pages", [cr])
items_cr = []
for pg in pages3:
    items_cr.extend(pg.get("message", {}).get("items", []))
print(f"\nCrossref items: {len(items_cr)}")
for it in items_cr[:8]:
    title = it.get("title", [""])[0] if it.get("title") else ""
    print(f"  - {title[:90]}")

# --- PubChem: what's actually IN the pugview files ---
pubchem_files = sorted(p.glob("pubchem_compound_*.json"))
print(f"\nPubChem files: {len(pubchem_files)}")
sample = json.loads(pubchem_files[0].read_text(encoding="utf-8", errors="ignore"))
# Find top-level TOC sections
toc_headings = []
def get_toc(obj, depth=0):
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
