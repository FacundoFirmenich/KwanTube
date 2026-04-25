"""Curate institutional public-data artifacts into analysis-ready CSV tables.

Axes extracted:
  A) Structural-spectral (RCSB): Wilson B, resolution, method, heterogeneity
     → informs bath coupling strength eta via B-factor proxy
  B) Literature (OpenAlex/Crossref/EuropePMC): study provenance
  C) Ligand (PubChem PUG-View): UV/Vis peaks, wavenumbers, spectral features
     → informs detectability / spectroscopic layer
"""

from __future__ import annotations

import csv
import json
import math
import re
import statistics
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _iter_latest_raw_json() -> Iterable[Path]:
    root = Path("data") / "raw" / "public_api"
    if not root.exists():
        return []
    runs = sorted([p for p in root.iterdir() if p.is_dir()])
    if not runs:
        return []
    return runs[-1]


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8", errors="replace"))


def _safe_float(val: Any, default: float = float("nan")) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def _safe_get(d: Dict[str, Any], path: List[str], default: Any = "") -> Any:
    cur: Any = d
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def _iter_strings(obj: Any) -> Iterable[str]:
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from _iter_strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _iter_strings(v)


def _iter_payload_pages(payload: Dict[str, Any]) -> Iterable[Dict[str, Any]]:
    """Yield one or multiple pages from stored payload."""
    pages = payload.get("pages")
    if isinstance(pages, list) and pages:
        for p in pages:
            if isinstance(p, dict):
                yield p
        return
    yield payload


def _write_progress(step: str, payload: Dict[str, str]) -> None:
    analysis_dir = Path("analysis")
    analysis_dir.mkdir(parents=True, exist_ok=True)
    progress_path = analysis_dir / "_progress_curate_compact.json"
    blob: Dict[str, str] = {
        "script": "curate_compact.py",
        "step": step,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    blob.update(payload)
    progress_path.write_text(str(blob), encoding="utf-8")


# ---------------------------------------------------------------------------
# Technique inference
# ---------------------------------------------------------------------------

def _infer_technique_tag(title: str) -> str:
    t = title.lower()
    mapping = [
        ("2d-ir", "2d-ir"),
        ("thz", "thz"),
        ("terahertz", "thz"),
        ("uv", "uv-vis"),
        ("fluorescence", "fluorescence"),
        ("nmr", "nmr"),
        ("raman", "raman"),
        ("infrared", "ir"),
        ("ftir", "ir"),
        ("spectroscopy", "spectroscopy"),
    ]
    for needle, tag in mapping:
        if needle in t:
            return tag
    return ""


def _is_relevant_study(title: str) -> bool:
    t = (title or "").lower()
    has_subject = ("microtubule" in t) or ("tubulin" in t)
    has_signal = any(
        k in t
        for k in [
            "spectroscopy", "2d-ir", "thz", "terahertz",
            "fluorescence", "uv", "nmr", "raman", "ftir", "infrared",
        ]
    )
    excludes = [
        "review for", "joint public review", "author comment",
        "recommendation:", "decision:",
    ]
    return has_subject and has_signal and not any(e in t for e in excludes)


# ---------------------------------------------------------------------------
# Axis A: RCSB structural metrics
# ---------------------------------------------------------------------------

def _extract_rcsb_structural_row(
    entry_id: str, payload: Dict[str, Any], source_file: str
) -> Dict[str, str]:
    """Extract structural quality metrics from one RCSB core entry payload.

    Returns a row with:
        - entry_id, method, resolution_angstrom
        - wilson_b_estimate (proxy for thermal/dynamic disorder → informs eta)
        - rfree, rfactor (quality indicators)
        - em_resolution (for cryo-EM entries)
        - experimental_method_count
    """
    # Experimental method
    exptl = payload.get("exptl", [])
    method = str(exptl[0].get("method", "")) if exptl else ""

    # Resolution
    res_combined = _safe_get(payload, ["rcsb_entry_info", "resolution_combined"], [])
    resolution = str(res_combined[0]) if res_combined else ""

    # Wilson B (thermal disorder proxy) — from validation report
    vrpt_diff = payload.get("pdbx_vrpt_summary_diffraction", [])
    wilson_b = ""
    rfree = ""
    rfactor = ""
    if isinstance(vrpt_diff, list) and vrpt_diff:
        vd = vrpt_diff[0]
        wb = vd.get("Wilson_B_estimate")
        wilson_b = str(wb) if wb is not None else ""
        rfree_val = vd.get("DCC_Rfree")
        rfree = str(rfree_val) if rfree_val is not None else ""
        r_val = vd.get("DCC_R")
        rfactor = str(r_val) if r_val is not None else ""

    # R-factors from refine block (if vrpt not available)
    if not rfree:
        refine = payload.get("refine", [])
        if refine:
            rfree = str(refine[0].get("ls_R_factor_R_free", ""))
            rfactor = str(refine[0].get("ls_R_factor_obs", ""))

    # EM resolution
    em_res = ""
    em_3d = payload.get("em_3d_reconstruction", [])
    if isinstance(em_3d, list) and em_3d:
        em_res = str(em_3d[0].get("resolution", ""))

    # Deposited atom count (proxy for system size)
    deposited_atoms = str(
        _safe_get(payload, ["rcsb_entry_info", "deposited_atom_count"], "")
    )

    # B-factor type flag
    bf_type = ""
    if isinstance(vrpt_diff, list) and vrpt_diff:
        bf_type = str(vrpt_diff[0].get("B_factor_type", ""))

    return {
        "dataset_id": f"rcsb::{entry_id}",
        "source_file": source_file,
        "entry_id": entry_id,
        "method": method,
        "resolution_angstrom": resolution,
        "em_resolution_angstrom": em_res,
        "wilson_b_estimate": wilson_b,
        "b_factor_type": bf_type,
        "rfree": rfree,
        "rfactor": rfactor,
        "deposited_atoms": deposited_atoms,
    }


def curate_structures(
    rows: List[Dict[str, str]], payload: Dict[str, Any], source_file: str
) -> None:
    """Curate RCSB search results (identifier list only)."""
    for item in payload.get("result_set", []):
        rows.append(
            {
                "dataset_id": f"rcsb::{item.get('identifier', '')}",
                "source_file": source_file,
                "entry_id": str(item.get("identifier", "")),
                "method": "",
                "resolution_angstrom": "",
                "em_resolution_angstrom": "",
                "wilson_b_estimate": "",
                "b_factor_type": "",
                "rfree": "",
                "rfactor": "",
                "deposited_atoms": "",
            }
        )


def curate_structures_core(
    rows: List[Dict[str, str]], payload: Dict[str, Any], source_file: str
) -> None:
    """Curate RCSB core entries (full payload with quality metrics)."""
    for item in payload.get("records", []):
        entry_id = str(item.get("entry_id", ""))
        if item.get("status") != "ok":
            continue
        core = item.get("payload", {})
        row = _extract_rcsb_structural_row(entry_id, core, source_file)
        rows.append(row)


# ---------------------------------------------------------------------------
# Axis A summary: bath parameter proxy from B-factors
# ---------------------------------------------------------------------------

def compute_structural_summary(structures: List[Dict[str, str]]) -> Dict[str, str]:
    """Compute structural disorder summary to constrain bath coupling eta.

    Physical interpretation:
        Wilson_B ~ 8*pi^2 * <u^2>  where <u^2> is mean-square displacement.
        For Ohmic bath: eta ~ gamma * M * omega_c / (pi * hbar)
        B-factor provides an empirical upper bound on <u^2> which enters
        the spectral density normalization.

    Returns:
        Dict with keys: n_bfactor, bfactor_median, bfactor_iqr, bfactor_mean,
                        n_resolution, resolution_median, n_xray, n_em, n_nmr,
                        eta_proxy_low, eta_proxy_high (dimensionless)
    """
    wb_vals = []
    res_vals = []
    methods: Dict[str, int] = {}
    for row in structures:
        wb = _safe_float(row.get("wilson_b_estimate", ""))
        if not math.isnan(wb) and wb > 0:
            wb_vals.append(wb)
        res = _safe_float(row.get("resolution_angstrom", ""))
        if not math.isnan(res) and 0 < res < 100:
            res_vals.append(res)
        m = row.get("method", "").upper()
        if m:
            methods[m] = methods.get(m, 0) + 1

    summary: Dict[str, str] = {}
    if wb_vals:
        wb_sorted = sorted(wb_vals)
        n = len(wb_vals)
        q1 = wb_sorted[n // 4]
        q3 = wb_sorted[3 * n // 4]
        med = statistics.median(wb_vals)
        mn = statistics.mean(wb_vals)
        # eta proxy: B-factor in Å² → <u²> in Å² → normalize to FMO reference
        # FMO: B ~ 12 Å² → eta_FMO ~ 0.15 (literature)
        # Tubulin: B_median → eta_proxy = eta_FMO * (B_median / B_FMO)
        # This is a linear-response estimate, not a derived MD value.
        eta_fmo_ref = 0.15
        b_fmo_ref = 12.0
        eta_proxy = eta_fmo_ref * (med / b_fmo_ref)
        # Bound by [0.1, 1.0] since that is the literature range for proteins
        eta_low = max(0.1, round(eta_proxy * 0.5, 3))
        eta_high = min(1.0, round(eta_proxy * 2.0, 3))
        summary.update({
            "n_bfactor": str(n),
            "bfactor_mean": f"{mn:.2f}",
            "bfactor_median": f"{med:.2f}",
            "bfactor_q1": f"{q1:.2f}",
            "bfactor_q3": f"{q3:.2f}",
            "bfactor_iqr": f"{q3 - q1:.2f}",
            "eta_proxy_low": str(eta_low),
            "eta_proxy_high": str(eta_high),
            "eta_proxy_central": f"{min(1.0, max(0.1, round(eta_proxy, 3))):.3f}",
            "eta_proxy_method": "BFactor_linear_response_FMO_rescaling",
        })
    else:
        summary.update({
            "n_bfactor": "0",
            "bfactor_mean": "",
            "bfactor_median": "",
            "bfactor_q1": "",
            "bfactor_q3": "",
            "bfactor_iqr": "",
            "eta_proxy_low": "",
            "eta_proxy_high": "",
            "eta_proxy_central": "",
            "eta_proxy_method": "",
        })

    if res_vals:
        summary["n_resolution"] = str(len(res_vals))
        summary["resolution_median"] = f"{statistics.median(res_vals):.2f}"
        summary["resolution_mean"] = f"{statistics.mean(res_vals):.2f}"
    else:
        summary["n_resolution"] = "0"
        summary["resolution_median"] = ""
        summary["resolution_mean"] = ""

    for tag, key in [
        ("X-RAY DIFFRACTION", "n_xray"),
        ("ELECTRON MICROSCOPY", "n_em"),
        ("SOLUTION NMR", "n_nmr"),
    ]:
        summary[key] = str(methods.get(tag, 0))

    return summary


# ---------------------------------------------------------------------------
# Axis B: Literature studies
# ---------------------------------------------------------------------------

def curate_studies(
    rows: List[Dict[str, str]],
    payload: Dict[str, Any],
    source_file: str,
    source_tag: str,
) -> None:
    for page in _iter_payload_pages(payload):
        if source_tag == "openalex":
            works = page.get("results", [])
            for w in works:
                doi = str(w.get("doi", "")).replace("https://doi.org/", "")
                title = str(w.get("title", ""))
                if _is_relevant_study(title):
                    rows.append(
                        {
                            "dataset_id": f"openalex::{w.get('id', '')}",
                            "source_file": source_file,
                            "doi": doi,
                            "title": title,
                            "year": str(w.get("publication_year", "")),
                            "technique_tag": _infer_technique_tag(title),
                        }
                    )
        elif source_tag == "crossref":
            works = _safe_get(page, ["message", "items"], [])
            for w in works:
                doi = str(w.get("DOI", ""))
                title_list = w.get("title", [])
                title = str(title_list[0]) if isinstance(title_list, list) and title_list else ""
                if _is_relevant_study(title):
                    rows.append(
                        {
                            "dataset_id": f"crossref::{doi}",
                            "source_file": source_file,
                            "doi": doi,
                            "title": title,
                            "year": str(
                                _safe_get(w, ["issued", "date-parts"], [[""]])[0][0]
                            ),
                            "technique_tag": _infer_technique_tag(title),
                        }
                    )
        elif source_tag == "europepmc":
            works = _safe_get(page, ["resultList", "result"], [])
            for w in works:
                title = str(w.get("title", ""))
                if not _is_relevant_study(title):
                    continue
                doi = str(w.get("doi", ""))
                rows.append(
                    {
                        "dataset_id": f"europepmc::{w.get('id', '')}",
                        "source_file": source_file,
                        "doi": doi,
                        "title": title,
                        "year": str(w.get("pubYear", "")),
                        "technique_tag": _infer_technique_tag(title),
                    }
                )


# ---------------------------------------------------------------------------
# Axis C: PubChem ligand spectral features
# ---------------------------------------------------------------------------

def _extract_spectral_quantities(
    text: str,
) -> Tuple[List[float], str, str, str, str, str, str]:
    """Extract UV/Vis peaks and related spectral data from free text."""
    wl_all = [float(x) for x in re.findall(r"(\d+(?:\.\d+)?)\s*nm", text, flags=re.IGNORECASE)]
    wl = [x for x in wl_all if 150.0 <= x <= 1200.0]

    # Convert wavenumbers to nm
    wn_matches = re.findall(r"(\d+(?:\.\d+)?)\s*cm-?1", text, flags=re.IGNORECASE)
    for wn in wn_matches:
        try:
            val = float(wn)
            if val > 0:
                wl.append(round(1e7 / val, 2))
        except (ValueError, ZeroDivisionError):
            continue

    fwhm_match = re.search(
        r"(?:FWHM|bandwidth|width)[^0-9]*(\d+(?:\.\d+)?)\s*(?:nm|cm-?1)", text, flags=re.IGNORECASE
    )
    snr_match = re.search(r"(?:SNR|S/N|ratio)[^0-9]*(\d+(?:\.\d+)?)", text, flags=re.IGNORECASE)
    bic_match = re.search(
        r"(?:Δ?\s*BIC|delta\s*BIC|Bayes)[^\-0-9]*(-?\d+(?:\.\d+)?)", text, flags=re.IGNORECASE
    )
    peak_match = re.search(
        r"(?:peak|max(?:imum)?\.?(?:\s+at)?|λmax|absorption)[^0-9]*(\d+(?:\.\d+)?)\s*nm",
        text,
        flags=re.IGNORECASE,
    )
    temp_match = re.search(r"(\d+(?:\.\d+)?)\s*(?:K|C|deg\s*C)", text, flags=re.IGNORECASE)

    fwhm = fwhm_match.group(1) if fwhm_match else ""
    snr = snr_match.group(1) if snr_match else ""
    delta_bic = bic_match.group(1) if bic_match else ""
    peak_nm = peak_match.group(1) if peak_match else ""
    temp = temp_match.group(1) if temp_match else ""
    baseline_note = (
        "derivative_or_reported"
        if ("baseline" in text.lower() or "derivative" in text.lower())
        else ""
    )
    return wl, fwhm, snr, temp, delta_bic, peak_nm, baseline_note


def curate_spectral_pubchem(
    rows: List[Dict[str, str]],
    payload: Dict[str, Any],
    source_file: str,
    dataset_id: str,
) -> None:
    """Extract UV/Vis and spectral features from PubChem PUG-View payload."""
    for text in _iter_strings(payload):
        wl_all = [
            float(x)
            for x in re.findall(r"(\d+(?:\.\d+)?)\s*nm", text, flags=re.IGNORECASE)
        ]
        wl = [x for x in wl_all if 150.0 <= x <= 1200.0]

        peak_matches = re.findall(
            r"(?:peak|max(?:imum)?\.?(?:\s+at)?|λmax|absorption)[^0-9]*(\d+(?:\.\d+)?)\s*nm",
            text,
            flags=re.IGNORECASE,
        )
        fwhm_match = re.search(
            r"(?:FWHM|bandwidth|width)[^0-9]*(\d+(?:\.\d+)?)\s*(?:nm|cm-?1)",
            text,
            flags=re.IGNORECASE,
        )
        snr_match = re.search(r"(?:SNR|S/N|ratio)[^0-9]*(\d+(?:\.\d+)?)", text, flags=re.IGNORECASE)
        bic_match = re.search(
            r"(?:Δ?\s*BIC|delta\s*BIC|Bayes)[^\-0-9]*(-?\d+(?:\.\d+)?)",
            text,
            flags=re.IGNORECASE,
        )
        temp_match = re.search(
            r"(\d+(?:\.\d+)?)\s*(?:K|C|deg\s*C)", text, flags=re.IGNORECASE
        )

        fwhm = fwhm_match.group(1) if fwhm_match else ""
        snr = snr_match.group(1) if snr_match else ""
        delta_bic = bic_match.group(1) if bic_match else ""
        temp = temp_match.group(1) if temp_match else ""
        baseline_note = (
            "derivative_or_reported"
            if ("baseline" in text.lower() or "derivative" in text.lower())
            else ""
        )

        for p_nm in (peak_matches if peak_matches else [""]):
            delta = ""
            if len(wl) >= 2:
                delta_val = max(wl) - min(wl)
                if 0.5 < delta_val < 300.0:
                    delta = f"{delta_val:.6g}"

            if not (delta or fwhm or snr or delta_bic or p_nm):
                continue

            rows.append(
                {
                    "dataset_id": dataset_id,
                    "delta_lambda_nm": delta,
                    "peak_nm": p_nm,
                    "fwhm_nm": fwhm,
                    "snr": snr,
                    "delta_bic": delta_bic,
                    "temperature_k": temp,
                    "baseline_note": baseline_note,
                    "prep_label": "pubchem_pugview",
                }
            )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    analysis_dir = Path("analysis")
    analysis_dir.mkdir(parents=True, exist_ok=True)
    print("[curate_compact] START")
    _write_progress("start", {})

    latest = _iter_latest_raw_json()
    if not isinstance(latest, Path):
        print("[curate_compact] ERROR: no raw run directory found")
        return

    all_json = sorted([p for p in latest.glob("*.json") if p.name != "fetch_summary.json"])
    total_json = len(all_json)
    print(f"[curate_compact] run directory: {latest.name}, files: {total_json}")

    structures: List[Dict[str, str]] = []
    studies: List[Dict[str, str]] = []
    spectral: List[Dict[str, str]] = []

    for idx, json_path in enumerate(all_json, start=1):
        payload = _read_json(json_path)
        name = json_path.stem

        if name == "rcsb_tubulin_core_entries":
            curate_structures_core(structures, payload, str(json_path))
        elif name.startswith("rcsb_"):
            curate_structures(structures, payload, str(json_path))
        elif name.startswith("openalex_"):
            curate_studies(studies, payload, str(json_path), "openalex")
        elif name.startswith("crossref_"):
            curate_studies(studies, payload, str(json_path), "crossref")
        elif name.startswith("europepmc_"):
            curate_studies(studies, payload, str(json_path), "europepmc")
        elif name.startswith("pubchem_compound_"):
            cid = name.replace("pubchem_compound_", "").replace("_pugview", "")
            curate_spectral_pubchem(spectral, payload, str(json_path), f"pubchem::{cid}")

        if idx == 1 or idx % 10 == 0 or idx == total_json:
            print(f"[curate_compact] progress {idx}/{total_json}")
            _write_progress(
                "curation",
                {
                    "current": str(idx),
                    "total": str(total_json),
                    "structures": str(len(structures)),
                    "studies": str(len(studies)),
                    "spectral": str(len(spectral)),
                },
            )

    # Deduplicate structures: prefer core entries (have wilson_b) over search results
    dedup_struct: Dict[str, Dict[str, str]] = {}
    for row in structures:
        key = row["dataset_id"]
        if key not in dedup_struct:
            dedup_struct[key] = row
            continue
        prev = dedup_struct[key]
        prev_score = (
            int(bool(prev.get("method")))
            + int(bool(prev.get("resolution_angstrom")))
            + int(bool(prev.get("wilson_b_estimate")))
        )
        cur_score = (
            int(bool(row.get("method")))
            + int(bool(row.get("resolution_angstrom")))
            + int(bool(row.get("wilson_b_estimate")))
        )
        if cur_score >= prev_score:
            dedup_struct[key] = row
    structures = list(dedup_struct.values())

    # -----------------------------------------------------------------------
    # Write structures_compact.csv (full structural metrics)
    # -----------------------------------------------------------------------
    struct_fields = [
        "dataset_id", "source_file", "entry_id", "method",
        "resolution_angstrom", "em_resolution_angstrom",
        "wilson_b_estimate", "b_factor_type",
        "rfree", "rfactor", "deposited_atoms",
    ]
    with (analysis_dir / "structures_compact.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=struct_fields)
        w.writeheader()
        for row in structures:
            w.writerow({k: row.get(k, "") for k in struct_fields})

    # -----------------------------------------------------------------------
    # Write bath_params_proxy.csv (B-factor → eta proxy)
    # Interpretation note is embedded in CSV as comment rows
    # -----------------------------------------------------------------------
    structural_summary = compute_structural_summary(structures)
    bath_params_path = analysis_dir / "bath_params_proxy.csv"
    bath_fields = list(structural_summary.keys())
    with bath_params_path.open("w", newline="", encoding="utf-8") as f:
        w2 = csv.DictWriter(f, fieldnames=bath_fields)
        w2.writeheader()
        w2.writerow(structural_summary)
    print(f"[curate_compact] bath_params_proxy: {structural_summary}")

    # -----------------------------------------------------------------------
    # Write studies_compact.csv
    # -----------------------------------------------------------------------
    with (analysis_dir / "studies_compact.csv").open("w", newline="", encoding="utf-8") as f:
        w3 = csv.DictWriter(
            f, fieldnames=["dataset_id", "source_file", "doi", "title", "year", "technique_tag"]
        )
        w3.writeheader()
        w3.writerows(studies)

    # -----------------------------------------------------------------------
    # Write spectral_compact.csv
    # -----------------------------------------------------------------------
    spectral_fields = [
        "dataset_id", "delta_lambda_nm", "peak_nm", "fwhm_nm",
        "snr", "delta_bic", "temperature_k", "baseline_note", "prep_label",
    ]
    with (analysis_dir / "spectral_compact.csv").open("w", newline="", encoding="utf-8") as f:
        w4 = csv.DictWriter(f, fieldnames=spectral_fields)
        w4.writeheader()
        w4.writerows(spectral)

    print(
        f"[curate_compact] structures={len(structures)} "
        f"studies={len(studies)} spectral={len(spectral)}"
    )
    (analysis_dir / "_done_curate_compact.flag").write_text(
        datetime.now(timezone.utc).isoformat(), encoding="utf-8"
    )
    _write_progress(
        "done",
        {
            "structures": str(len(structures)),
            "studies": str(len(studies)),
            "spectral": str(len(spectral)),
        },
    )
    print("[curate_compact] END")


if __name__ == "__main__":
    main()
