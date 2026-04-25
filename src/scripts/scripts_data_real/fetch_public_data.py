"""Fetch public institutional metadata for the manuscript reproducibility layer.

Production-oriented behavior:
- retry with exponential backoff
- paginated acquisition for major literature endpoints
- structured progress and run summary
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass(frozen=True)
class EndpointSpec:
    """Definition of one API request in the compact pipeline."""

    source: str
    name: str
    url: str


def _jdump(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    _jdump(path, payload)


def _http_get_json(url: str, timeout: int = 30, retries: int = 3, backoff_sec: float = 1.5) -> Dict[str, Any]:
    """Perform a GET request and decode JSON response with retry/backoff."""
    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        req = urllib.request.Request(
            url,
            headers={
                "User-Agent": "biofisicaquantiqaCLINE/production-public-pipeline",
                "Accept": "application/json",
            },
            method="GET",
        )
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = response.read().decode("utf-8")
            return json.loads(payload)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
            last_exc = exc
            if attempt < retries:
                sleep_s = backoff_sec * (2 ** (attempt - 1))
                print(f"[fetch_public_data] retry {attempt}/{retries - 1} after error: {exc}")
                time.sleep(sleep_s)
    if last_exc is None:
        raise RuntimeError("Unexpected fetch failure without captured exception")
    raise last_exc


def _write_progress(step: str, payload: Dict[str, Any]) -> None:
    """Persist lightweight heartbeat so UI users can see liveness."""
    analysis_dir = Path("analysis")
    analysis_dir.mkdir(parents=True, exist_ok=True)
    progress_path = analysis_dir / "_progress_fetch_public_data.json"
    blob = {
        "script": "fetch_public_data.py",
        "step": step,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    blob.update(payload)
    _jdump(progress_path, blob)


def build_endpoints(max_records: int) -> List[EndpointSpec]:
    """Construct singleton endpoints for non-paginated sources."""
    rcsb_search_payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_phrase",
                        "value": "tubulin",
                    },
                }
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_records}},
    }
    encoded = urllib.parse.quote(json.dumps(rcsb_search_payload, separators=(",", ":")))

    # Famous tubulin ligands (inhibitors/stabilizers)
    pubchem_cids = [
        6167, 36314, 5344, 5978, 10657, 5352103, 4122, 32248, 448799, 457813,
        148124, 11354606, 9854073, 10074128, 10321287, 5281862, 9949641,
        11351021, 53379371, 11414799, 6675804, 139369407, 104865, 3033830,
        5353431, 247034, 444695, 126941, 104850, 6440397, 6918290, 6918641
    ]

    # Curated list of 200 high-relevance tubulin inhibitor CIDs from PubChem (extracted via subagent)
    curated_cids = [
        172966289, 168475901, 11545920, 166642397, 168476159, 155804415, 162659980, 169494478, 145985142, 162674242,
        155813536, 6167, 162668817, 155816728, 168510236, 168510240, 168510440, 11609821, 15432611, 24749353,
        53262872, 167312487, 53262871, 162624503, 162665763, 170835989, 170836205, 60148288, 171475511, 66929931,
        168355527, 172966232, 155810671, 170835998, 171037358, 165177911, 166450626, 165412767, 171038014, 155812326,
        169494241, 171393604, 169450481, 168510334, 171713995, 168475882, 168475950, 168476054, 168510481, 165177913,
        165412768, 165413512, 169089875, 171346274, 168679661, 168679749, 156620384, 134868071, 131953485, 169450094,
        169450210, 169450726, 169450802, 163409053, 163409085, 171714110, 162716997, 137645287, 139369407, 137319628,
        165412665, 170836007, 170836279, 163322403, 166176963, 146469200, 125432004, 156019401, 171714317, 177838642,
        9843752, 131974320, 163408760, 163409072, 163409084, 163409092, 166177025, 166177033, 156696051, 172419027,
        169450095, 170453172, 139703315, 163408966, 168873437, 172638531, 172638543, 172677024, 172677108, 166642475,
        135417086, 134131455, 2359994, 5978, 5569, 6438581, 154632611, 3000319, 657289, 9962911,
        10607, 6446789, 13342, 429016, 5388993, 220401, 3717450, 429017, 5388992, 148124,
        28780, 11343137, 5388983, 5281828, 6324671, 184492, 233481, 148123, 2874761, 443593,
        11351021, 447865, 448013, 852489, 5353889, 43116, 43264, 241903, 2082, 4030,
        4122, 10838895, 3143, 241902, 6710780, 5672, 249332, 6477151, 448799, 40839,
        5352092, 10079877, 6437358, 335929, 4232041, 99957, 5354063, 547450, 9832825, 11488110,
        11957714, 11607738, 5365290, 5311497, 10125800, 11972518, 5468287, 29393, 15597809, 9810929,
        5488895, 16750137, 13200033, 6675804, 6711270, 20055492, 16051941, 7404291, 5915307, 11643449,
        17756418, 15560162, 46926355, 24721032, 45055483, 23643564, 21584522, 45073403, 44408087, 21864436,
        44417235, 45356262, 45358014, 5281829, 6305
    ]
    pubchem_cids = sorted(list(set(pubchem_cids + curated_cids)))
    print(f"[fetch_public_data] curated list merged: {len(pubchem_cids)} target compounds")

    return [
        EndpointSpec(
            source="rcsb",
            name="rcsb_tubulin_search",
            url=f"https://search.rcsb.org/rcsbsearch/v2/query?json={encoded}",
        ),
        *[
            EndpointSpec(
                source="pubchem",
                name=f"pubchem_compound_{cid}_pugview",
                url=f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON",
            )
            for cid in pubchem_cids
        ],
    ]


def _fetch_paginated_literature(max_records: int, max_pages: int, out_dir: Path, ts: str) -> List[Dict[str, Any]]:
    """Acquire paginated literature metadata from OpenAlex/Crossref/EuropePMC."""
    query_q = urllib.parse.quote("microtubule spectroscopy quantum")
    query_pmc = urllib.parse.quote('(\"microtubule\" OR \"tubulin\") AND (\"spectroscopy\" OR \"2D-IR\" OR \"THz\" OR \"fluorescence\")')
    summary: List[Dict[str, Any]] = []

    # OpenAlex
    try:
        oa_page_size = min(max_records, 200)
        openalex_pages: List[Dict[str, Any]] = []
        for page in range(1, max_pages + 1):
            url = f"https://api.openalex.org/works?search={query_q}&per-page={oa_page_size}&page={page}"
            print(f"[fetch_public_data] [openalex] page {page}/{max_pages}")
            data = _http_get_json(url)
            openalex_pages.append(data)
            _write_progress("fetch_openalex", {"run_id": ts, "page": page, "max_pages": max_pages})
            if not data.get("results"):
                break
        oa_file = out_dir / "openalex_microtubule_quantum.json"
        _write_json(oa_file, {"pages": openalex_pages, "timestamp_utc": ts})
        summary.append({"source": "openalex", "name": "openalex_microtubule_quantum", "status": "ok", "file": str(oa_file), "timestamp_utc": ts})
    except Exception as exc:
        print(f"[fetch_public_data] [openalex] failed: {exc}")
        summary.append({"source": "openalex", "name": "openalex_microtubule_quantum", "status": "error", "error": str(exc), "timestamp_utc": ts})

    # Crossref: offset pagination (capped at 10k by API)
    try:
        cr_page_size = min(max_records, 1000)
        crossref_pages: List[Dict[str, Any]] = []
        cr_max_pages = min(max_pages, 10) # offset 10000 is often the limit
        for page in range(0, cr_max_pages):
            offset = page * cr_page_size
            url = f"https://api.crossref.org/works?query={urllib.parse.quote('microtubule spectroscopy')}&rows={cr_page_size}&offset={offset}"
            print(f"[fetch_public_data] [crossref] page {page + 1}/{cr_max_pages} offset={offset}")
            data = _http_get_json(url)
            crossref_pages.append(data)
            _write_progress("fetch_crossref", {"run_id": ts, "page": page + 1, "max_pages": cr_max_pages})
            items = (((data.get("message") or {}).get("items")) or [])
            if not items:
                break
        cr_file = out_dir / "crossref_microtubule_spectroscopy.json"
        _write_json(cr_file, {"pages": crossref_pages, "timestamp_utc": ts})
        summary.append({"source": "crossref", "name": "crossref_microtubule_spectroscopy", "status": "ok", "file": str(cr_file), "timestamp_utc": ts})
    except Exception as exc:
        print(f"[fetch_public_data] [crossref] failed: {exc}")
        summary.append({"source": "crossref", "name": "crossref_microtubule_spectroscopy", "status": "error", "error": str(exc), "timestamp_utc": ts})

    # EuropePMC: page pagination
    try:
        ep_page_size = min(max_records, 1000)
        ep_pages: List[Dict[str, Any]] = []
        for page in range(1, max_pages + 1):
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={query_pmc}&page={page}&pageSize={ep_page_size}&format=json"
            print(f"[fetch_public_data] [europepmc] page {page}/{max_pages}")
            data = _http_get_json(url)
            ep_pages.append(data)
            _write_progress("fetch_europepmc", {"run_id": ts, "page": page, "max_pages": max_pages})
            results = (((data.get("resultList") or {}).get("result")) or [])
            if not results:
                break
        ep_file = out_dir / "europepmc_microtubule_spectroscopy.json"
        _write_json(ep_file, {"pages": ep_pages, "timestamp_utc": ts})
        summary.append({"source": "europepmc", "name": "europepmc_microtubule_spectroscopy", "status": "ok", "file": str(ep_file), "timestamp_utc": ts})
    except Exception as exc:
        print(f"[fetch_public_data] [europepmc] failed: {exc}")
        summary.append({"source": "europepmc", "name": "europepmc_microtubule_spectroscopy", "status": "error", "error": str(exc), "timestamp_utc": ts})

    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch compact public metadata for manuscript traceability.")
    parser.add_argument("--max-records", type=int, default=1000, help="Max records per endpoint.")
    parser.add_argument("--max-pages", type=int, default=15, help="Max pages for paginated literature sources.")
    parser.add_argument("--threads", type=int, default=16, help="Concurrent fetch threads.")
    args = parser.parse_args()

    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out_dir = Path("data") / "raw" / "public_api" / ts
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[fetch_public_data] START run={ts} max_records={args.max_records} threads={args.threads}")
    _write_progress("start", {"run_id": ts, "max_records": args.max_records})

    summary: List[Dict[str, Any]] = []

    # Paginated literature sources first.
    try:
        summary.extend(_fetch_paginated_literature(args.max_records, args.max_pages, out_dir, ts))
    except Exception as exc:  # pylint: disable=broad-except
        print(f"[fetch_public_data] paginated literature fetch failed: {exc}")
        summary.append({"source": "literature", "name": "paginated_fetch", "status": "error", "error": str(exc), "timestamp_utc": ts})

    endpoints = build_endpoints(args.max_records)
    total_endpoints = len(endpoints)
    
    print(f"[fetch_public_data] executing parallel endpoints: n={total_endpoints}")
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_spec = {executor.submit(_http_get_json, spec.url): spec for spec in endpoints}
        for idx, future in enumerate(as_completed(future_to_spec), start=1):
            spec = future_to_spec[future]
            item: Dict[str, Any] = {
                "source": spec.source,
                "name": spec.name,
                "url": spec.url,
                "timestamp_utc": ts,
                "status": "ok",
            }
            try:
                data = future.result()
                target = out_dir / f"{spec.name}.json"
                _write_json(target, data)
                item["file"] = str(target)
            except Exception as exc:
                item["status"] = "error"
                item["error"] = str(exc)
            
            summary.append(item)
            if idx == 1 or idx % 10 == 0 or idx == total_endpoints:
                print(f"[fetch_public_data] endpoint progress {idx}/{total_endpoints}")
                _write_progress(
                    "fetch_endpoints",
                    {
                        "run_id": ts,
                        "current": idx,
                        "total": total_endpoints,
                        "name": spec.name,
                        "status": item["status"],
                    },
                )

    # RCSB core-entry enrichment for structure-level fields used downstream.
    rcsb_search_file = out_dir / "rcsb_tubulin_search.json"
    if rcsb_search_file.exists():
        try:
            payload = json.loads(rcsb_search_file.read_text(encoding="utf-8"))
            identifiers = [x.get("identifier", "") for x in payload.get("result_set", []) if x.get("identifier")]
            core_records: List[Dict[str, Any]] = []
            total_ids = len(identifiers)
            print(f"[fetch_public_data] enriching RCSB core entries (parallel): n={total_ids}")
            
            with ThreadPoolExecutor(max_workers=args.threads) as executor:
                future_to_id = {executor.submit(_http_get_json, f"https://data.rcsb.org/rest/v1/core/entry/{eid}"): eid for eid in identifiers}
                for idx, future in enumerate(as_completed(future_to_id), start=1):
                    entry_id = future_to_id[future]
                    try:
                        core_records.append({
                            "entry_id": entry_id, 
                            "payload": future.result(), 
                            "url": f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}", 
                            "status": "ok"
                        })
                    except Exception as exc:
                        core_records.append({
                            "entry_id": entry_id, 
                            "url": f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}", 
                            "status": "error", 
                            "error": str(exc)
                        })
                    
                    if idx == 1 or idx % 50 == 0 or idx == total_ids:
                        print(f"[fetch_public_data] core enrichment progress {idx}/{total_ids}")
                        _write_progress(
                            "fetch_rcsb_core",
                            {
                                "run_id": ts,
                                "current": idx,
                                "total": total_ids,
                            },
                        )

            core_file = out_dir / "rcsb_tubulin_core_entries.json"
            _write_json(core_file, {"records": core_records, "timestamp_utc": ts})
            summary.append(
                {
                    "source": "rcsb",
                    "name": "rcsb_tubulin_core_entries",
                    "url": "https://data.rcsb.org/rest/v1/core/entry/{entry_id}",
                    "timestamp_utc": ts,
                    "status": "ok",
                    "file": str(core_file),
                }
            )
        except json.JSONDecodeError:
            pass

    _write_json(out_dir / "fetch_summary.json", {"summary": summary, "timestamp_utc": ts})
    done_flag = Path("analysis") / "_done_fetch_public_data.flag"
    done_flag.write_text(ts, encoding="utf-8")
    _write_progress("done", {"run_id": ts, "outputs_dir": str(out_dir), "items": len(summary)})
    print(f"[fetch_public_data] wrote outputs to: {out_dir}")
    print("[fetch_public_data] END")


if __name__ == "__main__":
    main()
