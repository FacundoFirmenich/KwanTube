"""Build compact provenance registry for public-data ingestion artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

# ---------------------------------------------------------------------------
# Path resolution
# WORKSPACE_ROOT is 4 levels up: KwanTube/src/scripts/data/ → … → biofisicaquantiqaCLINE/
# ---------------------------------------------------------------------------
_HERE          = Path(__file__).resolve()
KWANTUBE_ROOT  = _HERE.parents[3]   # …/KwanTube/
WORKSPACE_ROOT = _HERE.parents[4]   # …/biofisicaquantiqaCLINE/
PROJECT_ROOT   = KWANTUBE_ROOT      # alias for legacy downstream references

if str(KWANTUBE_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(KWANTUBE_ROOT / "src"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_progress(step: str, payload: Dict[str, str]) -> None:
    progress_dir = KWANTUBE_ROOT / "outputs_data" / "raw_json" / "progress"
    progress_dir.mkdir(parents=True, exist_ok=True)
    p = progress_dir / "_progress_build_registry.json"
    blob: Dict[str, str] = {
        "script": "build_registry.py",
        "step": step,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    blob.update(payload)
    p.write_text(json.dumps(blob, indent=2), encoding="utf-8")


def main() -> None:
    # Search for public_api snapshots in canonical companion-data dir first,
    # then fall back to KwanTube-internal outputs written by fetch_public_data.py.
    raw_root: Path | None = None
    for candidate in (
        WORKSPACE_ROOT / "data_downloaded_public_repos" / "raw" / "public_api",
        KWANTUBE_ROOT  / "outputs_data" / "raw_data" / "public_api",
    ):
        if candidate.exists():
            raw_root = candidate
            break

    out_dir = KWANTUBE_ROOT / "outputs_data" / "raw_csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, str]] = []
    print("[build_registry] START")
    _write_progress("start", {})

    if raw_root is None:
        print("[build_registry] no raw public_api directory found")
        _write_progress("done", {"rows": "0", "output": ""})
        return

    print(f"[build_registry] using raw_root: {raw_root}")
    run_dirs = sorted([p for p in raw_root.iterdir() if p.is_dir()])
    total_runs = len(run_dirs)
    for i, run_dir in enumerate(run_dirs, start=1):
        summary_path = run_dir / "fetch_summary.json"
        if not summary_path.exists():
            continue
        loaded = json.loads(summary_path.read_text(encoding="utf-8"))
        summary = loaded.get("summary", loaded) if isinstance(loaded, dict) else loaded
        if not isinstance(summary, list):
            summary = []
        print(f"[build_registry] run {i}/{total_runs}: {run_dir.name} items={len(summary)}")
        for item in summary:
            file_path = Path(item.get("file", "")) if item.get("file") else None
            checksum = sha256_file(file_path) if file_path and file_path.exists() else ""
            rows.append(
                {
                    "run_id": run_dir.name,
                    "source": str(item.get("source", "")),
                    "name": str(item.get("name", "")),
                    "url": str(item.get("url", "")),
                    "status": str(item.get("status", "")),
                    "file": str(item.get("file", "")),
                    "sha256": checksum,
                    "timestamp_utc": str(item.get("timestamp_utc", "")),
                }
            )
        _write_progress("processing_runs", {"current": str(i), "total": str(total_runs), "rows": str(len(rows))})

    out_csv = out_dir / "data_registry.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["run_id", "source", "name", "url", "status", "file", "sha256", "timestamp_utc"],
        )
        w.writeheader()
        w.writerows(rows)

    print(f"[build_registry] wrote {len(rows)} rows to {out_csv}")
    (out_dir / "_done_build_registry.flag").write_text(datetime.now(timezone.utc).isoformat(), encoding="utf-8")
    _write_progress("done", {"rows": str(len(rows)), "output": str(out_csv)})
    print("[build_registry] END")


if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    main()
