"""Build compact provenance registry for public-data ingestion artifacts."""

from __future__ import annotations
import sys

# Boilerplate para resolver importaciones y rutas desde la raiz del proyecto
PROJECT_ROOT = Path(__file__).resolve().parents[3] # retrocede desde src/scripts/scripts_data_real/ a la raiz
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_progress(step: str, payload: Dict[str, str]) -> None:
    progress_dir = PROJECT_ROOT / "outputs_data" / "raw_json"
    progress_dir.mkdir(parents=True, exist_ok=True)
    p = progress_dir / "_progress_build_registry.json"
    blob: Dict[str, str] = {
        "script": "build_registry.py",
        "step": step,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    blob.update(payload)
    p.write_text(str(blob), encoding="utf-8")


def main() -> None:
    raw_root = PROJECT_ROOT / "data" / "raw" / "public_api"
    out_dir = PROJECT_ROOT / "outputs_data" / "raw_csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, str]] = []
    print("[build_registry] START")
    _write_progress("start", {})

    if not raw_root.exists():
        print("[build_registry] no raw public_api directory found")
        return

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
    main()
