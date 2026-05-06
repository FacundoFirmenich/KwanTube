"""Run the compact pipeline end-to-end with one VS Code Run click.

This runner uses subprocess to call each script independently, which:
  - Avoids import-path conflicts between scripts in different folders.
  - Properly isolates each stage's environment.
  - Works correctly after the public_data / analysis folder reorganization.

Script locations after reorganization (v3.5.1+):
  public_data/  : fetch_public_data, curate_compact, build_registry, run_comparative_panels
  analysis/     : compute_detectability_metrics, export_claim_traceability
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
_HERE          = Path(__file__).resolve()
KWANTUBE_ROOT  = _HERE.parents[3]   # …/KwanTube/
WORKSPACE_ROOT = _HERE.parents[4]   # …/biofisicaquantiqaCLINE/
PYTHON         = sys.executable

SCRIPTS_DIR    = KWANTUBE_ROOT / "src" / "scripts"


def _ts() -> str:
    return datetime.now(timezone.utc).isoformat()


def main() -> None:
    analysis_dir = KWANTUBE_ROOT / "outputs_data" / "raw_json" / "structural"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    print(f"[run_pipeline_vscode] START {_ts()}")

    run_report = {
        "script": "run_pipeline_vscode.py",
        "started_at": _ts(),
        "stages": [],
        "status": "running",
        "project_root": str(KWANTUBE_ROOT),
    }

    # NOTE: fetch_public_data is intentionally excluded — data is already present
    # in data_downloaded_public_repos/. Run it manually only when refreshing data.
    stages: list[tuple[str, Path]] = [
        ("curate_compact",              SCRIPTS_DIR / "public_data" / "curate_compact.py"),
        ("build_registry",              SCRIPTS_DIR / "public_data" / "build_registry.py"),
        ("compute_detectability_metrics", SCRIPTS_DIR / "analysis"    / "compute_detectability_metrics.py"),
        ("run_comparative_panels",      SCRIPTS_DIR / "public_data" / "run_comparative_panels.py"),
        ("export_claim_traceability",   SCRIPTS_DIR / "analysis"    / "export_claim_traceability.py"),
    ]

    n = len(stages)
    for idx, (name, script_path) in enumerate(stages, start=1):
        print(f"[run_pipeline_vscode] [{idx}/{n}] START {name}  {_ts()}")
        stage_started = _ts()
        try:
            result = subprocess.run(
                [PYTHON, str(script_path)],
                check=True,
                capture_output=False,
                cwd=str(KWANTUBE_ROOT),
            )
            run_report["stages"].append({
                "name": name,
                "status": "ok",
                "started_at": stage_started,
                "ended_at": _ts(),
            })
            print(f"[run_pipeline_vscode] [{idx}/{n}] END   {name}  {_ts()}")
        except subprocess.CalledProcessError as exc:
            run_report["stages"].append({
                "name": name,
                "status": "error",
                "started_at": stage_started,
                "ended_at": _ts(),
                "returncode": exc.returncode,
            })
            run_report["status"] = "error"
            run_report["finished_at"] = _ts()
            (analysis_dir / "pipeline_run_report.json").write_text(
                json.dumps(run_report, ensure_ascii=False, indent=2), encoding="utf-8"
            )
            print(f"[run_pipeline_vscode] [{idx}/{n}] FAIL  {name}: returncode={exc.returncode}")
            raise

    (analysis_dir / "_done_run_pipeline_vscode.flag").write_text(_ts(), encoding="utf-8")
    run_report["status"] = "ok"
    run_report["finished_at"] = _ts()
    (analysis_dir / "pipeline_run_report.json").write_text(
        json.dumps(run_report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"[run_pipeline_vscode] END {_ts()}")


if __name__ == "__main__":
    main()
