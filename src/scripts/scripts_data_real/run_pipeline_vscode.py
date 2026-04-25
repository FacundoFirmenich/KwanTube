"""Run the compact pipeline end-to-end with one VS Code Run click.

This runner is intentionally simple: it imports and executes each script module's
``main()`` in deterministic order, while emitting clear stage logs.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import json

import build_registry
import compute_detectability_metrics
import curate_compact
import export_claim_traceability
import fetch_public_data
import run_comparative_panels


def _ts() -> str:
    return datetime.now(timezone.utc).isoformat()


def main() -> None:
    analysis_dir = Path("analysis")
    analysis_dir.mkdir(parents=True, exist_ok=True)
    print(f"[run_pipeline_vscode] START {_ts()}")

    run_report = {
        "script": "run_pipeline_vscode.py",
        "started_at": _ts(),
        "stages": [],
        "status": "running",
    }

    stages = [
        ("fetch_public_data", fetch_public_data.main),
        ("curate_compact", curate_compact.main),
        ("compute_detectability_metrics", compute_detectability_metrics.main),
        ("run_comparative_panels", run_comparative_panels.main),
        ("build_registry", build_registry.main),
        ("export_claim_traceability", export_claim_traceability.main),
    ]

    for idx, (name, fn) in enumerate(stages, start=1):
        print(f"[run_pipeline_vscode] [{idx}/{len(stages)}] START {name} {_ts()}")
        stage_started = _ts()
        try:
            fn()
            run_report["stages"].append({"name": name, "status": "ok", "started_at": stage_started, "ended_at": _ts()})
            print(f"[run_pipeline_vscode] [{idx}/{len(stages)}] END   {name} {_ts()}")
        except Exception as exc:  # pylint: disable=broad-except
            run_report["stages"].append(
                {
                    "name": name,
                    "status": "error",
                    "started_at": stage_started,
                    "ended_at": _ts(),
                    "error": str(exc),
                }
            )
            run_report["status"] = "error"
            run_report["finished_at"] = _ts()
            (analysis_dir / "pipeline_run_report.json").write_text(json.dumps(run_report, ensure_ascii=False, indent=2), encoding="utf-8")
            print(f"[run_pipeline_vscode] [{idx}/{len(stages)}] FAIL  {name}: {exc}")
            raise

    (analysis_dir / "_done_run_pipeline_vscode.flag").write_text(_ts(), encoding="utf-8")
    run_report["status"] = "ok"
    run_report["finished_at"] = _ts()
    (analysis_dir / "pipeline_run_report.json").write_text(json.dumps(run_report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[run_pipeline_vscode] END {_ts()}")


if __name__ == "__main__":
    main()
