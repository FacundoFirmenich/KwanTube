"""Validate KwanTube generated outputs for required integrity metadata."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.validate_integrity import validate_tree


def main() -> int:
    """Validate outputs_data checksums and write a compact JSON report."""

    parser = argparse.ArgumentParser(description="Validate KwanTube outputs")
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT / "outputs_data")
    parser.add_argument("--allow-missing-sidecar", action="store_true")
    parser.add_argument("--report", type=Path, default=PROJECT_ROOT / "outputs_data" / "raw_json" / "structural" / "outputs_validation_report.json")
    args = parser.parse_args()

    records = validate_tree(args.root, require_sidecar=not args.allow_missing_sidecar)
    payload = [record.to_dict() for record in records]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    bad = [r for r in records if r.status in {"missing", "missing_sidecar", "mismatch"}]
    print(f"validated={len(records)} bad={len(bad)} report={args.report}")
    return 1 if bad else 0


if __name__ == "__main__":
    from pathlib import Path as _RunAuditPath
    import sys as _run_audit_sys
    for _run_audit_parent in _RunAuditPath(__file__).resolve().parents:
        if (_run_audit_parent / "qmc_mt" / "run_audit.py").exists():
            _run_audit_sys.path.insert(0, str(_run_audit_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    raise SystemExit(main())