"""seal_outputs.py — Regenerate all SHA-256 sidecars for KwanTube artifacts.

Administrative tool: synchronizes the integrity ledger after any authorized
structural changes (Hamiltonian updates, directory reorganization, etc.).

Targets ONLY the file types that validate_integrity.py checks:
  *.npz  and  *.pkl

Flags:
  --dry-run    Print what would be sealed without writing.
  --report     Path to write a JSON summary of the sealing run.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Resolve PROJECT_ROOT: seal_outputs.py lives at src/scripts/validation/
PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.validate_integrity import (
    compute_file_sha256,
    iter_artifacts,
    sidecar_path,
    verify_file,
    write_sidecar,
)


OUTPUTS_ROOT = PROJECT_ROOT / "outputs_data"


def seal_all(root: Path, dry_run: bool = False) -> dict:
    """Write or overwrite sidecars for every .npz and .pkl artifact under root.

    Returns a summary dict suitable for JSON serialization.
    """
    sealed: list[str] = []
    skipped: list[str] = []
    errors: list[dict] = []

    artifacts = list(iter_artifacts(root))
    total = len(artifacts)
    print(f"[seal_outputs] Scanning {total} artifacts under {root}")

    for idx, path in enumerate(artifacts, start=1):
        try:
            pre = verify_file(path, require_sidecar=True)
            if pre.status == "ok":
                skipped.append(str(path))
                if idx % 10 == 0 or idx == total:
                    print(f"  [{idx}/{total}] OK (no change needed)  {path.name}")
                continue

            # Needs sealing: missing_sidecar or mismatch
            if dry_run:
                print(f"  [{idx}/{total}] DRY-RUN would seal  {path.name}  (status={pre.status})")
                sealed.append(str(path))
                continue

            rec = write_sidecar(path, overwrite=True)
            sealed.append(str(path))
            print(f"  [{idx}/{total}] SEALED  {path.name}  sha256={rec.sha256[:16]}...")

        except Exception as exc:
            errors.append({"path": str(path), "error": repr(exc)})
            print(f"  [{idx}/{total}] ERROR  {path.name}  {exc}", file=sys.stderr)

    print(f"\n[seal_outputs] Done.  sealed={len(sealed)}  skipped={len(skipped)}  errors={len(errors)}")
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "dry_run": dry_run,
        "root": str(root),
        "total_artifacts": total,
        "sealed": sealed,
        "skipped": skipped,
        "errors": errors,
    }


def validate_final(root: Path) -> tuple[int, int]:
    """Run a final validation pass and return (validated_count, bad_count)."""
    records = [verify_file(p) for p in iter_artifacts(root)]
    bad = [r for r in records if r.status in {"missing", "missing_sidecar", "mismatch"}]
    print(f"\n[seal_outputs] Final validation: validated={len(records)}  bad={len(bad)}")
    if bad:
        for r in bad:
            print(f"  [BAD]  {r.status}  {r.path}")
    else:
        print("  All artifacts CLEAN.")
    return len(records), len(bad)


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Seal KwanTube output artifacts")
    parser.add_argument("--root", type=Path, default=OUTPUTS_ROOT,
                        help=f"Root to scan (default: {OUTPUTS_ROOT})")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would be sealed without writing")
    parser.add_argument("--report", type=Path, default=None,
                        help="Optional path to write a JSON seal report")
    args = parser.parse_args()

    summary = seal_all(args.root, dry_run=args.dry_run)

    if not args.dry_run:
        validated, bad = validate_final(args.root)
        summary["final_validated"] = validated
        summary["final_bad"] = bad

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(f"\n[seal_outputs] Report written to: {args.report}")

    return 0 if summary.get("final_bad", 0) == 0 and not summary["errors"] else 1


if __name__ == "__main__":
    for _parent in Path(__file__).resolve().parents:
        if (_parent / "qmc_mt" / "run_audit.py").exists():
            sys.path.insert(0, str(_parent))
            break
    from qmc_mt.run_audit import install_run_audit as _install_run_audit
    _install_run_audit(__file__)
    raise SystemExit(main())
