"""validate_integrity.py — SHA-256 integrity utilities for KwanTube artifacts.

Patch v3.5.1 requires that generated binary artifacts can be checked against
sidecar ``.sha256`` files. This module provides small, dependency-free helpers
for creating and validating those sidecars for ``.npz`` and ``.pkl`` outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class IntegrityRecord:
    """Validation result for one artifact."""

    path: str
    sidecar: str
    sha256: str
    status: str
    message: str

    def to_dict(self) -> dict[str, str]:
        """Return a JSON-serializable representation."""

        return asdict(self)


def compute_file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Compute SHA-256 digest for a file without loading it fully into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sidecar_path(path: Path) -> Path:
    """Return canonical sidecar path, e.g. ``H_1JFF.npz.sha256``."""

    return path.with_suffix(path.suffix + ".sha256")


def write_sidecar(path: Path, overwrite: bool = False) -> IntegrityRecord:
    """Create a SHA-256 sidecar for ``path``.

    Parameters
    ----------
    path:
        Artifact to hash.
    overwrite:
        If False, an existing sidecar is left untouched and reported.
    """

    if not path.exists():
        raise FileNotFoundError(path)
    sidecar = sidecar_path(path)
    sha = compute_file_sha256(path)
    if sidecar.exists() and not overwrite:
        existing = sidecar.read_text(encoding="utf-8").strip().split()[0]
        status = "exists_ok" if existing == sha else "exists_mismatch"
        return IntegrityRecord(str(path), str(sidecar), sha, status, "sidecar already exists")
    sidecar.write_text(f"{sha}  {path.name}\n", encoding="utf-8")
    return IntegrityRecord(str(path), str(sidecar), sha, "written", "sidecar written")


def verify_file(path: Path, require_sidecar: bool = True) -> IntegrityRecord:
    """Verify one artifact against its SHA-256 sidecar."""

    sidecar = sidecar_path(path)
    if not path.exists():
        return IntegrityRecord(str(path), str(sidecar), "", "missing", "artifact missing")
    sha = compute_file_sha256(path)
    if not sidecar.exists():
        status = "missing_sidecar" if require_sidecar else "unchecked"
        return IntegrityRecord(str(path), str(sidecar), sha, status, "sidecar missing")
    expected = sidecar.read_text(encoding="utf-8").strip().split()[0]
    if expected == sha:
        return IntegrityRecord(str(path), str(sidecar), sha, "ok", "checksum verified")
    return IntegrityRecord(str(path), str(sidecar), sha, "mismatch", f"expected {expected}")


def is_temporary_artifact(path: Path) -> bool:
    """Return True for crash/WAL temporary files that should not require sidecars."""

    name = path.name.lower()
    return name.startswith("tmp") or name.startswith(".") or ".tmp" in name or name.endswith(".tmp")


def iter_artifacts(root: Path, patterns: Iterable[str] = ("*.npz", "*.pkl"), include_temporary: bool = False) -> Iterable[Path]:
    """Yield artifacts under ``root`` matching the configured patterns."""

    for pattern in patterns:
        for path in sorted(root.rglob(pattern)):
            if include_temporary or not is_temporary_artifact(path):
                yield path


def validate_tree(root: Path, require_sidecar: bool = True, include_temporary: bool = False) -> list[IntegrityRecord]:
    """Validate every supported artifact under ``root``."""

    return [verify_file(path, require_sidecar=require_sidecar) for path in iter_artifacts(root, include_temporary=include_temporary)]


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for checksum generation and validation."""

    parser = argparse.ArgumentParser(description="Validate KwanTube artifact checksums")
    parser.add_argument("root", nargs="?", default=str(PROJECT_ROOT / "outputs_data"))
    parser.add_argument("--write-missing", action="store_true", help="write missing sidecars")
    parser.add_argument("--overwrite", action="store_true", help="overwrite sidecars when writing")
    parser.add_argument("--allow-missing-sidecar", action="store_true")
    parser.add_argument("--include-temporary", action="store_true", help="include WAL/temp artifacts such as tmp*.npz and .*.npz")
    parser.add_argument("--strict", action="store_true", help="strict mode: include temporary artifacts and require all sidecars")
    parser.add_argument("--json", type=Path, default=None, help="optional JSON report path")
    args = parser.parse_args(argv)

    root = Path(args.root)
    include_temporary = bool(args.include_temporary or args.strict)
    require_sidecar = bool(args.strict or not args.allow_missing_sidecar)
    if args.write_missing:
        records = [write_sidecar(path, overwrite=args.overwrite) for path in iter_artifacts(root, include_temporary=include_temporary)]
    else:
        records = validate_tree(root, require_sidecar=require_sidecar, include_temporary=include_temporary)

    payload = [record.to_dict() for record in records]
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    bad = [r for r in records if r.status in {"missing", "missing_sidecar", "mismatch", "exists_mismatch"}]
    for record in records:
        print(f"[{record.status}] {record.path}")
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
    raise SystemExit(main(sys.argv[1:]))