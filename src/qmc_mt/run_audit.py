"""Runtime audit helpers for executable KwanTube scripts.

This module provides a lightweight execution-memory layer for scripts run from
``KwanTube/src``. It records verifiable UTC timestamps, mirrors console output
into an append-only text log, and hashes output artifacts generated during a
run.
"""

from __future__ import annotations

import atexit
import builtins
import hashlib
import json
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


AUDIT_SUFFIXES = {".json", ".npz", ".pkl", ".gz", ".txt", ".md"}


def utc_now_iso() -> str:
    """Return a timezone-aware UTC timestamp in ISO-8601 format."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def compute_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Compute the SHA-256 digest of ``path`` without loading it all at once."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _project_root(script_path: Path) -> Path:
    for parent in [script_path.parent, *script_path.parents]:
        if (parent / "outputs_data").exists() and (parent / "src").exists():
            return parent
        if parent.name == "src" and parent.parent.exists():
            return parent.parent
    return script_path.resolve().parents[2]


def _safe_console_text(value: Any) -> str:
    text = str(value)
    encoding = sys.stdout.encoding or "utf-8"
    return text.encode(encoding, errors="backslashreplace").decode(encoding, errors="replace")


def _iter_outputs(root: Path) -> Iterable[Path]:
    outputs_root = root / "outputs_data"
    if not outputs_root.exists():
        return []
    return (
        path
        for path in outputs_root.rglob("*")
        if path.is_file() and (path.suffix in AUDIT_SUFFIXES or "raw_txt+md" in str(path.parent))
    )


class RunAudit:
    """Append-only execution logger with output hashing."""

    def __init__(self, script_file: str | Path) -> None:
        self.script_path = Path(script_file).resolve()
        self.script_name = self.script_path.name
        self.root = _project_root(self.script_path)
        self.started_monotonic = time.monotonic()
        self.started_wall = time.time()
        self.started_at = utc_now_iso()
        self.log_path = self.root / "outputs_data" / "raw_txt+md" / "logs" / "execution_memory.log.txt"
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        self._original_print = builtins.print
        self._closed = False

    def install(self) -> "RunAudit":
        """Install timestamped, UTF-8-safe print mirroring and exit hook."""
        self.log_event("start", status="running", message=f"script={self.script_path}")
        builtins.print = self.print  # type: ignore[assignment]
        atexit.register(self.close)
        self.print(f"[RUN_AUDIT] START script={self.script_name} timestamp_utc={self.started_at}")
        return self

    def print(self, *args: Any, **kwargs: Any) -> None:
        """Timestamp console output and mirror it to the execution-memory log."""
        sep = kwargs.get("sep", " ")
        end = kwargs.get("end", "\n")
        file_obj = kwargs.get("file", sys.stdout)
        flush = bool(kwargs.get("flush", False))
        safe_args = [_safe_console_text(arg) for arg in args]
        message = sep.join(safe_args)
        prefix = f"[{utc_now_iso()}] [{self.script_name}]"
        if message:
            lines = [f"{prefix} {part}" if part else prefix for part in message.split("\n")]
        else:
            lines = [prefix]
        rendered = "\n".join(lines)
        self._original_print(rendered, end=end, file=file_obj, flush=flush)
        try:
            with self.log_path.open("a", encoding="utf-8", errors="backslashreplace") as handle:
                handle.write(rendered + ("" if end == "" else "\n"))
        except Exception:
            self._original_print(
                f"[{utc_now_iso()}] [run_audit] WARNING unable to write log {self.log_path}",
                file=sys.stderr,
            )

    def log_event(self, event: str, status: str = "ok", **payload: Any) -> None:
        """Append a structured event line to the execution-memory log."""
        record = {
            "timestamp_utc": utc_now_iso(),
            "script": self.script_name,
            "event": event,
            "status": status,
            **payload,
        }
        with self.log_path.open("a", encoding="utf-8", errors="backslashreplace") as handle:
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True, default=str) + "\n")

    def _changed_outputs(self) -> list[Path]:
        changed: list[Path] = []
        for path in _iter_outputs(self.root):
            try:
                if path.stat().st_mtime >= self.started_wall - 1.0 and path != self.log_path:
                    changed.append(path)
            except OSError:
                continue
        return sorted(changed)

    def close(self) -> None:
        """Record final status and SHA-256 hashes for outputs touched by this run."""
        if self._closed:
            return
        self._closed = True
        status = "ok"
        exc_type, exc_value, exc_tb = sys.exc_info()
        if exc_type is not None:
            status = "error"
            self.log_event(
                "exception",
                status="error",
                error="".join(traceback.format_exception(exc_type, exc_value, exc_tb)),
            )
        for output in self._changed_outputs():
            try:
                sha = compute_sha256(output)
                self.log_event(
                    "output_sha256",
                    output_path=str(output),
                    sha256=sha,
                    size_bytes=output.stat().st_size,
                )
            except Exception as exc:  # pragma: no cover - logging must not mask script result.
                self.log_event("output_sha256", status="error", output_path=str(output), error=repr(exc))
        elapsed = time.monotonic() - self.started_monotonic
        self.log_event("end", status=status, elapsed_s=round(elapsed, 6))
        self._original_print(
            f"[{utc_now_iso()}] [{self.script_name}] [RUN_AUDIT] END status={status} elapsed_s={elapsed:.3f} log={self.log_path}",
            flush=True,
        )


def install_run_audit(script_file: str | Path) -> RunAudit:
    """Install execution auditing for a script and return the audit object."""
    return RunAudit(script_file).install()
