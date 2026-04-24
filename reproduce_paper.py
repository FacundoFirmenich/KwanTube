#!/usr/bin/env python3
"""
Compatibility wrapper for the canonical reproduction entrypoint.

Canonical command:
    python reproduce_paper_results.py [--fast] [--full-roc]

This wrapper is intentionally thin to preserve backwards compatibility for
existing automation that still calls ``reproduce_paper.py``.
"""

from reproduce_paper_results import main


if __name__ == "__main__":
    raise SystemExit(main())
