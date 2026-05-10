# KwanTube Makefile — convenience targets for all platforms
# Requires GNU Make (Linux/macOS) or WSL/Git Bash on Windows.
#
# Usage examples:
#   make reproduce      # Full Phase 6-9 pipeline
#   make validate       # Run validation checks only
#   make figures        # Generate manuscript figures only
#   make seal           # Regenerate all SHA-256 sidecars
#   make test           # Run full test suite
#   make docker-run     # Run full pipeline inside Docker
#   make clean          # Remove generated figures (NOT data artifacts)

PYTHON ?= python
PYTEST ?= pytest

.PHONY: all reproduce validate figures seal test docker-build docker-run clean help

all: reproduce

# ─────────────────────────────────────────────────────────────
# FULL PIPELINE (Phases 6–9)
# ─────────────────────────────────────────────────────────────
reproduce:
	@echo "[KwanTube] Phase 6 — Bayesian HEOM hierarchy..."
	$(PYTHON) src/scripts/analysis/bayesian_heom_hierarchy_v2.py
	$(PYTHON) src/qmc_mt/sensitivity.py
	$(PYTHON) src/qmc_mt/sensitivity_priors.py
	$(PYTHON) src/qmc_mt/sbc_report.py
	@echo "[KwanTube] Phase 8 — Canonical validation..."
	$(PYTHON) src/scripts/validation/heom_pade_convergence.py
	$(PYTHON) src/scripts/validation/audit_lineage.py
	$(PYTHON) src/scripts/validation/reproduce_paper_results.py --mode paper
	@echo "[KwanTube] Phase 9 — Figures..."
	$(PYTHON) src/scripts/analysis/assemble_master_results.py
	$(PYTHON) src/scripts/analysis/extract_heom_production_figure.py
	$(PYTHON) src/scripts/figures/generate_paper_figures.py
	$(PYTHON) src/scripts/figures/extract_vector_figure.py
	@echo "[KwanTube] Phase 10 — Integrity seal..."
	$(PYTHON) src/scripts/validation/seal_outputs.py
	$(PYTHON) src/scripts/validation/validate_outputs.py
	@echo "[KwanTube] Done. All outputs in outputs_data/."

# ─────────────────────────────────────────────────────────────
# INDIVIDUAL TARGETS
# ─────────────────────────────────────────────────────────────
validate:
	$(PYTHON) src/scripts/validation/reproduce_paper_results.py --mode paper

figures:
	$(PYTHON) src/scripts/figures/generate_paper_figures.py
	$(PYTHON) src/scripts/analysis/extract_heom_production_figure.py

seal:
	$(PYTHON) src/scripts/validation/seal_outputs.py
	$(PYTHON) src/scripts/validation/validate_outputs.py

# ─────────────────────────────────────────────────────────────
# DATA PIPELINE (Phases 1–7)
# ─────────────────────────────────────────────────────────────
data:
	$(PYTHON) src/scripts/data/fetch_public_data.py
	$(PYTHON) src/scripts/data/curate_compact.py
	$(PYTHON) src/scripts/data/build_registry.py
	$(PYTHON) src/qmc_mt/pdb_tubulin_analysis.py
	$(PYTHON) src/scripts/data/build_hamiltonian.py
	$(PYTHON) src/scripts/data/compute_detectability_metrics.py
	$(PYTHON) src/scripts/data/run_comparative_panels.py
	$(PYTHON) src/scripts/data/export_claim_traceability.py

# ─────────────────────────────────────────────────────────────
# TESTS
# ─────────────────────────────────────────────────────────────
test:
	$(PYTEST) tests/ -v
	$(PYTHON) src/qmc_mt/test_ns_consistency.py

# ─────────────────────────────────────────────────────────────
# DOCKER
# ─────────────────────────────────────────────────────────────
docker-build:
	docker build -t kwantube:3.5.1.1 .

docker-run: docker-build
	docker compose up

# ─────────────────────────────────────────────────────────────
# INSTALL
# ─────────────────────────────────────────────────────────────
install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

# ─────────────────────────────────────────────────────────────
# CLEAN (safe — never deletes data artifacts)
# ─────────────────────────────────────────────────────────────
clean:
	find outputs_data/figures_final -name "*.png" -delete 2>/dev/null || true
	find outputs_data/figures_final -name "*.pdf" -delete 2>/dev/null || true
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true

# ─────────────────────────────────────────────────────────────
# HELP
# ─────────────────────────────────────────────────────────────
help:
	@echo "KwanTube v3.5.1.1 — Available targets:"
	@echo ""
	@echo "  make reproduce    Full pipeline (Phases 6-9) — recommended"
	@echo "  make data         Data acquisition pipeline (Phases 1-7)"
	@echo "  make validate     Run canonical validation only"
	@echo "  make figures      Generate manuscript figures only"
	@echo "  make seal         Regenerate all SHA-256 integrity sidecars"
	@echo "  make test         Run full test suite (pytest + NS consistency)"
	@echo "  make docker-run   Run full pipeline inside Docker"
	@echo "  make install      Install package in editable mode"
	@echo "  make clean        Remove generated figures (data artifacts preserved)"
	@echo ""
	@echo "  For the full pipeline map: see PIPELINE_MAP.md"
