# KwanTube — Dockerfile
# One-shot reproducibility image for KwanTube v3.5.1.1
#
# Usage:
#   docker build -t kwantube .
#   docker run --rm -v "$(pwd)/outputs_data:/app/outputs_data" kwantube
#
# This image reproduces Phases 6–9 of the KwanTube pipeline
# (Bayesian HEOM hierarchy → validation → figures → integrity seal).
# Pre-computed HEOM production trajectories are included in the repo;
# the expensive Phase 4 run is NOT re-executed.

FROM python:3.11-slim

LABEL maintainer="f.firmenich@cedesur.org"
LABEL version="3.5.1.1"
LABEL description="KwanTube: reproducible quantum coherence pipeline for tubulin"

# System dependencies for matplotlib/scipy
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    fonts-liberation \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy only what is needed for installation first (cache-friendly layer)
COPY requirements.txt pyproject.toml ./
COPY src/ ./src/

# Install the package and its dependencies
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt \
    && pip install --no-cache-dir -e .

# Copy the rest of the repository
COPY . .

# Default command: full Phases 6–9 reproduction
# Override with: docker run kwantube python src/scripts/validation/reproduce_paper_results.py
CMD ["bash", "-c", "\
    echo '[KwanTube Docker] Phase 6 — Bayesian HEOM hierarchy...' && \
    python src/scripts/analysis/bayesian_heom_hierarchy_v2.py && \
    echo '[KwanTube Docker] Phase 6 — Sensitivity analysis...' && \
    python src/qmc_mt/sensitivity.py && \
    python src/qmc_mt/sbc_report.py && \
    echo '[KwanTube Docker] Phase 8 — Canonical validation...' && \
    python src/scripts/validation/reproduce_paper_results.py --mode paper && \
    echo '[KwanTube Docker] Phase 9 — Figures...' && \
    python src/scripts/figures/generate_paper_figures.py && \
    python src/scripts/analysis/extract_heom_production_figure.py && \
    echo '[KwanTube Docker] Phase 10 — Integrity seal...' && \
    python src/scripts/validation/seal_outputs.py && \
    python src/scripts/validation/validate_outputs.py && \
    echo '[KwanTube Docker] Done. Outputs in outputs_data/.' \
"]
