#!/usr/bin/env bash
set -euo pipefail

echo "[KwanTube] Starting one-shot reproducibility workflow..."

cd "$(dirname "$0")"

if command -v python3 >/dev/null 2>&1; then
    PYTHON="python3"
elif command -v python >/dev/null 2>&1; then
    PYTHON="python"
else
    echo "[KwanTube] ERROR: Python is not installed or not in PATH."
    read -n 1 -s -r -p "Press any key to close..."
    exit 1
fi

echo "[KwanTube] Using: $($PYTHON --version)"

if [ ! -d ".venv" ]; then
    echo "[KwanTube] Creating virtual environment..."
    "$PYTHON" -m venv .venv
fi

source .venv/bin/activate

echo "[KwanTube] Installing dependencies..."
python -m pip install --upgrade pip
pip install -r requirements.txt

echo "[KwanTube] Running Bayesian HEOM hierarchy v2 analysis..."
python src/scripts/analysis/bayesian_heom_hierarchy_v2.py

echo "[KwanTube] Reproducing manuscript-level results..."
python src/scripts/validation/reproduce_paper_results.py --mode paper

echo "[KwanTube] Generating manuscript figures..."
python src/scripts/figures/generate_paper_figures.py
python src/scripts/analysis/extract_heom_production_figure.py

echo "[KwanTube] Sealing integrity ledger..."
python src/scripts/validation/seal_outputs.py
python src/scripts/validation/validate_outputs.py

echo "[KwanTube] Done."
echo "[KwanTube] Outputs: outputs_data/figures_final/ | validation: outputs_data/raw_json/structural/"

read -n 1 -s -r -p "Press any key to close..."
echo
