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
python scripts/bayesian_heom_hierarchy_v2.py src/heom_bayes_input_current.csv --output-dir heom_bayes_out_v2 --draws 20000

echo "[KwanTube] Reproducing manuscript-level results..."
python scripts/reproduce_paper_results.py --mode paper

echo "[KwanTube] Generating manuscript figures..."
python scripts/generate_paper_figures.py

echo "[KwanTube] Done."
echo "[KwanTube] Outputs: figures_final/ and validation_report.json."

read -n 1 -s -r -p "Press any key to close..."
echo
