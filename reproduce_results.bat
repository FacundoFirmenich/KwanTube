@echo off
setlocal enabledelayedexpansion

echo [KwanTube] Starting one-shot reproducibility workflow...

cd /d "%~dp0"

where py >nul 2>nul
if %errorlevel%==0 (
    set PY=py -3
) else (
    where python >nul 2>nul
    if %errorlevel%==0 (
        set PY=python
    ) else (
        echo [KwanTube] ERROR: Python was not found in PATH.
        pause
        exit /b 1
    )
)

echo [KwanTube] Using Python:
%PY% --version

if not exist ".venv" (
    echo [KwanTube] Creating virtual environment...
    %PY% -m venv .venv
)

call ".venv\Scripts\activate.bat"

echo [KwanTube] Installing dependencies...
python -m pip install --upgrade pip
pip install -r requirements.txt

echo [KwanTube] Running Bayesian HEOM hierarchy v2 analysis...
python src\scripts\analysis\bayesian_heom_hierarchy_v2.py

echo [KwanTube] Reproducing manuscript-level results...
python src\scripts\validation\reproduce_paper_results.py --mode paper

echo [KwanTube] Generating manuscript figures...
python src\scripts\figures\generate_paper_figures.py
python src\scripts\analysis\extract_heom_production_figure.py

echo [KwanTube] Sealing integrity ledger...
python src\scripts\validation\seal_outputs.py
python src\scripts\validation\validate_outputs.py

echo [KwanTube] Done.
echo [KwanTube] Outputs: outputs_data\figures_final\ ^| validation: outputs_data\raw_json\structural\
pause
