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

echo [KwanTube] Reproducing manuscript-level results...
python scripts\reproduce_paper_results.py --mode paper

echo [KwanTube] Done.
echo [KwanTube] Outputs: figures_final/ and validation_report.json.
pause
