"""Path resolution for KwanTube v3.5.1."""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUTS_DATA_DIR = PROJECT_ROOT / "outputs_data"
