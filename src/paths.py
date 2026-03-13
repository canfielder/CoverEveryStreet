"""Project path definitions."""

from pathlib import Path

ROOT_DIR = Path(__file__).parent.parent  # src/ -> project root
DATA_DIR = ROOT_DIR / "data"
NETWORK_CACHE_DIR = DATA_DIR / "networks"
SESSION_DIR = DATA_DIR / "sessions"
GPX_SAMPLE_DIR = DATA_DIR / "gps_tracking_files" / "test_case"

# Ensure data dirs exist at import time
NETWORK_CACHE_DIR.mkdir(parents=True, exist_ok=True)
SESSION_DIR.mkdir(parents=True, exist_ok=True)
