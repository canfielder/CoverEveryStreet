"""Project path definitions."""

from pathlib import Path

ROOT_DIR        = Path(__file__).parent.parent  # src/ -> project root
DATA_DIR        = ROOT_DIR / "data"
NETWORK_CACHE_DIR = DATA_DIR / "networks"
INBOX_DIR       = DATA_DIR / "gpx_inbox"
PROCESSED_DIR   = DATA_DIR / "gpx_processed"
PUBLISHED_DIR   = DATA_DIR / "published"
DB_PATH         = DATA_DIR / "walks.db"
GPX_SAMPLE_DIR  = DATA_DIR / "gps_tracking_files" / "test_case"

# Ensure required dirs exist at import time
for _d in (NETWORK_CACHE_DIR, INBOX_DIR, PROCESSED_DIR, PUBLISHED_DIR):
    _d.mkdir(parents=True, exist_ok=True)
