"""
Batch import Garmin FIT files into the GPX inbox.

Converts .fit files to .gpx using fit2gpx, then drops them into
data/gpx_inbox/ so the Streamlit app can process them.

Prerequisites:
    pip install fit2gpx       (or: uv add fit2gpx)

Usage:
    # Import all .fit files from a Garmin export folder
    uv run python -m scripts.import_batch --src ~/Downloads/garmin_export/

    # Import a single file
    uv run python -m scripts.import_batch --src ~/Downloads/2024-01-15.fit

The script is idempotent — already-imported files (same name in inbox or
processed) are skipped automatically.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from pathlib import Path

# Ensure project root is on sys.path when run as a script
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.paths import INBOX_DIR, PROCESSED_DIR

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def _already_imported(stem: str) -> bool:
    """Return True if a .gpx with this stem already exists in inbox or processed."""
    return (INBOX_DIR / f"{stem}.gpx").exists() or (PROCESSED_DIR / f"{stem}.gpx").exists()


def convert_fit(fit_path: Path, dest_dir: Path) -> Path | None:
    """Convert a single .fit file to .gpx using fit2gpx. Returns output path or None."""
    try:
        from fit2gpx import Converter
    except ImportError:
        logger.error("fit2gpx not installed. Run: uv add fit2gpx")
        sys.exit(1)

    stem = fit_path.stem
    out_path = dest_dir / f"{stem}.gpx"

    conv = Converter()
    _, gpx = conv.fit_to_gpx(f_in=str(fit_path))
    with open(out_path, "w") as fh:
        fh.write(gpx.to_xml())
    return out_path


def import_gpx(gpx_path: Path) -> Path | None:
    """Copy a .gpx file directly into the inbox."""
    dest = INBOX_DIR / gpx_path.name
    if dest.exists():
        logger.info("Skipping (already in inbox): %s", gpx_path.name)
        return None
    shutil.copy2(gpx_path, dest)
    return dest


def run(src: Path) -> None:
    INBOX_DIR.mkdir(parents=True, exist_ok=True)

    files: list[Path] = []
    if src.is_dir():
        files = sorted(src.glob("**/*.fit")) + sorted(src.glob("**/*.gpx"))
    elif src.is_file():
        files = [src]
    else:
        logger.error("Source not found: %s", src)
        sys.exit(1)

    if not files:
        logger.info("No .fit or .gpx files found in: %s", src)
        return

    converted, skipped, failed = 0, 0, 0
    for f in files:
        stem = f.stem
        if _already_imported(stem):
            logger.debug("Already imported, skipping: %s", f.name)
            skipped += 1
            continue

        if f.suffix.lower() == ".fit":
            logger.info("Converting: %s", f.name)
            try:
                out = convert_fit(f, INBOX_DIR)
                if out:
                    logger.info("  → %s", out.name)
                    converted += 1
            except Exception as exc:
                logger.warning("Failed to convert %s: %s", f.name, exc)
                failed += 1
        elif f.suffix.lower() == ".gpx":
            result = import_gpx(f)
            if result:
                logger.info("Copied: %s", result.name)
                converted += 1
            else:
                skipped += 1

    logger.info(
        "Done — %d imported, %d skipped (already exist), %d failed",
        converted, skipped, failed,
    )
    if converted:
        logger.info("Inbox now has %d pending file(s). Launch the app to process them.", converted)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", required=True, type=Path, help="Source file or directory")
    args = parser.parse_args()
    run(args.src)
