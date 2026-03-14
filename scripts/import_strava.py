"""
Import activities from a Strava bulk data export into the GPX inbox.

Reads directly from the export zip (no full extraction to disk).
Handles three file formats found in Strava exports:
  .gpx       — copied as-is
  .gpx.gz    — decompressed then copied
  .fit.gz    — decompressed then converted via fit2gpx

Activity metadata (name, type, date) is read from activities.csv inside
the zip and embedded into the output GPX so the app picks it up correctly.

Only Walk, Run, and Hike activity types are imported by default.
Other types (Weight Training, Yoga, etc.) are skipped — they have no
GPS track relevant to street coverage.

Prerequisites for .fit.gz files:
    uv add fit2gpx

Usage:
    # Import all walks/runs/hikes
    uv run python -m scripts.import_strava --zip data/garmin_exports/2026-03-13__strava_bulk_export.zip

    # Also include Ride activities
    uv run python -m scripts.import_strava --zip ... --types Walk Run Hike Ride

    # Dry run — show what would be imported without writing files
    uv run python -m scripts.import_strava --zip ... --dry-run
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import logging
import re
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.paths import INBOX_DIR, PROCESSED_DIR

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Strava activity type → normalized internal type
STRAVA_TYPE_MAP: dict[str, str] = {
    "walk":         "walk",
    "hike":         "walk",
    "run":          "run",
    "ride":         "cycle",
    "e-bike ride":  "cycle",
}

DEFAULT_IMPORT_TYPES = {"Walk", "Run", "Hike"}


def _safe_filename(name: str, date_str: str, activity_id: str) -> str:
    """Build a readable, filesystem-safe GPX filename."""
    safe = re.sub(r"[^\w\s-]", "", name).strip()
    safe = re.sub(r"\s+", "_", safe)[:60]
    return f"{date_str[:10]}_{safe}_{activity_id}.gpx"


def _already_imported(stem: str) -> bool:
    return (INBOX_DIR / f"{stem}.gpx").exists() or (PROCESSED_DIR / f"{stem}.gpx").exists()


def _inject_gpx_metadata(gpx_bytes: bytes, activity_name: str, activity_type: str) -> bytes:
    """
    Inject <name> and <type> tags into the first <trk> element of a GPX file.
    This ensures the app picks up the Strava activity name and type correctly
    rather than relying solely on pace-based inference.
    """
    text = gpx_bytes.decode("utf-8", errors="replace")

    norm_type = STRAVA_TYPE_MAP.get(activity_type.lower(), "walk")

    def _replace_trk(m):
        inner = m.group(1)
        # Remove any existing name/type tags
        inner = re.sub(r"<name>.*?</name>", "", inner, flags=re.DOTALL)
        inner = re.sub(r"<type>.*?</type>", "", inner, flags=re.DOTALL)
        return f"<trk><name>{activity_name}</name><type>{norm_type}</type>{inner}</trk>"

    text = re.sub(r"<trk>(.*?)</trk>", _replace_trk, text, flags=re.DOTALL)
    return text.encode("utf-8")


def _convert_fit_to_gpx(fit_bytes: bytes, activity_name: str, activity_type: str) -> bytes | None:
    """Convert raw FIT bytes to GPX bytes using fit2gpx. Returns None on failure."""
    try:
        from fit2gpx import Converter
    except ImportError:
        logger.warning("fit2gpx not installed — skipping .fit.gz files. Run: uv add fit2gpx")
        return None

    import tempfile, os
    tmp_in = tmp_out = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".fit", delete=False) as fin:
            fin.write(fit_bytes)
            tmp_in = fin.name
        with tempfile.NamedTemporaryFile(suffix=".gpx", delete=False) as fout:
            tmp_out = fout.name

        conv = Converter()
        conv.fit_to_gpx(f_in=tmp_in, f_out=tmp_out)
        gpx_bytes = Path(tmp_out).read_bytes()
        return _inject_gpx_metadata(gpx_bytes, activity_name, activity_type)
    except Exception as exc:
        logger.debug("fit2gpx conversion failed: %s", exc)
        return None
    finally:
        for p in (tmp_in, tmp_out):
            if p and os.path.exists(p):
                os.unlink(p)


def run(
    zip_path: Path,
    import_types: set[str] = DEFAULT_IMPORT_TYPES,
    dry_run: bool = False,
    limit: int | None = None,
) -> None:
    INBOX_DIR.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()

        # ── Load activities.csv ───────────────────────────────────────────────
        if "activities.csv" not in names:
            logger.error("activities.csv not found in zip — is this a Strava export?")
            sys.exit(1)

        csv_bytes = zf.read("activities.csv").decode("utf-8-sig")  # handle BOM
        reader = csv.DictReader(io.StringIO(csv_bytes))
        metadata: dict[str, dict] = {}
        for row in reader:
            fname = row.get("Filename", "").strip()
            if fname:
                metadata[fname] = {
                    "id":    row.get("Activity ID", ""),
                    "name":  row.get("Activity Name", "Unnamed"),
                    "type":  row.get("Activity Type", ""),
                    "date":  row.get("Activity Date", ""),
                }

        # ── Count what we'll process ──────────────────────────────────────────
        relevant = [
            (fname, meta)
            for fname, meta in metadata.items()
            if meta["type"] in import_types and fname in names
        ]

        total_matching = len(relevant)
        if limit and limit < total_matching:
            import random
            random.shuffle(relevant)
            relevant = relevant[:limit]
            logger.info(
                "Sampling %d of %d matching activities (types: %s)",
                limit, total_matching, sorted(import_types),
            )
        else:
            logger.info(
                "Found %d activities matching types %s (out of %d total)",
                total_matching, sorted(import_types), len(metadata),
            )

        # ── Process each activity ─────────────────────────────────────────────
        imported = skipped_dup = skipped_fmt = skipped_err = 0

        for i, (fname, meta) in enumerate(relevant, 1):
            activity_id = meta["id"]
            activity_name = meta["name"]
            activity_type = meta["type"]
            date_str = meta["date"][:10].replace(", ", "-").replace(" ", "-")

            out_stem = _safe_filename(activity_name, date_str, activity_id).removesuffix(".gpx")
            out_name = out_stem + ".gpx"

            if _already_imported(out_stem):
                skipped_dup += 1
                continue

            if i % 100 == 0 or i == len(relevant):
                logger.info("Processing %d/%d…", i, len(relevant))

            if dry_run:
                logger.info("  [dry-run] Would import: %s", out_name)
                imported += 1
                continue

            # Read and decompress the source file
            raw = zf.read(fname)
            if fname.endswith(".gz"):
                try:
                    raw = gzip.decompress(raw)
                except Exception as exc:
                    logger.warning("Failed to decompress %s: %s", fname, exc)
                    skipped_err += 1
                    continue

            # Convert to GPX
            if fname.endswith(".fit.gz") or fname.endswith(".fit"):
                gpx_bytes = _convert_fit_to_gpx(raw, activity_name, activity_type)
                if gpx_bytes is None:
                    skipped_fmt += 1
                    continue
            else:
                # Already GPX — inject metadata
                gpx_bytes = _inject_gpx_metadata(raw, activity_name, activity_type)

            dest = INBOX_DIR / out_name
            dest.write_bytes(gpx_bytes)
            imported += 1

    logger.info(
        "Done — %d imported, %d skipped (already exist), "
        "%d skipped (unsupported format), %d errors",
        imported, skipped_dup, skipped_fmt, skipped_err,
    )
    if skipped_fmt > 0:
        logger.info(
            "  %d .fit.gz files skipped. Install fit2gpx to convert them: uv add fit2gpx",
            skipped_fmt,
        )
    if imported > 0 and not dry_run:
        logger.info("Inbox now has files ready to process. Launch the app to review them.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zip", required=True, type=Path, help="Path to Strava export zip")
    parser.add_argument(
        "--types", nargs="+", default=sorted(DEFAULT_IMPORT_TYPES),
        help="Activity types to import (default: Walk Run Hike)",
    )
    parser.add_argument("--dry-run", action="store_true", help="Show what would be imported without writing files")
    parser.add_argument("--limit", type=int, default=None, help="Import a random sample of N activities")
    args = parser.parse_args()
    run(zip_path=args.zip, import_types=set(args.types), dry_run=args.dry_run, limit=args.limit)
