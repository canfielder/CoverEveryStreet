"""
GPX inbox scanning and ingestion pipeline.

Scans data/gpx_inbox/ for .gpx files, deduplicates against the database
by content_hash and start_time, and returns pending files ready for user
review before being written to the database.

Dedup logic (two-layer):
  1. content_hash  — exact duplicate bytes (same file uploaded twice)
  2. start_time    — same activity re-exported with a different name
     (e.g. activity renamed in Garmin then re-exported)

Files that pass dedup become "pending" — they are parsed but NOT yet
written to the database until the user confirms companion / type metadata
in the UI and presses Process.

After processing, files are moved from gpx_inbox/ → gpx_processed/.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from src import database as db
from src.gpx_parser import parse_gpx_file
from src.paths import INBOX_DIR, PROCESSED_DIR

logger = logging.getLogger(__name__)


# ── Public API ────────────────────────────────────────────────────────────────

def scan_inbox() -> list[dict]:
    """
    Scan INBOX_DIR for .gpx files and return parsed dicts for new files.

    Each returned dict is a parsed track dict (from gpx_parser) plus an
    extra key:
        ``inbox_path`` (Path) — absolute path to the source file

    Files are skipped (with a debug log) if:
      - content_hash already exists in the database
      - start_time already exists in the database

    Returns:
        List of pending track dicts, sorted chronologically by start_time.
    """
    gpx_files = sorted(INBOX_DIR.glob("*.gpx"))
    if not gpx_files:
        logger.info("Inbox is empty: %s", INBOX_DIR)
        return []

    known_hashes = db.get_all_content_hashes()
    known_start_times = db.get_all_start_times()

    pending: list[dict] = []
    for path in gpx_files:
        raw = path.read_bytes()
        parsed = parse_gpx_file(raw, filename=path.name)
        if parsed is None:
            logger.warning("Could not parse inbox file: %s", path.name)
            continue

        if parsed["content_hash"] in known_hashes:
            logger.debug("Skipping duplicate (hash match): %s", path.name)
            continue

        start_iso = parsed["start_time"].isoformat() if parsed["start_time"] else None
        if start_iso and start_iso in known_start_times:
            logger.debug("Skipping duplicate (start_time match): %s", path.name)
            continue

        parsed["inbox_path"] = path
        pending.append(parsed)

    pending.sort(key=lambda t: t["start_time"] or "")
    logger.info("Inbox scan: %d new / %d total files", len(pending), len(gpx_files))
    return pending


def move_to_processed(inbox_path: Path) -> Path:
    """
    Move a file from gpx_inbox/ to gpx_processed/.

    If a file with the same name already exists in gpx_processed/, the
    new file is renamed with a numeric suffix to avoid collision.

    Returns the destination Path.
    """
    dest = PROCESSED_DIR / inbox_path.name
    if dest.exists():
        stem, suffix = inbox_path.stem, inbox_path.suffix
        i = 1
        while dest.exists():
            dest = PROCESSED_DIR / f"{stem}_{i}{suffix}"
            i += 1

    shutil.move(str(inbox_path), str(dest))
    logger.info("Moved %s → %s", inbox_path.name, dest.name)
    return dest


def ingest_activity(
    track: dict,
    companions: list[str],
    activity_type_override: str | None,
    network_gdf,
    snap_distance_m: float,
) -> int:
    """
    Persist a single pending activity to the database.

    Steps:
      1. Insert activity row (filename, metadata).
      2. Match track to blocks, compute coverage ratios.
      3. Bulk-insert block_walk rows (skipping any block_exclusions).
      4. Move file from inbox → processed.

    Args:
        track:                  Parsed track dict (must include ``inbox_path``).
        companions:             List of companion names selected by the user.
        activity_type_override: If the user changed the activity type in the UI,
                                use this value; otherwise use track["activity_type"].
        network_gdf:            Current network GeoDataFrame (with block_id column).
        snap_distance_m:        Snap distance used for matching.

    Returns:
        The new activity_id (int).
    """
    from src.track_matcher import match_track_to_blocks

    activity_type = activity_type_override or track["activity_type"]

    activity_id = db.insert_activity(
        filename=track["filename"],
        content_hash=track["content_hash"],
        start_time=track["start_time"],
        end_time=track["end_time"],
        activity_type=activity_type,
        companions=companions,
        raw_gpx_type=track["raw_gpx_type"],
        track_name=track["track_name"],
    )

    block_pairs = match_track_to_blocks(network_gdf, track, activity_id, snap_distance_m)
    db.insert_block_walks(block_pairs)

    move_to_processed(track["inbox_path"])
    logger.info(
        "Ingested activity %d (%s, %d companions, %d blocks matched)",
        activity_id, activity_type, len(companions), len(block_pairs),
    )
    return activity_id
