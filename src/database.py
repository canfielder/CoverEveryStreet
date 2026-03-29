"""
SQLite persistence layer.

All reads and writes to walks.db go through this module.
Schema: activities, block_walks, block_exclusions.
"""

from __future__ import annotations

import json
import logging
import sqlite3
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

from src.paths import DB_PATH

logger = logging.getLogger(__name__)

# ── Schema ────────────────────────────────────────────────────────────────────

_SCHEMA_VERSION = 1  # increment when adding to _MIGRATIONS

_MIGRATIONS: list[tuple[int, str]] = [
    # (target_version, SQL statement)
    # Example for future use:
    # (2, "ALTER TABLE activities ADD COLUMN notes TEXT"),
]

_SCHEMA = """
CREATE TABLE IF NOT EXISTS activities (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    filename      TEXT    NOT NULL,
    content_hash  TEXT    NOT NULL,
    start_time    TEXT    NOT NULL UNIQUE,   -- ISO-8601, true dedup key
    end_time      TEXT,
    activity_type TEXT    NOT NULL,          -- normalized: walk / run / cycle …
    companions    TEXT    NOT NULL DEFAULT '[]',  -- JSON array of strings
    raw_gpx_type  TEXT,                      -- original <type> tag from GPX
    track_name    TEXT,
    processed_at  TEXT    NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS block_walks (
    block_id        TEXT    NOT NULL,
    activity_id     INTEGER NOT NULL REFERENCES activities(id) ON DELETE CASCADE,
    coverage_ratio  REAL    NOT NULL,        -- 0.0–1.0 auto; 1.0 manual
    source          TEXT    NOT NULL,        -- "auto" | "manual"
    PRIMARY KEY (block_id, activity_id)
);

CREATE TABLE IF NOT EXISTS block_exclusions (
    block_id    TEXT    NOT NULL,
    activity_id INTEGER NOT NULL REFERENCES activities(id) ON DELETE CASCADE,
    PRIMARY KEY (block_id, activity_id)
);

CREATE INDEX IF NOT EXISTS idx_block_walks_block   ON block_walks (block_id);
CREATE INDEX IF NOT EXISTS idx_block_walks_activity ON block_walks (activity_id);
CREATE INDEX IF NOT EXISTS idx_activities_start    ON activities (start_time);
"""


# ── Connection management ─────────────────────────────────────────────────────

def _connect(path: Path = DB_PATH) -> sqlite3.Connection:
    conn = sqlite3.connect(path, detect_types=sqlite3.PARSE_DECLTYPES)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA journal_mode = WAL")
    return conn


@contextmanager
def _db(path: Path = DB_PATH):
    conn = _connect(path)
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def init_db(path: Path = DB_PATH) -> None:
    """Create tables if they don't exist and apply pending migrations. Safe to call on every app start."""
    with _db(path) as conn:
        conn.executescript(_SCHEMA)
        _apply_migrations(conn)
    logger.info("Database initialised: %s", path)


def _apply_migrations(conn: sqlite3.Connection) -> None:
    """Apply any pending migrations inside an existing _db() context."""
    conn.execute(
        "CREATE TABLE IF NOT EXISTS schema_version (version INTEGER NOT NULL)"
    )
    row = conn.execute("SELECT version FROM schema_version").fetchone()
    current = row["version"] if row else 0
    if current == 0:
        conn.execute("INSERT INTO schema_version (version) VALUES (0)")

    for version, sql in _MIGRATIONS:
        if version > current:
            conn.execute(sql)
            conn.execute("UPDATE schema_version SET version = ?", (version,))
            logger.info("Applied migration to schema version %d", version)


# ── Activities ────────────────────────────────────────────────────────────────

def get_all_content_hashes(path: Path = DB_PATH) -> set[str]:
    with _db(path) as conn:
        rows = conn.execute("SELECT content_hash FROM activities").fetchall()
    return {r["content_hash"] for r in rows}


def get_all_start_times(path: Path = DB_PATH) -> set[str]:
    with _db(path) as conn:
        rows = conn.execute("SELECT start_time FROM activities").fetchall()
    return {r["start_time"] for r in rows}


def get_activity_by_start_time(start_time: str, path: Path = DB_PATH) -> sqlite3.Row | None:
    with _db(path) as conn:
        return conn.execute(
            "SELECT * FROM activities WHERE start_time = ?", (start_time,)
        ).fetchone()


def insert_activity(
    filename: str,
    content_hash: str,
    start_time: datetime,
    end_time: datetime | None,
    activity_type: str,
    companions: list[str],
    raw_gpx_type: str | None,
    track_name: str | None,
    path: Path = DB_PATH,
) -> int:
    """Insert a new activity row and return its id."""
    with _db(path) as conn:
        cur = conn.execute(
            """INSERT INTO activities
               (filename, content_hash, start_time, end_time,
                activity_type, companions, raw_gpx_type, track_name)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                filename,
                content_hash,
                start_time.isoformat(),
                end_time.isoformat() if end_time else None,
                activity_type,
                json.dumps(companions),
                raw_gpx_type,
                track_name,
            ),
        )
        return cur.lastrowid


def update_activity(
    activity_id: int,
    activity_type: str,
    companions: list[str],
    path: Path = DB_PATH,
) -> None:
    """Update editable fields on an existing activity."""
    with _db(path) as conn:
        conn.execute(
            "UPDATE activities SET activity_type=?, companions=? WHERE id=?",
            (activity_type, json.dumps(companions), activity_id),
        )


def get_all_activities(path: Path = DB_PATH) -> list[dict]:
    """Return all activities sorted by start_time descending."""
    with _db(path) as conn:
        rows = conn.execute(
            "SELECT * FROM activities ORDER BY start_time DESC"
        ).fetchall()
    return [_activity_row_to_dict(r) for r in rows]


def _activity_row_to_dict(row: sqlite3.Row) -> dict:
    d = dict(row)
    d["companions"] = json.loads(d["companions"] or "[]")
    return d


# ── Block walks ───────────────────────────────────────────────────────────────

def insert_block_walks(
    block_activity_pairs: list[tuple[str, int, float, str]],
    path: Path = DB_PATH,
) -> None:
    """
    Bulk-insert (block_id, activity_id, coverage_ratio, source) rows.
    Skips pairs that exist in block_exclusions.
    Ignores duplicate (block_id, activity_id) pairs.
    """
    if not block_activity_pairs:
        return
    with _db(path) as conn:
        exclusions = {
            (r["block_id"], r["activity_id"])
            for r in conn.execute("SELECT block_id, activity_id FROM block_exclusions").fetchall()
        }
        to_insert = [
            (bid, aid, ratio, src)
            for bid, aid, ratio, src in block_activity_pairs
            if (bid, aid) not in exclusions
        ]
        conn.executemany(
            """INSERT OR IGNORE INTO block_walks
               (block_id, activity_id, coverage_ratio, source)
               VALUES (?, ?, ?, ?)""",
            to_insert,
        )
        logger.info("Inserted %d block_walk rows (%d skipped by exclusions)",
                    len(to_insert), len(block_activity_pairs) - len(to_insert))


def insert_manual_block_walk(
    block_id: str, activity_id: int, path: Path = DB_PATH
) -> None:
    """Manually assign a block to an activity (coverage_ratio=1.0, source=manual)."""
    with _db(path) as conn:
        conn.execute(
            """INSERT OR REPLACE INTO block_walks
               (block_id, activity_id, coverage_ratio, source)
               VALUES (?, ?, 1.0, 'manual')""",
            (block_id, activity_id),
        )
        # Clear any exclusion that might exist for this pair
        conn.execute(
            "DELETE FROM block_exclusions WHERE block_id=? AND activity_id=?",
            (block_id, activity_id),
        )


def remove_block_walk(
    block_id: str, activity_id: int, source: str, path: Path = DB_PATH
) -> None:
    """
    Remove a block assignment.
    - Auto removals → also write to block_exclusions (prevent re-match on reprocess)
    - Manual removals → just delete from block_walks
    """
    with _db(path) as conn:
        conn.execute(
            "DELETE FROM block_walks WHERE block_id=? AND activity_id=?",
            (block_id, activity_id),
        )
        if source == "auto":
            conn.execute(
                "INSERT OR IGNORE INTO block_exclusions (block_id, activity_id) VALUES (?, ?)",
                (block_id, activity_id),
            )


def get_walked_blocks(
    start_date: datetime | None = None,
    end_date: datetime | None = None,
    coverage_threshold: float = 0.80,
    path: Path = DB_PATH,
) -> dict[str, dict]:
    """
    Return walked blocks in the date range as {block_id: {activity_type, companions, walk_date}}.
    For each block, only the earliest qualifying activity counts.
    """
    params: list = [coverage_threshold]
    date_filter = ""
    if start_date:
        date_filter += " AND a.start_time >= ?"
        params.append(start_date.isoformat())
    if end_date:
        date_filter += " AND a.start_time <= ?"
        params.append(end_date.isoformat())

    query = f"""
        SELECT
            bw.block_id,
            a.activity_type,
            a.companions,
            a.start_time AS walk_date,
            bw.source
        FROM block_walks bw
        JOIN activities a ON a.id = bw.activity_id
        WHERE bw.coverage_ratio >= ?
        {date_filter}
        ORDER BY bw.block_id, a.start_time ASC
    """
    with _db(path) as conn:
        rows = conn.execute(query, params).fetchall()

    # Keep only the earliest activity per block
    result: dict[str, dict] = {}
    for row in rows:
        bid = row["block_id"]
        if bid not in result:
            result[bid] = {
                "activity_type": row["activity_type"],
                "companions":    json.loads(row["companions"] or "[]"),
                "walk_date":     row["walk_date"],
                "source":        row["source"],
            }
    return result


def get_block_assignments(block_id: str, path: Path = DB_PATH) -> list[dict]:
    """Return all activity assignments for a block (for the assignment panel)."""
    with _db(path) as conn:
        rows = conn.execute(
            """SELECT bw.activity_id, bw.coverage_ratio, bw.source,
                      a.filename, a.track_name, a.start_time, a.activity_type, a.companions
               FROM block_walks bw
               JOIN activities a ON a.id = bw.activity_id
               WHERE bw.block_id = ?
               ORDER BY a.start_time""",
            (block_id,),
        ).fetchall()
    result = []
    for r in rows:
        d = dict(r)
        d["companions"] = json.loads(d["companions"] or "[]")
        result.append(d)
    return result


# ── Network change reconciliation ─────────────────────────────────────────────

def remove_auto_walks_for_blocks(block_ids: list[str], path: Path = DB_PATH) -> int:
    """
    Delete auto block_walk rows for blocks that no longer exist after a
    network refresh. Returns count of deleted rows.
    """
    if not block_ids:
        return 0
    placeholders = ",".join("?" * len(block_ids))
    with _db(path) as conn:
        cur = conn.execute(
            f"DELETE FROM block_walks WHERE block_id IN ({placeholders}) AND source='auto'",
            block_ids,
        )
        return cur.rowcount


def flag_manual_walks_for_review(block_ids: list[str], path: Path = DB_PATH) -> list[dict]:
    """
    Return manual block_walk rows for gone blocks so the UI can alert the user.
    Does NOT delete them — the UI handles that after user confirmation.
    """
    if not block_ids:
        return []
    placeholders = ",".join("?" * len(block_ids))
    with _db(path) as conn:
        rows = conn.execute(
            f"""SELECT bw.block_id, a.track_name, a.start_time
                FROM block_walks bw
                JOIN activities a ON a.id = bw.activity_id
                WHERE bw.block_id IN ({placeholders}) AND bw.source='manual'""",
            block_ids,
        ).fetchall()
    return [dict(r) for r in rows]
