"""
Publish a map snapshot for the read-only Streamlit Cloud viewer.

Reads the current database state + the cached network for a given area,
joins walked-block metadata onto the network GDF, and writes:

    data/published/map_data.gpkg   — GeoPackage with walked status + metadata
    data/published/meta.json       — timestamp + summary stats

Run this after processing new activities:
    uv run python -m scripts.publish_snapshot --place "Charlotte, North Carolina"

The output files should be committed to git so Streamlit Cloud can read them.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Ensure project root is on sys.path when run as a script
sys.path.insert(0, str(Path(__file__).parent.parent))

import src.database as db
from src.network import _cache_key, get_or_fetch_network
from src.paths import PUBLISHED_DIR


def publish(
    place_name: str | None = None,
    bbox: tuple | None = None,
    coverage_threshold: float = 0.80,
) -> None:
    if not place_name and not bbox:
        raise ValueError("Provide --place or --bbox")

    print(f"Loading network for: {place_name or bbox}")
    network_gdf, _ = get_or_fetch_network(place_name=place_name, bbox=bbox)

    print("Fetching walked blocks from database…")
    walked_blocks = db.get_walked_blocks(coverage_threshold=coverage_threshold)

    gdf = network_gdf.copy()
    gdf["_walked"] = gdf["block_id"].isin(walked_blocks)
    gdf["_walk_date"] = gdf["block_id"].map(
        lambda bid: walked_blocks[bid]["walk_date"] if bid in walked_blocks else None
    )
    gdf["_activity_type"] = gdf["block_id"].map(
        lambda bid: walked_blocks[bid]["activity_type"] if bid in walked_blocks else None
    )
    gdf["_companions"] = gdf["block_id"].map(
        lambda bid: ", ".join(walked_blocks[bid]["companions"]) if bid in walked_blocks else ""
    )

    PUBLISHED_DIR.mkdir(parents=True, exist_ok=True)
    history_dir = PUBLISHED_DIR / "history"
    history_dir.mkdir(exist_ok=True)

    now = datetime.now(tz=timezone.utc)
    timestamp = now.strftime("%Y%m%d_%H%M%S")

    out_path = PUBLISHED_DIR / "map_data.gpkg"
    gdf.to_file(out_path, driver="GPKG")
    print(f"Wrote: {out_path}")

    # Timestamped archive copy for local rollback
    archive_path = history_dir / f"map_data_{timestamp}.gpkg"
    gdf.to_file(archive_path, driver="GPKG")
    print(f"Archived: {archive_path}")

    n_walked = gdf["_walked"].sum()
    total = len(gdf["block_id"].unique())
    meta = {
        "published_at": now.strftime("%Y-%m-%d %H:%M UTC"),
        "n_walked_blocks": int(n_walked),
        "total_blocks": int(total),
        "pct_walked": round(n_walked / total * 100, 1) if total else 0.0,
        "coverage_threshold": coverage_threshold,
    }
    meta_path = PUBLISHED_DIR / "meta.json"
    meta_path.write_text(json.dumps(meta, indent=2))

    # Timestamped archive copy of meta
    archive_meta = history_dir / f"meta_{timestamp}.json"
    archive_meta.write_text(json.dumps(meta, indent=2))

    print(f"Wrote: {meta_path}")
    print(f"Summary: {n_walked}/{total} blocks walked ({meta['pct_walked']}%)")
    print(f"\nTo roll back: cp {archive_path} {out_path}")


def _parse_bbox(s: str) -> tuple:
    parts = [float(x) for x in s.split(",")]
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("bbox must be 'north,south,east,west'")
    return tuple(parts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--place", help="Geocodable place name")
    parser.add_argument("--bbox", type=_parse_bbox, help="north,south,east,west")
    parser.add_argument("--threshold", type=float, default=0.80, help="Coverage threshold (default 0.80)")
    args = parser.parse_args()
    publish(place_name=args.place, bbox=args.bbox, coverage_threshold=args.threshold)
