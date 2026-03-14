"""
Street network fetching and caching.

Replaces the entire R pipeline (generate_network_00 → 03):
  - osmdata API query
  - dodgr edge table construction
  - Segment_Labeling.R recursive walk algorithm
  - Parallel_Edges.R detection

osmnx with simplify=True produces the same result — LineStrings between
intersections and dead ends — in a single function call.
"""

from __future__ import annotations

import hashlib
import logging
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import osmnx as ox

from config import (
    EXCLUDE_HIGHWAY_TYPES,
    METRIC_CRS,
    OSM_NETWORK_TYPE,
    WGS84_CRS,
)
from src.paths import NETWORK_CACHE_DIR

logger = logging.getLogger(__name__)


# ── Cache helpers ─────────────────────────────────────────────────────────────

def _cache_key(place_name: str | None, bbox: tuple[float, float, float, float] | None) -> str:
    """Deterministic cache key from place name or bounding box."""
    if place_name:
        raw = place_name.lower().strip()
    elif bbox:
        # Round to 3 decimal places (~111m precision) for stable keys
        raw = "_".join(f"{v:.3f}" for v in bbox)
    else:
        raise ValueError("Either place_name or bbox must be provided.")
    return hashlib.md5(raw.encode()).hexdigest()[:12]


def _cache_path(key: str) -> Path:
    return NETWORK_CACHE_DIR / f"network_{key}.gpkg"


def _cache_age_days(key: str) -> float | None:
    path = _cache_path(key)
    if not path.exists():
        return None
    mtime = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)
    return (datetime.now(tz=timezone.utc) - mtime).total_seconds() / 86400


def cache_exists(key: str) -> bool:
    return _cache_path(key).exists()


def invalidate_cache(key: str) -> None:
    path = _cache_path(key)
    if path.exists():
        path.unlink()
        logger.info("Invalidated network cache: %s", key)


def list_cached_networks() -> list[dict]:
    """Return metadata for all cached networks."""
    results = []
    for f in NETWORK_CACHE_DIR.glob("network_*.gpkg"):
        age = _cache_age_days(f.stem.removeprefix("network_"))
        results.append({"file": f.name, "age_days": round(age, 1) if age else None})
    return results


# ── Network fetching ──────────────────────────────────────────────────────────

def _build_network_gdf(G) -> gpd.GeoDataFrame:
    """
    Convert an osmnx graph to a clean GeoDataFrame of street segments.

    Each row is one simplified segment (LineString between two intersections
    or dead ends) — equivalent to the R project's segment_id concept.
    """
    _, edges = ox.graph_to_gdfs(G)
    gdf = edges.reset_index()

    # Filter excluded highway types.
    # Special case: service roads that have a name are real streets (some named
    # residential streets are tagged highway=service in OSM) — keep them.
    # Only filter unnamed service roads (driveways, parking aisles, alleys).
    import pandas as pd

    def _is_excluded(row) -> bool:
        hw = row["highway"]
        hw_types = hw if isinstance(hw, list) else [hw]
        if not all(t in EXCLUDE_HIGHWAY_TYPES for t in hw_types):
            return False
        if all(t == "service" for t in hw_types):
            name = row["name"]
            if isinstance(name, list):
                name = name[0] if name else None
            # pd.notna required — np.nan is truthy in Python so `if name` alone fails
            has_name = pd.notna(name) and bool(str(name).strip())
            return not has_name  # keep named, exclude unnamed
        return True

    mask_keep = ~gdf.apply(_is_excluded, axis=1)
    gdf = gdf[mask_keep].copy()

    # Stable segment ID from OSM node IDs (reproducible across fetches)
    gdf["segment_id"] = gdf["u"].astype(str) + "_" + gdf["v"].astype(str)

    # Drop reverse-direction duplicates for bidirectional streets.
    # Sort (u, v) so both directions share the same key, keep first occurrence.
    gdf["_dedup_key"] = gdf.apply(
        lambda r: "_".join(sorted([str(r["u"]), str(r["v"])])), axis=1
    )
    gdf = gdf.drop_duplicates(subset="_dedup_key").drop(columns="_dedup_key")

    # Compute length in meters using projected CRS
    gdf_proj = gdf.to_crs(METRIC_CRS)
    gdf["length_m"] = gdf_proj.geometry.length

    # Keep useful columns only
    keep_cols = ["segment_id", "u", "v", "name", "highway", "length_m", "geometry"]
    gdf = gdf[[c for c in keep_cols if c in gdf.columns]].copy()
    gdf = gdf.set_crs(WGS84_CRS, allow_override=True)

    logger.info("Network built: %d segments, %.1f km total", len(gdf), gdf["length_m"].sum() / 1000)

    # Group edges into blocks (named street sections between cross streets)
    from src.block_builder import build_blocks
    gdf = build_blocks(gdf)

    return gdf


def fetch_network(
    place_name: str | None = None,
    bbox: tuple[float, float, float, float] | None = None,
) -> gpd.GeoDataFrame:
    """
    Fetch and simplify the street network from OSM.

    Args:
        place_name: Geocodable place string, e.g. "Charlotte, North Carolina".
        bbox: (north, south, east, west) bounding box in WGS84 degrees.

    Returns:
        GeoDataFrame with one row per street segment.
    """
    ox.settings.log_console = False
    ox.settings.use_cache = True

    logger.info("Fetching OSM network for: %s %s", place_name, bbox)

    if place_name:
        G = ox.graph_from_place(
            place_name,
            network_type=OSM_NETWORK_TYPE,
            simplify=True,
        )
    elif bbox:
        north, south, east, west = bbox
        G = ox.graph_from_bbox(
            north=north, south=south, east=east, west=west,
            network_type=OSM_NETWORK_TYPE,
            simplify=True,
        )
    else:
        raise ValueError("Either place_name or bbox must be provided.")

    return _build_network_gdf(G)


def get_or_fetch_network(
    place_name: str | None = None,
    bbox: tuple[float, float, float, float] | None = None,
    force_refresh: bool = False,
) -> tuple[gpd.GeoDataFrame, str]:
    """
    Load the network from disk cache if fresh, otherwise fetch from OSM.

    Returns:
        (GeoDataFrame, cache_key)
    """
    key = _cache_key(place_name, bbox)
    path = _cache_path(key)

    if force_refresh:
        invalidate_cache(key)

    if cache_exists(key):
        logger.info("Loading network from cache: %s", path)
        gdf = gpd.read_file(path)
        return gdf, key

    logger.info("Cache miss — fetching from OSM")
    gdf = fetch_network(place_name=place_name, bbox=bbox)
    gdf.to_file(path, driver="GPKG")
    logger.info("Network cached to: %s", path)
    return gdf, key
