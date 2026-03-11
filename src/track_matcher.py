"""
GPS track-to-street matching.

Completes the logic left unfinished in the R project's
eda_paint_lines_w_gps_tracking.Rmd — using a spatial index (STRtree)
instead of a brute-force loop for performance.
"""

from __future__ import annotations

import logging
from datetime import datetime

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import MultiLineString

from config import DEFAULT_SNAP_DISTANCE_M, METRIC_CRS, WGS84_CRS

logger = logging.getLogger(__name__)


def match_tracks_to_segments(
    network_gdf: gpd.GeoDataFrame,
    tracks: list[dict],
    snap_distance_m: float = DEFAULT_SNAP_DISTANCE_M,
) -> gpd.GeoDataFrame:
    """
    Mark street segments as walked based on proximity to GPS tracks.

    Algorithm:
      1. Project network and tracks to EPSG:3857 (metres).
      2. Build a Shapely STRtree spatial index over all segment geometries.
      3. For each track, buffer the track line by snap_distance_m.
      4. Query the STRtree for candidate segments that intersect the buffer.
      5. Run an exact distance check on candidates.
      6. Record the earliest walk date and cumulative walk count per segment.

    Args:
        network_gdf:    Street network from core.network (WGS84 CRS).
        tracks:         Parsed track dicts from core.gpx_parser.
        snap_distance_m: Maximum metres from track to segment to count as walked.

    Returns:
        network_gdf with new columns:
            walked      (bool)
            walk_date   (datetime | None) — date of first walk covering this segment
            walk_count  (int)             — number of tracks covering this segment
    """
    if not tracks:
        network_gdf = network_gdf.copy()
        network_gdf["walked"] = False
        network_gdf["walk_date"] = None
        network_gdf["walk_count"] = 0
        return network_gdf

    # Project network to metric CRS for accurate distance calculations
    net_proj = network_gdf.to_crs(METRIC_CRS)
    geoms = net_proj.geometry.values

    # Build spatial index once across all segments
    tree = shapely.STRtree(geoms)

    n = len(network_gdf)
    walk_count = np.zeros(n, dtype=int)
    walk_dates: list[datetime | None] = [None] * n

    for track in tracks:
        track_line = track["track_line"]
        start_time = track.get("start_time")

        # Project track geometry to metric CRS
        track_proj = (
            gpd.GeoSeries([track_line], crs=WGS84_CRS)
            .to_crs(METRIC_CRS)
            .iloc[0]
        )

        # Buffer the track to get the search envelope
        track_buffer = track_proj.buffer(snap_distance_m)

        # STRtree query returns integer positional indices into geoms array
        candidate_idxs = tree.query(track_buffer, predicate="intersects")

        # Exact distance check on candidates (intersects predicate is already
        # exact for polygons vs lines, but explicit check makes the threshold clear)
        for idx in candidate_idxs:
            if geoms[idx].distance(track_proj) <= snap_distance_m:
                walk_count[idx] += 1
                if start_time is not None:
                    if walk_dates[idx] is None or start_time < walk_dates[idx]:
                        walk_dates[idx] = start_time

    result = network_gdf.copy()
    result["walked"] = walk_count > 0
    result["walk_count"] = walk_count
    result["walk_date"] = walk_dates

    walked_n = int(result["walked"].sum())
    walked_km = result.loc[result["walked"], "length_m"].sum() / 1000
    logger.info(
        "Matched %d segments (%.1f km) from %d track(s) at snap=%.0fm",
        walked_n, walked_km, len(tracks), snap_distance_m,
    )
    return result


def merge_walk_results(
    existing: gpd.GeoDataFrame,
    new_match: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Merge a newly matched result into an existing walked network.

    Used when the user uploads additional GPX files incrementally —
    previously walked segments are preserved.
    """
    merged = existing.copy()

    new_walked_mask = new_match["walked"].values
    merged.loc[new_walked_mask, "walked"] = True
    merged["walk_count"] = (
        merged["walk_count"].fillna(0).astype(int)
        + new_match["walk_count"].fillna(0).astype(int)
    )

    # Keep the earliest walk date across both results
    for idx in np.where(new_walked_mask)[0]:
        new_date = new_match.iloc[idx]["walk_date"]
        old_date = merged.iloc[idx]["walk_date"]
        if new_date is not None:
            if old_date is None or new_date < old_date:
                merged.iloc[idx, merged.columns.get_loc("walk_date")] = new_date

    return merged
