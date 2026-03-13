"""
GPS track-to-block matching.

Matches GPS tracks to street blocks (not individual edges). A block is
marked as walked only if the fraction of the block's total length that
falls within snap_distance_m of the track meets coverage_threshold.
When a block is walked, ALL of its edges are marked walked.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from datetime import datetime

import geopandas as gpd
import numpy as np
import shapely

from config import DEFAULT_COVERAGE_THRESHOLD, DEFAULT_SNAP_DISTANCE_M, METRIC_CRS, WGS84_CRS

logger = logging.getLogger(__name__)


def match_tracks_to_segments(
    network_gdf: gpd.GeoDataFrame,
    tracks: list[dict],
    snap_distance_m: float = DEFAULT_SNAP_DISTANCE_M,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
) -> gpd.GeoDataFrame:
    """
    Mark street edges as walked based on block-level GPS coverage.

    Algorithm:
      1. Project network and tracks to EPSG:3857 (metres).
      2. Build a Shapely STRtree spatial index over all edge geometries.
      3. For each track, buffer the track line by snap_distance_m.
      4. Query the STRtree for candidate edges that intersect the buffer.
      5. For each candidate edge, compute the length of the edge that
         falls inside the buffer and accumulate it per block.
      6. For each block where covered_length / block_length >= threshold,
         mark ALL edges in that block as walked.

    Args:
        network_gdf:        Edge GDF with block_id and block_length_m columns.
        tracks:             Parsed track dicts from src.gpx_parser.
        snap_distance_m:    Buffer radius in metres around the GPS track.
        coverage_threshold: Fraction of block length (0–1) required to count
                            as walked. Default 0.95 (95%).

    Returns:
        network_gdf with new columns:
            walked      (bool)
            walk_date   (datetime | None)
            walk_count  (int)
    """
    if not tracks:
        result = network_gdf.copy()
        result["walked"] = False
        result["walk_date"] = None
        result["walk_count"] = 0
        return result

    # Project to metric CRS for accurate distance calculations
    net_proj = network_gdf.to_crs(METRIC_CRS)
    geoms = net_proj.geometry.values
    block_ids = network_gdf["block_id"].values
    block_lengths = network_gdf["block_length_m"].values

    # Build positional index: block_id → list of row indices into geoms
    block_to_indices: dict[str, list[int]] = defaultdict(list)
    for i, bid in enumerate(block_ids):
        block_to_indices[bid].append(i)

    # Build spatial index once for all edges
    tree = shapely.STRtree(geoms)

    n = len(network_gdf)
    walk_count = np.zeros(n, dtype=int)
    walk_dates: list[datetime | None] = [None] * n

    for track in tracks:
        start_time = track.get("start_time")

        # Project track to metric CRS
        track_proj = (
            gpd.GeoSeries([track["track_line"]], crs=WGS84_CRS)
            .to_crs(METRIC_CRS)
            .iloc[0]
        )

        # Buffer = "everywhere I walked"
        track_buffer = track_proj.buffer(snap_distance_m)

        # Accumulate covered length per block from candidate edges
        block_covered: dict[str, float] = defaultdict(float)
        candidate_idxs = tree.query(track_buffer, predicate="intersects")

        for idx in candidate_idxs:
            seg = geoms[idx]
            if seg.length == 0:
                continue
            covered = seg.intersection(track_buffer).length
            if covered > 0:
                block_covered[block_ids[idx]] += covered

        # Mark blocks that meet the coverage threshold
        for bid, covered_length in block_covered.items():
            indices = block_to_indices[bid]
            # All edges in a block share the same block_length_m
            total_length = block_lengths[indices[0]]
            if total_length > 0 and covered_length / total_length >= coverage_threshold:
                for idx in indices:
                    walk_count[idx] += 1
                    if start_time is not None:
                        if walk_dates[idx] is None or start_time < walk_dates[idx]:
                            walk_dates[idx] = start_time

    result = network_gdf.copy()
    result["walked"] = walk_count > 0
    result["walk_count"] = walk_count
    result["walk_date"] = walk_dates

    walked_blocks = len({bid for i, bid in enumerate(block_ids) if walk_count[i] > 0})
    walked_km = result.loc[result["walked"], "length_m"].sum() / 1000
    logger.info(
        "Walked %d blocks (%.1f km) from %d track(s) | snap=%.0fm threshold=%.0f%%",
        walked_blocks, walked_km, len(tracks), snap_distance_m, coverage_threshold * 100,
    )
    return result
