"""
GPS track-to-block matching.

Matches a GPS track to street blocks. For each block, computes the
coverage ratio (0.0–1.0): fraction of the block's total length that
falls within snap_distance_m of the track.

Returns (block_id, activity_id, coverage_ratio, "auto") tuples ready
for bulk-insert into block_walks.
"""

from __future__ import annotations

import logging
from collections import defaultdict

import geopandas as gpd
import shapely

from config import DEFAULT_SNAP_DISTANCE_M, METRIC_CRS, WGS84_CRS

logger = logging.getLogger(__name__)


def match_track_to_blocks(
    network_gdf: gpd.GeoDataFrame,
    track: dict,
    activity_id: int,
    snap_distance_m: float = DEFAULT_SNAP_DISTANCE_M,
) -> list[tuple[str, int, float, str]]:
    """
    Match a single GPS track against all blocks in the network.

    Args:
        network_gdf:     Edge GDF with block_id and block_length_m columns.
        track:           Parsed track dict from src.gpx_parser.
        activity_id:     DB id for the activity owning this track.
        snap_distance_m: Buffer radius in metres around the GPS track.

    Returns:
        List of (block_id, activity_id, coverage_ratio, "auto") tuples
        for every block where coverage_ratio > 0.
    """
    # Project network to metric CRS
    net_proj = network_gdf.to_crs(METRIC_CRS)
    geoms = net_proj.geometry.values
    block_ids = network_gdf["block_id"].values
    block_lengths = network_gdf["block_length_m"].values

    # block_id → representative block_length_m (all edges in a block share the same value)
    block_length_map: dict[str, float] = {}
    for i, bid in enumerate(block_ids):
        if bid not in block_length_map:
            block_length_map[bid] = float(block_lengths[i])

    # Spatial index over all edge geometries
    tree = shapely.STRtree(geoms)

    # Project track
    track_proj = (
        gpd.GeoSeries([track["track_line"]], crs=WGS84_CRS)
        .to_crs(METRIC_CRS)
        .iloc[0]
    )
    track_buffer = track_proj.buffer(snap_distance_m)

    # Accumulate covered length per block
    block_covered: dict[str, float] = defaultdict(float)
    for idx in tree.query(track_buffer, predicate="intersects"):
        seg = geoms[idx]
        if seg.length == 0:
            continue
        covered = seg.intersection(track_buffer).length
        if covered > 0:
            block_covered[block_ids[idx]] += covered

    # Build result tuples — include ALL blocks with any coverage > 0
    results: list[tuple[str, int, float, str]] = []
    for bid, covered_length in block_covered.items():
        total_length = block_length_map.get(bid, 0.0)
        if total_length <= 0:
            continue
        ratio = min(covered_length / total_length, 1.0)
        results.append((bid, activity_id, ratio, "auto"))

    logger.debug(
        "Activity %d: matched %d blocks (snap=%.0fm)",
        activity_id, len(results), snap_distance_m,
    )
    return results


def match_tracks_to_blocks(
    network_gdf: gpd.GeoDataFrame,
    tracks_with_ids: list[tuple[dict, int]],
    snap_distance_m: float = DEFAULT_SNAP_DISTANCE_M,
) -> list[tuple[str, int, float, str]]:
    """
    Match multiple GPS tracks against the network.

    Args:
        network_gdf:      Edge GDF with block_id and block_length_m columns.
        tracks_with_ids:  List of (track_dict, activity_id) pairs.
        snap_distance_m:  Buffer radius in metres.

    Returns:
        Flat list of (block_id, activity_id, coverage_ratio, "auto") tuples
        for all tracks combined.
    """
    all_results: list[tuple[str, int, float, str]] = []
    for track, activity_id in tracks_with_ids:
        all_results.extend(
            match_track_to_blocks(network_gdf, track, activity_id, snap_distance_m)
        )

    walked_blocks = len({r[0] for r in all_results})
    logger.info(
        "Matched %d track(s) → %d block assignments | snap=%.0fm",
        len(tracks_with_ids), walked_blocks, snap_distance_m,
    )
    return all_results
