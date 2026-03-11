"""
GPX file parsing.

Converts raw .gpx bytes (Garmin Connect / Strava format) into structured
track dictionaries containing a Shapely geometry and metadata.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

import gpxpy
import gpxpy.gpx
from shapely.geometry import LineString, MultiLineString, Point

logger = logging.getLogger(__name__)


def parse_gpx_file(file_bytes: bytes, filename: str = "") -> dict | None:
    """
    Parse a single GPX file into a structured track dictionary.

    Args:
        file_bytes: Raw bytes of the .gpx file.
        filename: Original filename, used as fallback track name.

    Returns:
        Dict with keys:
            track_name    : str
            activity_type : str | None
            start_time    : datetime | None
            end_time      : datetime | None
            point_count   : int
            bbox          : (min_lon, min_lat, max_lon, max_lat)
            track_line    : MultiLineString — full path, gap-aware
            points        : list of (lon, lat) tuples
        Returns None if the file contains no usable track points.
    """
    try:
        gpx = gpxpy.parse(file_bytes.decode("utf-8") if isinstance(file_bytes, bytes) else file_bytes)
    except Exception as exc:
        logger.warning("Failed to parse GPX file '%s': %s", filename, exc)
        return None

    # Collect segments as separate LineStrings to preserve GPS dropout gaps
    linestrings: list[LineString] = []
    all_points: list[tuple[float, float]] = []
    all_times: list[datetime] = []

    for track in gpx.tracks:
        for segment in track.segments:
            coords = [
                (pt.longitude, pt.latitude)
                for pt in segment.points
                if pt.longitude is not None and pt.latitude is not None
            ]
            if len(coords) >= 2:
                linestrings.append(LineString(coords))
                all_points.extend(coords)

            for pt in segment.points:
                if pt.time:
                    all_times.append(pt.time)

    # Also handle waypoints-only GPX files (less common)
    for wpt in gpx.waypoints:
        if wpt.longitude is not None and wpt.latitude is not None:
            all_points.append((wpt.longitude, wpt.latitude))

    if not all_points:
        logger.warning("No track points found in '%s'", filename)
        return None

    if not linestrings:
        logger.warning("Not enough consecutive points to form a line in '%s'", filename)
        return None

    track_line = MultiLineString(linestrings) if len(linestrings) > 1 else linestrings[0]

    lons = [p[0] for p in all_points]
    lats = [p[1] for p in all_points]

    # Prefer track name from GPX metadata; fall back to filename
    track_name = filename
    activity_type = None
    if gpx.tracks:
        first_track = gpx.tracks[0]
        if first_track.name:
            track_name = first_track.name
        if first_track.type:
            activity_type = first_track.type

    return {
        "track_name": track_name,
        "filename": filename,
        "activity_type": activity_type,
        "start_time": min(all_times) if all_times else None,
        "end_time": max(all_times) if all_times else None,
        "point_count": len(all_points),
        "bbox": (min(lons), min(lats), max(lons), max(lats)),
        "track_line": track_line,
        "points": all_points,
    }


def parse_multiple_gpx(files: list) -> list[dict]:
    """
    Parse multiple GPX files concurrently.

    Args:
        files: List of Streamlit UploadedFile objects or (bytes, filename) tuples.

    Returns:
        List of parsed track dicts (failed files are silently skipped).
    """
    def _read(f) -> tuple[bytes, str]:
        if hasattr(f, "read"):
            # Streamlit UploadedFile
            data = f.read()
            name = getattr(f, "name", "unknown.gpx")
            # Reset for potential re-reads
            if hasattr(f, "seek"):
                f.seek(0)
            return data, name
        # (bytes, filename) tuple
        return f[0], f[1]

    results: list[dict] = []
    with ThreadPoolExecutor(max_workers=min(len(files), 8)) as executor:
        futures = {executor.submit(_read, f): f for f in files}
        for future in as_completed(futures):
            raw_bytes, name = future.result()
            parsed = parse_gpx_file(raw_bytes, filename=name)
            if parsed:
                results.append(parsed)

    # Sort chronologically by start_time so history displays correctly
    results.sort(key=lambda t: t["start_time"] or datetime.min)
    logger.info("Parsed %d / %d GPX files successfully", len(results), len(files))
    return results


def combined_bbox(tracks: list[dict]) -> tuple[float, float, float, float]:
    """
    Compute the bounding box spanning all provided tracks.

    Returns:
        (min_lon, min_lat, max_lon, max_lat)
    """
    min_lons, min_lats, max_lons, max_lats = zip(*[t["bbox"] for t in tracks])
    return min(min_lons), min(min_lats), max(max_lons), max(max_lats)


def bbox_to_osmnx(bbox: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
    """
    Convert (min_lon, min_lat, max_lon, max_lat) to osmnx (north, south, east, west).
    Adds a small padding buffer so edge streets are fully included.
    """
    min_lon, min_lat, max_lon, max_lat = bbox
    pad = 0.005  # ~500m padding at mid-latitudes
    return (
        max_lat + pad,   # north
        min_lat - pad,   # south
        max_lon + pad,   # east
        min_lon - pad,   # west
    )
