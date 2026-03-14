"""
GPX file parsing.

Converts raw .gpx bytes (Garmin Connect / Strava format) into structured
track dictionaries containing a Shapely geometry and metadata.
Activity type is normalized via ACTIVITY_TYPE_MAP; pace-based fallback
is used when the GPX <type> tag is absent.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import gpxpy
import gpxpy.gpx
from shapely.geometry import LineString, MultiLineString

from config import ACTIVITY_TYPE_MAP, PACE_RUN_THRESHOLD_KPH

logger = logging.getLogger(__name__)


def _normalize_activity_type(raw: str | None, tracks: list) -> tuple[str, str | None]:
    """
    Return (normalized_type, raw_type).

    Priority:
      1. GPX <type> tag → look up in ACTIVITY_TYPE_MAP
      2. Pace calculated from track timestamps
      3. Default to "walk"
    """
    raw_type = raw.strip().lower() if raw else None

    if raw_type and raw_type in ACTIVITY_TYPE_MAP:
        return ACTIVITY_TYPE_MAP[raw_type], raw

    if raw_type:
        # Unknown type string — still store raw, default to walk
        logger.debug("Unknown activity type '%s', defaulting to 'walk'", raw)
        return "walk", raw

    # Pace-based fallback
    avg_kph = _average_speed_kph(tracks)
    if avg_kph is not None:
        normalized = "run" if avg_kph > PACE_RUN_THRESHOLD_KPH else "walk"
        logger.debug("No GPX type tag; pace %.1f kph → '%s'", avg_kph, normalized)
        return normalized, None

    return "walk", None


def _average_speed_kph(tracks: list) -> float | None:
    """Compute average speed in km/h from track point timestamps and coords."""
    from math import radians, cos, sin, asin, sqrt

    def haversine(lon1, lat1, lon2, lat2) -> float:
        R = 6371.0
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
        dlon, dlat = lon2 - lon1, lat2 - lat1
        a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
        return 2 * R * asin(sqrt(a))

    total_km, total_sec = 0.0, 0.0
    for track in tracks:
        for segment in track.segments:
            pts = [p for p in segment.points if p.time and p.longitude and p.latitude]
            for i in range(1, len(pts)):
                p0, p1 = pts[i - 1], pts[i]
                total_km += haversine(p0.longitude, p0.latitude, p1.longitude, p1.latitude)
                total_sec += (p1.time - p0.time).total_seconds()

    if total_sec < 60:
        return None
    return (total_km / total_sec) * 3600


def parse_gpx_file(file_bytes: bytes, filename: str = "") -> dict | None:
    """
    Parse a single GPX file into a structured track dictionary.

    Returns dict with keys:
        track_name      : str
        filename        : str
        content_hash    : str  — MD5 hex of raw file bytes (inbox dedup)
        activity_type   : str  — normalized ("walk", "run", …)
        raw_gpx_type    : str | None  — original <type> tag
        start_time      : datetime | None (timezone-aware)
        end_time        : datetime | None (timezone-aware)
        point_count     : int
        bbox            : (min_lon, min_lat, max_lon, max_lat)
        track_line      : LineString | MultiLineString
        points          : list of (lon, lat) tuples
    Returns None if the file contains no usable track points.
    """
    import hashlib
    raw_bytes = file_bytes if isinstance(file_bytes, bytes) else file_bytes.encode()
    content_hash = hashlib.md5(raw_bytes).hexdigest()

    try:
        gpx = gpxpy.parse(raw_bytes.decode("utf-8"))
    except Exception as exc:
        logger.warning("Failed to parse GPX file '%s': %s", filename, exc)
        return None

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
                    t = pt.time
                    if t.tzinfo is None:
                        t = t.replace(tzinfo=timezone.utc)
                    all_times.append(t)

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

    # Track name
    track_name = filename
    raw_type = None
    if gpx.tracks:
        first_track = gpx.tracks[0]
        if first_track.name:
            track_name = first_track.name
        raw_type = first_track.type or None

    activity_type, raw_gpx_type = _normalize_activity_type(raw_type, gpx.tracks)

    return {
        "track_name":   track_name,
        "filename":     filename,
        "content_hash": content_hash,
        "activity_type":  activity_type,
        "raw_gpx_type":   raw_gpx_type,
        "start_time":   min(all_times) if all_times else None,
        "end_time":     max(all_times) if all_times else None,
        "point_count":  len(all_points),
        "bbox":         (min(lons), min(lats), max(lons), max(lats)),
        "track_line":   track_line,
        "points":       all_points,
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
