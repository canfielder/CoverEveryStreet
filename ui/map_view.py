"""
Folium map construction.

Renders street blocks colored by activity category (solo walk, solo run,
walk+Kelsey, …), unwalked blocks in gray, and optional GPX track overlays.
"""

from __future__ import annotations

import folium
import geopandas as gpd
import pandas as pd

from config import (
    ACTIVITY_COLORS,
    COLOR_FALLBACK,
    COLOR_GPX_TRACK,
    COLOR_UNWALKED,
    MAP_TILES,
    MAP_ZOOM_START,
    WGS84_CRS,
)

_HIGHLIGHT = {"weight": 6, "color": "#facc15", "opacity": 1.0}


def build_map(
    network_gdf: gpd.GeoDataFrame,
    walked_blocks: dict[str, dict],
    tracks: list[dict] | None = None,
) -> folium.Map:
    """
    Build a Folium map showing block coverage.

    Args:
        network_gdf:   Edge GDF with block_id, name, block_length_m columns.
        walked_blocks: {block_id: {activity_type, companions, walk_date, source}}
                       from database.get_walked_blocks().
        tracks:        Parsed track dicts (optional) — drawn as orange lines.

    Returns:
        Configured folium.Map ready for st_folium().
    """
    gdf = network_gdf.to_crs(WGS84_CRS) if network_gdf.crs.to_epsg() != 4326 else network_gdf

    centroid = gdf.geometry.unary_union.centroid
    m = folium.Map(
        location=[centroid.y, centroid.x],
        zoom_start=MAP_ZOOM_START,
        tiles=MAP_TILES,
    )

    # Attach walk metadata to each edge row
    gdf = gdf.copy()
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

    # Round display columns
    for col in ("length_m", "block_length_m"):
        if col in gdf.columns:
            gdf[col] = gdf[col].round(1)

    # ── Unwalked layer ────────────────────────────────────────────────────────
    unwalked = gdf[~gdf["_walked"]]
    if not unwalked.empty:
        _fields = [c for c in ["block_id", "name", "block_length_m", "highway"] if c in unwalked.columns]
        _aliases = {"block_id": "Block ID", "name": "Street", "block_length_m": "Block length (m)", "highway": "Type"}
        folium.GeoJson(
            unwalked.__geo_interface__,
            name="Unwalked Streets",
            show=True,
            style_function=lambda _: {"color": COLOR_UNWALKED, "weight": 2, "opacity": 0.6},
            highlight_function=lambda _: _HIGHLIGHT,
            tooltip=folium.GeoJsonTooltip(
                fields=_fields,
                aliases=[_aliases.get(f, f) for f in _fields],
            ),
            popup=folium.GeoJsonPopup(
                fields=_fields,
                aliases=[_aliases.get(f, f) for f in _fields],
            ),
        ).add_to(m)

    # ── Walked layers (one FeatureGroup per color category) ───────────────────
    walked = gdf[gdf["_walked"]].copy()
    if not walked.empty:
        walked["_color"] = walked.apply(_edge_color, axis=1)

        # Group by color to minimize number of GeoJson layers
        for color, group in walked.groupby("_color"):
            category = _color_label(color)
            _fields = [c for c in ["block_id", "name", "block_length_m", "_walk_date", "_activity_type", "_companions"] if c in group.columns]
            _aliases = {
                "block_id": "Block ID", "name": "Street",
                "block_length_m": "Block length (m)",
                "_walk_date": "First walked", "_activity_type": "Activity",
                "_companions": "Companions",
            }
            folium.GeoJson(
                group.__geo_interface__,
                name=category,
                style_function=lambda _, c=color: {"color": c, "weight": 4, "opacity": 0.9},
                highlight_function=lambda _: _HIGHLIGHT,
                tooltip=folium.GeoJsonTooltip(
                    fields=_fields,
                    aliases=[_aliases.get(f, f) for f in _fields],
                ),
                popup=folium.GeoJsonPopup(
                    fields=_fields,
                    aliases=[_aliases.get(f, f) for f in _fields],
                ),
            ).add_to(m)

    # ── GPX track overlays ────────────────────────────────────────────────────
    if tracks:
        track_group = folium.FeatureGroup(name="GPS Tracks", show=True)
        for track in tracks:
            line = track["track_line"]
            label = track.get("track_name", "GPS Track")
            if line.geom_type == "LineString":
                folium.PolyLine(
                    [(y, x) for x, y in line.coords],
                    color=COLOR_GPX_TRACK, weight=3, opacity=0.7, tooltip=label,
                ).add_to(track_group)
            elif line.geom_type == "MultiLineString":
                for seg in line.geoms:
                    folium.PolyLine(
                        [(y, x) for x, y in seg.coords],
                        color=COLOR_GPX_TRACK, weight=3, opacity=0.7, tooltip=label,
                    ).add_to(track_group)
        track_group.add_to(m)

    folium.LayerControl(collapsed=False).add_to(m)
    return m


_MINI_TRACK_COLORS = ["#f97316", "#3b82f6", "#a855f7", "#ec4899", "#22c55e"]


def build_mini_map(
    block_id: str,
    network_gdf: gpd.GeoDataFrame,
    tracks: list[dict] | None = None,
) -> folium.Map:
    """Mini Folium map zoomed to a specific block with GPS track overlays."""
    gdf = network_gdf.to_crs(WGS84_CRS) if network_gdf.crs.to_epsg() != 4326 else network_gdf

    block_edges = gdf[gdf["block_id"] == block_id]
    if block_edges.empty:
        block_edges = gdf  # fallback: show all

    bounds = block_edges.geometry.total_bounds  # (minx, miny, maxx, maxy)
    pad_x = (bounds[2] - bounds[0]) * 0.3 or 0.002
    pad_y = (bounds[3] - bounds[1]) * 0.3 or 0.002
    sw = [bounds[1] - pad_y, bounds[0] - pad_x]
    ne = [bounds[3] + pad_y, bounds[2] + pad_x]

    centroid = block_edges.geometry.unary_union.centroid
    m = folium.Map(
        location=[centroid.y, centroid.x],
        tiles=MAP_TILES,
        zoom_start=17,
    )
    m.fit_bounds([sw, ne])

    # All network edges as thin gray context
    folium.GeoJson(
        gdf.__geo_interface__,
        style_function=lambda _: {"color": COLOR_UNWALKED, "weight": 1.5, "opacity": 0.5},
    ).add_to(m)

    # Selected block highlighted yellow
    folium.GeoJson(
        block_edges.__geo_interface__,
        style_function=lambda _: {"color": "#facc15", "weight": 6, "opacity": 1.0},
    ).add_to(m)

    # Track overlays
    if tracks:
        for i, track in enumerate(tracks):
            color = _MINI_TRACK_COLORS[i % len(_MINI_TRACK_COLORS)]
            line = track["track_line"]
            label = track.get("track_name", "GPS Track")
            if line.geom_type == "LineString":
                folium.PolyLine(
                    [(y, x) for x, y in line.coords],
                    color=color, weight=3, opacity=0.85, tooltip=label,
                ).add_to(m)
            elif line.geom_type == "MultiLineString":
                for seg in line.geoms:
                    folium.PolyLine(
                        [(y, x) for x, y in seg.coords],
                        color=color, weight=3, opacity=0.85, tooltip=label,
                    ).add_to(m)

    return m


# ── Color helpers ─────────────────────────────────────────────────────────────

def _edge_color(row) -> str:
    """Return hex color for a walked edge based on activity type + companions."""
    atype = row.get("_activity_type") or "walk"
    companions_str = row.get("_companions") or ""
    companions = frozenset(c.strip().lower() for c in companions_str.split(",") if c.strip())
    key = (atype, companions)
    if key in ACTIVITY_COLORS:
        return ACTIVITY_COLORS[key]
    # Try without companions as fallback
    if (atype, frozenset()) in ACTIVITY_COLORS:
        return ACTIVITY_COLORS[(atype, frozenset())]
    return COLOR_FALLBACK


def _color_label(color: str) -> str:
    """Human-readable layer name for a given hex color."""
    reverse = {v: k for k, v in ACTIVITY_COLORS.items()}
    if color not in reverse:
        return "Walked Streets"
    atype, companions = reverse[color]
    label = atype.capitalize()
    if companions:
        label += " + " + ", ".join(sorted(c.capitalize() for c in companions))
    return label
