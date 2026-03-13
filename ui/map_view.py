"""
Folium map construction.

Renders walked (green) and unwalked (gray) street segments with
optional GPX track overlays and street name tooltips.
"""

from __future__ import annotations

import folium
import geopandas as gpd
import pandas as pd

from config import COLOR_GPX_TRACK, COLOR_UNWALKED, COLOR_WALKED, MAP_TILES, MAP_ZOOM_START, WGS84_CRS


def build_map(
    network_gdf: gpd.GeoDataFrame,
    tracks: list[dict] | None = None,
) -> folium.Map:
    """
    Build a Folium map displaying street coverage.

    Args:
        network_gdf: Street network with 'walked', 'walk_date', 'walk_count' columns.
        tracks:      Parsed track dicts (optional) — drawn as orange lines.

    Returns:
        Configured folium.Map ready for st_folium().
    """
    # Ensure WGS84 for display
    gdf = network_gdf.to_crs(WGS84_CRS) if network_gdf.crs.to_epsg() != 4326 else network_gdf

    # Center on network centroid
    centroid = gdf.geometry.unary_union.centroid
    m = folium.Map(
        location=[centroid.y, centroid.x],
        zoom_start=MAP_ZOOM_START,
        tiles=MAP_TILES,
    )

    # Prepare display columns shared by both layers
    def _prep_gdf(df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        df = df.copy()
        if "walk_date" in df.columns:
            df["walk_date"] = df["walk_date"].apply(
                lambda x: x.strftime("%Y-%m-%d") if pd.notna(x) else ""
            )
        for col in ("length_m", "block_length_m"):
            if col in df.columns:
                df[col] = df[col].round(1)
        return df

    _HIGHLIGHT = {"weight": 6, "color": "#facc15", "opacity": 1.0}  # yellow on hover/click

    # ── Unwalked streets layer ────────────────────────────────────────────────
    unwalked = _prep_gdf(gdf[~gdf["walked"]])
    if not unwalked.empty:
        _unwalked_fields = [c for c in ["block_id", "name", "block_length_m", "highway"] if c in unwalked.columns]
        _unwalked_aliases = {"block_id": "Block ID", "name": "Street", "block_length_m": "Block length (m)", "highway": "Type"}
        folium.GeoJson(
            unwalked.__geo_interface__,
            name="Unwalked Streets",
            show=True,
            style_function=lambda _: {"color": COLOR_UNWALKED, "weight": 2, "opacity": 0.6},
            highlight_function=lambda _: _HIGHLIGHT,
            tooltip=folium.GeoJsonTooltip(
                fields=_unwalked_fields,
                aliases=[_unwalked_aliases.get(f, f) for f in _unwalked_fields],
            ),
            popup=folium.GeoJsonPopup(
                fields=_unwalked_fields,
                aliases=[_unwalked_aliases.get(f, f) for f in _unwalked_fields],
            ),
        ).add_to(m)

    # ── Walked streets layer ──────────────────────────────────────────────────
    walked = _prep_gdf(gdf[gdf["walked"]])
    if not walked.empty:
        _walked_fields = [c for c in ["block_id", "name", "block_length_m", "walk_date", "walk_count"] if c in walked.columns]
        _walked_aliases = {
            "block_id": "Block ID", "name": "Street", "block_length_m": "Block length (m)",
            "walk_date": "First walked", "walk_count": "Times walked",
        }
        folium.GeoJson(
            walked.__geo_interface__,
            name="Walked Streets",
            style_function=lambda _: {"color": COLOR_WALKED, "weight": 4, "opacity": 0.9},
            highlight_function=lambda _: _HIGHLIGHT,
            tooltip=folium.GeoJsonTooltip(
                fields=_walked_fields,
                aliases=[_walked_aliases.get(f, f) for f in _walked_fields],
            ),
            popup=folium.GeoJsonPopup(
                fields=_walked_fields,
                aliases=[_walked_aliases.get(f, f) for f in _walked_fields],
            ),
        ).add_to(m)

    # ── GPX track overlays ────────────────────────────────────────────────────
    if tracks:
        track_group = folium.FeatureGroup(name="GPS Tracks", show=True)
        for track in tracks:
            line = track["track_line"]
            # Handle both LineString and MultiLineString
            if line.geom_type == "LineString":
                coords = [(y, x) for x, y in line.coords]
                folium.PolyLine(
                    coords,
                    color=COLOR_GPX_TRACK,
                    weight=3,
                    opacity=0.7,
                    tooltip=track.get("track_name", "GPS Track"),
                ).add_to(track_group)
            elif line.geom_type == "MultiLineString":
                for segment in line.geoms:
                    coords = [(y, x) for x, y in segment.coords]
                    folium.PolyLine(
                        coords,
                        color=COLOR_GPX_TRACK,
                        weight=3,
                        opacity=0.7,
                        tooltip=track.get("track_name", "GPS Track"),
                    ).add_to(track_group)
        track_group.add_to(m)

    folium.LayerControl(collapsed=False).add_to(m)
    return m
