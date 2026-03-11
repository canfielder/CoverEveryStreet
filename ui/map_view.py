"""
Folium map construction.

Renders walked (green) and unwalked (gray) street segments with
optional GPX track overlays and street name tooltips.
"""

from __future__ import annotations

import folium
import geopandas as gpd

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

    # ── Unwalked streets layer ────────────────────────────────────────────────
    unwalked = gdf[~gdf["walked"]]
    if not unwalked.empty:
        folium.GeoJson(
            unwalked.__geo_interface__,
            name="Unwalked Streets",
            style_function=lambda _: {
                "color": COLOR_UNWALKED,
                "weight": 2,
                "opacity": 0.6,
            },
        ).add_to(m)

    # ── Walked streets layer ──────────────────────────────────────────────────
    walked = gdf[gdf["walked"]]
    if not walked.empty:
        tooltip_fields = [c for c in ["name", "highway", "walk_date", "walk_count", "length_m"] if c in walked.columns]
        tooltip_aliases = {
            "name": "Street",
            "highway": "Type",
            "walk_date": "First walked",
            "walk_count": "Times walked",
            "length_m": "Length (m)",
        }

        folium.GeoJson(
            walked.__geo_interface__,
            name="Walked Streets",
            style_function=lambda _: {
                "color": COLOR_WALKED,
                "weight": 4,
                "opacity": 0.9,
            },
            tooltip=folium.GeoJsonTooltip(
                fields=tooltip_fields,
                aliases=[tooltip_aliases.get(f, f) for f in tooltip_fields],
                localize=True,
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
