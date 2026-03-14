"""
Walk Every Street — read-only viewer for Streamlit Cloud.

Reads a pre-published GeoPackage snapshot from data/published/map_data.gpkg
and renders the coverage map. No database, no file uploads, no processing.

The published snapshot is committed to git by the local app after each
processing run (see scripts/publish_snapshot.py).
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

# Ensure project root is on sys.path so src.*, config resolve correctly
sys.path.insert(0, str(Path(__file__).parent.parent))

import folium
import geopandas as gpd
import streamlit as st
from streamlit_folium import st_folium

from config import (
    ACTIVITY_COLORS,
    COLOR_FALLBACK,
    COLOR_UNWALKED,
    MAP_TILES,
    MAP_ZOOM_START,
    WGS84_CRS,
)
from src.paths import PUBLISHED_DIR

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

SNAPSHOT_PATH = PUBLISHED_DIR / "map_data.gpkg"
META_PATH = PUBLISHED_DIR / "meta.json"

st.set_page_config(
    page_title="Walk Every Street — Map",
    page_icon="🗺️",
    layout="wide",
    initial_sidebar_state="collapsed",
)

_HIGHLIGHT = {"weight": 6, "color": "#facc15", "opacity": 1.0}


@st.cache_data(show_spinner=False)
def _load_snapshot():
    if not SNAPSHOT_PATH.exists():
        return None, None
    gdf = gpd.read_file(SNAPSHOT_PATH)
    meta = json.loads(META_PATH.read_text()) if META_PATH.exists() else {}
    return gdf, meta


def main() -> None:
    st.title("🗺️ Walk Every Street")

    gdf, meta = _load_snapshot()

    if gdf is None:
        st.warning(
            "No published snapshot found. Run the local app and publish a snapshot first.",
            icon="⚠️",
        )
        return

    published_at = meta.get("published_at", "unknown")
    n_walked = int(gdf["_walked"].sum()) if "_walked" in gdf.columns else 0
    total = len(gdf)
    pct = n_walked / total * 100 if total else 0.0

    col1, col2, col3 = st.columns(3)
    col1.metric("Blocks Walked", f"{n_walked:,}", f"{pct:.1f}%")
    col2.metric("Total Blocks", f"{total:,}")
    col3.caption(f"Published: {published_at}")
    st.progress(min(pct / 100, 1.0))

    st.divider()

    gdf_wgs = gdf.to_crs(WGS84_CRS) if gdf.crs and gdf.crs.to_epsg() != 4326 else gdf
    centroid = gdf_wgs.geometry.unary_union.centroid
    m = folium.Map(
        location=[centroid.y, centroid.x],
        zoom_start=MAP_ZOOM_START,
        tiles=MAP_TILES,
    )

    # Round display columns
    for col in ("length_m", "block_length_m"):
        if col in gdf_wgs.columns:
            gdf_wgs = gdf_wgs.copy()
            gdf_wgs[col] = gdf_wgs[col].round(1)

    # Unwalked
    unwalked = gdf_wgs[~gdf_wgs["_walked"]] if "_walked" in gdf_wgs.columns else gdf_wgs
    if not unwalked.empty:
        _fields = [c for c in ["block_id", "name", "block_length_m"] if c in unwalked.columns]
        folium.GeoJson(
            unwalked.__geo_interface__,
            name="Unwalked Streets",
            style_function=lambda _: {"color": COLOR_UNWALKED, "weight": 2, "opacity": 0.6},
            highlight_function=lambda _: _HIGHLIGHT,
            tooltip=folium.GeoJsonTooltip(fields=_fields),
        ).add_to(m)

    # Walked — color by category
    walked = gdf_wgs[gdf_wgs["_walked"]] if "_walked" in gdf_wgs.columns else gpd.GeoDataFrame()
    if not walked.empty:
        walked = walked.copy()
        walked["_color"] = walked.apply(_edge_color, axis=1)
        for color, group in walked.groupby("_color"):
            _all = ["block_id", "name", "block_length_m", "_walk_date", "_activity_type"]
            _fields = [c for c in _all if c in group.columns]
            folium.GeoJson(
                group.__geo_interface__,
                name=_color_label(color),
                style_function=lambda _, c=color: {"color": c, "weight": 4, "opacity": 0.9},
                highlight_function=lambda _: _HIGHLIGHT,
                tooltip=folium.GeoJsonTooltip(fields=_fields),
            ).add_to(m)

    folium.LayerControl(collapsed=False).add_to(m)
    st_folium(m, use_container_width=True, height=700, returned_objects=[])


def _edge_color(row) -> str:
    atype = row.get("_activity_type") or "walk"
    companions_str = row.get("_companions") or ""
    companions = frozenset(c.strip().lower() for c in companions_str.split(",") if c.strip())
    key = (atype, companions)
    return ACTIVITY_COLORS.get(key, ACTIVITY_COLORS.get((atype, frozenset()), COLOR_FALLBACK))


def _color_label(color: str) -> str:
    reverse = {v: k for k, v in ACTIVITY_COLORS.items()}
    if color not in reverse:
        return "Walked"
    atype, companions = reverse[color]
    label = atype.capitalize()
    if companions:
        label += " + " + ", ".join(sorted(c.capitalize() for c in companions))
    return label


if __name__ == "__main__":
    main()
