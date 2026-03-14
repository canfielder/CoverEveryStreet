"""
Walk Every Street — Editor app.

Full-featured local tool for processing GPX activities and assigning them
to street blocks. Not for public deployment.

Data flow:
  1. Init database (creates tables if needed).
  2. Sidebar: area selection + settings.
  3. Load street network from disk cache or fetch fresh from OSM.
  4. Scan GPX inbox for new files (cached in session_state — only re-scans on
     page load or after processing, not on every map click).
  5. After processing, fetch walked blocks from DB (date-filtered).
  6. Render single-color Folium map + Activities tab.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

# Ensure project root is on sys.path so src.*, ui.*, config resolve correctly
sys.path.insert(0, str(Path(__file__).parent.parent))
from datetime import datetime, timezone

import geopandas as gpd
import shapely
import streamlit as st
from streamlit_folium import st_folium

import src.database as db
from src.inbox import scan_inbox
from src.network import _cache_age_days, get_or_fetch_network
from ui.block_inspector import render_block_inspector
from ui.map_view import build_map
from ui.pending_panel import render_pending_panel
from ui.sidebar import render_sidebar
from ui.stats_panel import render_activity_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

st.set_page_config(
    page_title="Walk Every Street — Editor",
    page_icon="✏️",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ── Cached network loader ──────────────────────────────────────────────────────

@st.cache_data(show_spinner=False)
def _load_network(place_name: str | None, bbox: tuple | None, force_refresh: bool):
    gdf, key = get_or_fetch_network(
        place_name=place_name,
        bbox=bbox,
        force_refresh=force_refresh,
    )
    return gdf, key


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    db.init_db()

    cache_info = st.session_state.get("network_cache_info")
    cfg = render_sidebar(network_cache_info=cache_info)

    place_name = cfg["place_name"]
    bbox = cfg["bbox"]

    if not place_name and not bbox:
        st.error("Please specify an area in the sidebar.")
        return

    # ── Load network ───────────────────────────────────────────────────────────
    area_label = place_name or "bounding box"
    with st.spinner(f"Loading street network for {area_label}…"):
        try:
            network_gdf, cache_key = _load_network(place_name, bbox, cfg["force_refresh"])
        except Exception as exc:
            st.error(f"Failed to load street network: {exc}")
            logger.exception("Network load failed")
            return

    age = _cache_age_days(cache_key)
    n_unique_blocks = network_gdf["block_id"].nunique() if "block_id" in network_gdf.columns else 0
    st.session_state["network_cache_info"] = {
        "age_days": age,
        "block_count": n_unique_blocks,
    }

    # ── Inbox scan (cached — skip on map click reruns) ─────────────────────────
    if "pending" not in st.session_state:
        st.session_state["pending"] = scan_inbox()
    pending = st.session_state["pending"]

    if pending:
        render_pending_panel(
            pending=pending,
            network_gdf=network_gdf,
            snap_distance_m=cfg["snap_distance"],
        )
        st.divider()

    # ── Fetch walked blocks from DB ────────────────────────────────────────────
    start_dt = (
        datetime(cfg["start_date"].year, cfg["start_date"].month, cfg["start_date"].day,
                 tzinfo=timezone.utc)
        if cfg.get("start_date") else None
    )
    end_dt = (
        datetime(cfg["end_date"].year, cfg["end_date"].month, cfg["end_date"].day,
                 23, 59, 59, tzinfo=timezone.utc)
        if cfg.get("end_date") else None
    )

    walked_blocks = db.get_walked_blocks(
        start_date=start_dt,
        end_date=end_dt,
        coverage_threshold=cfg["coverage_threshold"],
    )
    activities = db.get_all_activities()

    # ── Tabs: Map / Activities ─────────────────────────────────────────────────
    tab_map, tab_activities = st.tabs(["Map", "Activities"])

    with tab_map:
        _map_key = (id(network_gdf), frozenset(walked_blocks.keys()))
        if st.session_state.get("_map_cache_key") != _map_key:
            st.session_state["_folium_map"] = build_map(
                network_gdf, walked_blocks=walked_blocks, color_by_activity=False
            )
            st.session_state["_map_cache_key"] = _map_key
        folium_map = st.session_state["_folium_map"]

        col_map, col_right = st.columns([3, 2])

        with col_map:
            map_data = st_folium(
                folium_map,
                use_container_width=True,
                height=680,
                returned_objects=["last_object_clicked"],
                key="main_map",
            )

        clicked = map_data.get("last_object_clicked") if map_data else None
        if clicked and isinstance(clicked, dict):
            lat, lng = clicked.get("lat"), clicked.get("lng")
            if lat is not None and (lat, lng) != st.session_state.get("_last_click"):
                st.session_state["_last_click"] = (lat, lng)
                found_block = _find_nearest_block(lat, lng, network_gdf)
                if found_block:
                    st.session_state["selected_block_id"] = found_block

        selected_block_id = st.session_state.get("selected_block_id")
        with col_right:
            if selected_block_id:
                render_block_inspector(
                    selected_block_id, network_gdf, activities, walked_blocks,
                    coverage_threshold=cfg["coverage_threshold"],
                )
            else:
                st.caption("Click a block on the map to inspect it.")

    with tab_activities:
        render_activity_table(activities)


def _get_block_tree(network_gdf: gpd.GeoDataFrame):
    from config import METRIC_CRS

    tree_id = id(network_gdf)
    if st.session_state.get("_tree_net_id") != tree_id:
        net_proj = network_gdf.to_crs(METRIC_CRS)
        st.session_state["_block_tree"] = shapely.STRtree(net_proj.geometry.values)
        st.session_state["_block_tree_net"] = net_proj
        st.session_state["_tree_net_id"] = tree_id
    return st.session_state["_block_tree"], st.session_state["_block_tree_net"]


def _find_nearest_block(lat: float, lng: float, network_gdf: gpd.GeoDataFrame) -> str | None:
    from shapely.geometry import Point

    from config import METRIC_CRS, WGS84_CRS

    click_pt = gpd.GeoSeries([Point(lng, lat)], crs=WGS84_CRS).to_crs(METRIC_CRS).iloc[0]
    tree, net_proj = _get_block_tree(network_gdf)
    idx = tree.nearest(click_pt)
    if net_proj.geometry.iloc[idx].distance(click_pt) > 150:
        return None
    return network_gdf.iloc[idx]["block_id"]


if __name__ == "__main__":
    main()
