"""
Walk Every Street — Streamlit app entry point.

Data flow:
  1. Sidebar collects area + GPX files + settings.
  2. Street network loaded from disk cache or fetched fresh from OSM.
  3. GPX files parsed into track geometries.
  4. Tracks matched to street segments via spatial index.
  5. Results merged into persisted WalkState and saved to disk.
  6. Folium map and stats panel rendered.
"""

from __future__ import annotations

import logging

import geopandas as gpd
import streamlit as st
from streamlit_folium import st_folium

from src.gpx_parser import bbox_to_osmnx, combined_bbox, parse_multiple_gpx
from src.network import _cache_key, cache_exists, get_or_fetch_network
from src.track_matcher import match_tracks_to_segments
from src.walk_tracker import WalkState
from ui.map_view import build_map
from ui.sidebar import render_sidebar
from ui.stats_panel import render_stats

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

st.set_page_config(
    page_title="Walk Every Street",
    page_icon="🗺️",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ── Cached network loader ─────────────────────────────────────────────────────

@st.cache_data(ttl=86400 * 30, show_spinner=False)
def _load_network(place_name: str | None, bbox: tuple | None, force_refresh: bool):
    """Streamlit-cached wrapper around get_or_fetch_network."""
    gdf, key = get_or_fetch_network(
        place_name=place_name,
        bbox=bbox,
        force_refresh=force_refresh,
    )
    return gdf, key


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    cfg = render_sidebar()

    # ── Landing state (nothing processed yet) ────────────────────────────────
    if not cfg["process"] and "network_gdf" not in st.session_state:
        st.markdown("## 🗺️ Walk Every Street")
        st.markdown(
            "Track your progress walking every street in your city. "
            "Configure the area and upload your GPX files in the sidebar, "
            "then click **▶ Process Walks**."
        )
        st.info(
            "**Getting started:**\n"
            "1. Enter a city name (or bounding box) in the sidebar.\n"
            "2. Upload one or more `.gpx` files exported from Garmin Connect or Strava.\n"
            "3. Click **▶ Process Walks**.\n\n"
            "Your progress is saved between sessions — you can upload new files "
            "at any time to accumulate coverage.",
            icon="👟",
        )
        return

    tracks: list[dict] = []
    network_gdf = None
    cache_key = None

    with st.status("Processing...", expanded=True) as status:

        # ── Step 1: Parse GPX files ───────────────────────────────────────────
        if cfg["uploaded_files"]:
            n = len(cfg["uploaded_files"])
            st.write(f"📂 Parsing {n} GPX file(s)...")
            tracks = parse_multiple_gpx(cfg["uploaded_files"])
            if not tracks:
                status.update(label="No valid GPX data found.", state="error")
                st.warning("No valid GPX data found in the uploaded files.")
                return
            st.write(f"✅ Parsed {len(tracks)} track(s).")

        # ── Step 2: Determine area ────────────────────────────────────────────
        place_name = cfg["place_name"]
        bbox = cfg["bbox"]

        if cfg["area_method"] == "Auto from GPX":
            if not tracks:
                status.update(label="No GPX files uploaded.", state="error")
                st.error("Please upload at least one GPX file to use Auto from GPX mode.")
                return
            raw_bbox = combined_bbox(tracks)
            bbox = bbox_to_osmnx(raw_bbox)
            place_name = None

        if not place_name and not bbox:
            status.update(label="No area specified.", state="error")
            st.error("Please specify an area in the sidebar.")
            return

        # ── Step 3: Load network ──────────────────────────────────────────────
        area_label = place_name or "bounding box"
        ck = _cache_key(place_name, bbox)
        if cache_exists(ck) and not cfg["force_refresh"]:
            st.write(f"🌍 Loading street network for {area_label} from cache...")
        else:
            st.write(
                f"🌍 Fetching street network for {area_label} from OSM "
                f"and building blocks — **first run only, takes 1–3 min**..."
            )
        try:
            network_gdf, cache_key = _load_network(place_name, bbox, cfg["force_refresh"])
        except Exception as exc:
            status.update(label="Failed to load street network.", state="error")
            st.error(f"Failed to load street network: {exc}")
            logger.exception("Network load failed")
            return
        st.write(f"✅ Network loaded ({len(network_gdf):,} blocks).")

        st.session_state["network_gdf"] = network_gdf
        st.session_state["cache_key"] = cache_key

        # ── Step 4: Load or reset persisted walk state ────────────────────────
        st.write("💾 Loading saved walk history...")
        walk_state = WalkState.load(cache_key)
        if cfg["reset_walks"]:
            walk_state.reset()
            st.write("🗑️ Walk history cleared.")
        else:
            st.write(f"✅ Walk history loaded ({len(walk_state.walked_block_ids):,} blocks previously walked).")

        # ── Step 5: Match new tracks to segments ──────────────────────────────
        if tracks:
            st.write(f"📍 Matching {len(tracks)} track(s) to street segments...")
            matched_gdf = match_tracks_to_segments(
                network_gdf, tracks,
                snap_distance_m=cfg["snap_distance"],
                coverage_threshold=cfg["coverage_threshold"],
            )
            newly_walked = int(matched_gdf["walked"].sum())
            walk_state.merge_from_gdf(matched_gdf)
            walk_state.save()
            st.write(f"✅ Matched {newly_walked:,} segments from new uploads.")

        # ── Step 6: Apply persisted state to network ──────────────────────────
        st.write("🗺️ Building map...")
        display_gdf = walk_state.apply_to_network(network_gdf)
        summary = walk_state.summary(display_gdf)

        status.update(
            label=f"Done — {summary['pct_walked']:.1f}% of streets walked.",
            state="complete",
            expanded=False,
        )

    render_stats(summary, tracks=tracks if tracks else None)

    st.divider()

    col_map, col_info = st.columns([4, 1])

    with col_info:
        st.markdown("**Legend**")
        st.markdown(f"🟢 Walked &nbsp;&nbsp; 🔘 Unwalked")
        if tracks:
            st.markdown("🟠 GPS track")
        st.caption(
            f"Network: {summary['total_segments']:,} segments · "
            f"{summary['total_km']:.0f} km"
        )
        if not tracks:
            st.caption("Upload GPX files and re-process to update your coverage.")

    with col_map:
        folium_map = build_map(display_gdf, tracks=tracks if tracks else None)
        st_folium(folium_map, use_container_width=True, height=680, returned_objects=[])


if __name__ == "__main__":
    main()
