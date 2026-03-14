"""Streamlit sidebar — area selection and settings."""

from __future__ import annotations

from datetime import date, timedelta

import streamlit as st

from config import DEFAULT_COVERAGE_THRESHOLD, DEFAULT_SNAP_DISTANCE_M


def render_sidebar(network_cache_info: dict | None = None) -> dict:
    """
    Render all sidebar controls and return a config dict.

    Args:
        network_cache_info: Optional dict with keys 'age_days' (float) and
                            'block_count' (int) describing the loaded network.

    Returns dict with keys:
        area_method         : "Place Name" | "Bounding Box"
        place_name          : str | None
        bbox                : (north, south, east, west) | None
        snap_distance       : float (metres)
        coverage_threshold  : float (0–1)
        force_refresh       : bool  — re-download network from OSM
        start_date          : date | None
        end_date            : date | None
    """
    with st.sidebar:
        st.title("🗺️ Walk Every Street")
        st.caption("Drop GPX files in the inbox folder, then reload the app to process them.")

        st.divider()

        # ── 1. Area selection ─────────────────────────────────────────────────
        st.subheader("1. Select Area")
        area_method = st.radio(
            "How to define the area",
            ["Place Name", "Bounding Box"],
            index=0,
            label_visibility="collapsed",
        )

        place_name: str | None = None
        bbox: tuple | None = None

        if area_method == "Place Name":
            place_name = st.text_input(
                "City or place",
                value="Charlotte, North Carolina",
                placeholder="e.g. Asheville, North Carolina",
            )
        else:
            st.caption("Enter coordinates in decimal degrees (WGS84).")
            col_n, col_s = st.columns(2)
            col_e, col_w = st.columns(2)
            north = col_n.number_input("North", value=35.35, format="%.4f", step=0.01)
            south = col_s.number_input("South", value=35.10, format="%.4f", step=0.01)
            east  = col_e.number_input("East",  value=-80.65, format="%.4f", step=0.01)
            west  = col_w.number_input("West",  value=-80.95, format="%.4f", step=0.01)
            bbox = (north, south, east, west)

        # Network cache status
        if network_cache_info:
            age = network_cache_info.get("age_days")
            n_blocks = network_cache_info.get("block_count", 0)
            age_str = f"{age:.0f}d old" if age is not None else "unknown age"
            st.caption(f"Loaded network: {n_blocks:,} blocks · {age_str}")

        st.divider()

        # ── 2. Date filter ────────────────────────────────────────────────────
        st.subheader("2. Date Filter")
        st.caption("Show only activities within this date range.")

        use_date_filter = st.checkbox("Enable date filter", value=False)
        start_date: date | None = None
        end_date: date | None = None

        if use_date_filter:
            today = date.today()
            col_s, col_e = st.columns(2)
            start_date = col_s.date_input("From", value=today - timedelta(days=365))
            end_date   = col_e.date_input("To",   value=today)

        st.divider()

        # ── 3. Settings ───────────────────────────────────────────────────────
        st.subheader("3. Settings")
        snap_distance = st.slider(
            "GPS snap distance (metres)",
            min_value=5,
            max_value=50,
            value=DEFAULT_SNAP_DISTANCE_M,
            step=5,
            help=(
                "Street segments within this distance of your GPS track "
                "will be marked as walked."
            ),
        )

        coverage_threshold = st.slider(
            "Coverage required (%)",
            min_value=50,
            max_value=100,
            value=int(DEFAULT_COVERAGE_THRESHOLD * 100),
            step=5,
            help="Fraction of a block that must be covered to count as walked.",
        )

        with st.expander("Advanced"):
            force_refresh = st.checkbox(
                "Re-download street network from OSM",
                value=False,
                help="Clears the local cache and fetches a fresh copy of the street network.",
            )

    return {
        "area_method":          area_method,
        "place_name":           place_name,
        "bbox":                 bbox,
        "snap_distance":        snap_distance,
        "coverage_threshold":   coverage_threshold / 100,
        "force_refresh":        force_refresh,
        "start_date":           start_date,
        "end_date":             end_date,
    }
