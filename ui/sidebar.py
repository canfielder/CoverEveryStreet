"""Streamlit sidebar — area selection, GPX upload, and settings."""

from __future__ import annotations

import streamlit as st

from config import DEFAULT_COVERAGE_THRESHOLD, DEFAULT_SNAP_DISTANCE_M


def render_sidebar() -> dict:
    """
    Render all sidebar controls and return a config dict.

    Returns a dict consumed by app.py with keys:
        area_method     : "Place Name" | "Bounding Box" | "Auto from GPX"
        place_name      : str | None
        bbox            : (north, south, east, west) | None  — osmnx format
        uploaded_files  : list of UploadedFile
        snap_distance   : float (metres)
        force_refresh   : bool
        process         : bool  — True when the Process button is clicked
        reset_walks     : bool
    """
    with st.sidebar:
        st.title("🗺️ Walk Every Street")
        st.caption("Upload your GPX files and see how much of your city you've covered.")

        st.divider()

        # ── 1. Area selection ─────────────────────────────────────────────────
        st.subheader("1. Select Area")
        area_method = st.radio(
            "How to define the area",
            ["Place Name", "Bounding Box", "Auto from GPX"],
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

        elif area_method == "Bounding Box":
            st.caption("Enter coordinates in decimal degrees (WGS84).")
            col_n, col_s = st.columns(2)
            col_e, col_w = st.columns(2)
            north = col_n.number_input("North", value=35.35, format="%.4f", step=0.01)
            south = col_s.number_input("South", value=35.10, format="%.4f", step=0.01)
            east = col_e.number_input("East", value=-80.65, format="%.4f", step=0.01)
            west = col_w.number_input("West", value=-80.95, format="%.4f", step=0.01)
            bbox = (north, south, east, west)

        else:  # Auto from GPX
            st.info("The area will be detected from your uploaded GPX files.", icon="📍")

        st.divider()

        # ── 2. GPX upload ─────────────────────────────────────────────────────
        st.subheader("2. Upload GPX Files")
        uploaded_files = st.file_uploader(
            "Drag and drop or click to browse",
            type=["gpx"],
            accept_multiple_files=True,
            label_visibility="collapsed",
        )
        if uploaded_files:
            st.caption(f"📂 {len(uploaded_files)} file(s) selected")

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
                "will be marked as walked. Increase if segments are being missed; "
                "decrease to reduce false positives."
            ),
        )

        coverage_threshold = st.slider(
            "Segment coverage required (%)",
            min_value=50,
            max_value=100,
            value=int(DEFAULT_COVERAGE_THRESHOLD * 100),
            step=5,
            help=(
                "How much of a street segment must fall within the snap distance "
                "to count as fully walked. 95% means you need to have covered "
                "nearly the whole block — lower this if segments are being missed."
            ),
        )

        with st.expander("Advanced"):
            force_refresh = st.checkbox(
                "Re-download street network from OSM",
                value=False,
                help="Clears the local cache and fetches a fresh copy of the street network.",
            )
            reset_walks = st.checkbox(
                "Reset all walk history for this area",
                value=False,
                help="Permanently clears your saved progress for the selected area.",
            )

        st.divider()

        # ── 4. Process button ─────────────────────────────────────────────────
        process = st.button(
            "▶ Process Walks",
            type="primary",
            use_container_width=True,
            disabled=(
                area_method == "Auto from GPX" and not uploaded_files
            ),
        )

    return {
        "area_method": area_method,
        "place_name": place_name,
        "bbox": bbox,
        "uploaded_files": uploaded_files or [],
        "snap_distance": snap_distance,
        "coverage_threshold": coverage_threshold / 100,
        "force_refresh": force_refresh,
        "reset_walks": reset_walks,
        "process": process,
    }
