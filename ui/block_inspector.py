"""
Block Inspector panel.

Renders in the right column when a block is selected:
  - Header with block name / length and close button
  - List of existing activity assignments with Remove buttons
  - Add-assignment selectbox + Assign button
  - Mini verification map showing the block + its GPS tracks
"""

from __future__ import annotations

import re
from pathlib import Path

import geopandas as gpd
import streamlit as st
from streamlit_folium import st_folium

import src.database as db
from src.gpx_parser import parse_gpx_file
from src.paths import PROCESSED_DIR
from ui.map_view import build_mini_map


def render_block_inspector(
    block_id: str,
    network_gdf: gpd.GeoDataFrame,
    activities: list[dict],
    walked_blocks: dict[str, dict],
    coverage_threshold: float = 0.80,
) -> None:
    """Render the block inspector in the right column."""

    # ── Header ────────────────────────────────────────────────────────────────
    block_edges = network_gdf[network_gdf["block_id"] == block_id]
    block_name = ""
    block_len = None
    if not block_edges.empty:
        names = block_edges["name"].dropna().unique() if "name" in block_edges.columns else []
        block_name = names[0] if len(names) > 0 else block_id
        if "block_length_m" in block_edges.columns:
            block_len = block_edges["block_length_m"].iloc[0]

    col_title, col_close = st.columns([5, 1])
    with col_title:
        st.markdown("**Block Inspector**")
        st.caption(f"{block_name}  ·  {block_id}" + (f"  ·  {block_len:.0f} m" if block_len else ""))
    with col_close:
        if st.button("✕", key="inspector_close", help="Close inspector"):
            st.session_state["selected_block_id"] = None
            st.rerun()  # fragment-only: UI change only

    # ── Existing assignments ──────────────────────────────────────────────────
    all_assignments = db.get_block_assignments(block_id)
    assignments = [a for a in all_assignments if a["coverage_ratio"] >= coverage_threshold]
    partial_matches = sorted(
        [a for a in all_assignments if a["coverage_ratio"] < coverage_threshold],
        key=lambda a: a["coverage_ratio"],
        reverse=True,
    )

    if assignments:
        st.markdown("**Assignments**")
        for asgn in assignments:
            _render_assignment_row(block_id, asgn)
    elif partial_matches:
        st.markdown("**Partial matches**")
        st.caption(f"Below {coverage_threshold * 100:.0f}% threshold — click Assign to confirm.")
        for asgn in partial_matches:
            _render_partial_match_row(block_id, asgn)
    else:
        st.caption("No activities matched this block.")

    st.divider()

    # ── Add assignment ────────────────────────────────────────────────────────
    if activities:
        st.markdown("**Add assignment**")
        options = {
            _activity_label(a): a["id"] for a in activities
        }
        chosen_label = st.selectbox(
            "Activity",
            options=list(options.keys()),
            key=f"inspector_add_{block_id}",
            label_visibility="collapsed",
        )
        if st.button("Assign", key=f"inspector_assign_{block_id}"):
            db.insert_manual_block_walk(block_id, options[chosen_label])
            st.rerun(scope="app")  # walked_blocks changed → full rerun
    else:
        st.caption("No activities available.")

    st.divider()

    # ── Mini verification map ─────────────────────────────────────────────────
    st.markdown("**Verification map**")
    _mini_key = (block_id, id(network_gdf), tuple(a["activity_id"] for a in assignments))
    if st.session_state.get("_mini_map_key") != _mini_key:
        tracks = _load_tracks_for_assignments(assignments)
        st.session_state["_mini_map"] = build_mini_map(block_id, network_gdf, tracks=tracks or None)
        st.session_state["_mini_map_key"] = _mini_key
    mini_map = st.session_state["_mini_map"]
    st_folium(mini_map, use_container_width=True, height=300, returned_objects=[], key=f"mini_map_{block_id}")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _render_assignment_row(block_id: str, asgn: dict) -> None:
    date_str = (asgn["start_time"] or "")[:10]
    name = asgn.get("track_name") or asgn.get("filename") or f"Activity {asgn['activity_id']}"
    source = asgn.get("source", "auto")
    coverage = asgn.get("coverage_ratio", 0.0)
    badge = "🤖 auto" if source == "auto" else "✋ manual"

    col_info, col_btn = st.columns([4, 1])
    with col_info:
        st.markdown(
            f"**{name}**  \n"
            f"{date_str}  ·  {badge}  ·  {coverage * 100:.0f}% coverage"
        )
    with col_btn:
        if st.button("Remove", key=f"remove_{block_id}_{asgn['activity_id']}"):
            db.remove_block_walk(block_id, asgn["activity_id"], source)
            st.rerun(scope="app")  # walked_blocks changed → full rerun


def _render_partial_match_row(block_id: str, asgn: dict) -> None:
    date_str = (asgn["start_time"] or "")[:10]
    name = asgn.get("track_name") or asgn.get("filename") or f"Activity {asgn['activity_id']}"
    coverage = asgn.get("coverage_ratio", 0.0)

    col_info, col_btn = st.columns([4, 1])
    with col_info:
        st.markdown(
            f"**{name}**  \n"
            f"{date_str}  ·  {coverage * 100:.0f}% coverage"
        )
    with col_btn:
        if st.button("Assign", key=f"partial_{block_id}_{asgn['activity_id']}"):
            db.insert_manual_block_walk(block_id, asgn["activity_id"])
            st.rerun(scope="app")  # walked_blocks changed → full rerun


def _activity_label(a: dict) -> str:
    date_str = (a.get("start_time") or "")[:10]
    name = a.get("track_name") or a.get("filename") or f"#{a['id']}"
    atype = a.get("activity_type", "walk")
    return f"{date_str} — {name} ({atype})"


def _load_tracks_for_assignments(assignments: list[dict]) -> list[dict]:
    """Parse GPX files for each assignment; silently skip missing files."""
    tracks: list[dict] = []
    for asgn in assignments:
        filename = asgn.get("filename")
        if not filename:
            continue
        gpx_path = _resolve_gpx_path(PROCESSED_DIR, filename)
        if gpx_path is None:
            st.caption(f"GPX file not found: {filename}")
            continue
        try:
            raw = gpx_path.read_bytes()
            parsed = parse_gpx_file(raw, filename=gpx_path.name)
            if parsed:
                tracks.append(parsed)
        except Exception:
            st.caption(f"Could not read GPX: {filename}")
    return tracks


def _resolve_gpx_path(processed_dir: Path, filename: str) -> Path | None:
    """Find a GPX file by name, handling numeric collision suffixes."""
    direct = processed_dir / filename
    if direct.exists():
        return direct
    stem = Path(filename).stem
    suffix = Path(filename).suffix
    pattern = re.compile(rf"^{re.escape(stem)}(_\d+)?{re.escape(suffix)}$")
    candidates = sorted(processed_dir.glob(f"{stem}*{suffix}"))
    matches = [p for p in candidates if pattern.match(p.name)]
    return matches[0] if matches else None
