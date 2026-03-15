"""Batch block assignment panel."""

from __future__ import annotations

import geopandas as gpd
import streamlit as st

import src.database as db


def render_batch_panel(
    selected_block_ids: set[str],
    network_gdf: gpd.GeoDataFrame,
    activities: list[dict],
) -> None:
    """Render the batch assignment panel in the right column."""

    col_title, col_clear = st.columns([3, 1])
    with col_title:
        n = len(selected_block_ids)
        st.markdown(f"**Batch Assign** — {n} block{'s' if n != 1 else ''}")
    with col_clear:
        if st.button("Clear", key="batch_clear", disabled=n == 0):
            st.session_state["batch_selected_blocks"] = set()
            st.rerun()

    if not selected_block_ids:
        st.caption("Click blocks on the map to select them.")
    else:
        _render_selected_list(selected_block_ids, network_gdf)

    st.divider()

    if not activities:
        st.caption("No activities available.")
        return

    options = {_activity_label(a): a["id"] for a in activities}
    chosen_label = st.selectbox(
        "Assign to activity",
        options=list(options.keys()),
        key="batch_activity_select",
        disabled=not selected_block_ids,
    )

    if st.button(
        f"Assign {len(selected_block_ids)} block{'s' if len(selected_block_ids) != 1 else ''}",
        key="batch_assign",
        type="primary",
        disabled=not selected_block_ids,
    ):
        activity_id = options[chosen_label]
        for bid in selected_block_ids:
            db.insert_manual_block_walk(bid, activity_id)
        st.session_state["batch_selected_blocks"] = set()
        st.rerun(scope="app")  # walked_blocks changed → full rerun


# ── Helpers ───────────────────────────────────────────────────────────────────

def _render_selected_list(selected_block_ids: set[str], network_gdf: gpd.GeoDataFrame) -> None:
    for bid in sorted(selected_block_ids):
        block_edges = network_gdf[network_gdf["block_id"] == bid]
        if not block_edges.empty and "name" in block_edges.columns:
            names = block_edges["name"].dropna().unique()
            street_name = names[0] if len(names) > 0 else "—"
        else:
            street_name = "—"

        col_name, col_btn = st.columns([4, 1])
        with col_name:
            st.markdown(f"☑ **{street_name}**  \n{bid}")
        with col_btn:
            if st.button("✕", key=f"batch_remove_{bid}", help="Remove from selection"):
                st.session_state["batch_selected_blocks"].discard(bid)
                st.rerun()


def _activity_label(a: dict) -> str:
    date_str = (a.get("start_time") or "")[:10]
    name = a.get("track_name") or a.get("filename") or f"#{a['id']}"
    atype = a.get("activity_type", "walk")
    return f"{date_str} — {name} ({atype})"
