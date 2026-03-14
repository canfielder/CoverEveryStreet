"""Statistics panel — coverage metrics and activity history."""

from __future__ import annotations

import pandas as pd
import streamlit as st


def render_metrics(
    network_gdf,
    walked_blocks: dict[str, dict],
    activities: list[dict],
) -> None:
    """Headline metrics + progress bar only."""
    blocks = (
        network_gdf[["block_id", "block_length_m"]]
        .drop_duplicates("block_id")
        .copy()
    )
    total_blocks = len(blocks)
    total_km = blocks["block_length_m"].sum() / 1000

    walked_ids = set(walked_blocks.keys())
    walked_mask = blocks["block_id"].isin(walked_ids)
    n_walked = walked_mask.sum()
    walked_km = blocks.loc[walked_mask, "block_length_m"].sum() / 1000
    unwalked_km = total_km - walked_km
    pct = (n_walked / total_blocks * 100) if total_blocks else 0.0

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Blocks Walked", f"{n_walked:,}", f"{pct:.1f}% of {total_blocks:,}")
    c2.metric("Distance Walked", f"{walked_km:,.1f} km", f"{unwalked_km:,.1f} km remaining")
    c3.metric("Total Activities", f"{len(activities):,}")
    c4.metric("Network Size", f"{total_km:,.1f} km")

    st.progress(min(pct / 100, 1.0), text=f"**{pct:.1f}%** of blocks walked")


def render_detail_tables(
    network_gdf,
    walked_blocks: dict[str, dict],
    activities: list[dict],
) -> None:
    """By category breakdown + Activity history tables."""
    blocks = (
        network_gdf[["block_id", "block_length_m"]]
        .drop_duplicates("block_id")
        .copy()
    )

    # ── Category breakdown ────────────────────────────────────────────────────
    if walked_blocks:
        st.subheader("By category")
        category_counts: dict[str, int] = {}
        category_km: dict[str, float] = {}

        for bid, info in walked_blocks.items():
            atype = info.get("activity_type", "walk")
            companions = frozenset(
                c.strip().lower() for c in info.get("companions", [])
            )
            key = (atype, companions)
            label = _category_label(key)
            category_counts[label] = category_counts.get(label, 0) + 1
            block_km = blocks.loc[blocks["block_id"] == bid, "block_length_m"].sum() / 1000
            category_km[label] = category_km.get(label, 0.0) + block_km

        cat_rows = [
            {"Category": lbl, "Blocks": cnt, "Distance (km)": f"{category_km.get(lbl, 0):.1f}"}
            for lbl, cnt in sorted(category_counts.items(), key=lambda x: -x[1])
        ]
        st.dataframe(pd.DataFrame(cat_rows), use_container_width=True, hide_index=True)

    # ── Activity history ──────────────────────────────────────────────────────
    render_activity_table(activities)


def render_activity_table(activities: list[dict]) -> None:
    """Activity history table — used by the editor Activities tab."""
    if activities:
        st.subheader("Activity history")
        rows = []
        for a in activities:
            companions = a.get("companions") or []
            rows.append({
                "Date": a.get("start_time", "")[:10],
                "Name": a.get("track_name") or a.get("filename") or "—",
                "Type": (a.get("activity_type") or "—").capitalize(),
                "Companions": ", ".join(companions) if companions else "Solo",
                "File": a.get("filename", "—"),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)


def render_stats(
    network_gdf,
    walked_blocks: dict[str, dict],
    activities: list[dict],
) -> None:
    """Backwards-compat wrapper — renders metrics + detail tables."""
    render_metrics(network_gdf, walked_blocks, activities)
    render_detail_tables(network_gdf, walked_blocks, activities)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _category_label(key: tuple) -> str:
    atype, companions = key
    label = atype.capitalize()
    if companions:
        label += " + " + ", ".join(sorted(c.capitalize() for c in companions))
    return label
