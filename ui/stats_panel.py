"""Statistics panel — coverage metrics and walk history."""

from __future__ import annotations

import pandas as pd
import streamlit as st


def render_stats(summary: dict, tracks: list[dict] | None = None) -> None:
    """
    Render the coverage statistics panel.

    Args:
        summary: Dict from WalkState.summary() with keys:
                    walked_segments, total_segments, pct_walked,
                    walked_km, total_km, unwalked_km.
        tracks:  Parsed track dicts (optional) for the history table.
    """
    pct_blocks = summary["pct_walked"]
    pct_streets = summary.get("pct_streets", 0.0)

    # ── Headline metrics ──────────────────────────────────────────────────────
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric(
        "Streets Completed",
        f"{summary.get('walked_streets', 0):,}",
        f"{pct_streets:.1f}% of {summary.get('total_streets', 0):,}",
    )
    c2.metric(
        "Blocks Walked",
        f"{summary['walked_segments']:,}",
        f"{pct_blocks:.1f}% of {summary['total_segments']:,}",
    )
    c3.metric(
        "Distance Walked",
        f"{summary['walked_km']:,.1f} km",
        f"{summary['unwalked_km']:,.1f} km remaining",
    )
    c4.metric(
        "Block Coverage",
        f"{pct_blocks:.1f}%",
    )
    c5.metric(
        "Network Size",
        f"{summary['total_km']:,.1f} km",
    )

    # ── Progress bars ─────────────────────────────────────────────────────────
    st.progress(
        min(pct_streets / 100, 1.0),
        text=f"**{pct_streets:.1f}%** of streets fully completed",
    )
    st.progress(
        min(pct_blocks / 100, 1.0),
        text=f"**{pct_blocks:.1f}%** of blocks walked",
    )

    # ── Walk history table ────────────────────────────────────────────────────
    if tracks:
        st.subheader("Recent uploads")
        rows = []
        for t in tracks:
            start = t.get("start_time")
            rows.append({
                "File": t.get("filename", t.get("track_name", "—")),
                "Activity": t.get("activity_type") or "—",
                "Date": start.strftime("%Y-%m-%d") if start else "—",
                "Points": f"{t.get('point_count', 0):,}",
            })
        st.dataframe(
            pd.DataFrame(rows),
            use_container_width=True,
            hide_index=True,
        )
