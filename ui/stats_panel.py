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
    pct = summary["pct_walked"]

    # ── Headline metrics ──────────────────────────────────────────────────────
    c1, c2, c3, c4 = st.columns(4)
    c1.metric(
        "Streets Walked",
        f"{summary['walked_segments']:,}",
        f"of {summary['total_segments']:,} total",
    )
    c2.metric(
        "Distance Walked",
        f"{summary['walked_km']:.1f} km",
        f"{summary['unwalked_km']:.1f} km remaining",
    )
    c3.metric(
        "Coverage",
        f"{pct:.1f}%",
    )
    c4.metric(
        "Network Size",
        f"{summary['total_km']:.1f} km",
    )

    # ── Progress bar ──────────────────────────────────────────────────────────
    st.progress(
        min(pct / 100, 1.0),
        text=f"**{pct:.1f}%** of streets walked",
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
