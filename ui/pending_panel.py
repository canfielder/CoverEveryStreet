"""
Pending activities panel.

Shows GPX files found in the inbox that haven't been ingested yet.
For each file the user can:
  - Set the activity type (walk / run / cycle)
  - Tag companions
  - Confirm (Process) or skip

When the user clicks "Process All", all confirmed activities are
ingested into the database in one pass.
"""

from __future__ import annotations

import streamlit as st

from config import KNOWN_COMPANIONS


def render_pending_panel(
    pending: list[dict],
    network_gdf,
    snap_distance_m: float,
) -> int:
    """
    Render the pending activities review UI.

    Args:
        pending:         List of parsed track dicts from inbox.scan_inbox().
        network_gdf:     Current network GDF (needed for block matching).
        snap_distance_m: Snap distance used for block matching.

    Returns:
        Number of activities successfully ingested (0 if nothing processed).
    """
    if not pending:
        return 0

    st.subheader(f"📥 {len(pending)} new activity{'s' if len(pending) != 1 else ''} in inbox")
    st.caption("Review each activity, then click **Process All** to save.")

    # Per-activity metadata editors
    activity_configs: list[dict] = []
    for i, track in enumerate(pending):
        with st.expander(
            _activity_label(track),
            expanded=(i == 0),
        ):
            col_type, col_comp = st.columns(2)

            chosen_type = col_type.selectbox(
                "Activity type",
                options=["walk", "run", "cycle"],
                index=["walk", "run", "cycle"].index(
                    track.get("activity_type", "walk")
                    if track.get("activity_type") in ("walk", "run", "cycle")
                    else "walk"
                ),
                key=f"type_{i}",
            )

            chosen_companions = col_comp.multiselect(
                "Companions",
                options=KNOWN_COMPANIONS,
                default=[],
                key=f"companions_{i}",
            )

            activity_configs.append({
                "track": track,
                "activity_type": chosen_type,
                "companions": chosen_companions,
            })

    st.divider()
    col_btn, col_info = st.columns([1, 3])
    process = col_btn.button(
        "▶ Process All",
        type="primary",
        use_container_width=True,
    )
    col_info.caption(
        f"This will match {len(pending)} activity/activities to street blocks "
        "and move files to the processed folder."
    )

    if not process:
        return 0

    # Ingest all pending activities
    from src.inbox import ingest_activity

    ingested = 0
    progress = st.progress(0, text="Processing activities…")
    for step, cfg in enumerate(activity_configs):
        track = cfg["track"]
        try:
            ingest_activity(
                track=track,
                companions=cfg["companions"],
                activity_type_override=cfg["activity_type"],
                network_gdf=network_gdf,
                snap_distance_m=snap_distance_m,
            )
            ingested += 1
        except Exception as exc:
            st.warning(f"Failed to ingest {track['filename']}: {exc}")
        progress.progress(
            (step + 1) / len(activity_configs),
            text=f"Processed {step + 1} / {len(activity_configs)}",
        )

    progress.empty()
    if ingested:
        st.success(f"✅ Ingested {ingested} activity/activities. Refreshing map…")
        st.rerun()

    return ingested


# ── Helpers ───────────────────────────────────────────────────────────────────

def _activity_label(track: dict) -> str:
    start = track.get("start_time")
    date_str = start.strftime("%Y-%m-%d") if start else "Unknown date"
    name = track.get("track_name") or track.get("filename") or "Unnamed"
    atype = track.get("activity_type", "walk").capitalize()
    return f"{date_str}  ·  {atype}  ·  {name}"
