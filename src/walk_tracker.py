"""
Persistent walk state management.

Saves which segments have been walked so progress accumulates across
Streamlit sessions — new GPX uploads add to existing history rather
than replacing it.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import pandas as pd

from config import WGS84_CRS
from paths import SESSION_DIR

logger = logging.getLogger(__name__)


@dataclass
class WalkState:
    """Persistent record of walked segments for one city/area."""

    cache_key: str
    walked_ids: set[str] = field(default_factory=set)
    # segment_id → earliest walk datetime
    walk_dates: dict[str, datetime] = field(default_factory=dict)
    # segment_id → cumulative walk count
    walk_counts: dict[str, int] = field(default_factory=dict)

    # ── Persistence ──────────────────────────────────────────────────────────

    @property
    def _path(self) -> Path:
        return SESSION_DIR / f"walks_{self.cache_key}.parquet"

    def save(self) -> None:
        """Serialize walk state to Parquet."""
        if not self.walked_ids:
            return
        records = [
            {
                "segment_id": sid,
                "walk_date": self.walk_dates.get(sid),
                "walk_count": self.walk_counts.get(sid, 1),
            }
            for sid in self.walked_ids
        ]
        df = pd.DataFrame(records)
        df.to_parquet(self._path, index=False)
        logger.info("Walk state saved: %d segments → %s", len(self.walked_ids), self._path)

    @classmethod
    def load(cls, cache_key: str) -> "WalkState":
        """Load walk state from disk, or return an empty state if none exists."""
        state = cls(cache_key=cache_key)
        path = SESSION_DIR / f"walks_{cache_key}.parquet"
        if not path.exists():
            return state
        try:
            df = pd.read_parquet(path)
            state.walked_ids = set(df["segment_id"].tolist())
            state.walk_dates = dict(zip(df["segment_id"], df["walk_date"]))
            state.walk_counts = dict(zip(df["segment_id"], df["walk_count"]))
            logger.info("Walk state loaded: %d segments from %s", len(state.walked_ids), path)
        except Exception as exc:
            logger.warning("Could not load walk state from %s: %s", path, exc)
        return state

    # ── Merging ───────────────────────────────────────────────────────────────

    def merge_from_gdf(self, matched_gdf: gpd.GeoDataFrame) -> None:
        """
        Incorporate newly matched segments into this walk state.

        matched_gdf must have columns: segment_id, walked, walk_date, walk_count.
        """
        walked = matched_gdf[matched_gdf["walked"]]
        for _, row in walked.iterrows():
            sid = row["segment_id"]
            self.walked_ids.add(sid)

            new_date = row.get("walk_date")
            if new_date is not None:
                old_date = self.walk_dates.get(sid)
                if old_date is None or new_date < old_date:
                    self.walk_dates[sid] = new_date

            self.walk_counts[sid] = (
                self.walk_counts.get(sid, 0) + int(row.get("walk_count", 1))
            )

    # ── Applying to network ───────────────────────────────────────────────────

    def apply_to_network(self, network_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Return a copy of network_gdf with walked/walk_date/walk_count
        columns populated from persisted state.
        """
        gdf = network_gdf.copy()
        gdf["walked"] = gdf["segment_id"].isin(self.walked_ids)
        gdf["walk_date"] = gdf["segment_id"].map(self.walk_dates)
        gdf["walk_count"] = gdf["segment_id"].map(self.walk_counts).fillna(0).astype(int)
        return gdf

    # ── Summary ───────────────────────────────────────────────────────────────

    def summary(self, network_gdf: gpd.GeoDataFrame) -> dict:
        """Compute coverage statistics against the full network."""
        total = len(network_gdf)
        walked_mask = network_gdf["segment_id"].isin(self.walked_ids)
        walked = int(walked_mask.sum())
        total_km = network_gdf["length_m"].sum() / 1000
        walked_km = network_gdf.loc[walked_mask, "length_m"].sum() / 1000
        pct = 100 * walked / total if total > 0 else 0.0
        return {
            "total_segments": total,
            "walked_segments": walked,
            "unwalked_segments": total - walked,
            "pct_walked": round(pct, 2),
            "total_km": round(total_km, 2),
            "walked_km": round(walked_km, 2),
            "unwalked_km": round(total_km - walked_km, 2),
        }

    def reset(self) -> None:
        """Clear all walk history and delete the persisted file."""
        self.walked_ids.clear()
        self.walk_dates.clear()
        self.walk_counts.clear()
        if self._path.exists():
            self._path.unlink()
            logger.info("Walk state reset for cache key: %s", self.cache_key)
