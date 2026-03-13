"""
Persistent walk state management.

Stores which BLOCKS have been walked (not individual edges). A segment
is displayed as walked if its block_id is in the persisted walked set.
Progress accumulates across sessions — new GPX uploads add to history.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import pandas as pd

from src.paths import SESSION_DIR

logger = logging.getLogger(__name__)


@dataclass
class WalkState:
    """Persistent record of walked blocks for one city/area."""

    cache_key: str
    walked_block_ids: set[str] = field(default_factory=set)
    # block_id → earliest walk datetime
    walk_dates: dict[str, datetime] = field(default_factory=dict)
    # block_id → cumulative walk count
    walk_counts: dict[str, int] = field(default_factory=dict)

    # ── Persistence ──────────────────────────────────────────────────────────

    @property
    def _path(self) -> Path:
        return SESSION_DIR / f"walks_{self.cache_key}.csv"

    def save(self) -> None:
        """Serialize walk state to CSV."""
        if not self.walked_block_ids:
            return
        records = [
            {
                "block_id": bid,
                "walk_date": self.walk_dates.get(bid),
                "walk_count": self.walk_counts.get(bid, 1),
            }
            for bid in self.walked_block_ids
        ]
        pd.DataFrame(records).to_csv(self._path, index=False)
        logger.info("Walk state saved: %d blocks → %s", len(self.walked_block_ids), self._path)

    @classmethod
    def load(cls, cache_key: str) -> "WalkState":
        """Load walk state from disk, or return empty state if none exists."""
        state = cls(cache_key=cache_key)
        path = SESSION_DIR / f"walks_{cache_key}.csv"
        if not path.exists():
            return state
        try:
            df = pd.read_csv(path, parse_dates=["walk_date"])
            state.walked_block_ids = set(df["block_id"].tolist())
            state.walk_dates = dict(zip(df["block_id"], df["walk_date"]))
            state.walk_counts = dict(zip(df["block_id"], df["walk_count"]))
            logger.info("Walk state loaded: %d blocks from %s", len(state.walked_block_ids), path)
        except Exception as exc:
            logger.warning("Could not load walk state from %s: %s", path, exc)
        return state

    # ── Merging ───────────────────────────────────────────────────────────────

    def merge_from_gdf(self, matched_gdf: gpd.GeoDataFrame) -> None:
        """
        Incorporate newly walked blocks from a matched GDF into this state.
        matched_gdf must have columns: block_id, walked, walk_date, walk_count.
        """
        walked = matched_gdf[matched_gdf["walked"]].drop_duplicates("block_id")
        for _, row in walked.iterrows():
            bid = row["block_id"]
            self.walked_block_ids.add(bid)

            new_date = row.get("walk_date")
            if pd.notna(new_date) and new_date is not None:
                old_date = self.walk_dates.get(bid)
                if old_date is None or new_date < old_date:
                    self.walk_dates[bid] = new_date

            self.walk_counts[bid] = (
                self.walk_counts.get(bid, 0) + int(row.get("walk_count", 1))
            )

    # ── Applying to network ───────────────────────────────────────────────────

    def apply_to_network(self, network_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Return a copy of network_gdf with walked/walk_date/walk_count columns
        populated from persisted block state.
        """
        gdf = network_gdf.copy()
        gdf["walked"] = gdf["block_id"].isin(self.walked_block_ids)
        gdf["walk_date"] = gdf["block_id"].map(self.walk_dates)
        gdf["walk_count"] = gdf["block_id"].map(self.walk_counts).fillna(0).astype(int)
        return gdf

    # ── Summary ───────────────────────────────────────────────────────────────

    def summary(self, network_gdf: gpd.GeoDataFrame) -> dict:
        """Compute block-level coverage statistics."""
        # Count unique blocks
        all_blocks = network_gdf["block_id"].unique()
        total_blocks = len(all_blocks)
        walked_blocks = sum(1 for b in all_blocks if b in self.walked_block_ids)

        # Length uses the edge-level data (block_length_m per block via first edge)
        block_lengths = (
            network_gdf.drop_duplicates("block_id")
            .set_index("block_id")["block_length_m"]
        )
        total_km = block_lengths.sum() / 1000
        walked_km = block_lengths[block_lengths.index.isin(self.walked_block_ids)].sum() / 1000
        pct = 100 * walked_blocks / total_blocks if total_blocks > 0 else 0.0

        # Street-level stats: a street is fully walked when ALL its blocks are walked
        unique_blocks = network_gdf.drop_duplicates("block_id")[["name", "block_id"]]
        street_blocks = unique_blocks.groupby("name")["block_id"].apply(set)
        total_streets = len(street_blocks)
        walked_streets = sum(
            1 for blocks in street_blocks if blocks.issubset(self.walked_block_ids)
        )
        pct_streets = 100 * walked_streets / total_streets if total_streets > 0 else 0.0

        return {
            "total_segments": total_blocks,
            "walked_segments": walked_blocks,
            "unwalked_segments": total_blocks - walked_blocks,
            "pct_walked": round(pct, 2),
            "total_km": round(total_km, 2),
            "walked_km": round(walked_km, 2),
            "unwalked_km": round(total_km - walked_km, 2),
            "total_streets": total_streets,
            "walked_streets": walked_streets,
            "pct_streets": round(pct_streets, 2),
        }

    def reset(self) -> None:
        """Clear all walk history and delete the persisted file."""
        self.walked_block_ids.clear()
        self.walk_dates.clear()
        self.walk_counts.clear()
        if self._path.exists():
            self._path.unlink()
            logger.info("Walk state reset for cache key: %s", self.cache_key)
