"""
Block identification.

Groups osmnx simplified edges into "blocks" — sections of named street
between cross streets. Unnamed edges are dropped. Parallel edges of the
same name (divided roads) share a block_id so that walking either side
counts as covering the block.

Algorithm overview:
  1. Drop unnamed edges.
  2. Build a graph from remaining named edges.
  3. Identify "boundary nodes" — nodes where two or more differently-named
     streets meet. These are the cross-street intersections that delimit blocks.
  4. Use the virtual-node technique: at each boundary node, give every
     incident edge its own copy of that node. This disconnects same-named
     edges that share a cross-street node without affecting edges that share
     a mid-block node. Each connected component of the resulting per-name
     subgraph is one block.
  5. Merge parallel same-name blocks within PARALLEL_MERGE_DISTANCE_M metres
     (the two sides of a divided road).
  6. Add block_length_m = sum of edge lengths for all edges in a block.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from itertools import combinations

import geopandas as gpd
import networkx as nx
import shapely

from config import METRIC_CRS

logger = logging.getLogger(__name__)

# Blocks with the same street name whose geometries are within this distance
# (metres) are considered the two sides of a divided road and merged.
PARALLEL_MERGE_DISTANCE_M = 25


def build_blocks(network_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Assign block_id to each edge and remove unnamed edges.

    Args:
        network_gdf: Edge GDF from src.network with columns:
                     segment_id, u, v, name, highway, length_m, geometry.

    Returns:
        GDF with unnamed edges removed and two new columns:
            block_id       (str)   — shared by all edges in the same block
            block_length_m (float) — total length of all edges in the block
    """
    # ── Normalize name column ─────────────────────────────────────────────────
    # OSM name can be a list (e.g. ["Main St", "Route 1"]) — take the first.
    def _normalize_name(val) -> str | None:
        if isinstance(val, list):
            val = val[0] if val else None
        if val is None:
            return None
        return str(val).strip() or None

    gdf = network_gdf.copy()
    gdf["name"] = gdf["name"].apply(_normalize_name)

    # ── Find boundary nodes using the FULL network ────────────────────────────
    # Build a graph from ALL edges (named + unnamed) so that unnamed streets
    # (e.g. "35th Street" missing a name tag in OSM) still act as block
    # boundaries for the named streets that cross them.
    # Unnamed edges get name="" so they count as a distinct name at any node
    # they share with a named street.
    G_full = nx.Graph()
    G_full.add_edges_from(
        (row.u, row.v, {"name": row.name if row.name is not None else ""})
        for row in gdf.itertuples(index=False)
    )
    node_incident_names_full: dict[int, set[str]] = {
        node: {G_full[node][nb]["name"] for nb in G_full.neighbors(node)}
        for node in G_full.nodes()
    }
    boundary_nodes: set[int] = {
        node for node, names in node_incident_names_full.items()
        if len(names) > 1
    }

    # ── Drop unnamed (after boundary detection) ───────────────────────────────
    named = gdf[gdf["name"].notna()].copy()
    n_dropped = len(gdf) - len(named)
    logger.info("Named edges: %d (dropped %d unnamed)", len(named), n_dropped)
    gdf = named

    if gdf.empty:
        return gdf

    # ── Assign block IDs via virtual-node technique ───────────────────────────
    # For each street name, build a subgraph where boundary nodes are replaced
    # by per-edge virtual copies. This splits the graph at cross-street nodes
    # so each connected component is exactly one block.
    edge_to_block: dict[str, str] = {}
    block_counter = 0

    for name, name_gdf in gdf.groupby("name"):
        EG = nx.Graph()
        EG.add_edges_from(
            (
                (row.u, row.segment_id) if row.u in boundary_nodes else row.u,
                (row.v, row.segment_id) if row.v in boundary_nodes else row.v,
                {"segment_id": row.segment_id},
            )
            for row in name_gdf.itertuples(index=False)
        )

        for component in nx.connected_components(EG):
            block_id = f"b{block_counter:07d}"
            block_counter += 1
            for u, v, data in EG.subgraph(component).edges(data=True):
                edge_to_block[data["segment_id"]] = block_id

    gdf["block_id"] = gdf["segment_id"].map(edge_to_block)
    gdf = gdf.dropna(subset=["block_id"]).copy()

    # ── Merge parallel same-name blocks (divided roads) ───────────────────────
    gdf = _merge_parallel_blocks(gdf)

    # ── Add block_length_m ────────────────────────────────────────────────────
    block_lengths = gdf.groupby("block_id")["length_m"].sum()
    gdf["block_length_m"] = gdf["block_id"].map(block_lengths)

    n_blocks = gdf["block_id"].nunique()
    logger.info("Built %d blocks from %d named edges", n_blocks, len(gdf))
    return gdf


def _merge_parallel_blocks(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Merge same-name blocks that are geometrically close — the two sides of
    a divided road. Uses a spatial index for efficiency.
    """
    proj = gdf.to_crs(METRIC_CRS)

    # Dissolve edges by block_id to get one geometry per block
    block_geoms: dict[str, object] = proj.dissolve(by="block_id")["geometry"].to_dict()
    block_names: dict[str, str] = (
        gdf.drop_duplicates("block_id").set_index("block_id")["name"].to_dict()
    )

    # Group blocks by street name for scoped comparison
    name_to_blocks: dict[str, list[str]] = defaultdict(list)
    for bid, name in block_names.items():
        name_to_blocks[name].append(bid)

    # Union-find for merging
    parent: dict[str, str] = {bid: bid for bid in block_geoms}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    merge_count = 0
    for name, block_ids in name_to_blocks.items():
        if len(block_ids) < 2:
            continue

        geom_list = [block_geoms[bid] for bid in block_ids]
        tree = shapely.STRtree(geom_list)

        for i, bid1 in enumerate(block_ids):
            # Buffer the block geometry and query for close neighbours
            search_area = geom_list[i].buffer(PARALLEL_MERGE_DISTANCE_M)
            candidates = tree.query(search_area)
            for j in candidates:
                if j <= i:
                    continue
                bid2 = block_ids[j]
                dist = geom_list[i].distance(geom_list[j])
                if 0 < dist <= PARALLEL_MERGE_DISTANCE_M:
                    # Also require centroids to be close — prevents corner-proximity
                    # false merges on divided boulevards where the endpoint of block N
                    # on the north side is ~road-width from the endpoint of block N+1
                    # on the south side. Centroids of truly parallel blocks are
                    # ~road-width apart; corner-adjacent blocks are ~block-length apart.
                    centroid_dist = geom_list[i].centroid.distance(geom_list[j].centroid)
                    if centroid_dist <= PARALLEL_MERGE_DISTANCE_M:
                        union(bid1, bid2)
                        merge_count += 1

    if merge_count:
        gdf = gdf.copy()
        gdf["block_id"] = gdf["block_id"].map(lambda x: find(x))
        logger.info("Merged %d parallel block pairs (divided roads)", merge_count)

    return gdf
