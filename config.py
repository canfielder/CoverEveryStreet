"""App-wide configuration constants."""

# ── Coordinate Reference Systems ──────────────────────────────────────────────
WGS84_CRS = "EPSG:4326"
METRIC_CRS = "EPSG:3857"  # Web Mercator — for metre-based distance calculations

# ── Network Fetching ──────────────────────────────────────────────────────────
OSM_NETWORK_TYPE = "drive"

EXCLUDE_HIGHWAY_TYPES = {
    # Highways — not walkable
    "motorway",
    "motorway_link",
    "trunk",
    "trunk_link",
    # Unnamed service roads (driveways, parking aisles) handled separately in network.py
    "service",
}

NETWORK_CACHE_TTL_DAYS = 30

# ── Track Matching ────────────────────────────────────────────────────────────
DEFAULT_SNAP_DISTANCE_M = 20

# Minimum fraction of a segment's length that must fall within snap_distance_m
# of a GPS track for the segment to count as "walked". 0.95 = 95%.
DEFAULT_COVERAGE_THRESHOLD = 0.95

# ── Map Display ───────────────────────────────────────────────────────────────
COLOR_WALKED = "#22c55e"
COLOR_UNWALKED = "#94a3b8"
COLOR_GPX_TRACK = "#f97316"

MAP_ZOOM_START = 14
MAP_TILES = "CartoDB positron"
