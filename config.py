"""App-wide configuration constants."""

# ── Coordinate Reference Systems ──────────────────────────────────────────────
WGS84_CRS = "EPSG:4326"
METRIC_CRS = "EPSG:3857"  # Web Mercator — for metre-based distance calculations

# ── Network Fetching ──────────────────────────────────────────────────────────
OSM_NETWORK_TYPE = "walk"

EXCLUDE_HIGHWAY_TYPES = {
    "service",
    "cycleway",
    "steps",
    "bridleway",
    "track",
    "path",
    "corridor",
    "elevator",
    "escalator",
}

NETWORK_CACHE_TTL_DAYS = 30

# ── Track Matching ────────────────────────────────────────────────────────────
DEFAULT_SNAP_DISTANCE_M = 20

# ── Map Display ───────────────────────────────────────────────────────────────
COLOR_WALKED = "#22c55e"
COLOR_UNWALKED = "#94a3b8"
COLOR_GPX_TRACK = "#f97316"

MAP_ZOOM_START = 14
MAP_TILES = "CartoDB positron"
