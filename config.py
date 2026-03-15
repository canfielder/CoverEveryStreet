"""App-wide configuration constants."""

# ── Coordinate Reference Systems ──────────────────────────────────────────────
WGS84_CRS  = "EPSG:4326"
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

# ── Track Matching ────────────────────────────────────────────────────────────
DEFAULT_SNAP_DISTANCE_M    = 20
DEFAULT_COVERAGE_THRESHOLD = 0.80  # applied at view time, not write time

# ── Activity Type Normalization ───────────────────────────────────────────────
# Maps raw GPX <type> strings to internal normalized types.
# Add new entries here as needed — no schema changes required.
ACTIVITY_TYPE_MAP: dict[str, str] = {
    "running":  "run",
    "run":      "run",
    "walking":  "walk",
    "walk":     "walk",
    "hiking":   "walk",
    "hike":     "walk",
    "cycling":  "cycle",
    "biking":   "cycle",
}

# Speed threshold for pace-based fallback when GPX has no <type> tag
PACE_RUN_THRESHOLD_KPH = 9.7  # > this → "run", else "walk"

# ── Activity Colors ───────────────────────────────────────────────────────────
# Key: (activity_type, frozenset(companions))
# Add new combinations here to extend color coding — no code changes needed.
ACTIVITY_COLORS: dict[tuple, str] = {
    ("walk", frozenset()):             "#22c55e",  # solo walk   — green
    ("run",  frozenset()):             "#3b82f6",  # solo run    — blue
    ("walk", frozenset(["kelsey"])):   "#a855f7",  # walk+Kelsey — purple
    ("run",  frozenset(["kelsey"])):   "#ec4899",  # run+Kelsey  — pink
}

COLOR_UNWALKED  = "#94a3b8"  # gray
COLOR_FALLBACK  = "#f97316"  # orange — any category not explicitly mapped
COLOR_GPX_TRACK = "#f97316"

# ── Map Display ───────────────────────────────────────────────────────────────
MAP_ZOOM_START = 14
MAP_TILES      = "CartoDB positron"

# ── Known Companions ─────────────────────────────────────────────────────────
# Drives the companion multi-select in the pending panel UI.
# Add names here to make them available as options.
KNOWN_COMPANIONS: list[str] = ["Kelsey"]
