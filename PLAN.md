# Walk Every Street — Refactor Plan

This document is the implementation reference for the Python Streamlit refactor. Work through sections in order. Each section is self-contained; cross-references are noted inline.

---

## Section 1: Data Cleanup

### Files and Folders to Delete

Remove the following from `data/`:

| Path | Reason |
|---|---|
| `data/archive/` | Old R `.RDS` files — no longer needed |
| `data/googledrive_transfer/` | Old R transfer artifacts |
| `data/open_streets_map_tiles/` | Unused tile cache |
| `data/intersection_nodes_2021-08-24.RDS` | R artifact |
| `data/segment_mapping_no_parallel_seperation_2021-08-24.RDS` | R artifact |
| `data/segmented_street_network_parallel_cases_labeled_10km.RDS` | R artifact |
| `data/zip_code_charlotte.csv` | R artifact |
| `data/networks/` contents | Regenerable from OSM (delete files, keep folder or recreate) |
| `data/sessions/` contents | Replaced by SQLite |

### Keep

- `data/gps_tracking_files/test_case/` — test GPX files (`walk_solo.gpx`, `walk_kelsey.gpx`, `run_solo.gpx`)

### New Folder Structure After Cleanup

```
data/
  gpx_inbox/              # drop new GPX files here — auto-detected on app start
  gpx_processed/          # files moved here after successful processing
  networks/               # .gpkg network cache (gitignored, regenerable)
  published/
    map_data.gpkg         # pre-rendered export for viewer app (git-tracked)
  walks.db                # SQLite database (gitignored — local only)
  gps_tracking_files/
    test_case/            # test GPX files — keep
```

---

## Section 2: Critical Issues to Fix Before Building

### 2.1 Stable Block IDs (CRITICAL)

**Problem:** Current block IDs (`b0000001`, `b0000002`, ...) are positional counters assigned during network construction. If the OSM network is force-refreshed, block IDs are reassigned from scratch. The SQLite database becomes orphaned — all rows in `block_walks` point to block IDs that no longer exist in the network.

**Fix:** Derive `block_id` as a deterministic hash of the block's constituent OSM segment IDs.

Algorithm:
1. For each block component, collect the sorted list of OSM segment IDs (e.g., `["1234_5678", "5678_9012"]`)
2. Join with `|` separator: `"1234_5678|5678_9012"`
3. `block_id = "b_" + md5("1234_5678|5678_9012")[:12]`

Result: `block_id` is stable as long as the underlying OSM node IDs don't change. OSM node IDs are stable for existing roads; road splits or merges are a rare edge case (see Section 10).

**Implementation:** Change `src/block_builder.py` to collect `segment_ids` per block component and hash them instead of incrementing a counter. After this change, force a network cache rebuild by clearing `data/networks/`.

### 2.2 Store Raw Coverage Ratio (Not Boolean)

**Problem:** The snap distance threshold is currently applied at GPX processing time. Blocks are marked walked or not-walked at write time. If the user wants to change the threshold (e.g., tighten it after reviewing false positives), all GPX files must be re-processed.

**Fix:** Store the raw `coverage_ratio` (0.0–1.0) per block per activity in `block_walks`. The threshold becomes a view-time filter. Users can retroactively adjust the threshold without touching the database.

- `coverage_ratio = 1.0` for manual assignments
- `coverage_ratio = <float>` for auto-matched blocks (fraction of block length covered by GPS track)

### 2.3 Activity Type Normalization

**Problem:** Different GPS platforms write different strings to the GPX `<type>` tag. Garmin writes `"running"`, `"walking"`, `"hiking"`. Strava writes different values. The internal color-coding system requires a normalized type.

**Fix:** Add a normalization layer in `src/gpx_parser.py`. Store both the raw string and the normalized type.

Normalization map — store in `config.py`:

```python
ACTIVITY_TYPE_MAP = {
    "running": "run",
    "run":     "run",
    "walking": "walk",
    "walk":    "walk",
    "hiking":  "walk",
    "hike":    "walk",
    # extensible — add new mappings here as needed
}
```

Pace-based fallback (when no `<type>` tag is present or the value is unmapped):
- Compute average speed from track start/end time and total distance
- Average speed > 9.7 km/h (6 mph) → `"run"`
- Otherwise → `"walk"`

### 2.4 `start_time` as the Dedup Key

**Problem:** Content hash alone is insufficient for deduplication. If a user exports a GPX from Garmin, renames the activity, and re-exports, the content hash changes but it is the same activity. Double-processing would write duplicate `block_walks` rows.

**Fix:** Use `start_time` as the `UNIQUE` constraint on the `activities` table. This is the true dedup key — two exports of the same activity have the same GPS start timestamp.

Content hash is still stored and used for inbox file dedup — detecting if a specific file has been dropped into `gpx_inbox/` before (prevent scanning the same file twice).

---

## Section 3: SQLite Schema

File: `data/walks.db`

```sql
CREATE TABLE activities (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    filename      TEXT NOT NULL,
    content_hash  TEXT NOT NULL,              -- inbox dedup: same file dropped twice
    start_time    TIMESTAMP NOT NULL UNIQUE,  -- true activity dedup key
    end_time      TIMESTAMP,
    activity_type TEXT NOT NULL,              -- normalized: "walk", "run", etc.
    companions    TEXT NOT NULL DEFAULT '[]', -- JSON array: [] or ["kelsey"]
    raw_gpx_type  TEXT,                       -- original string from GPX <type> tag
    track_name    TEXT,                       -- GPX activity name field
    processed_at  TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE block_walks (
    block_id        TEXT    NOT NULL,
    activity_id     INTEGER NOT NULL REFERENCES activities(id) ON DELETE CASCADE,
    coverage_ratio  REAL    NOT NULL,   -- 0.0–1.0 for auto; 1.0 for manual
    source          TEXT    NOT NULL,   -- "auto" | "manual"
    PRIMARY KEY (block_id, activity_id)
);

CREATE TABLE block_exclusions (
    block_id    TEXT    NOT NULL,
    activity_id INTEGER NOT NULL REFERENCES activities(id) ON DELETE CASCADE,
    PRIMARY KEY (block_id, activity_id)
    -- Populated when user removes an auto-matched block
    -- Prevents re-processing the same activity from re-adding the false positive
);
```

### Key Query Patterns

**Map render — get walked blocks with color info:**
For each block, find the earliest activity within the current date range, then look up `activity_type` and `companions` to derive color.

```sql
SELECT bw.block_id, a.activity_type, a.companions, MIN(a.start_time) AS walk_date
FROM block_walks bw
JOIN activities a ON a.id = bw.activity_id
WHERE a.start_time >= :start AND a.start_time <= :end
  AND bw.coverage_ratio >= :threshold
GROUP BY bw.block_id;
```

**Inbox dedup check:**
```sql
SELECT content_hash FROM activities;
SELECT start_time   FROM activities;
```

**Re-process guard — skip excluded blocks:**
```sql
INSERT INTO block_walks (block_id, activity_id, coverage_ratio, source)
SELECT :block_id, :activity_id, :ratio, 'auto'
WHERE NOT EXISTS (
    SELECT 1 FROM block_exclusions
    WHERE block_id = :block_id AND activity_id = :activity_id
);
```

**Removing an auto-matched block:**
```sql
DELETE FROM block_walks     WHERE block_id = :bid AND activity_id = :aid;
INSERT INTO block_exclusions (block_id, activity_id) VALUES (:bid, :aid);
```

**Removing a manual assignment:**
```sql
DELETE FROM block_walks WHERE block_id = :bid AND activity_id = :aid;
-- No exclusion row needed for manual assignments
```

---

## Section 4: Activity Categories and Colors

Stored in `config.py`. The key is `(activity_type, frozenset(companions))`. Adding a new category requires only adding a line to this dict — no code changes elsewhere.

```python
ACTIVITY_COLORS = {
    ("walk", frozenset()):             "#22c55e",  # solo walk — green
    ("run",  frozenset()):             "#3b82f6",  # solo run — blue
    ("walk", frozenset(["kelsey"])):   "#a855f7",  # walk with Kelsey — purple
}
COLOR_UNWALKED = "#94a3b8"  # slate gray — unwalked blocks
COLOR_FALLBACK = "#f97316"  # orange — walked but no matching category entry
```

**Color derivation at render time:**
1. For each walked block, retrieve the earliest activity's `activity_type` and `companions` (parsed from JSON)
2. Build lookup key: `(activity_type, frozenset(companions))`
3. Look up in `ACTIVITY_COLORS`; if missing, use `COLOR_FALLBACK`
4. Unwalked blocks always get `COLOR_UNWALKED`

---

## Section 5: New and Changed Files

### New Files

| File | Purpose |
|---|---|
| `src/database.py` | SQLite layer — all DB reads/writes. Schema creation, all query functions. |
| `src/inbox.py` | Scan `gpx_inbox/`, detect unprocessed files via content hash and `start_time` dedup. |
| `ui/pending_panel.py` | Pending activities UI — companion tagging, type override, process button. |
| `viewer.py` | Read-only Streamlit Cloud app (~50 lines). |
| `scripts/import_batch.py` | Garmin bulk FIT export → inbox. Separate script; requires `fit2gpx`. |

### Changed Files

| File | Change |
|---|---|
| `src/block_builder.py` | Stable hash-based block IDs (see Section 2.1) |
| `src/track_matcher.py` | Return `coverage_ratio` floats instead of booleans (see Section 2.2) |
| `src/gpx_parser.py` | Activity type normalization + pace-based fallback (see Section 2.3) |
| `src/paths.py` | Add `INBOX_DIR`, `PROCESSED_DIR`, `PUBLISHED_DIR`, `DB_PATH` |
| `src/walk_tracker.py` | **Delete** — replaced entirely by `src/database.py` |
| `config.py` | Add `ACTIVITY_TYPE_MAP`, `ACTIVITY_COLORS`; remove old color constants |
| `app.py` | New flow: inbox scan → pending panel → map |
| `ui/map_view.py` | Multi-color rendering by category + assign mode toggle |
| `ui/sidebar.py` | Add date range filter; remove file upload widget (now inbox-based) |
| `ui/stats_panel.py` | Category breakdown stats (km by category) |
| `.gitignore` | Add `walks.db`; carve out `data/published/` |

---

## Section 6: Application Architecture

### Local App (`app.py`) — Full Featured

**Flow on every app start:**

1. Scan `data/gpx_inbox/` for `.gpx` files
2. For each file:
   - Compute content hash → check against `activities.content_hash`
   - Parse `start_time` → check against `activities.start_time`
3. Files not in DB → shown in **Pending Activities** panel
4. Files matching an existing `start_time` but different `content_hash` → shown with warning: `"Re-export of [activity name] — update or skip?"`
5. User tags companions (multi-select, extensible list), confirms or overrides `activity_type`
6. Click **Process** → track matcher runs → `block_walks` rows written → file moved to `gpx_processed/`
7. Map renders using date filter from sidebar

**Map Interaction — Two Modes:**

- **View mode** (default): hover tooltip shows `block_id`, street name, block length, activity info
- **Assign mode** (toggle button above map): clicking a block opens the assignment panel

**Assignment Panel (Assign mode):**

- Shows selected block: `block_id`, street name, current assignment if any
- Dropdown: assign to an activity (sorted by date, most recent first)
- "Remove assignment" button if block is currently assigned:
  - Removing an auto-matched block → `DELETE FROM block_walks` + `INSERT INTO block_exclusions`
  - Removing a manual assignment → `DELETE FROM block_walks` only
- Confirm button → writes to DB, map re-renders

**Activities Panel (sidebar expander):**

- Table of all processed activities
- Inline edit: `activity_type`, `companions`
- Save → `UPDATE activities` row → map re-colors (no GPX re-processing required)

**Publish Button:**

1. Join network blocks with current walk state (all-time or respecting current date filter — user chooses)
2. Export to `data/published/map_data.gpkg` with columns:
   - `block_id`, `geometry`, `color`, `street_name`, `block_length_m`, `walked`, `activity_type`, `companions`, `walk_date`
3. User then commits and pushes manually:
   ```
   git add data/published/map_data.gpkg
   git commit -m "publish map update"
   git push
   ```

### Viewer App (`viewer.py`) — Streamlit Cloud

Approximately 50 lines. No sidebar controls for processing, no assignment panel, no DB access.

**Logic:**
1. Load `data/published/map_data.gpkg` via `geopandas.read_file()`
2. Optional: date range slider (enabled if `walk_date` column is present — always include it in the published file)
3. Filter rows by date range
4. Render Folium map using pre-computed `color` column
5. Show minimal stats (total km walked)

**Deployment:** Streamlit Cloud with `viewer.py` as the entry point. The `data/published/` directory is git-tracked so Streamlit Cloud can read it directly from the repo.

---

## Section 7: Gitignore Updates

Add or update the following in `.gitignore`:

```gitignore
# SQLite database — local only, contains personal walk history
data/walks.db

# Generated/cached network data — regenerable from OSM
data/networks/*.gpkg
cache/

# Processed GPX files — already in DB, no need to track
data/gpx_processed/

# Published snapshot is git-tracked for Streamlit Cloud viewer
# data/published/ is intentionally NOT gitignored
```

Note: `data/gpx_inbox/` should also not be gitignored, but its contents are typically transient. Consider adding `data/gpx_inbox/*.gpx` if you don't want inbox files tracked.

---

## Section 8: Batch Import Script (`scripts/import_batch.py`)

For importing historical Garmin bulk export (FIT files). Kept separate from the main app — requires `fit2gpx`, which is not in `pyproject.toml` and should not be added (it pulls in heavy dependencies not needed for normal operation).

**Usage:**
```bash
uv run python scripts/import_batch.py \
  --source ~/Downloads/garmin_export.zip \
  --before 2026-01-01
```

**Steps:**
1. Unzip the archive to a temp directory
2. Locate FIT files in `DI_CONNECT/DI-Connect-Fitness/`
3. Parse `summarizedActivities.json` for activity metadata (name, type, date)
4. Filter to activities with date `< --before`
5. Convert each FIT file to GPX via `fit2gpx`
6. For each converted GPX, parse `start_time` and compute content hash
7. Check `start_time` against existing DB records — skip duplicates
8. Copy new GPX files to `data/gpx_inbox/`
9. Print summary:
   ```
   Added to inbox:      47 files
   Skipped (duplicate): 12 files
   Skipped (after cutoff): 3 files
   ```

**After running:** Open the local app to review pending activities (companion tagging, type confirmation) and click Process.

---

## Section 9: Build Order

Execute steps in this order. The app should remain functional (or at least launchable) at each step.

| Step | Task | Notes |
|---|---|---|
| 1 | **Data cleanup** | Delete R artifacts; create new folder structure (`gpx_inbox/`, `gpx_processed/`, `published/`) |
| 2 | **`src/paths.py`** | Add `INBOX_DIR`, `PROCESSED_DIR`, `PUBLISHED_DIR`, `DB_PATH` |
| 3 | **`config.py`** | Add `ACTIVITY_TYPE_MAP`, `ACTIVITY_COLORS`; keep old color constants temporarily to avoid breaking existing imports |
| 4 | **`src/block_builder.py`** | Stable hash-based block IDs. After merging: clear `data/networks/` to force cache rebuild on next app launch |
| 5 | **`src/database.py`** | Create SQLite layer: schema creation on first connect, all query functions (`get_walked_blocks`, `insert_activity`, `insert_block_walk`, `add_exclusion`, etc.) |
| 6 | **`src/gpx_parser.py`** | Activity type normalization using `ACTIVITY_TYPE_MAP`; pace-based fallback |
| 7 | **`src/track_matcher.py`** | Return `coverage_ratio` floats instead of booleans |
| 8 | **`src/inbox.py`** | Inbox scanner: content hash + `start_time` dedup logic; returns lists of new files and re-export candidates |
| 9 | **`ui/pending_panel.py`** | Pending activities UI: companion multi-select, type override dropdown, Process button |
| 10 | **`app.py`** | Wire new flow together: inbox scan → pending panel → map |
| 11 | **`ui/map_view.py`** | Multi-color block rendering by category; assign mode toggle |
| 12 | **`ui/sidebar.py`** | Date range filter; remove file upload widget |
| 13 | **`ui/stats_panel.py`** | Category breakdown: km walked per category |
| 14 | **Delete `src/walk_tracker.py`** | Only after `database.py` is wired in and app is verified working |
| 15 | **`viewer.py`** | Read-only Streamlit Cloud app |
| 16 | **`.gitignore` updates** | Add `walks.db`; carve out `data/published/` |
| 17 | **`scripts/import_batch.py`** | Garmin batch import script |

---

## Section 10: Open Questions and Future Work

**Strava API integration**
OAuth2 flow to fetch activities automatically from Strava, rather than requiring manual GPX export. Out of scope for this refactor. The inbox model accommodates it — Strava activities would be fetched and written to `gpx_inbox/` by a separate sync script.

**Multi-area support**
The current DB schema works for one geographic area (Charlotte). Future enhancement: add `area_id` to `block_walks` and `activities` to track multiple cities. Would require a corresponding `areas` table and area selection in the app. Accept the single-area limitation for now.

**Network refresh strategy**

When the OSM network is force-refreshed, only re-process the blocks that actually changed — not all GPX files against all blocks. The unchanged blocks keep their walk history completely untouched.

**Network refresh is always manual.** The app never auto-refreshes. It shows the network cache age ("Last fetched: 47 days ago") in the sidebar so the user can decide when a refresh is warranted. The existing `NETWORK_CACHE_TTL_DAYS` constant should be removed — the cache never auto-expires.

When the user presses **Refresh Network**:
1. Before refresh: snapshot current set of block IDs from the network cache
2. Fetch fresh network from OSM, rebuild blocks
3. After refresh: compute the new set of block IDs
4. Derive three sets:
   - **Unchanged** (`old ∩ new`) — walk history stays, zero re-processing
   - **New** (`new − old`) — run all GPX files against only these blocks, scoped to their geographic bounding box (only GPX tracks that overlap the affected area are candidates)
   - **Gone** (`old − new`) — delete their `block_walks` auto entries; flag any manual entries as "needs review" in the UI
5. Show a summary: "Refresh complete — N new blocks matched, M blocks removed, X manual assignments need review"
6. GPX files are the source of truth — they live permanently in `gpx_processed/` for exactly this reason. Never delete them.

Manual assignments for gone blocks cannot be re-derived from GPX (that was the point — the GPX didn't match). These must be flagged and the user re-applies the small number of corrections. Manual exclusions for gone blocks can simply be dropped (the block no longer exists).

This makes a typical network refresh (one new road splits a few blocks) very cheap: only a handful of new blocks need matching, only GPX files geographically near the change are candidates.

**Viewer date filter**
The published GeoPackage includes a `walk_date` column. The viewer app should expose a date range slider so viewers can see the map as it looked at any point in time. Always include `walk_date` in the publish export.

**Stats by category**
`ui/stats_panel.py` should show kilometers walked broken down by category — solo walk, solo run, walk with Kelsey, etc. — using the same `ACTIVITY_COLORS` dict keys as the map rendering.

**`fit2gpx` in CI**
The batch import script requires `fit2gpx`. If you ever want to run it in a CI environment, it would need to be a separate optional dependency group in `pyproject.toml` (e.g., `[project.optional-dependencies] batch = ["fit2gpx"]`). For now, run it manually with `uv run --with fit2gpx`.
