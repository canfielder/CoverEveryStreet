.PHONY: help install install-dev sync lock launch lint test clean cache-clear

# ── Variables ─────────────────────────────────────────────────────────────────
APP       := app.py
VENV      := .venv
PYTHON    := $(VENV)/bin/python
UV        := uv

# ── Default target ────────────────────────────────────────────────────────────
help:
	@echo "🗺️  Walk Every Street — available commands"
	@echo ""
	@echo "  make install       📦 Create venv and install dependencies (prod)"
	@echo "  make install-dev   🔧 Create venv and install dependencies (dev + prod)"
	@echo "  make sync          🔄 Sync venv to match pyproject.toml (uv sync)"
	@echo "  make lock          🔒 Regenerate uv.lock from pyproject.toml"
	@echo "  make launch        🚀 Launch the Streamlit app"
	@echo "  make lint          🧹 Lint and format with ruff"
	@echo "  make test          🧪 Run tests with pytest"
	@echo "  make cache-clear   🗑️  Delete cached street networks and session data"
	@echo "  make clean         💣 Remove venv and Python cache files"
	@echo ""

# ── Environment setup ─────────────────────────────────────────────────────────
install:
	@echo "📦 Creating virtual environment and installing dependencies..."
	$(UV) venv --python 3.11 $(VENV)
	$(UV) pip install --python $(PYTHON) -r requirements.txt
	@echo "✅ Installation complete. Run 'make launch' to start the app."

install-dev:
	@echo "🔧 Installing dependencies (dev + prod) into virtual environment..."
	$(UV) venv --python 3.11 $(VENV)
	$(UV) pip install --python $(PYTHON) -e ".[dev]"
	@echo "✅ Dev installation complete."

sync:
	@echo "🔄 Syncing virtual environment to pyproject.toml..."
	$(UV) sync
	@echo "✅ Environment synced."

lock:
	@echo "🔒 Regenerating uv.lock from pyproject.toml..."
	$(UV) lock
	@echo "✅ Lock file updated."

# ── App ───────────────────────────────────────────────────────────────────────
launch:
	@echo "🚀 Launching Walk Every Street Streamlit app..."
	@echo "   Open http://localhost:8501 in your browser"
	$(UV) run streamlit run $(APP)

# ── Code quality ──────────────────────────────────────────────────────────────
lint:
	@echo "🧹 Running ruff linter and formatter..."
	$(UV) run ruff check --fix .
	$(UV) run ruff format .
	@echo "✅ Linting complete."

test:
	@echo "🧪 Running test suite..."
	$(UV) run pytest -v --cov=core --cov=ui --cov-report=term-missing
	@echo "✅ Tests complete."

# ── Data management ───────────────────────────────────────────────────────────
cache-clear:
	@echo "🗑️  Clearing cached street networks and session data..."
	rm -rf data/networks/*.gpkg
	rm -rf data/sessions/*.csv
	@echo "✅ Cache cleared. Networks will be re-fetched on next launch."

# ── Cleanup ───────────────────────────────────────────────────────────────────
clean:
	@echo "💣 Removing virtual environment and Python cache files..."
	rm -rf $(VENV)
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	@echo "✅ Clean complete. Run 'make install' to start fresh."
