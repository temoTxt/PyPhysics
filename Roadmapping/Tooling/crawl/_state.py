"""Per-crawler 'last run' cursor tracking.

State lives in Roadmapping/Tooling/triage/_state.json (gitignored).
Each crawler stores its own last-run timestamp / cursor under a top-level
key matching the crawler name ('s2', 'arxiv', 'crossref', 'zotero').
"""

import json
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
STATE_PATH = REPO_ROOT / "Roadmapping" / "Tooling" / "triage" / "_state.json"


def load_state() -> dict:
    if not STATE_PATH.exists():
        return {}
    try:
        return json.loads(STATE_PATH.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def save_state(state: dict) -> None:
    STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    STATE_PATH.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")


def get_cursor(source: str, default: str | None = None) -> str | None:
    """Return the last-run cursor (ISO 8601 UTC datetime string) for a source."""
    return load_state().get(source, {}).get("last_run", default)


def set_cursor(source: str, when: str | None = None) -> None:
    """Update the last-run cursor for a source. Defaults to now (UTC)."""
    state = load_state()
    state.setdefault(source, {})["last_run"] = when or datetime.now(timezone.utc).isoformat()
    save_state(state)


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).isoformat()
