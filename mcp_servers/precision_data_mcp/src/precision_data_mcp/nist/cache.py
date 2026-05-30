"""On-disk JSON cache for NIST ASD lookups.

Cache directory resolution order:
  1. ``NIST_MCP_CACHE_DIR`` env var, if set.
  2. Platform default ``~/.cache/nist_mcp/``.

Entries are keyed by a stable hash of ``CACHE_SCHEMA_VERSION`` + the normalized
query, so bumping ``CACHE_SCHEMA_VERSION`` after a parser/schema change
invalidates old entries instead of returning stale-format payloads. Each entry
stores a fetch timestamp; pass ``ttl`` (seconds) or ``refresh=True`` to bypass a
cached entry. Clear the cache by deleting the directory.
"""

import hashlib
import json
import os
import time
from functools import lru_cache
from pathlib import Path

CACHE_SCHEMA_VERSION = "1"


@lru_cache(maxsize=1)
def cache_dir() -> Path:
    env = os.environ.get("NIST_MCP_CACHE_DIR")
    base = Path(env).expanduser() if env else Path.home() / ".cache" / "nist_mcp"
    base = base.resolve()
    base.mkdir(parents=True, exist_ok=True)
    return base


def _path(query: str) -> Path:
    raw = f"{CACHE_SCHEMA_VERSION}\x00{query}".encode("utf-8")
    return cache_dir() / f"{hashlib.sha256(raw).hexdigest()}.json"


def cache_get(query: str, ttl: float | None = None) -> dict | None:
    """Return the cached payload for ``query``, or None if absent/expired/corrupt."""
    p = _path(query)
    if not p.exists():
        return None
    try:
        entry = json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None
    if ttl is not None and (time.time() - entry.get("_fetched_at", 0)) > ttl:
        return None
    return entry.get("payload")


def cache_put(query: str, payload: dict) -> None:
    """Write ``payload`` for ``query``, tagged with the schema version + timestamp."""
    entry = {
        "_schema": CACHE_SCHEMA_VERSION,
        "_fetched_at": time.time(),
        "payload": payload,
    }
    try:
        _path(query).write_text(json.dumps(entry), encoding="utf-8")
    except OSError:
        pass
