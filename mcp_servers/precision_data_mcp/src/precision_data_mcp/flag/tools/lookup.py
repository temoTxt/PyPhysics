"""FLAG-namespace lookup implementations.

Backed by ``precision_data_mcp/flag/data.json``, seeded from the FLAG 2024
report's N_f=2+1+1 averages.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

_DATA_PATH = Path(__file__).resolve().parent.parent / "data.json"


def _load() -> dict:
    with _DATA_PATH.open(encoding="utf-8") as f:
        return json.load(f)


def _cache_key(tool_name: str, args: dict) -> str:
    payload = json.dumps({"tool": tool_name, "args": args}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def _stamp(record: dict, *, tool_name: str, args: dict, source_revision: dict) -> dict:
    out = dict(record)
    out["retrieved_at"] = source_revision.get("retrieved_at", datetime.now(timezone.utc).isoformat())
    out["cache_key"] = _cache_key(tool_name, args)
    out["source_revision"] = {"flag_edition": source_revision.get("flag_edition", "unknown")}
    return out


def get_quantity(quantity_id: str) -> dict:
    db = _load()
    src_rev = db["$source_revision"]
    qs = db["quantities"]
    if quantity_id not in qs:
        raise KeyError(f"unknown FLAG quantity: {quantity_id!r}; known: {sorted(qs.keys())}")
    return _stamp(qs[quantity_id], tool_name="flag.get_quantity", args={"quantity_id": quantity_id}, source_revision=src_rev)


def list_quantities() -> list[dict]:
    db = _load()
    return [{"quantity_id": k, "description": v.get("description", "")} for k, v in db["quantities"].items()]
