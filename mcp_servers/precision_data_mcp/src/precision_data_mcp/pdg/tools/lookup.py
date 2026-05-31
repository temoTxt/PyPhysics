"""PDG-namespace lookup implementations.

Backed by ``precision_data_mcp/pdg/data.json``, a hand-curated seed of PDG 2024
Review of Particle Physics values for the most-cited particles. Adding particles
is a mechanical follow-up PR per the seed's `$source_revision` block.

Per the umbrella REFRESH_POLICY.md, every returned record carries::

    {value, uncertainty, unit, source, retrieved_at, cache_key}

and multi-value disagreements (neutron lifetime puzzle, muon g-2 lattice-vs-Fermilab
tension) return *all* values with method labels rather than a single averaged value.
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
    """Decorate a raw seed record with retrieved_at + cache_key + source_revision."""
    out = dict(record)
    out["retrieved_at"] = source_revision.get("retrieved_at", datetime.now(timezone.utc).isoformat())
    out["cache_key"] = _cache_key(tool_name, args)
    out["source_revision"] = {"pdg_edition": source_revision.get("pdg_edition", "unknown")}
    return out


def _stamp_field(value_record: dict | None, *, tool_name: str, args: dict, source_revision: dict) -> dict | None:
    """Stamp a quantity sub-record (e.g. mass, lifetime). Handles single-value and multi-value cases."""
    if value_record is None:
        return None
    if "values" in value_record:
        return {
            "values": [_stamp(v, tool_name=tool_name, args=args, source_revision=source_revision) for v in value_record["values"]],
            "_note": value_record.get("_note"),
        }
    return _stamp(value_record, tool_name=tool_name, args=args, source_revision=source_revision)


def _resolve_particle_id(particle_id: str | int) -> tuple[str, dict]:
    """Accept a PDG MCID (int or stringified int) or a particle name/symbol; return (mcid_str, particle_dict)."""
    db = _load()
    particles = db["particles"]

    mcid_str = str(particle_id)
    if mcid_str in particles:
        return mcid_str, particles[mcid_str]

    needle = mcid_str.lower()
    for mcid, info in particles.items():
        if info["name"].lower() == needle or info["symbol"].lower() == needle:
            return mcid, info

    raise KeyError(f"unknown PDG particle: {particle_id!r}; known: {sorted(particles.keys())}")


def get_particle(particle_id: str | int) -> dict:
    """Return all PDG properties for a single particle (mass, lifetime, charge, spin).

    Multi-value fields (e.g. neutron lifetime under the bottle-vs-beam puzzle) are
    returned as ``{"values": [...]}`` per the umbrella's disagreement-representation
    rule.
    """
    db = _load()
    src_rev = db["$source_revision"]
    mcid, info = _resolve_particle_id(particle_id)

    args = {"particle_id": particle_id}
    return {
        "particle_id": mcid,
        "name": info["name"],
        "symbol": info["symbol"],
        "mass": _stamp_field(info.get("mass"), tool_name="pdg.get_particle", args=args, source_revision=src_rev),
        "lifetime": _stamp_field(info.get("lifetime"), tool_name="pdg.get_particle", args=args, source_revision=src_rev),
        "charge": _stamp_field(info.get("charge"), tool_name="pdg.get_particle", args=args, source_revision=src_rev),
        "spin": _stamp_field(info.get("spin"), tool_name="pdg.get_particle", args=args, source_revision=src_rev),
    }


def get_anomaly(particle_id: str | int) -> dict:
    """Return the anomalous magnetic moment record for a particle.

    The muon a_mu entry surfaces the canonical disagreement (Fermilab E989 +
    BNL E821 experimental value vs Theory Initiative 2020 data-driven SM
    prediction vs BMW 2021 lattice-QCD calculation) as a ``values`` list per the
    umbrella's disagreement-representation rule.
    """
    db = _load()
    src_rev = db["$source_revision"]
    mcid, _ = _resolve_particle_id(particle_id)

    anomalies = db.get("anomalies", {})
    if mcid not in anomalies:
        return {"particle_id": mcid, "status": "no_anomaly_data", "note": "Only leptons have curated anomalous-moment data; add via follow-up PR."}

    raw = anomalies[mcid]
    args = {"particle_id": particle_id}
    return {
        "particle_id": mcid,
        "particle": raw["particle"],
        "values": [_stamp(v, tool_name="pdg.get_anomaly", args=args, source_revision=src_rev) for v in raw["values"]],
        "_note": raw.get("_note"),
    }


def list_particles() -> list[dict]:
    """Enumerate seeded particles for agent discovery."""
    db = _load()
    out = []
    for mcid, info in db["particles"].items():
        out.append({"particle_id": mcid, "name": info["name"], "symbol": info["symbol"]})
    return out
