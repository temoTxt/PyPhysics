"""QED-namespace lookup implementations.

Backed by ``precision_data_mcp/qed/data.json``, a hand-curated seed of bound-state
QED precision observables (Lamb shifts, hyperfine, bound-electron g-factors,
high-precision transitions) for the species the Bethe-Salpeter (#50) and
hydrogenic-Z-scan (#82) campaigns actively cite.

Per the umbrella REFRESH_POLICY.md, every returned record carries::

    {value, uncertainty, unit, source, retrieved_at, cache_key}
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from precision_data_mcp.safety import apply_safety_contract

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
    out["source_revision"] = {"data_revision": source_revision.get("data_revision", "unknown")}
    return apply_safety_contract(out)


def _resolve_species(species: str) -> tuple[str, dict]:
    """Accept canonical species ID (``H``, ``He II``, ``Li III``, ``Si XIV``, ``muonic_H``, ``muonium``, ``positronium``)
    or any documented alias; return (canonical_id, species_dict).
    """
    db = _load()
    species_map = db["species"]

    if species in species_map:
        return species, species_map[species]

    needle = species.lower()
    for canonical, info in species_map.items():
        if canonical.lower() == needle:
            return canonical, info
        for alias in info.get("_aliases", []):
            if alias.lower() == needle:
                return canonical, info

    raise KeyError(f"unknown QED species: {species!r}; known: {sorted(species_map.keys())}")


def get_lamb_shift(species: str, transition: str, method: str | None = None) -> dict:
    """Lamb-shift value for a given species + transition (e.g. ``H``, ``"2S1/2-2P1/2"``).

    Multi-source entries (e.g. H 2S1/2-2P1/2 has both a CODATA-2018 adjusted value
    and the direct Lundeen-Pipkin 1981 RF measurement) follow the umbrella's
    disagreement-representation rule: when ``method=None`` returns all values as
    a list; otherwise returns the single matching value. Use the ``method`` arg
    or filter on ``safe_for_model_verification`` to pick the experimental anchor.
    """
    db = _load()
    src_rev = db["$source_revision"]
    canonical, info = _resolve_species(species)
    shifts = info.get("lamb_shifts", {})
    if transition not in shifts:
        raise KeyError(f"no Lamb-shift data for species={canonical!r}, transition={transition!r}; known transitions: {sorted(shifts.keys())}")
    entry = shifts[transition]
    args_base = {"species": species, "transition": transition, "method": method}

    if "values" in entry:
        candidates = entry["values"]
        if method is not None:
            for rec in candidates:
                if rec.get("method") == method:
                    return _stamp(rec, tool_name="qed.get_lamb_shift", args=args_base, source_revision=src_rev)
            raise KeyError(f"no Lamb-shift method {method!r} for species={canonical!r}, transition={transition!r}; known methods: {[r.get('method') for r in candidates]}")
        return {
            "species": canonical,
            "transition": transition,
            "values": [_stamp(r, tool_name="qed.get_lamb_shift", args={"species": species, "transition": transition, "method": r.get("method")}, source_revision=src_rev) for r in candidates],
            "_note": entry.get("_note"),
        }

    return _stamp(entry, tool_name="qed.get_lamb_shift", args=args_base, source_revision=src_rev)


def get_hyperfine(species: str, level: str) -> dict:
    """Hyperfine-splitting value (e.g. ``H``, ``"1s2S1/2"``; ``positronium``, ``"1s_ortho_para"``)."""
    db = _load()
    src_rev = db["$source_revision"]
    canonical, info = _resolve_species(species)
    hfs = info.get("hyperfine", {})
    if level not in hfs:
        raise KeyError(f"no hyperfine data for species={canonical!r}, level={level!r}; known levels: {sorted(hfs.keys())}")
    args = {"species": species, "level": level}
    return _stamp(hfs[level], tool_name="qed.get_hyperfine", args=args, source_revision=src_rev)


def get_g_factor(species: str, electron_state: str) -> dict:
    """Bound-electron g-factor for a hydrogenic-ion electron state (e.g. ``Si XIV``, ``"1s2S1/2"``)."""
    db = _load()
    src_rev = db["$source_revision"]
    canonical, info = _resolve_species(species)
    gs = info.get("g_factor_bound_electron", {})
    if electron_state not in gs:
        raise KeyError(f"no g-factor data for species={canonical!r}, electron_state={electron_state!r}; known states: {sorted(gs.keys())}")
    args = {"species": species, "electron_state": electron_state}
    return _stamp(gs[electron_state], tool_name="qed.get_g_factor", args=args, source_revision=src_rev)


def get_transition_precision(species: str, transition: str) -> dict:
    """High-precision transition frequency (e.g. ``H``, ``"1S-2S"``)."""
    db = _load()
    src_rev = db["$source_revision"]
    canonical, info = _resolve_species(species)
    trs = info.get("transitions", {})
    if transition not in trs:
        raise KeyError(f"no transition data for species={canonical!r}, transition={transition!r}; known: {sorted(trs.keys())}")
    args = {"species": species, "transition": transition}
    return _stamp(trs[transition], tool_name="qed.get_transition_precision", args=args, source_revision=src_rev)


def list_species() -> list[dict]:
    db = _load()
    return [
        {"species": canonical, "description": info.get("description", ""), "aliases": info.get("_aliases", [])}
        for canonical, info in db["species"].items()
    ]


def list_observables(species: str) -> dict:
    """For a given species, list which observables are curated."""
    canonical, info = _resolve_species(species)
    return {
        "species": canonical,
        "lamb_shifts": sorted(info.get("lamb_shifts", {}).keys()),
        "hyperfine": sorted(info.get("hyperfine", {}).keys()),
        "g_factor_bound_electron": sorted(info.get("g_factor_bound_electron", {}).keys()),
        "transitions": sorted(info.get("transitions", {}).keys()),
    }
