"""Nuclear-namespace lookup implementations.

Backed by ``precision_data_mcp/nuclear/data.json``. The proton-charge-radius
puzzle is the load-bearing disagreement-representation demonstrator: calling
``nuclear.get_charge_radius("proton")`` returns *all* values with method labels,
never an averaged "best" value, per umbrella issue #92 §Resolved-decisions #5.
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
    out["source_revision"] = {"data_revision": source_revision.get("data_revision", "unknown")}
    return out


def _resolve_nucleus(nucleus: str) -> tuple[str, dict]:
    db = _load()
    nuclei = db["nuclei"]

    if nucleus in nuclei:
        return nucleus, nuclei[nucleus]

    needle = nucleus.lower()
    for canonical, info in nuclei.items():
        if canonical.lower() == needle:
            return canonical, info
        for alias in info.get("_aliases", []):
            if alias.lower() == needle:
                return canonical, info

    raise KeyError(f"unknown nucleus: {nucleus!r}; known: {sorted(nuclei.keys())}")


def get_charge_radius(nucleus: str, method: str | None = None) -> dict | list[dict]:
    """Return charge-radius value(s) for a nucleus.

    If ``method`` is None and the nucleus has multi-value disagreement (the
    proton-radius puzzle), returns ALL values as a list — never an averaged
    single value, per umbrella §Resolved-decisions #5. If ``method`` is
    specified, returns the single matching value.

    Args:
        nucleus: nucleus ID or alias (``p``, ``proton``, ``d``, ``deuteron``, ``n``, ``neutron``).
        method: optional method label to select a single value
            (``muonic_hydrogen_spectroscopy``, ``electron_scattering``,
            ``1S-3S_spectroscopy``, ``PRad_electron_scattering``,
            ``CODATA_2018_pre_muonic_average``, ...).
    """
    db = _load()
    src_rev = db["$source_revision"]
    canonical, info = _resolve_nucleus(nucleus)

    cr = info.get("charge_radius")
    if cr is None:
        raise KeyError(f"no charge-radius data for nucleus={canonical!r}")

    args_base = {"nucleus": nucleus, "method": method}

    if "values" in cr:
        candidates = cr["values"]
        if method is not None:
            for rec in candidates:
                if rec.get("method") == method:
                    return _stamp(rec, tool_name="nuclear.get_charge_radius", args=args_base, source_revision=src_rev)
            raise KeyError(f"no charge-radius measurement for nucleus={canonical!r}, method={method!r}; known methods: {[r.get('method') for r in candidates]}")
        return {
            "nucleus": canonical,
            "values": [_stamp(r, tool_name="nuclear.get_charge_radius", args={"nucleus": nucleus, "method": r.get("method")}, source_revision=src_rev) for r in candidates],
            "_note": cr.get("_note"),
        }

    if cr.get("value") is None:
        return {**_stamp(cr, tool_name="nuclear.get_charge_radius", args=args_base, source_revision=src_rev), "nucleus": canonical}
    return _stamp(cr, tool_name="nuclear.get_charge_radius", args=args_base, source_revision=src_rev)


def get_magnetic_moment(nucleus: str) -> dict:
    db = _load()
    src_rev = db["$source_revision"]
    canonical, info = _resolve_nucleus(nucleus)
    mm = info.get("magnetic_moment")
    if mm is None:
        raise KeyError(f"no magnetic-moment data for nucleus={canonical!r}")
    return _stamp(mm, tool_name="nuclear.get_magnetic_moment", args={"nucleus": nucleus}, source_revision=src_rev)


def list_nuclei() -> list[dict]:
    db = _load()
    return [{"nucleus": canonical, "description": info.get("description", ""), "aliases": info.get("_aliases", [])} for canonical, info in db["nuclei"].items()]
