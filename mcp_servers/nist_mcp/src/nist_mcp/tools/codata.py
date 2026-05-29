"""CODATA fundamental physical constants via ``scipy.constants``.

scipy ships a specific CODATA release bundled with its version (scipy 1.11-1.13
-> CODATA 2018; scipy >= 1.15 -> CODATA 2022). The campaign cites CODATA-2018,
so ``get_constant`` reports the release it actually used in the ``source`` field
rather than letting a silent version bump change the numbers.
"""

import scipy
import scipy.constants as sc

_CUU_URL = "https://physics.nist.gov/cuu/Constants/"


def codata_release() -> str:
    """Best-effort label for the CODATA release scipy bundled (version included)."""
    ver = scipy.__version__
    try:
        major, minor = (int(x) for x in ver.split(".")[:2])
    except ValueError:
        return f"CODATA (scipy {ver})"
    if (major, minor) >= (1, 15):
        year = "2022"
    elif (major, minor) >= (1, 6):
        year = "2018"
    else:
        year = "2014"
    return f"CODATA {year} (scipy {ver})"


def _record(key: str) -> dict:
    value, unit, uncertainty = sc.physical_constants[key]
    return {
        "quantity": key,
        "value": value,
        "uncertainty": uncertainty,
        "unit": unit or "dimensionless",
        "source": codata_release(),
        "reference": _CUU_URL,
    }


def get_constant(name: str) -> dict:
    """Look up a CODATA constant by name.

    Returns one stable, ``status``-discriminated shape so a client never has to
    branch on record-vs-list:
      - ``{"status": "ok", "match": <record>}`` for a unique hit.
      - ``{"status": "ambiguous", "candidates": [keys...]}`` for >1 substring hit.
      - ``{"status": "not_found", "candidates": []}`` otherwise.
    """
    pcs = sc.physical_constants

    if name in pcs:
        return {"status": "ok", "match": _record(name)}

    lowered = {k.lower(): k for k in pcs}
    if name.lower() in lowered:
        return {"status": "ok", "match": _record(lowered[name.lower()])}

    candidates = sorted(k for k in pcs if name.lower() in k.lower())
    if len(candidates) == 1:
        return {"status": "ok", "match": _record(candidates[0])}
    if candidates:
        return {"status": "ambiguous", "candidates": candidates}
    return {"status": "not_found", "candidates": []}
