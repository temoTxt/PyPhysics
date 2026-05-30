"""NIST Atomic Spectra Database (ASD) retrieval.

Source-of-record: https://www.nist.gov/pml/atomic-spectra-database
The landing page links to the Lines and Levels query forms, which resolve to the
CGI endpoints below. Values are returned in ASD's NATIVE units (levels in cm^-1;
line wavelengths in nm, with an air/vacuum flag derived from the column name) —
there is no silent conversion. Converting to MHz/eV is the caller's job.

The CGI parameter sets were confirmed against live H I queries on 2026-05-29.
This is a scraped form, not a stable API: the committed fixtures under
tests/fixtures/ pin the response layout so a NIST form change surfaces as a test
failure rather than silent garbage. Query parameters are passed via requests'
``params=`` (URL-encoded), never string-interpolated into the URL.
"""

import urllib.parse

import requests

from precision_data_mcp.nist.cache import cache_get, cache_put

LANDING_PAGE = "https://www.nist.gov/pml/atomic-spectra-database"
LEVELS_URL = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl"
LINES_URL = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"
USER_AGENT = "PyPhysics-nist-mcp/0.1 (research; github.com/temoTxt/PyPhysics)"
TIMEOUT = 30


# ───── query parameter builders ─────

def _levels_params(species: str) -> dict:
    return {
        "de": 0, "spectrum": species, "units": 0, "format": 3, "output": 0,
        "page_size": 15, "multiplet_ordered": 0, "conf_out": "on",
        "term_out": "on", "level_out": "on", "unc_out": "on", "j_out": "on",
        "temp": "", "submit": "Retrieve Data",
    }


def _lines_params(species: str) -> dict:
    return {
        "spectra": species, "limits_type": 0, "low_w": "", "upp_w": "", "unit": 1,
        "submit": "Retrieve Data", "de": 0, "format": 3, "line_out": 0,
        "en_unit": 0, "output": 0, "bibrefs": 1, "show_obs_wl": 1,
        "show_calc_wl": 1, "unc_out": 1, "order_out": 0, "max_low_enrg": "",
        "show_av": 2, "max_upp_enrg": "", "tsb_value": 0, "min_str": "",
        "A_out": 0, "intens_out": "on", "max_str": "", "allowed_out": 1,
        "forbid_out": 1, "min_accur": "", "min_intens": "", "conf_out": "on",
        "term_out": "on", "enrg_out": "on", "J_out": "on", "g_out": "on",
    }


# ───── tab-delimited parsing ─────

def _split(text: str) -> list[list[str]]:
    rows = []
    for line in text.splitlines():
        if line.strip() == "":
            continue
        rows.append([c.strip().strip('"') for c in line.split("\t")])
    return rows


def _parse_header(text: str) -> list[str]:
    rows = _split(text)
    return rows[0] if rows else []


def _parse_tsv(text: str) -> list[dict]:
    """Return data rows as dicts keyed by the header columns."""
    rows = _split(text)
    if len(rows) < 2:
        return []
    header = rows[0]
    return [dict(zip(header, r)) for r in rows[1:]]


def _to_float(s: str | None) -> float | None:
    if s is None:
        return None
    s = s.strip().strip("[]()").strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


# ───── record mapping ─────

def _level_record(row: dict) -> dict:
    conf = row.get("Configuration", "")
    term = row.get("Term", "")
    j = row.get("J", "")
    label = " ".join(p for p in (conf, term, j) if p)
    return {
        "quantity": f"energy level {label}".strip(),
        "value": _to_float(row.get("Level (cm-1)")),
        "uncertainty": _to_float(row.get("Uncertainty (cm-1)")),
        "unit": "cm^-1",
        "source": "NIST ASD (levels)",
        "reference": LANDING_PAGE,
    }


def _wl_columns(header: list[str]) -> tuple[str | None, str | None]:
    obs = next((h for h in header if h.startswith("obs_wl")), None)
    ritz = next((h for h in header if h.startswith("ritz_wl")), None)
    return obs, ritz


def _medium(col: str | None) -> str:
    if col and "vac" in col:
        return "vacuum"
    if col and "air" in col:
        return "air"
    return "unknown"


def _transition_record(row: dict, obs_col, ritz_col, ref_col) -> dict:
    obs = _to_float(row.get(obs_col)) if obs_col else None
    ritz = _to_float(row.get(ritz_col)) if ritz_col else None
    if obs is not None:
        value, unc, kind, medium = obs, _to_float(row.get("unc_obs_wl")), "observed", _medium(obs_col)
    else:
        value, unc, kind, medium = ritz, _to_float(row.get("unc_ritz_wl")), "ritz", _medium(ritz_col)
    conf_i, term_i = row.get("conf_i", ""), row.get("term_i", "")
    conf_k, term_k = row.get("conf_k", ""), row.get("term_k", "")
    label = f"{conf_i} {term_i} - {conf_k} {term_k}".strip()
    ref = row.get(ref_col) if ref_col else ""
    return {
        "quantity": f"{kind} wavelength {label}".strip(),
        "value": value,
        "uncertainty": unc,
        "unit": f"nm ({medium})",
        "source": "NIST ASD (lines)",
        "reference": ref or LANDING_PAGE,
    }


def _transition_records(text: str) -> list[dict]:
    header = _parse_header(text)
    obs_col, ritz_col = _wl_columns(header)
    ref_col = "line_ref" if "line_ref" in header else ("tp_ref" if "tp_ref" in header else None)
    recs = [_transition_record(r, obs_col, ritz_col, ref_col) for r in _parse_tsv(text)]
    return [r for r in recs if r["value"] is not None]


# ───── networked retrieval (cached) ─────

def _fetch(url: str, params: dict, refresh: bool, ttl: float | None) -> tuple[dict | None, str | None]:
    query = url + "?" + urllib.parse.urlencode(sorted(params.items()))
    if not refresh:
        cached = cache_get(query, ttl=ttl)
        if cached is not None:
            return cached, None
    try:
        resp = requests.get(url, params=params, headers={"User-Agent": USER_AGENT}, timeout=TIMEOUT)
        resp.raise_for_status()
    except requests.RequestException as exc:
        return None, str(exc)
    payload = {"text": resp.text}
    cache_put(query, payload)
    return payload, None


def get_levels(species: str, refresh: bool = False, ttl: float | None = None) -> dict:
    """ASD energy levels for ``species`` (e.g. 'H I', 'He I', 'Li III'). Native unit: cm^-1."""
    payload, err = _fetch(LEVELS_URL, _levels_params(species), refresh, ttl)
    base = {"source": "NIST ASD (levels)", "query": {"species": species, "endpoint": LEVELS_URL}}
    if err:
        return {**base, "results": [], "error": err}
    results = [r for r in (_level_record(r) for r in _parse_tsv(payload["text"])) if r["value"] is not None]
    return {**base, "results": results}


def get_transitions(species: str, refresh: bool = False, ttl: float | None = None) -> dict:
    """ASD transitions for ``species``. Native unit: nm, with an air/vacuum flag in ``unit``."""
    payload, err = _fetch(LINES_URL, _lines_params(species), refresh, ttl)
    base = {"source": "NIST ASD (lines)", "query": {"species": species, "endpoint": LINES_URL}}
    if err:
        return {**base, "results": [], "error": err}
    return {**base, "results": _transition_records(payload["text"])}
