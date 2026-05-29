"""ASD parser tests against committed fixtures (offline) + a live-gated test."""

import os
from pathlib import Path

import pytest

from nist_mcp.tools import asd

FIX = Path(__file__).parent / "fixtures"
SCHEMA = {"quantity", "value", "uncertainty", "unit", "source", "reference"}


def test_parse_levels_fixture_native_units_and_uncertainty():
    text = (FIX / "asd_levels_sample.tsv").read_text(encoding="utf-8")
    recs = [r for r in (asd._level_record(r) for r in asd._parse_tsv(text)) if r["value"] is not None]
    assert recs
    assert all(SCHEMA <= set(r) for r in recs)
    assert all(r["unit"] == "cm^-1" for r in recs)          # native unit, no conversion
    # the H I 2p 2P* 3/2 level is present, with an uncertainty
    assert any(abs(r["value"] - 82259.2850014) < 1e-3 for r in recs)
    assert all(r["uncertainty"] is not None for r in recs[:3])


def test_parse_lines_fixture_air_vacuum_flag():
    text = (FIX / "asd_lines_sample.tsv").read_text(encoding="utf-8")
    header = asd._parse_header(text)
    obs_col, ritz_col = asd._wl_columns(header)
    assert asd._medium(obs_col or ritz_col) in {"vacuum", "air"}
    recs = asd._transition_records(text)
    assert recs
    assert all(SCHEMA <= set(r) for r in recs)
    assert all("nm" in r["unit"] for r in recs)


@pytest.mark.skipif(os.environ.get("NIST_MCP_LIVE") != "1", reason="live network test; set NIST_MCP_LIVE=1")
def test_live_levels_h_native_units():
    out = asd.get_levels("H I", refresh=True)
    assert out["results"], "live ASD returned no levels"
    assert all(r["unit"] == "cm^-1" for r in out["results"])
    # ASD coverage is asserted here, NOT the MHz splittings (ASD does not supply those).
