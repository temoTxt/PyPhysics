"""Offline tests for get_constant (no network)."""

from nist_mcp.tools import codata

SCHEMA = {"quantity", "value", "uncertainty", "unit", "source", "reference"}


def test_electron_g_factor_ok():
    res = codata.get_constant("electron g factor")
    assert res["status"] == "ok"
    rec = res["match"]
    assert SCHEMA <= set(rec)
    assert rec["value"] != 0
    assert "CODATA" in rec["source"]
    assert rec["uncertainty"] is not None


def test_source_reports_codata_release():
    rec = codata.get_constant("Rydberg constant")["match"]
    assert "CODATA" in rec["source"]
    assert "scipy" in rec["source"]


def test_ambiguous_returns_stable_shape():
    res = codata.get_constant("mass")
    assert res["status"] == "ambiguous"
    assert len(res["candidates"]) > 1
    assert "match" not in res


def test_not_found():
    res = codata.get_constant("definitely not a real constant zzz")
    assert res["status"] == "not_found"
    assert res["candidates"] == []
