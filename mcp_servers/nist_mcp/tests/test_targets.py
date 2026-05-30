"""Offline test guarding the scope boundary: the precision splittings are curated
and cited, never silently re-routed through ASD."""

from nist_mcp.tools import targets

SCHEMA = {"quantity", "value", "uncertainty", "unit", "source", "reference"}


def test_targets_curated_and_cited():
    items = targets.list_dirac_targets()
    assert len(items) >= 4
    for t in items:
        assert SCHEMA <= set(t)
        assert t["reference"], "every curated target must carry a citation"
        assert "curated" in t["source"].lower()
        assert t["unit"] == "MHz"
        assert t["value"] is not None


def test_positronium_and_muonium_present():
    quantities = " ".join(t["quantity"].lower() for t in targets.list_dirac_targets())
    assert "positronium" in quantities
    assert "muonium" in quantities
