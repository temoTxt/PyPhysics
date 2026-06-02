"""Offline test guarding the scope boundary: the precision splittings are curated
and cited, never silently re-routed through ASD."""

from precision_data_mcp.nist.tools import targets

SCHEMA = {"quantity", "value", "uncertainty", "unit", "source", "reference", "value_class", "safe_for_model_verification"}


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


def test_lamb_shift_target_flagged_as_codata_adjusted():
    """The H 2S-2P Lamb shift transcribed here is the CODATA-2018 value, which combines
    experiment + theory. Must NOT be safe_for_model_verification."""
    items = targets.list_dirac_targets()
    lamb = next(t for t in items if "Lamb shift" in t["quantity"])
    assert lamb["value_class"] == "codata_adjusted"
    assert lamb["safe_for_model_verification"] is False
    assert "warning" in lamb


def test_direct_measurements_flagged_as_experimental():
    """The 21-cm maser, He fine structure, positronium and muonium hyperfine values
    are direct experimental measurements — must be safe_for_model_verification."""
    items = targets.list_dirac_targets()
    direct_measurements = [t for t in items if t["quantity"] not in ("hydrogen 2S_1/2 - 2P_1/2 Lamb shift",)]
    for t in direct_measurements:
        assert t["value_class"] == "experimental", f"{t['quantity']} should be experimental, got {t['value_class']}"
        assert t["safe_for_model_verification"] is True
        assert "warning" not in t
