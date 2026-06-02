"""Tests for the flag.* namespace (FLAG 2024 lattice-QCD averages)."""

import asyncio

import pytest

from precision_data_mcp.flag.tools import lookup
from precision_data_mcp.server import mcp


REQUIRED_FIELDS = {"value", "uncertainty", "unit", "source", "retrieved_at", "cache_key"}


def test_f_pi_smoke():
    r = lookup.get_quantity("f_pi")
    assert r["value"] == pytest.approx(130.2, rel=1e-3)
    assert r["unit"] == "MeV"
    assert r["source"] == "flag2024_averages"
    assert r["source_revision"]["flag_edition"] == "2024"
    assert REQUIRED_FIELDS.issubset(set(r.keys()))


def test_g_A_present():
    r = lookup.get_quantity("g_A")
    assert r["value"] == pytest.approx(1.246, rel=1e-3)


def test_list_quantities_includes_canonical():
    qs = {q["quantity_id"] for q in lookup.list_quantities()}
    assert {"f_pi", "f_K", "g_A", "sigma_piN", "f_K_over_f_pi"}.issubset(qs)


def test_unknown_raises():
    with pytest.raises(KeyError):
        lookup.get_quantity("xyzzy_quantity")


class TestServerRegistration:
    def test_nuclear_and_flag_tools_registered(self):
        tools = asyncio.run(mcp.list_tools())
        names = {t.name for t in tools}
        expected = {
            "nuclear.get_charge_radius",
            "nuclear.get_magnetic_moment",
            "nuclear.list_nuclei",
            "flag.get_quantity",
            "flag.list_quantities",
        }
        missing = expected - names
        assert not missing, f"missing nuclear.*/flag.* tools: {missing}"

    def test_earlier_namespaces_still_present(self):
        tools = asyncio.run(mcp.list_tools())
        names = {t.name for t in tools}
        assert "nist.get_constant" in names
        assert "pdg.get_particle" in names
        assert "qed.get_lamb_shift" in names
