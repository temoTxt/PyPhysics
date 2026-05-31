"""Tests for the qed.* namespace lookup tools."""

import pytest

from precision_data_mcp.qed.tools import lookup


REQUIRED_FIELDS = {"value", "uncertainty", "unit", "source", "retrieved_at", "cache_key"}


def _assert_record(rec):
    missing = REQUIRED_FIELDS - set(rec.keys())
    assert not missing, f"missing umbrella contract fields: {missing}"
    assert isinstance(rec["source"], str) and rec["source"], "source must be a non-empty cite_key"


class TestHydrogen:
    def test_lamb_shift_canonical_value(self):
        r = lookup.get_lamb_shift("H", "2S1/2-2P1/2")
        assert r["value"] == pytest.approx(1057.845, rel=1e-6)
        assert r["unit"] == "MHz"
        assert r["source"] == "codata2018_constants"
        _assert_record(r)

    def test_hyperfine_21cm(self):
        r = lookup.get_hyperfine("H", "1s2S1/2")
        assert r["value"] == pytest.approx(1420.405751768, rel=1e-12)
        assert r["unit"] == "MHz"

    def test_1s_2s_transition_hz(self):
        r = lookup.get_transition_precision("H", "1S-2S")
        assert r["value"] == 2466061413187035
        assert r["unit"] == "Hz"


class TestHydrogenicIons:
    def test_si13_g_factor(self):
        r = lookup.get_g_factor("Si XIV", "1s2S1/2")
        assert r["value"] == pytest.approx(1.9953395931, rel=1e-9)
        assert r["source"] == "sturm2013_si13_g_factor"

    def test_he_g_factor(self):
        r = lookup.get_g_factor("He II", "1s2S1/2")
        assert r["value"] == pytest.approx(1.9999857940, rel=1e-9)

    def test_li2_hyperfine(self):
        r = lookup.get_hyperfine("Li III", "1s2S1/2")
        assert r["value"] == pytest.approx(11890.018, rel=1e-6)


class TestExoticAtoms:
    def test_muonic_hydrogen_lamb_shift(self):
        r = lookup.get_lamb_shift("muonic_H", "2S1/2-2P3/2")
        assert r["source"] == "pohl2010_muonic_hydrogen"
        assert r["unit"] == "GHz"

    def test_muonium_hyperfine(self):
        r = lookup.get_hyperfine("muonium", "1s2S1/2")
        assert r["value"] == pytest.approx(4463.302776, rel=1e-9)

    def test_positronium_hyperfine(self):
        r = lookup.get_hyperfine("positronium", "1s_ortho_para")
        assert r["value"] == 203389


class TestAliasResolution:
    def test_hydrogen_alias(self):
        r = lookup.get_lamb_shift("hydrogen", "2S1/2-2P1/2")
        assert r["value"] == pytest.approx(1057.845, rel=1e-6)

    def test_he_plus_alias(self):
        r = lookup.get_g_factor("He+", "1s2S1/2")
        assert r["value"] == pytest.approx(1.9999857940, rel=1e-9)


class TestErrorPaths:
    def test_unknown_species(self):
        with pytest.raises(KeyError):
            lookup.get_lamb_shift("xenon", "1s2S1/2")

    def test_unknown_observable_for_known_species(self):
        with pytest.raises(KeyError):
            lookup.get_lamb_shift("H", "99S-100P")


class TestEnumeration:
    def test_list_species_includes_all_seven(self):
        ss = {s["species"] for s in lookup.list_species()}
        assert {"H", "He II", "Li III", "Si XIV", "muonic_H", "muonium", "positronium"} == ss

    def test_list_observables_for_H(self):
        o = lookup.list_observables("H")
        assert o["species"] == "H"
        assert "2S1/2-2P1/2" in o["lamb_shifts"]
        assert "1s2S1/2" in o["hyperfine"]
        assert "1S-2S" in o["transitions"]


class TestProvenance:
    def test_cache_key_deterministic_per_arg_tuple(self):
        a = lookup.get_lamb_shift("H", "2S1/2-2P1/2")
        b = lookup.get_lamb_shift("H", "2S1/2-2P1/2")
        assert a["cache_key"] == b["cache_key"]

    def test_cache_key_differs_by_args(self):
        a = lookup.get_lamb_shift("H", "2S1/2-2P1/2")
        b = lookup.get_hyperfine("H", "1s2S1/2")
        assert a["cache_key"] != b["cache_key"]
