"""Tests for the qed.* namespace lookup tools."""

import pytest

from precision_data_mcp.qed.tools import lookup


REQUIRED_FIELDS = {"value", "uncertainty", "unit", "source", "retrieved_at", "cache_key"}


def _assert_record(rec):
    missing = REQUIRED_FIELDS - set(rec.keys())
    assert not missing, f"missing umbrella contract fields: {missing}"
    assert isinstance(rec["source"], str) and rec["source"], "source must be a non-empty cite_key"


class TestHydrogen:
    def test_lamb_shift_multivalue_returns_both(self):
        """Post issue #108: H Lamb shift is multi-value (CODATA + direct Lundeen RF)."""
        r = lookup.get_lamb_shift("H", "2S1/2-2P1/2")
        assert "values" in r
        methods = {v.get("method") for v in r["values"]}
        assert "CODATA_2018_global_adjustment" in methods
        assert "direct_microwave_RF_between_2S_and_2P" in methods
        for v in r["values"]:
            _assert_record(v)
        assert all(v["value"] == pytest.approx(1057.845, rel=1e-6) for v in r["values"])

    def test_lamb_shift_method_selection_codata(self):
        r = lookup.get_lamb_shift("H", "2S1/2-2P1/2", method="CODATA_2018_global_adjustment")
        assert r["source"] == "codata2018_constants"
        _assert_record(r)

    def test_lamb_shift_method_selection_lundeen(self):
        r = lookup.get_lamb_shift("H", "2S1/2-2P1/2", method="direct_microwave_RF_between_2S_and_2P")
        assert r["source"] == "lundeen1981_2s2p_microwave"
        assert r["value_class"] == "experimental"
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
        # Multi-value entry; check first value matches
        r = lookup.get_lamb_shift("hydrogen", "2S1/2-2P1/2")
        assert r["values"][0]["value"] == pytest.approx(1057.845, rel=1e-6)

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
        a = lookup.get_lamb_shift("H", "2S1/2-2P1/2", method="CODATA_2018_global_adjustment")
        b = lookup.get_lamb_shift("H", "2S1/2-2P1/2", method="CODATA_2018_global_adjustment")
        assert a["cache_key"] == b["cache_key"]

    def test_cache_key_differs_by_args(self):
        a = lookup.get_lamb_shift("H", "2S1/2-2P1/2", method="CODATA_2018_global_adjustment")
        b = lookup.get_hyperfine("H", "1s2S1/2")
        assert a["cache_key"] != b["cache_key"]


class TestHydrogenPrecisionSpectroscopyExtension:
    """Issue #108: new H transitions seeded in qed/data.json's species.H.transitions."""

    def test_1s_3s_fleurbaey(self):
        r = lookup.get_transition_precision("H", "1S-3S")
        assert r["source"] == "fleurbaey2018_1s_3s_spectroscopy"
        assert r["value_class"] == "experimental"
        assert r["safe_for_model_verification"] is True
        _assert_record(r)

    def test_2s_4p_beyer(self):
        r = lookup.get_transition_precision("H", "2S-4P1/2")
        assert r["source"] == "beyer2017_2s4p"
        assert r["value_class"] == "experimental"

    def test_2s_8s_de_beauvoir(self):
        r = lookup.get_transition_precision("H", "2S-8S1/2")
        assert r["source"] == "debeauvoir1997_2s_8s_8d"
        assert r["value_class"] == "experimental"

    def test_h_d_isotope_shift_parthey(self):
        r = lookup.get_transition_precision("H", "1S-2S_H-D_isotope_shift")
        assert r["source"] == "parthey2010_h_d_isotope_shift"
        assert r["value_class"] == "experimental"

    def test_lundeen_lamb_shift_is_safe_alternative_to_codata(self):
        """Per issue #108: a Lundeen direct-RF Lamb shift is a safe-for-model-verification
        alternative to the codata_adjusted entry."""
        lundeen = lookup.get_lamb_shift("H", "2S1/2-2P1/2", method="direct_microwave_RF_between_2S_and_2P")
        assert lundeen["safe_for_model_verification"] is True
        codata = lookup.get_lamb_shift("H", "2S1/2-2P1/2", method="CODATA_2018_global_adjustment")
        assert codata["safe_for_model_verification"] is False
        # And their numerical values match (the codata uses lundeen as input)
        assert lundeen["value"] == codata["value"]
