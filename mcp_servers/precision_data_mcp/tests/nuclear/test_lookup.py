"""Tests for the nuclear.* namespace, with special attention to the proton-radius
puzzle as the umbrella's disagreement-representation demonstrator."""

import pytest

from precision_data_mcp.nuclear.tools import lookup


REQUIRED_FIELDS = {"value", "unit", "source", "retrieved_at", "cache_key"}


class TestProtonRadiusPuzzle:
    """This class is the load-bearing acceptance test for issue #98 — the
    proton-radius puzzle must surface as multiple values with method labels."""

    def test_returns_multiple_values_when_method_unspecified(self):
        r = lookup.get_charge_radius("proton")
        assert "values" in r, "proton charge radius must be multi-valued — that IS the puzzle"
        assert len(r["values"]) >= 3, f"expected >=3 puzzle values, got {len(r['values'])}"

    def test_all_puzzle_values_have_method_labels(self):
        r = lookup.get_charge_radius("proton")
        for v in r["values"]:
            assert v.get("method"), f"every puzzle value must carry a method label, missing on {v}"

    def test_canonical_puzzle_methods_present(self):
        r = lookup.get_charge_radius("proton")
        methods = {v["method"] for v in r["values"]}
        assert "muonic_hydrogen_spectroscopy" in methods, "Pohl 2010 muonic-H must be in the seed"
        assert "electron_scattering" in methods, "Bernauer 2010 A1 Mainz must be in the seed"

    def test_muonic_h_value_matches_pohl_2010(self):
        r = lookup.get_charge_radius("proton", method="muonic_hydrogen_spectroscopy")
        assert r["value"] == pytest.approx(0.84087, rel=1e-4)
        assert r["source"] == "pohl2010_muonic_hydrogen"

    def test_electron_scattering_value_matches_bernauer_2010(self):
        r = lookup.get_charge_radius("proton", method="electron_scattering")
        assert r["value"] == pytest.approx(0.879, rel=1e-2)
        assert r["source"] == "bernauer2010_a1_mainz"

    def test_disagreement_at_5_sigma(self):
        """The whole point of the puzzle: the values must disagree significantly."""
        all_values = lookup.get_charge_radius("proton")["values"]
        muonic = next(v for v in all_values if v["method"] == "muonic_hydrogen_spectroscopy")
        scattering = next(v for v in all_values if v["method"] == "electron_scattering")
        gap = abs(muonic["value"] - scattering["value"])
        combined_uncert = (muonic["uncertainty"]**2 + scattering["uncertainty"]**2) ** 0.5
        sigma = gap / combined_uncert
        assert sigma > 4.0, f"puzzle should be at >=4 sigma; got {sigma:.1f}"

    def test_unknown_method_raises(self):
        with pytest.raises(KeyError):
            lookup.get_charge_radius("proton", method="psychic_intuition")


class TestDeuteronAnalogue:
    def test_deuteron_charge_radius_also_multi_valued(self):
        r = lookup.get_charge_radius("deuteron")
        assert "values" in r
        assert len(r["values"]) >= 2


class TestNeutronSpecialCase:
    def test_neutron_has_mean_square_radius_not_simple_radius(self):
        r = lookup.get_charge_radius("neutron")
        assert r.get("value") is None  # neutron is electrically neutral; single radius not defined


class TestMagneticMoments:
    def test_proton_magnetic_moment(self):
        m = lookup.get_magnetic_moment("proton")
        assert m["value"] == pytest.approx(2.79284734463, rel=1e-9)
        assert m["unit"] == "mu_N"

    def test_neutron_magnetic_moment_is_negative(self):
        m = lookup.get_magnetic_moment("neutron")
        assert m["value"] < 0


class TestAliasResolution:
    def test_p_plus_alias(self):
        r = lookup.get_charge_radius("p+", method="muonic_hydrogen_spectroscopy")
        assert r["value"] == pytest.approx(0.84087, rel=1e-4)
