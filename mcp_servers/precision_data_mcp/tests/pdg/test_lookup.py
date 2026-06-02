"""Tests for the pdg.* namespace lookup tools.

Verifies (a) basic single-value retrieval, (b) cite-key discipline (every source
field is a `cite_key` string, never a free-form URL/DOI), (c) the umbrella's
disagreement-representation rule for the neutron-lifetime and muon-g-2 cases,
and (d) every returned record carries the umbrella's contract fields.
"""

import pytest

from precision_data_mcp.pdg.tools import lookup


REQUIRED_RECORD_FIELDS = {"value", "unit", "source", "retrieved_at", "cache_key"}


def _assert_record_has_contract_fields(record: dict, *, allow_none_value: bool = False):
    """A quantity record must carry the umbrella's six contract fields."""
    missing = REQUIRED_RECORD_FIELDS - set(record.keys())
    assert not missing, f"record missing umbrella contract fields: {missing}"
    if not allow_none_value:
        assert record["value"] is not None
    assert isinstance(record["source"], str) and record["source"], "source must be a non-empty cite_key string"


class TestSinglePerticleRetrieval:
    def test_electron_by_mcid(self):
        e = lookup.get_particle("11")
        assert e["name"] == "electron"
        assert e["symbol"] == "e-"
        _assert_record_has_contract_fields(e["mass"])
        assert e["mass"]["value"] == pytest.approx(0.51099895069, rel=1e-9)
        assert e["mass"]["unit"] == "MeV/c^2"

    def test_electron_by_name(self):
        e = lookup.get_particle("electron")
        assert e["particle_id"] == "11"

    def test_muon_by_symbol(self):
        mu = lookup.get_particle("mu-")
        assert mu["name"] == "muon"
        _assert_record_has_contract_fields(mu["lifetime"])
        assert mu["lifetime"]["value"] == pytest.approx(2.1969811e-6, rel=1e-6)

    def test_proton_mass_cites_codata_not_pdg(self):
        """Proton mass is more precisely known via CODATA than via the PDG average."""
        p = lookup.get_particle("proton")
        assert p["mass"]["source"] == "codata2018_constants"

    def test_unknown_raises(self):
        with pytest.raises(KeyError):
            lookup.get_particle("xyzzy")


class TestDisagreementRepresentation:
    """Per umbrella issue #92 §Resolved-decisions #5, multi-value disagreements
    return all values with method labels — never an averaged single value."""

    def test_neutron_lifetime_returns_both_bottle_and_beam(self):
        n = lookup.get_particle("neutron")
        lifetime = n["lifetime"]
        assert "values" in lifetime, "neutron lifetime must be multi-value (bottle vs beam puzzle)"
        methods = {v["method"] for v in lifetime["values"]}
        assert methods == {"ultracold_neutron_bottle", "beam_measurement"}
        for v in lifetime["values"]:
            _assert_record_has_contract_fields(v)

    def test_muon_g_minus_2_returns_all_three(self):
        anomaly = lookup.get_anomaly("muon")
        assert "values" in anomaly
        methods = {v["method"] for v in anomaly["values"]}
        assert "BNL_E821_plus_Fermilab_E989_combined" in methods
        assert "lattice_QCD_hadronic_vacuum_polarisation" in methods
        assert "data_driven_dispersive_HVP" in methods
        for v in anomaly["values"]:
            _assert_record_has_contract_fields(v)

    def test_electron_g_minus_2_returns_experiment_and_sm(self):
        anomaly = lookup.get_anomaly("electron")
        assert "values" in anomaly
        methods = {v["method"] for v in anomaly["values"]}
        assert "Penning_trap_measurement" in methods
        assert "Standard_Model_prediction" in methods


class TestCacheKeyAndProvenance:
    def test_cache_key_is_deterministic(self):
        a = lookup.get_particle("electron")
        b = lookup.get_particle("electron")
        assert a["mass"]["cache_key"] == b["mass"]["cache_key"]

    def test_cache_key_differs_per_particle(self):
        e = lookup.get_particle("electron")
        mu = lookup.get_particle("muon")
        assert e["mass"]["cache_key"] != mu["mass"]["cache_key"]

    def test_retrieved_at_present_on_every_record(self):
        e = lookup.get_particle("electron")
        for field in ("mass", "charge", "spin"):
            assert e[field]["retrieved_at"]

    def test_source_revision_present(self):
        e = lookup.get_particle("electron")
        assert e["mass"]["source_revision"]["pdg_edition"] == "2024"


class TestListParticles:
    def test_list_includes_canonical_leptons_and_hadrons(self):
        ps = lookup.list_particles()
        names = {p["name"] for p in ps}
        assert {"electron", "muon", "tau", "proton", "neutron", "Higgs boson"}.issubset(names)
