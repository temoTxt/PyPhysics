"""Cross-namespace safety contract tests.

Every value-bearing return across every namespace must carry value_class +
safe_for_model_verification + (when unsafe) a warning field. These tests
exercise the contract across qed.*, pdg.*, nuclear.*, flag.*, and the curated
nist.* sub-tools — the umbrella REFRESH_POLICY.md §1 guarantee 6.
"""

import pytest

from precision_data_mcp import safety
from precision_data_mcp.flag.tools import lookup as flag_lookup
from precision_data_mcp.nist.tools import codata as nist_codata
from precision_data_mcp.nist.tools import targets as nist_targets
from precision_data_mcp.nuclear.tools import lookup as nuclear_lookup
from precision_data_mcp.pdg.tools import lookup as pdg_lookup
from precision_data_mcp.qed.tools import lookup as qed_lookup


def _assert_contract(record: dict):
    assert "value_class" in record, f"missing value_class on {record}"
    assert record["value_class"] in safety.ALL_VALUE_CLASSES, f"unknown value_class: {record['value_class']}"
    assert "safe_for_model_verification" in record
    if record["safe_for_model_verification"] is False:
        assert "warning" in record, f"unsafe record must carry a warning: {record}"
    else:
        # safe records must not falsely advertise a warning
        assert "warning" not in record


class TestQedSafetyContract:
    def test_h_lamb_shift_multivalue_each_carries_contract(self):
        """Per issue #108: H Lamb shift is multi-value (CODATA + Lundeen direct RF)."""
        r = qed_lookup.get_lamb_shift("H", "2S1/2-2P1/2")
        for v in r["values"]:
            _assert_contract(v)
        # CODATA entry unsafe, Lundeen entry safe — both surfaced
        classes = {v["value_class"] for v in r["values"]}
        assert "codata_adjusted" in classes
        assert "experimental" in classes

    def test_h_hyperfine_safe_experimental(self):
        r = qed_lookup.get_hyperfine("H", "1s2S1/2")
        _assert_contract(r)
        assert r["safe_for_model_verification"] is True

    def test_h_1s_2s_safe_experimental(self):
        r = qed_lookup.get_transition_precision("H", "1S-2S")
        _assert_contract(r)
        assert r["safe_for_model_verification"] is True

    def test_muonic_h_safe_experimental(self):
        r = qed_lookup.get_lamb_shift("muonic_H", "2S1/2-2P3/2")
        _assert_contract(r)
        assert r["safe_for_model_verification"] is True

    def test_si13_g_factor_safe(self):
        r = qed_lookup.get_g_factor("Si XIV", "1s2S1/2")
        _assert_contract(r)
        assert r["safe_for_model_verification"] is True


class TestPdgSafetyContract:
    def test_muon_mass_safe_world_average(self):
        p = pdg_lookup.get_particle("muon")
        _assert_contract(p["mass"])
        assert p["mass"]["safe_for_model_verification"] is True

    def test_proton_mass_unsafe_codata(self):
        p = pdg_lookup.get_particle("proton")
        _assert_contract(p["mass"])
        assert p["mass"]["safe_for_model_verification"] is False

    def test_definitional_fields_safe(self):
        p = pdg_lookup.get_particle("electron")
        for field in ("charge", "spin"):
            _assert_contract(p[field])
            assert p[field]["safe_for_model_verification"] is True

    def test_muon_anomaly_mixed_classes(self):
        """a_mu values include experimental, lattice_qcd, and sm_prediction — different safety."""
        a = pdg_lookup.get_anomaly("muon")
        classes = {v["value_class"] for v in a["values"]}
        assert "experimental_world_average" in classes
        assert "theoretical_lattice_qcd" in classes
        assert "theoretical_sm_prediction" in classes
        # Each individual value must carry the contract
        for v in a["values"]:
            _assert_contract(v)

    def test_neutron_lifetime_multivalue_each_safe(self):
        n = pdg_lookup.get_particle("neutron")
        for v in n["lifetime"]["values"]:
            _assert_contract(v)


class TestNuclearSafetyContract:
    def test_proton_radius_puzzle_each_value_carries_contract(self):
        r = nuclear_lookup.get_charge_radius("proton")
        for v in r["values"]:
            _assert_contract(v)
        # Some are experimental (Pohl, Bernauer, Fleurbaey, PRad); CODATA pre-puzzle is codata_adjusted
        classes = {v["value_class"] for v in r["values"]}
        assert "experimental" in classes
        assert "codata_adjusted" in classes

    def test_proton_magnetic_moment_unsafe_codata(self):
        r = nuclear_lookup.get_magnetic_moment("proton")
        _assert_contract(r)
        assert r["safe_for_model_verification"] is False


class TestFlagSafetyContract:
    def test_f_pi_unsafe_lattice_qcd(self):
        r = flag_lookup.get_quantity("f_pi")
        _assert_contract(r)
        assert r["value_class"] == "theoretical_lattice_qcd"
        assert r["safe_for_model_verification"] is False

    def test_m_pi_charged_safe_experimental(self):
        """FLAG quotes the PDG experimental pion mass as input, not as a lattice calculation."""
        r = flag_lookup.get_quantity("m_pi_charged")
        _assert_contract(r)
        assert r["safe_for_model_verification"] is True


class TestNistSafetyContract:
    def test_codata_get_constant_unsafe(self):
        r = nist_codata.get_constant("electron g factor")["match"]
        _assert_contract(r)
        assert r["safe_for_model_verification"] is False

    def test_dirac_targets_mixed(self):
        items = nist_targets.list_dirac_targets()
        for t in items:
            _assert_contract(t)


class TestSafetyHelperItself:
    def test_missing_value_class_raises(self):
        with pytest.raises(ValueError, match="missing required value_class"):
            safety.apply_safety_contract({"value": 1.0, "uncertainty": 0.1})

    def test_unknown_value_class_raises(self):
        with pytest.raises(ValueError, match="unknown value_class"):
            safety.apply_safety_contract({"value": 1.0, "value_class": "made_up"})

    def test_safe_classes_get_no_warning(self):
        out = safety.apply_safety_contract({"value": 1, "value_class": "experimental"})
        assert out["safe_for_model_verification"] is True
        assert "warning" not in out

    def test_unsafe_classes_get_warning(self):
        out = safety.apply_safety_contract({"value": 1, "value_class": "codata_adjusted"})
        assert out["safe_for_model_verification"] is False
        assert "warning" in out
        assert "CODATA" in out["warning"]
