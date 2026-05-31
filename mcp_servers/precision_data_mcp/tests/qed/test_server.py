"""Smoke test: the umbrella server registers the expected qed.* tools."""

import asyncio

from precision_data_mcp.server import mcp


EXPECTED_QED_TOOLS = {
    "qed.get_lamb_shift",
    "qed.get_hyperfine",
    "qed.get_g_factor",
    "qed.get_transition_precision",
    "qed.list_species",
    "qed.list_observables",
}


def test_qed_namespace_tools_registered():
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    missing = EXPECTED_QED_TOOLS - names
    assert not missing, f"missing qed.* tools: {missing}"


def test_pdg_and_nist_still_present():
    """Regression: adding qed.* must not displace earlier namespaces."""
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    assert "nist.get_constant" in names
    assert "pdg.get_particle" in names
