"""Smoke test: the umbrella server registers the expected pdg.* tools."""

import asyncio

from precision_data_mcp.server import mcp


EXPECTED_PDG_TOOLS = {
    "pdg.get_particle",
    "pdg.get_anomaly",
    "pdg.list_particles",
}


def test_pdg_namespace_tools_registered():
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    missing = EXPECTED_PDG_TOOLS - names
    assert not missing, f"missing pdg.* tools: {missing}"


def test_nist_namespace_still_present_after_pdg_addition():
    """Regression: adding the pdg.* namespace must not displace the nist.* one."""
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    assert "nist.get_constant" in names
    assert "nist.get_levels" in names
