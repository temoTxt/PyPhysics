"""Deterministic smoke test: the umbrella server registers the four expected nist.* tools.

No stdio invocation (which can hang), so this is safe in CI. Tests against the umbrella
``precision_data_mcp.server`` rather than against a per-namespace server, since per the
[#99](https://github.com/temoTxt/PyPhysics/issues/99) migration there is one umbrella
runtime registering all namespaced tools.
"""

import asyncio

from precision_data_mcp.server import mcp


EXPECTED_NIST_TOOLS = {
    "nist.get_constant",
    "nist.get_levels",
    "nist.get_transitions",
    "nist.list_dirac_targets",
}


def test_nist_namespace_tools_registered():
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    missing = EXPECTED_NIST_TOOLS - names
    assert not missing, f"missing nist.* tools: {missing}"


def test_no_unprefixed_nist_tools():
    """Regression: the un-prefixed names from the original standalone nist_mcp must NOT
    survive the migration. If they appear here, the umbrella server is double-registering."""
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    forbidden = {"get_constant", "get_levels", "get_transitions", "list_dirac_targets"}
    leaked = forbidden & names
    assert not leaked, f"un-prefixed tool names leaked from the standalone nist_mcp: {leaked}"
