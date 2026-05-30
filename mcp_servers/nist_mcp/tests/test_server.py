"""Deterministic smoke test: the server registers exactly the four expected tools.

No stdio invocation (which can hang), so this is safe in CI.
"""

import asyncio

from nist_mcp.server import mcp


def test_four_tools_registered():
    tools = asyncio.run(mcp.list_tools())
    names = sorted(t.name for t in tools)
    assert names == ["get_constant", "get_levels", "get_transitions", "list_dirac_targets"]
