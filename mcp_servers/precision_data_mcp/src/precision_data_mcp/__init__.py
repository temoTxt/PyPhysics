"""precision_data_mcp — umbrella MCP server for precision-physics experimental data.

Exposes namespaced tools across PDG, NIST ASD + CODATA, bound-state QED,
nuclear / particle structure, FLAG lattice-QCD, and (deferred) astronomical /
cosmological observables. Every tool returns ``{value, uncertainty, unit, source,
retrieved_at, cache_key}`` where ``source`` is a ``cite_key`` resolving against
``Roadmapping/History/Bibliography/``.

Tracks issue [#92](https://github.com/temoTxt/PyPhysics/issues/92). See
``REFRESH_POLICY.md`` for the shared caching + refresh schema across namespaces.

Entry point:
    precision-data-mcp           # registered console script
    python -m precision_data_mcp # equivalent
"""

from precision_data_mcp.server import main

__all__ = ["main"]
