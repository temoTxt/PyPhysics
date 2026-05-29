"""NIST MCP server — exposes the NIST Atomic Spectra Database (ASD) and CODATA
fundamental physical constants as Claude-callable MCP tools, plus a curated,
individually-cited registry of the precision splittings used by the repo's
Dirac-equation exploration.

Entry point:
    nist-mcp                  # registered console script
    python -m nist_mcp        # equivalent
"""

from nist_mcp.server import main

__all__ = ["main"]
