"""PyPhysics MCP server — exposes the repo's bibliography, chapters, and
validators as Claude-callable MCP tools.

Entry point:
    pyphysics-mcp                  # registered console script
    python -m pyphysics_mcp        # equivalent
"""

from pyphysics_mcp.server import main

__all__ = ["main"]
