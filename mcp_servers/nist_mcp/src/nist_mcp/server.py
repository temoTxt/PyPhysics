"""FastMCP server registering the four NIST tools.

Run as:
    nist-mcp                  # via console_script entry point
    python -m nist_mcp        # equivalent

The server speaks MCP over stdio. Register it with Claude Code by adding to your
MCP config (``.mcp.json`` or ``~/.claude/mcp.json``):

    {
      "mcpServers": {
        "nist": {"command": "nist-mcp"}
      }
    }

Tools (4) and their data sources:
  get_constant        -> CODATA via scipy.constants (offline)
  get_levels          -> NIST ASD energy levels (cm^-1)
  get_transitions     -> NIST ASD transitions (nm, air/vacuum flagged)
  list_dirac_targets  -> curated, cited precision splittings (NOT from ASD)
"""

from mcp.server.fastmcp import FastMCP

from nist_mcp.tools import asd, codata, targets

mcp = FastMCP("nist")


@mcp.tool()
def get_constant(name: str) -> dict:
    """CODATA fundamental physical constant by name (offline, via scipy.constants).

    Returns a status-discriminated dict: {"status": "ok", "match": record} for a
    unique hit, {"status": "ambiguous", "candidates": [...]} for several substring
    hits, or {"status": "not_found", "candidates": []}. The record's ``source``
    names the CODATA release used (scipy bundles a specific release).

    Args:
        name: constant name, e.g. "electron g factor", "Rydberg constant".
    """
    return codata.get_constant(name)


@mcp.tool()
def get_levels(species: str, refresh: bool = False) -> dict:
    """NIST ASD energy levels for a species. Native unit: cm^-1 (no conversion).

    Source: NIST Atomic Spectra Database (https://www.nist.gov/pml/atomic-spectra-database).
    Note: for light elements (e.g. hydrogen) ASD levels are largely theoretical/Ritz
    values, not independent measurements.

    Args:
        species: spectrum label, e.g. "H I", "He I", "Li III".
        refresh: bypass the on-disk cache and re-fetch.
    """
    return asd.get_levels(species, refresh=refresh)


@mcp.tool()
def get_transitions(species: str, refresh: bool = False) -> dict:
    """NIST ASD transitions for a species. Native unit: nm, with an air/vacuum flag.

    Source: NIST Atomic Spectra Database (https://www.nist.gov/pml/atomic-spectra-database).

    Args:
        species: spectrum label, e.g. "H I", "He I", "Li III".
        refresh: bypass the on-disk cache and re-fetch.
    """
    return asd.get_transitions(species, refresh=refresh)


@mcp.tool()
def list_dirac_targets() -> list[dict]:
    """Curated, individually-cited precision splittings for the Dirac-equation exploration.

    These MHz/sub-MHz values (hydrogen fine structure, Lamb shift, 1S hyperfine,
    helium 3P fine structure) are NOT sourced from ASD — ASD does not carry them
    at this precision. Each record cites its provenance. The electron g-factor is
    delegated to get_constant (CODATA).
    """
    return targets.list_dirac_targets()


def main() -> None:
    """stdio MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
