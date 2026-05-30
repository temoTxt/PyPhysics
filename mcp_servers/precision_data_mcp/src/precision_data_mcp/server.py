"""FastMCP server registering the umbrella precision_data tool namespaces.

Run as:
    precision-data-mcp           # via console_script entry point
    python -m precision_data_mcp # equivalent

The server speaks MCP over stdio. Register it with Claude Code by adding to your
MCP config (``.mcp.json`` or ``~/.claude/mcp.json``):

    {
      "mcpServers": {
        "precision_data": {"command": "precision-data-mcp"}
      }
    }

This server is the umbrella tracked by issue [#92](https://github.com/temoTxt/PyPhysics/issues/92).
It exposes precision-physics experimental data via namespaced tools:

  nist.*    -> NIST ASD energy levels / transitions + CODATA constants
              (migrated from the standalone nist_mcp under issue #99)
  pdg.*     -> PDG particle masses, lifetimes, branching ratios, g-2 (issue #96)
  qed.*     -> Bound-state QED: Lamb shift, hyperfine, g-factor (issue #97)
  nuclear.* -> Charge radii (incl. proton-radius puzzle), magnetic moments (issue #98)
  flag.*    -> FLAG lattice-QCD averages (issue #98)
  astro.*   -> Astronomical / gravitational observables (deferred per umbrella #92)
  cosmo.*   -> Cosmological parameters (deferred per umbrella #92)

Every tool returns a dict carrying ``{value, uncertainty, unit, source, retrieved_at,
cache_key}`` where ``source`` is a ``cite_key`` resolving against
``Roadmapping/History/Bibliography/``. See ``REFRESH_POLICY.md`` for the shared
caching + refresh schema.
"""

from mcp.server.fastmcp import FastMCP

from precision_data_mcp.nist.tools import asd, codata, targets

mcp = FastMCP("precision_data")


# ---------------------------------------------------------------------------
# nist.* namespace (migrated from nist_mcp under issue #99; same contract)
# ---------------------------------------------------------------------------


@mcp.tool(name="nist.get_constant")
def nist_get_constant(name: str) -> dict:
    """CODATA fundamental physical constant by name (offline, via scipy.constants).

    Returns a status-discriminated dict: {"status": "ok", "match": record} for a
    unique hit, {"status": "ambiguous", "candidates": [...]} for several substring
    hits, or {"status": "not_found", "candidates": []}. The record's ``source``
    names the CODATA release used (scipy bundles a specific release).

    Args:
        name: constant name, e.g. "electron g factor", "Rydberg constant".
    """
    return codata.get_constant(name)


@mcp.tool(name="nist.get_levels")
def nist_get_levels(species: str, refresh: bool = False) -> dict:
    """NIST ASD energy levels for a species. Native unit: cm^-1 (no conversion).

    Source: NIST Atomic Spectra Database (https://www.nist.gov/pml/atomic-spectra-database).
    Note: for light elements (e.g. hydrogen) ASD levels are largely theoretical/Ritz
    values, not independent measurements.

    Args:
        species: spectrum label, e.g. "H I", "He I", "Li III".
        refresh: bypass the on-disk cache and re-fetch.
    """
    return asd.get_levels(species, refresh=refresh)


@mcp.tool(name="nist.get_transitions")
def nist_get_transitions(species: str, refresh: bool = False) -> dict:
    """NIST ASD transitions for a species. Native unit: nm, with an air/vacuum flag.

    Source: NIST Atomic Spectra Database (https://www.nist.gov/pml/atomic-spectra-database).

    Args:
        species: spectrum label, e.g. "H I", "He I", "Li III".
        refresh: bypass the on-disk cache and re-fetch.
    """
    return asd.get_transitions(species, refresh=refresh)


@mcp.tool(name="nist.list_dirac_targets")
def nist_list_dirac_targets() -> list[dict]:
    """Curated, individually-cited precision splittings for the Dirac-equation exploration.

    These MHz/sub-MHz values (hydrogen fine structure, Lamb shift, 1S hyperfine,
    helium 3P fine structure) are NOT sourced from ASD — ASD does not carry them
    at this precision. Each record cites its provenance. The electron g-factor is
    delegated to nist.get_constant (CODATA).

    Per the umbrella issue #92 §Resolved-decisions #5 (disagreement representation):
    when issue #97 lands the qed.* namespace these targets will be re-routed to the
    qed.* tools, since they are bound-state QED observables proper.
    """
    return targets.list_dirac_targets()


# ---------------------------------------------------------------------------
# Future namespaces — to be registered when their sub-issues land:
#   pdg.*     (issue #96)
#   qed.*     (issue #97)
#   nuclear.* + flag.*  (issue #98)
# ---------------------------------------------------------------------------


def main() -> None:
    """stdio MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
