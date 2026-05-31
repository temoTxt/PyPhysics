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
from precision_data_mcp.pdg.tools import lookup as pdg_lookup
from precision_data_mcp.qed.tools import lookup as qed_lookup
from precision_data_mcp.nuclear.tools import lookup as nuclear_lookup
from precision_data_mcp.flag.tools import lookup as flag_lookup

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
# pdg.* namespace (issue #96 — Particle Data Group reference values)
# ---------------------------------------------------------------------------


@mcp.tool(name="pdg.get_particle")
def pdg_get_particle(particle_id: str) -> dict:
    """PDG mass, lifetime, charge, and spin for a particle.

    Multi-value fields (e.g. the neutron lifetime under the bottle-vs-beam
    puzzle) return a ``{"values": [...]}`` list per the umbrella's
    disagreement-representation rule rather than an averaged single value.

    Args:
        particle_id: PDG MCID (e.g. "11" for electron, "13" for muon) or
            a particle name / symbol ("electron", "mu-", "proton", "p+").
    """
    return pdg_lookup.get_particle(particle_id)


@mcp.tool(name="pdg.get_anomaly")
def pdg_get_anomaly(particle_id: str) -> dict:
    """Anomalous magnetic moment for a particle.

    The muon a_mu entry surfaces the canonical Fermilab-E989 vs Theory-Initiative
    vs BMW-lattice-QCD disagreement explicitly as a ``values`` list.

    Args:
        particle_id: PDG MCID or particle name / symbol; currently only
            leptons (electron MCID 11, muon MCID 13) have curated anomaly data.
    """
    return pdg_lookup.get_anomaly(particle_id)


@mcp.tool(name="pdg.list_particles")
def pdg_list_particles() -> list[dict]:
    """Enumerate seeded particles with their MCIDs / names / symbols."""
    return pdg_lookup.list_particles()


# ---------------------------------------------------------------------------
# qed.* namespace (issue #97 — bound-state QED precision observables)
# ---------------------------------------------------------------------------


@mcp.tool(name="qed.get_lamb_shift")
def qed_get_lamb_shift(species: str, transition: str) -> dict:
    """Lamb-shift value for a given species + transition.

    Examples:
        qed.get_lamb_shift("H", "2S1/2-2P1/2")  # hydrogen 1057.845(9) MHz
        qed.get_lamb_shift("muonic_H", "2S1/2-2P3/2")  # the proton-radius-puzzle trigger
    """
    return qed_lookup.get_lamb_shift(species, transition)


@mcp.tool(name="qed.get_hyperfine")
def qed_get_hyperfine(species: str, level: str) -> dict:
    """Hyperfine-splitting value for a given species + level.

    Examples:
        qed.get_hyperfine("H", "1s2S1/2")              # 1420.405... MHz (21-cm)
        qed.get_hyperfine("muonium", "1s2S1/2")        # 4463.302... MHz
        qed.get_hyperfine("positronium", "1s_ortho_para")  # ~203389 MHz
    """
    return qed_lookup.get_hyperfine(species, level)


@mcp.tool(name="qed.get_g_factor")
def qed_get_g_factor(species: str, electron_state: str) -> dict:
    """Bound-electron g-factor for a hydrogenic-ion electron state.

    Examples:
        qed.get_g_factor("Si XIV", "1s2S1/2")  # Sturm et al. 2013 — issue #82 target
        qed.get_g_factor("He II", "1s2S1/2")   # Kohler et al. 2016
    """
    return qed_lookup.get_g_factor(species, electron_state)


@mcp.tool(name="qed.get_transition_precision")
def qed_get_transition_precision(species: str, transition: str) -> dict:
    """High-precision transition frequency.

    Examples:
        qed.get_transition_precision("H", "1S-2S")  # Hänsch et al. 2011 absolute frequency
    """
    return qed_lookup.get_transition_precision(species, transition)


@mcp.tool(name="qed.list_species")
def qed_list_species() -> list[dict]:
    """Enumerate seeded species (H, He II, Li III, Si XIV, muonic_H, muonium, positronium)."""
    return qed_lookup.list_species()


@mcp.tool(name="qed.list_observables")
def qed_list_observables(species: str) -> dict:
    """For a given species, list which observables are curated (Lamb shifts / hyperfine / g-factor / transitions)."""
    return qed_lookup.list_observables(species)


# ---------------------------------------------------------------------------
# nuclear.* namespace (issue #98 — nuclear / particle structure observables)
# ---------------------------------------------------------------------------


@mcp.tool(name="nuclear.get_charge_radius")
def nuclear_get_charge_radius(nucleus: str, method: str | None = None) -> dict:
    """Charge-radius value(s) for a nucleus.

    The proton-radius puzzle is the canonical demonstrator of the umbrella's
    disagreement-representation rule: ``nuclear.get_charge_radius("proton")``
    returns ALL values (muonic-H, electron-scattering, 1S-3S, PRad, CODATA-pre-puzzle)
    with method labels, never an averaged "best" value.

    Args:
        nucleus: ``p``/``proton``, ``n``/``neutron``, ``d``/``deuteron``, ...
        method: optional method label to select a single value when disambiguation is desired.
    """
    return nuclear_lookup.get_charge_radius(nucleus, method=method)


@mcp.tool(name="nuclear.get_magnetic_moment")
def nuclear_get_magnetic_moment(nucleus: str) -> dict:
    """Nuclear magnetic moment in units of nuclear magnetons."""
    return nuclear_lookup.get_magnetic_moment(nucleus)


@mcp.tool(name="nuclear.list_nuclei")
def nuclear_list_nuclei() -> list[dict]:
    """Enumerate seeded nuclei."""
    return nuclear_lookup.list_nuclei()


# ---------------------------------------------------------------------------
# flag.* namespace (issue #98 — FLAG lattice-QCD averages)
# ---------------------------------------------------------------------------


@mcp.tool(name="flag.get_quantity")
def flag_get_quantity(quantity_id: str) -> dict:
    """FLAG-averaged lattice-QCD quantity (e.g. ``f_pi``, ``g_A``, ``sigma_piN``)."""
    return flag_lookup.get_quantity(quantity_id)


@mcp.tool(name="flag.list_quantities")
def flag_list_quantities() -> list[dict]:
    """Enumerate seeded FLAG quantities."""
    return flag_lookup.list_quantities()


# ---------------------------------------------------------------------------
# Future namespaces — deferred per umbrella #92 §Deferred:
#   astro.*   — GPS clocks, Mercury perihelion, GP-B, Hulse-Taylor, LIGO
#   cosmo.*   — H_0, CMB, Planck params, BBN
# ---------------------------------------------------------------------------


def main() -> None:
    """stdio MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
