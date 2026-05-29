# Task 90: nist mcp

## Objective

Build a new MCP server, `nist_mcp`, that exposes NIST experimental reference data — atomic energy levels and transitions from the **NIST Atomic Spectra Database (ASD)** and **CODATA fundamental physical constants** — as Claude-callable tools, with every returned quantity carrying **value, uncertainty (precision), unit, and source/citation**. The server lives alongside the existing server at [mcp_servers/pyphysics_mcp/](../../mcp_servers/pyphysics_mcp/) and is registered as a second entry in [.mcp.json](../../.mcp.json).

**Source-of-record (per issue comment, MRK):** the NIST Atomic Spectra Database landing page <https://www.nist.gov/pml/atomic-spectra-database> is the canonical entry point for all ASD data to pull; the actual machine-readable queries go through the Lines and Levels query forms that page links to (the `lines1.pl` / `energy1.pl` CGI endpoints). Document this URL as the source-of-record in `asd.py` and the README.

**Primary consumer — the Dirac-equation exploration (per issue comment):** the first values this server must be able to retrieve are the ones that support the repo's Dirac-equation work, i.e. the Bethe–Salpeter precision-predictions campaign under [Roadmapping/Quantum_Mechanics/Bethe_Salpeter/](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/) and the dual-Dirac g-factor experimental design. Those documents currently hand-type experimental comparison values from CODATA-2018 / PDG-2024; this server makes them programmatically retrievable with their uncertainties.

## Background

The repo already ships one FastMCP server, `pyphysics_mcp`, which establishes the pattern this task mirrors:

- **Server wiring** — [mcp_servers/pyphysics_mcp/src/pyphysics_mcp/server.py:30](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/server.py#L30) constructs `mcp = FastMCP("pyphysics")` and registers each tool with an `@mcp.tool()`-decorated wrapper (lines 35-128) that delegates to a function in a `tools/` module. `main()` calls `mcp.run()` ([server.py:133](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/server.py#L133)).
- **Package entry points** — [src/pyphysics_mcp/__init__.py](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/__init__.py) re-exports `main`; [src/pyphysics_mcp/__main__.py](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/__main__.py) enables `python -m`. The console script is declared in [pyproject.toml:12-13](../../mcp_servers/pyphysics_mcp/pyproject.toml#L12-L13) as `pyphysics-mcp = "pyphysics_mcp:main"`, build backend `hatchling`, packages `src/pyphysics_mcp` ([pyproject.toml:19-20](../../mcp_servers/pyphysics_mcp/pyproject.toml#L19-L20)).
- **Path/config resolution** — [src/pyphysics_mcp/repo.py:23-42](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/repo.py#L23-L42) shows the env-var-then-walk-up pattern (`PYPHYSICS_REPO_PATH`) with an `@lru_cache` and a clear `RepoNotFound` error. `nist_mcp` reuses this shape for its cache directory.
- **Tool module style** — [src/pyphysics_mcp/tools/validation.py](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/tools/validation.py) shows tools returning a structured dict (`{status, returncode, stdout, stderr, ...}`); the other tool modules (`bibliography.py`, `chapters.py`, `animations.py`) return typed dicts/lists parsed in-process.
- **Registration** — [.mcp.json:2-6](../../.mcp.json#L2-L6) registers the existing server purely by console-script command (`"command": "pyphysics-mcp"`). A `nist` key is added the same way.
- **README pattern** — [mcp_servers/pyphysics_mcp/README.md](../../mcp_servers/pyphysics_mcp/README.md) documents install (`uv tool install ./mcp_servers/pyphysics_mcp`), Claude Code registration, the env-var resolution order, the architecture tree, and per-tool verification steps. `nist_mcp` gets an analogous README.

Why this is wanted (from the issue): derivations and verification cards in this repo (e.g. the g-factor / DRQM findings recorded in [CLAUDE.md](../../CLAUDE.md), "Three open findings flagged for author review") compare theory against experimental values that are currently looked up by hand. A NIST MCP lets those comparisons fetch the value **and** its stated uncertainty in one typed call, so error bars are correct and lookups are reproducible.

**Concrete Dirac-exploration target values.** The Bethe–Salpeter campaign README ([Roadmapping/Quantum_Mechanics/Bethe_Salpeter/README.md:3](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/README.md#L3)) states results are compared "against modern measured values (CODATA-2018 / PDG-2024)," and the cross-comparison table ([Bethe_Salpeter/10_CrossComparison.md:61-68](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md#L61-L68)) lists the exact hand-typed measurements the dual-Dirac predictions are checked against — these are the values the ASD tools should be able to source or corroborate:

- **Hydrogen `2P₃/₂–2P₁/₂` fine structure** ≈ 10,969.13(10) MHz — the §14 fine-structure pivot (PR C).
- **Hydrogen `2S₁/₂–2P₁/₂` Lamb shift** = 1,057.845(9) MHz — the §§19–21 headline (PR E).
- **Hydrogen 1S hyperfine (21-cm)** = 1,420.405,751,768(2) MHz — the §22 headline (PR F).
- **Helium `³P₀–³P₁` fine structure** = 29,616.952 MHz (PR I).
- **Electron g-factor** `g_s = -2.00231930…` (CODATA) — the dual-Dirac g-factor formula in [Roadmapping/Experimental_Design_Plan_Dual_Theory.md:167](../../Roadmapping/Experimental_Design_Plan_Dual_Theory.md#L167) and the DRQM `r_e` finding.

The fine-structure / Lamb-shift / hyperfine section documents that consume these values are [Bethe_Salpeter/03_FineStructure.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md), [05_LambShift.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/05_LambShift.md), and [06_Hyperfine.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/06_Hyperfine.md). The electron g-factor comes from the CODATA path (`get_constant`); the hydrogen and helium level energies / fine-structure splittings come from the ASD path (`get_levels` / `get_transitions`) — so the server's two data sources jointly cover the campaign's comparison column.

Resolving the issue's "API vs scraping" open question — concrete, grounded approaches:

- **CODATA constants** are available offline via `scipy.constants.physical_constants`, a dict mapping each constant name → `(value, unit, uncertainty)`. This is the canonical, no-network source for `get_constant` and already carries the uncertainty the issue asks for.
- **NIST ASD** is reached from the landing page <https://www.nist.gov/pml/atomic-spectra-database> (the documented source-of-record), which links to the server-side CGI query forms that return machine-readable output: the lines endpoint (`https://physics.nist.gov/cgi-bin/ASD/lines1.pl`) and the levels endpoint (`https://physics.nist.gov/cgi-bin/ASD/energy1.pl`). Both accept query parameters selecting a delimited (CSV / tab-delimited) output format and include observed/uncertainty columns. The implementer confirms exact parameter names/values against a live query during implementation and records both the landing-page URL and the chosen CGI parameters in the module docstring.

## Implementation Plan

1. **Scaffold the sub-project** `mcp_servers/nist_mcp/` mirroring `pyphysics_mcp`'s layout:
   ```
   mcp_servers/nist_mcp/
   ├── pyproject.toml
   ├── README.md
   └── src/nist_mcp/
       ├── __init__.py        # re-exports main
       ├── __main__.py        # python -m nist_mcp
       ├── server.py          # FastMCP("nist") + @mcp.tool() wrappers
       ├── cache.py           # cache-dir resolution + on-disk JSON cache
       └── tools/
           ├── __init__.py
           ├── codata.py      # get_constant
           └── asd.py         # get_levels, get_transitions
   ```

2. **`pyproject.toml`** — copy the `pyphysics_mcp` structure ([pyproject.toml](../../mcp_servers/pyphysics_mcp/pyproject.toml)): name `nist-mcp`, `requires-python = ">=3.10"`, build backend `hatchling`, `packages = ["src/nist_mcp"]`, console script `nist-mcp = "nist_mcp:main"`. Dependencies: `mcp>=1.0`, `scipy>=1.11` (CODATA), `requests>=2.31` (ASD HTTP). `__init__.py` / `__main__.py` re-export/call `nist_mcp.server.main`, matching the existing two-line files.

3. **`cache.py`** — resolve a cache directory in the order: `NIST_MCP_CACHE_DIR` env var → platform default (`~/.cache/nist_mcp/`), creating it if absent (`@lru_cache`, clear error on failure — same shape as `repo.py:repo_root`). Provide `cache_get(key) -> dict | None` and `cache_put(key, payload)` writing one JSON file per key (key = stable hash of the normalized query). ASD lookups consult the cache before any network call; CODATA needs no cache.

4. **`tools/codata.py`** — `get_constant(name: str) -> dict`:
   - Look `name` up in `scipy.constants.physical_constants`; if no exact key, do case-insensitive substring matching and, when ambiguous, return the candidate key list rather than guessing.
   - Return the shared record schema (see step 6): `{quantity, value, uncertainty, unit, source: "CODATA (scipy.constants)", reference}`. `uncertainty` is the scipy-reported standard uncertainty (may be `0.0` for exactly-defined SI constants — surface that honestly).

5. **`tools/asd.py`** — two functions reaching the ASD via the source-of-record landing page's Lines/Levels CGI forms, using `requests`, parsing the delimited response, and mapping value/uncertainty columns into the shared schema. Record the landing-page URL (<https://www.nist.gov/pml/atomic-spectra-database>) and the resolved CGI parameters in the module docstring.
   - `get_levels(species: str, ...) -> list[dict]` → energy levels for an ion/atom (e.g. `"H I"`, `"He I"`, `"Li III"`), each with energy value + uncertainty + unit + reference.
   - `get_transitions(species: str, ...) -> list[dict]` → transition wavelengths/frequencies with uncertainties.
   - Build the request URL from documented ASD parameters selecting delimited output; parse defensively (skip header/comment lines, tolerate empty results), and return `{"results": [...], "source": "NIST ASD", "query": {...}}`. On HTTP error or empty result, return `{"results": [], "error": ...}` rather than raising. Route every call through `cache.py`.

5a. **Dirac-exploration target set.** In `asd.py`, define a small `DIRAC_TARGETS` registry (a list of `{species, n, term/level labels, observable}` records) covering the Bethe–Salpeter campaign's hydrogen and helium comparison points — hydrogen `2P₃/₂`/`2P₁/₂`/`2S₁/₂` levels (fine structure + Lamb shift), hydrogen 1S (hyperfine), and helium `³P_J` levels — so the campaign documents can request exactly these without re-deriving query strings. Surface them either as documented example arguments to `get_levels`/`get_transitions` in the README or as an optional `list_dirac_targets()` helper. This registry seeds the live ASD test (step 10) and is the concrete bridge to [Bethe_Salpeter/10_CrossComparison.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md). The electron g-factor target is served by `get_constant` (CODATA), not ASD.

6. **Shared output schema** — every quantity-bearing record returned by any tool is a dict with keys: `quantity` (str), `value` (float), `uncertainty` (float | None), `unit` (str), `source` (str), `reference` (str | None). Document this once in the README and keep it consistent across `codata.py` and `asd.py`.

7. **`server.py`** — `mcp = FastMCP("nist")`; register `get_levels`, `get_transitions`, `get_constant` as thin `@mcp.tool()` wrappers delegating to the `tools/` functions (mirroring [pyphysics_mcp/server.py:35-128](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/server.py#L35-L128)); `def main(): mcp.run()`.

8. **Register the server** — add a `nist` entry to [.mcp.json](../../.mcp.json) alongside `pyphysics`:
   ```json
   {
     "mcpServers": {
       "pyphysics": { "command": "pyphysics-mcp" },
       "nist": { "command": "nist-mcp" }
     }
   }
   ```

9. **`README.md`** — document install (`uv tool install ./mcp_servers/nist_mcp`), Claude Code registration, the three tools + shared schema, the cache-dir resolution order, and per-tool verification examples (e.g. `get_constant("electron g factor")`, `get_levels("H I")`), mirroring [pyphysics_mcp/README.md](../../mcp_servers/pyphysics_mcp/README.md).

10. **Tests** — add `mcp_servers/nist_mcp/tests/`: an offline, deterministic test for `get_constant` (assert a known CODATA constant returns the expected value/uncertainty/unit and that the schema keys are present), and a network-gated test for the ASD tools (skipped unless an opt-in env var like `NIST_MCP_LIVE=1` is set, so CI stays offline) plus an offline parser test that feeds a captured ASD response fixture through the parsing function.

## Files to Modify

| File | Change |
|---|---|
| `mcp_servers/nist_mcp/pyproject.toml` | New. Package metadata, `nist-mcp` console script, deps `mcp`, `scipy`, `requests`; hatchling build. |
| `mcp_servers/nist_mcp/README.md` | New. Install / registration / tool docs / shared schema / verification examples; documents the ASD landing-page source-of-record and the Dirac-exploration target examples. |
| `mcp_servers/nist_mcp/src/nist_mcp/__init__.py` | New. Re-export `main` from `nist_mcp.server`. |
| `mcp_servers/nist_mcp/src/nist_mcp/__main__.py` | New. `python -m nist_mcp` entry calling `main`. |
| `mcp_servers/nist_mcp/src/nist_mcp/server.py` | New. `FastMCP("nist")` + `@mcp.tool()` wrappers for `get_levels`, `get_transitions`, `get_constant`; `main()`. |
| `mcp_servers/nist_mcp/src/nist_mcp/cache.py` | New. Cache-dir resolution (`NIST_MCP_CACHE_DIR` → default) + on-disk JSON `cache_get`/`cache_put`. |
| `mcp_servers/nist_mcp/src/nist_mcp/tools/__init__.py` | New. Package marker. |
| `mcp_servers/nist_mcp/src/nist_mcp/tools/codata.py` | New. `get_constant` over `scipy.constants.physical_constants`. |
| `mcp_servers/nist_mcp/src/nist_mcp/tools/asd.py` | New. `get_levels`, `get_transitions` via the ASD landing-page Lines/Levels CGI forms, parsing value+uncertainty; `DIRAC_TARGETS` registry (+ optional `list_dirac_targets`) for the Bethe–Salpeter comparison points. |
| `mcp_servers/nist_mcp/tests/test_codata.py` | New. Offline test for `get_constant` value/uncertainty/unit + schema. |
| `mcp_servers/nist_mcp/tests/test_asd.py` | New. Offline parser test (fixture) + network-gated live test (`NIST_MCP_LIVE=1`). |
| `.mcp.json` | Add a `nist` server entry (`"command": "nist-mcp"`) next to `pyphysics`. |

## Dependencies

- **New Python deps (declared in `mcp_servers/nist_mcp/pyproject.toml` only, not the root project):** `mcp>=1.0` (FastMCP, already used by `pyphysics_mcp`), `scipy>=1.11` (CODATA constants via `scipy.constants.physical_constants`), `requests>=2.31` (ASD HTTP queries). `pytest` as a dev/test dependency for the test files.
- **External service:** NIST Atomic Spectra Database CGI endpoints (`physics.nist.gov/cgi-bin/ASD/lines1.pl`, `.../energy1.pl`) — network access required for live ASD calls; CODATA path is fully offline.
- **Env vars:** `NIST_MCP_CACHE_DIR` (optional, overrides default cache location); `NIST_MCP_LIVE` (optional, opt-in for the live network test).

## Acceptance Criteria

- [ ] `uv tool install ./mcp_servers/nist_mcp` succeeds and installs a `nist-mcp` console script.
- [ ] `.mcp.json` contains both a `pyphysics` and a `nist` server entry, and remains valid JSON.
- [ ] `get_constant(name)` returns a record with `value`, `uncertainty`, `unit`, and `source` for a known CODATA constant (e.g. the electron g-factor), sourced from `scipy.constants.physical_constants`.
- [ ] `get_levels(species)` and `get_transitions(species)` return ASD records where each quantity carries `value`, `uncertainty`, `unit`, and `source` ("NIST ASD"), and return an empty-results record (not an exception) on no-match or network failure.
- [ ] Every tool's returned quantity conforms to the shared schema (`quantity, value, uncertainty, unit, source, reference`).
- [ ] ASD lookups are cached on disk and a repeated identical query is served from cache without a second network call.
- [ ] `asd.py` documents <https://www.nist.gov/pml/atomic-spectra-database> as the source-of-record and exposes a Dirac-exploration target set covering the Bethe–Salpeter hydrogen/helium comparison points.
- [ ] At least one Dirac-exploration value is retrievable end-to-end with its uncertainty: hydrogen `2P` fine-structure levels (or the `2P₃/₂–2P₁/₂` transition) via the ASD path, and the electron g-factor via `get_constant` — matching the comparison column in [Bethe_Salpeter/10_CrossComparison.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md).
- [ ] `mcp_servers/nist_mcp/README.md` documents install, Claude Code registration, the three tools, the shared schema, and cache-dir resolution.
- [ ] Offline tests pass without network; the live ASD test is skipped unless `NIST_MCP_LIVE=1`.

## Testing

Reviewer commands (run from the repo root):

```bash
# Install the new server.
uv tool install ./mcp_servers/nist_mcp
nist-mcp --help 2>&1 | head -3            # should not error (server expects stdio from a client)

# Offline test suite (no network).
uv run --project mcp_servers/nist_mcp pytest mcp_servers/nist_mcp/tests -q

# Opt-in live ASD check (hits physics.nist.gov).
NIST_MCP_LIVE=1 uv run --project mcp_servers/nist_mcp pytest mcp_servers/nist_mcp/tests/test_asd.py -q

# Validate .mcp.json is still well-formed and has both servers.
uv run python -c "import json; d=json.load(open('.mcp.json')); print(sorted(d['mcpServers']))"
```

Tests added:

- `mcp_servers/nist_mcp/tests/test_codata.py` — asserts `get_constant("electron g factor")` (or another stable CODATA key) returns the expected value/uncertainty/unit and that all shared-schema keys are present; asserts ambiguous/unknown names return candidate lists rather than raising.
- `mcp_servers/nist_mcp/tests/test_asd.py` — offline test feeding a captured ASD delimited-response fixture through the `asd.py` parser and asserting value+uncertainty extraction and schema conformance; plus a `NIST_MCP_LIVE`-gated test that performs a real `get_levels`/`get_transitions` call against the Dirac-exploration target set (hydrogen `2P` fine-structure levels and helium `³P` levels) and asserts non-empty, schema-conforming results with non-null uncertainties.

After registering with Claude Code (restart required), confirm the deferred tools list includes `mcp__nist__get_constant`, `mcp__nist__get_levels`, and `mcp__nist__get_transitions`.
