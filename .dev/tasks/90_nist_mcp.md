# Task 90: nist mcp

## Objective

Build a new MCP server, `nist_mcp`, that exposes NIST experimental reference data — atomic energy levels and transitions from the **NIST Atomic Spectra Database (ASD)** and **CODATA fundamental physical constants** — as Claude-callable tools, with every returned quantity carrying **value, uncertainty (precision), unit, and source/citation**. The server lives alongside the existing server at [mcp_servers/pyphysics_mcp/](../../mcp_servers/pyphysics_mcp/) and is registered as a second entry in [.mcp.json](../../.mcp.json).

**Source-of-record (per issue comment, MRK):** the NIST Atomic Spectra Database landing page <https://www.nist.gov/pml/atomic-spectra-database> is the canonical entry point for all ASD data to pull; the actual machine-readable queries go through the Lines and Levels query forms that page links to (the `lines1.pl` / `energy1.pl` CGI endpoints). Document this URL as the source-of-record in `asd.py` and the README.

**Primary consumer — the Dirac-equation exploration (per issue comment):** the first values this server must be able to retrieve are the ones that support the repo's Dirac-equation work, i.e. the Bethe–Salpeter precision-predictions campaign under [Roadmapping/Quantum_Mechanics/Bethe_Salpeter/](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/) and the dual-Dirac g-factor experimental design. Those documents currently hand-type experimental comparison values from CODATA-2018 / PDG-2024. This server makes those values programmatically retrievable with their uncertainties across **three honest retrieval paths**, because they do not all come from the same place: CODATA constants via `get_constant`; a curated, individually-cited registry (`list_dirac_targets`) for the MHz/sub-MHz precision splittings that the ASD general tables do **not** carry; and general ASD level/line lookups (`get_levels` / `get_transitions`) for the quantities ASD genuinely covers. The "Critical scope correction" in Background is load-bearing — the ASD path is **not** the source for the cross-comparison table's precision splittings.

## Background

The repo already ships one FastMCP server, `pyphysics_mcp`, which establishes the pattern this task mirrors:

- **Server wiring** — [mcp_servers/pyphysics_mcp/src/pyphysics_mcp/server.py:30](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/server.py#L30) constructs `mcp = FastMCP("pyphysics")` and registers each tool with an `@mcp.tool()`-decorated wrapper (lines 35-128) that delegates to a function in a `tools/` module. `main()` calls `mcp.run()` ([server.py:133](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/server.py#L133)).
- **Package entry points** — [src/pyphysics_mcp/__init__.py](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/__init__.py) re-exports `main`; [src/pyphysics_mcp/__main__.py](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/__main__.py) enables `python -m`. The console script is declared in [pyproject.toml:12-13](../../mcp_servers/pyphysics_mcp/pyproject.toml#L12-L13) as `pyphysics-mcp = "pyphysics_mcp:main"`, build backend `hatchling`, packages `src/pyphysics_mcp` ([pyproject.toml:19-20](../../mcp_servers/pyphysics_mcp/pyproject.toml#L19-L20)).
- **Path/config resolution** — [src/pyphysics_mcp/repo.py:23-42](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/repo.py#L23-L42) shows the env-var-then-walk-up pattern (`PYPHYSICS_REPO_PATH`) with an `@lru_cache` and a clear `RepoNotFound` error. `nist_mcp` reuses this shape for its cache directory.
- **Tool module style** — [src/pyphysics_mcp/tools/validation.py](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/tools/validation.py) shows tools returning a structured dict (`{status, returncode, stdout, stderr, ...}`); the other tool modules (`bibliography.py`, `chapters.py`, `animations.py`) return typed dicts/lists parsed in-process.
- **Registration** — [.mcp.json:2-6](../../.mcp.json#L2-L6) registers the existing server purely by console-script command (`"command": "pyphysics-mcp"`). A `nist` key is added the same way.
- **README pattern** — [mcp_servers/pyphysics_mcp/README.md](../../mcp_servers/pyphysics_mcp/README.md) documents install (`uv tool install ./mcp_servers/pyphysics_mcp`), Claude Code registration, the env-var resolution order, the architecture tree, and per-tool verification steps. `nist_mcp` gets an analogous README.

Why this is wanted (from the issue): derivations and verification cards in this repo (e.g. the g-factor / DRQM findings recorded in [CLAUDE.md](../../CLAUDE.md), "Three open findings flagged for author review") compare theory against experimental values that are currently looked up by hand. A NIST MCP lets those comparisons fetch the value **and** its stated uncertainty in one typed call, so error bars are correct and lookups are reproducible.

**Concrete Dirac-exploration target values.** The Bethe–Salpeter campaign README ([Roadmapping/Quantum_Mechanics/Bethe_Salpeter/README.md:3](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/README.md#L3)) states results are compared "against modern measured values (CODATA-2018 / PDG-2024)," and the cross-comparison table ([Bethe_Salpeter/10_CrossComparison.md:61-68](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md#L61-L68)) lists the exact hand-typed measurements the dual-Dirac predictions are checked against:

- **Hydrogen `2P₃/₂–2P₁/₂` fine structure** ≈ 10,969.13(10) MHz — the §14 fine-structure pivot (PR C).
- **Hydrogen `2S₁/₂–2P₁/₂` Lamb shift** = 1,057.845(9) MHz — the §§19–21 headline (PR E).
- **Hydrogen 1S hyperfine (21-cm)** = 1,420.405,751,768(2) MHz — the §22 headline (PR F).
- **Helium `³P₀–³P₁` fine structure** = 29,616.952 MHz (PR I).
- **Electron g-factor** `g_s = -2.00231930…` (CODATA) — the dual-Dirac g-factor formula in [Roadmapping/Experimental_Design_Plan_Dual_Theory.md:167](../../Roadmapping/Experimental_Design_Plan_Dual_Theory.md#L167) and the DRQM `r_e` finding.

**Critical scope correction — do not over-promise the ASD path.** The four splittings above are MHz / sub-MHz quantities from dedicated precision-spectroscopy experiments. The NIST ASD general level/line tables do **not** carry them at that precision, and for light elements (notably hydrogen) ASD's tabulated level energies are largely *theoretical/Ritz* values rather than independent measurements — so pulling them to "corroborate" a theory prediction would be circular. The retrieval responsibilities therefore split three ways, and the plan must not blur them:

- **Electron g-factor** → CODATA (`get_constant`).
- **The MHz/sub-MHz splittings** (Lamb shift, 1S hyperfine, H 2P fine structure, He ³P) → a **curated, individually-cited registry** (`list_dirac_targets`, see step 5a), *not* scraped from ASD. Each entry cites its own source (the value the campaign already hand-typed, now made retrievable with provenance).
- **ASD path** (`get_levels` / `get_transitions`) → the general-purpose retrieval layer for level energies (cm⁻¹) and transition wavelengths (nm) across species. It is genuinely useful for the broader hydrogenic / high-Z work (e.g. `Li III`) and as a structural cross-check, but it is **not** the source for the cross-comparison table's precision splittings.

The fine-structure / Lamb-shift / hyperfine section documents that consume these values are [Bethe_Salpeter/03_FineStructure.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md), [05_LambShift.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/05_LambShift.md), and [06_Hyperfine.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/06_Hyperfine.md).

Resolving the issue's "API vs scraping" open question — concrete, grounded approaches:

- **CODATA constants** are available offline via `scipy.constants.physical_constants`, a dict mapping each constant name → `(value, unit, uncertainty)`. This is the canonical, no-network source for `get_constant` and already carries the uncertainty the issue asks for. **Version caveat:** scipy ships whichever CODATA release its version bundles (scipy 1.11 → CODATA-2018; later scipy → CODATA-2022), and the campaign cites **CODATA-2018**. `get_constant` must therefore report the release actually used (`scipy.constants.value`/`scipy.constants.codata` exposes it) in the `source` field, and the dependency pin should target a scipy whose CODATA release matches the docs — a mismatch must be surfaced, never silently served.
- **NIST ASD** is reached from the landing page <https://www.nist.gov/pml/atomic-spectra-database> (the documented source-of-record), which links to the server-side CGI query forms that return machine-readable output: the lines endpoint (`https://physics.nist.gov/cgi-bin/ASD/lines1.pl`) and the levels endpoint (`https://physics.nist.gov/cgi-bin/ASD/energy1.pl`). Both accept query parameters selecting a delimited (CSV / tab-delimited) output format and include observed/uncertainty columns. Constraints the implementation must respect:
  - **Native units, no silent conversion.** Return ASD values verbatim — level energies in **cm⁻¹**, line wavelengths in **nm** with an explicit **air-vs-vacuum** flag (a classic ASD footgun) — and put the native unit in the `unit` field. Conversion to MHz/eV is the caller's responsibility, documented but not applied inside the tool.
  - **Capture provenance.** Map ASD's per-record bibliographic/reference codes into the `reference` field where present.
  - **Request hygiene.** Build query parameters via `requests` `params=` (URL-encoded, never string-interpolated into the URL), send a descriptive `User-Agent` identifying this repo, set a timeout, and bound retries — physics.nist.gov is a shared public service.
  - **Pin the brittle bits.** The exact CGI parameter set is scraped, not a stable API; record the landing-page URL and the resolved parameters in the module docstring and back the parser with a committed response fixture (step 10) so a NIST form change surfaces as a test failure, not silent garbage. (`astroquery.nist` was considered as an alternative but only wraps *lines*, pulls in the heavy `astropy` dependency, and does not expose levels — direct HTTP covers both uniformly.)

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
           ├── asd.py         # get_levels, get_transitions
           └── targets.py     # DIRAC_TARGETS curated+cited registry, list_dirac_targets
   ```

2. **`pyproject.toml`** — copy the `pyphysics_mcp` structure ([pyproject.toml](../../mcp_servers/pyphysics_mcp/pyproject.toml)): name `nist-mcp`, `requires-python = ">=3.10"`, build backend `hatchling`, `packages = ["src/nist_mcp"]`, console script `nist-mcp = "nist_mcp:main"`. Dependencies: `mcp>=1.0`, `scipy>=1.11` (CODATA), `requests>=2.31` (ASD HTTP). `__init__.py` / `__main__.py` re-export/call `nist_mcp.server.main`, matching the existing two-line files.

3. **`cache.py`** — resolve a cache directory in the order: `NIST_MCP_CACHE_DIR` env var → platform default (`~/.cache/nist_mcp/`), creating it if absent (`@lru_cache`, clear error on failure — same shape as `repo.py:repo_root`). Provide `cache_get(key) -> dict | None` and `cache_put(key, payload)` writing one JSON file per key (key = stable hash of the normalized query **plus a `CACHE_SCHEMA_VERSION` constant**, so a parser/schema change invalidates stale entries instead of returning old-format payloads). Store a fetch timestamp in each entry and honor an optional TTL / `refresh=True` to bypass a cached entry and re-fetch. Document that the cache is cleared by deleting the directory. ASD lookups consult the cache before any network call; CODATA needs no cache.

4. **`tools/codata.py`** — `get_constant(name: str) -> dict`:
   - Look `name` up in `scipy.constants.physical_constants`. Always return **one stable shape** with a `status` discriminator, so an MCP client never has to branch on record-vs-list: `{"status": "ok", "match": <record>}` for a unique hit; `{"status": "ambiguous", "candidates": [keys…]}` for multiple case-insensitive substring hits; `{"status": "not_found", "candidates": []}` otherwise. Never guess among ambiguous matches.
   - The `record` follows the shared schema (step 6). `source` names the **CODATA release actually used** (e.g. `"CODATA 2018 (scipy 1.11)"`) since scipy bundles whichever release its version ships — see the version caveat in Background. `uncertainty` is scipy's standard uncertainty (may be `0.0` for exactly-defined SI constants — surface that honestly rather than as missing data).

5. **`tools/asd.py`** — two functions reaching the ASD via the source-of-record landing page's Lines/Levels CGI forms, using `requests`, parsing the delimited response, and mapping value/uncertainty columns into the shared schema. Honor every constraint in the Background ASD bullet (native units cm⁻¹ / nm + air-vs-vacuum flag; no silent conversion; reference codes into `reference`; `params=` encoding; `User-Agent` + timeout + bounded retry; docstring records landing-page URL + resolved CGI params). Record the landing-page URL (<https://www.nist.gov/pml/atomic-spectra-database>) in the module docstring.
   - `get_levels(species: str, ...) -> list[dict]` → energy levels for an ion/atom (e.g. `"H I"`, `"He I"`, `"Li III"`), each with energy value + uncertainty + unit (cm⁻¹) + reference.
   - `get_transitions(species: str, ...) -> list[dict]` → transition wavelengths/frequencies with uncertainties (nm, air/vacuum flagged).
   - Build the request URL from documented ASD parameters selecting delimited output; parse defensively (skip header/comment lines, tolerate empty results), and return `{"results": [...], "source": "NIST ASD", "query": {...}}`. On HTTP error or empty result, return `{"results": [], "error": ...}` rather than raising. Route every call through `cache.py`.

5a. **`tools/targets.py` — curated, cited Dirac-exploration registry (the honest home for the precision splittings).** Per the Background scope correction, the MHz/sub-MHz splittings are **not** retrievable from ASD at the needed precision, so they live here as hand-curated records, each carrying `value + uncertainty + unit + an individual literature citation`: hydrogen `2P₃/₂–2P₁/₂` fine structure, `2S₁/₂–2P₁/₂` Lamb shift, 1S hyperfine, and helium `³P₀–³P₁` — sourced from the values the campaign already hand-typed in [Bethe_Salpeter/10_CrossComparison.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md), now made retrievable *with provenance*. Expose `list_dirac_targets() -> list[dict]` returning these in the shared schema, and **register it as a fourth MCP tool** (step 7). Each record MUST cite its source so the curated origin is never confused with an ASD/CODATA fetch (`source` like `"curated (PDG-2024)"`, `reference` the citation). The electron g-factor is delegated to `get_constant` (CODATA), not duplicated here.

6. **Shared output schema** — every quantity-bearing record returned by any tool is a dict with keys: `quantity` (str), `value` (float), `uncertainty` (float | None), `unit` (str), `source` (str), `reference` (str | None). Document this once in the README and keep it consistent across `codata.py` and `asd.py`.

7. **`server.py`** — `mcp = FastMCP("nist")`; register **four** tools — `get_levels`, `get_transitions`, `get_constant`, `list_dirac_targets` — as thin `@mcp.tool()` wrappers delegating to the `tools/` functions (mirroring [pyphysics_mcp/server.py:35-128](../../mcp_servers/pyphysics_mcp/src/pyphysics_mcp/server.py#L35-L128)); `def main(): mcp.run()`. Each wrapper's docstring states the data source so the source-of-truth (ASD vs CODATA vs curated) is visible to the calling agent.

8. **Register the server** — add a `nist` entry to [.mcp.json](../../.mcp.json) alongside `pyphysics`:
   ```json
   {
     "mcpServers": {
       "pyphysics": { "command": "pyphysics-mcp" },
       "nist": { "command": "nist-mcp" }
     }
   }
   ```

9. **`README.md`** — document install (`uv tool install ./mcp_servers/nist_mcp`), Claude Code registration, the **four** tools + shared schema, the cache-dir resolution order, the source-of-record URL, and the three-way retrieval split (ASD vs CODATA vs curated registry) so users do not mistake a curated splitting for an ASD fetch. Include per-tool verification examples (e.g. `get_constant("electron g factor")`, `get_levels("H I")`, `list_dirac_targets()`), mirroring [pyphysics_mcp/README.md](../../mcp_servers/pyphysics_mcp/README.md).

10. **Tests** — add `mcp_servers/nist_mcp/tests/`:
    - **Deterministic server smoke test** (replaces the unreliable `nist-mcp --help` check): import the package and assert the `FastMCP` app registers exactly the four expected tools — no stdio invocation, no hanging.
    - **`get_constant`** (offline): a known CODATA constant returns expected value/uncertainty/unit, all shared-schema keys present, and the `source` names a CODATA release; ambiguous/unknown names return the `status`-discriminated shape (never raise).
    - **`list_dirac_targets`** (offline): every curated record carries a non-empty `reference` and a `source` marked curated — guarding against the precision splittings being silently re-routed through ASD.
    - **ASD parser** (offline): feed a **committed response fixture** (`tests/fixtures/`) through the parser and assert value+uncertainty extraction, native-unit preservation, and the air/vacuum flag.
    - **Live ASD** (`NIST_MCP_LIVE=1`-gated, skipped by default so CI stays offline): a real `get_levels`/`get_transitions` call returns non-empty, schema-conforming records in native units. It asserts ASD coverage — **not** the MHz splittings, which ASD does not supply.

## Files to Modify

| File | Change |
|---|---|
| `mcp_servers/nist_mcp/pyproject.toml` | New. Package metadata, `nist-mcp` console script, deps `mcp`, `scipy`, `requests`; hatchling build. |
| `mcp_servers/nist_mcp/README.md` | New. Install / registration / tool docs / shared schema / verification examples; documents the ASD landing-page source-of-record and the Dirac-exploration target examples. |
| `mcp_servers/nist_mcp/src/nist_mcp/__init__.py` | New. Re-export `main` from `nist_mcp.server`. |
| `mcp_servers/nist_mcp/src/nist_mcp/__main__.py` | New. `python -m nist_mcp` entry calling `main`. |
| `mcp_servers/nist_mcp/src/nist_mcp/server.py` | New. `FastMCP("nist")` + `@mcp.tool()` wrappers for the four tools `get_levels`, `get_transitions`, `get_constant`, `list_dirac_targets`; `main()`. |
| `mcp_servers/nist_mcp/src/nist_mcp/cache.py` | New. Cache-dir resolution (`NIST_MCP_CACHE_DIR` → default) + on-disk JSON `cache_get`/`cache_put` with `CACHE_SCHEMA_VERSION` keying, timestamps, and TTL/`refresh` bypass. |
| `mcp_servers/nist_mcp/src/nist_mcp/tools/__init__.py` | New. Package marker. |
| `mcp_servers/nist_mcp/src/nist_mcp/tools/codata.py` | New. `get_constant` over `scipy.constants.physical_constants`; `status`-discriminated return; reports CODATA release in `source`. |
| `mcp_servers/nist_mcp/src/nist_mcp/tools/asd.py` | New. `get_levels`, `get_transitions` via the ASD landing-page Lines/Levels CGI forms; native units (cm⁻¹ / nm + air-vacuum flag), reference codes, `params=` encoding, `User-Agent`/timeout/retry. |
| `mcp_servers/nist_mcp/src/nist_mcp/tools/targets.py` | New. `DIRAC_TARGETS` curated+cited registry of the campaign's MHz/sub-MHz splittings + `list_dirac_targets`. |
| `mcp_servers/nist_mcp/tests/test_codata.py` | New. Offline `get_constant`: value/uncertainty/unit + schema + CODATA-release + `status` shapes. |
| `mcp_servers/nist_mcp/tests/test_asd.py` | New. Offline parser test (fixture) + network-gated live test (`NIST_MCP_LIVE=1`) asserting native-unit ASD coverage (not the splittings). |
| `mcp_servers/nist_mcp/tests/test_targets.py` | New. Offline test that every curated Dirac target carries a citation + curated `source`. |
| `mcp_servers/nist_mcp/tests/test_server.py` | New. Deterministic smoke test: the `FastMCP` app registers exactly the four expected tools. |
| `mcp_servers/nist_mcp/tests/fixtures/asd_levels_sample.tsv` | New. Captured ASD delimited response, committed so the parser is testable offline and NIST form changes surface as failures. |
| `.mcp.json` | Add a `nist` server entry (`"command": "nist-mcp"`) next to `pyphysics`. |

## Dependencies

- **New Python deps (declared in `mcp_servers/nist_mcp/pyproject.toml` only, not the root project):** `mcp>=1.0` (FastMCP, already used by `pyphysics_mcp`), `scipy` pinned to a release whose bundled CODATA matches the campaign's **CODATA-2018** (e.g. `scipy>=1.11,<1.14`; verify the bundled release at implementation time and record it), `requests>=2.31` (ASD HTTP queries). `pytest` as a dev/test dependency. Deliberately **not** depending on `astropy`/`astroquery` (see Background ASD bullet).
- **External service:** NIST Atomic Spectra Database CGI endpoints (`physics.nist.gov/cgi-bin/ASD/lines1.pl`, `.../energy1.pl`) — network access required for live ASD calls; CODATA path is fully offline.
- **Env vars:** `NIST_MCP_CACHE_DIR` (optional, overrides default cache location); `NIST_MCP_LIVE` (optional, opt-in for the live network test).

## Acceptance Criteria

- [ ] `uv tool install ./mcp_servers/nist_mcp` succeeds and installs a `nist-mcp` console script, and a deterministic smoke test confirms the `FastMCP` app registers exactly four tools (`get_levels`, `get_transitions`, `get_constant`, `list_dirac_targets`) — without relying on `--help`.
- [ ] `.mcp.json` contains both a `pyphysics` and a `nist` server entry, and remains valid JSON.
- [ ] `get_constant(name)` returns a single `status`-discriminated shape (`ok` / `ambiguous` / `not_found`); the `ok` record carries `value`, `uncertainty`, `unit`, and a `source` that names the **CODATA release used** (e.g. "CODATA 2018").
- [ ] `get_levels(species)` and `get_transitions(species)` return ASD records in **native units** (levels cm⁻¹; lines nm with an air/vacuum flag) with `value`, `uncertainty`, `source` ("NIST ASD"), and `reference` where available — and return an empty-results record (not an exception) on no-match or network failure. No silent unit conversion.
- [ ] Every tool's returned quantity conforms to the shared schema (`quantity, value, uncertainty, unit, source, reference`).
- [ ] ASD lookups are cached on disk keyed by query + `CACHE_SCHEMA_VERSION`; a repeated identical query is served from cache without a second network call, and `refresh=True` bypasses the cache.
- [ ] `list_dirac_targets()` returns the campaign's MHz/sub-MHz splittings (Lamb shift, 1S hyperfine, H 2P fine structure, He ³P) as curated records, **each with a non-empty citation and a `source` marked curated** — and the plan/README state explicitly that these are **not** sourced from ASD (which lacks them at this precision).
- [ ] End-to-end Dirac-exploration coverage is demonstrated honestly: the precision splittings via `list_dirac_targets()` (cited), the electron g-factor via `get_constant` (CODATA), and native-unit level/line data for a test species (e.g. `H I` or `Li III`) via the ASD path — matching the measurement column of [Bethe_Salpeter/10_CrossComparison.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md) by source, not by pretending ASD supplies the splittings.
- [ ] `mcp_servers/nist_mcp/README.md` documents install, Claude Code registration, the four tools, the shared schema, the source-of-record URL, the three-way retrieval split, and cache-dir resolution.
- [ ] Offline tests pass without network; the live ASD test is skipped unless `NIST_MCP_LIVE=1`.

## Testing

Reviewer commands (run from the repo root):

```bash
# Install the new server.
uv tool install ./mcp_servers/nist_mcp

# Deterministic registration check (no stdio, no hang) — asserts the four tools exist.
uv run --project mcp_servers/nist_mcp python -c "from nist_mcp.server import mcp; import asyncio; print(sorted(t.name for t in asyncio.run(mcp.list_tools())))"

# Offline test suite (no network).
uv run --project mcp_servers/nist_mcp pytest mcp_servers/nist_mcp/tests -q

# Opt-in live ASD check (hits physics.nist.gov).
NIST_MCP_LIVE=1 uv run --project mcp_servers/nist_mcp pytest mcp_servers/nist_mcp/tests/test_asd.py -q

# Validate .mcp.json is still well-formed and has both servers.
uv run python -c "import json; d=json.load(open('.mcp.json')); print(sorted(d['mcpServers']))"
```

Tests added:

- `mcp_servers/nist_mcp/tests/test_server.py` — deterministic smoke test asserting the `FastMCP` app registers exactly the four expected tools.
- `mcp_servers/nist_mcp/tests/test_codata.py` — asserts `get_constant("electron g factor")` (or another stable CODATA key) returns the expected value/uncertainty/unit, all shared-schema keys present, and a `source` naming the CODATA release; asserts ambiguous/unknown names return the `status`-discriminated shape rather than raising.
- `mcp_servers/nist_mcp/tests/test_targets.py` — asserts every `list_dirac_targets()` record carries a non-empty citation and a curated `source`, guarding the scope boundary (splittings are curated, not ASD-sourced).
- `mcp_servers/nist_mcp/tests/test_asd.py` — offline test feeding the committed `tests/fixtures/asd_levels_sample.tsv` through the `asd.py` parser and asserting value+uncertainty extraction, native-unit preservation, and the air/vacuum flag; plus a `NIST_MCP_LIVE`-gated test performing a real `get_levels`/`get_transitions` call for a test species (e.g. `H I` / `Li III`) and asserting non-empty, schema-conforming native-unit results — explicitly *not* asserting the MHz splittings.

After registering with Claude Code (restart required), confirm the deferred tools list includes `mcp__nist__get_constant`, `mcp__nist__get_levels`, `mcp__nist__get_transitions`, and `mcp__nist__list_dirac_targets`.
