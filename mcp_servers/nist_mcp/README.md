# nist-mcp

MCP server exposing NIST experimental reference data — the [NIST Atomic Spectra Database (ASD)](https://www.nist.gov/pml/atomic-spectra-database) and CODATA fundamental physical constants — as Claude-callable tools, with every returned quantity carrying **value, uncertainty, unit, and source**.

Companion to [`pyphysics-mcp`](../pyphysics_mcp/). Built for issue [#90](https://github.com/temoTxt/PyPhysics/issues/90); plan: [`.dev/tasks/90_nist_mcp.md`](../../.dev/tasks/90_nist_mcp.md).

## Source-of-record

All ASD data is pulled from the landing page <https://www.nist.gov/pml/atomic-spectra-database>, whose Lines/Levels query forms resolve to the CGI endpoints `energy1.pl` (levels) and `lines1.pl` (lines).

## The three-way retrieval split (read this first)

Not every value comes from the same place, and the server does not pretend otherwise:

| Quantity kind | Tool | Source |
|---|---|---|
| Fundamental constants (e.g. electron g-factor) | `get_constant` | CODATA via `scipy.constants` (offline) |
| Atomic energy levels (cm⁻¹), transitions (nm) | `get_levels`, `get_transitions` | NIST ASD |
| MHz/sub-MHz precision splittings (Lamb shift, hyperfine, fine structure) | `list_dirac_targets` | **Curated + individually cited** — *not* ASD |

The precision splittings are **not** retrievable from the ASD general tables at the needed precision, and for light elements (hydrogen) ASD level energies are largely *theoretical/Ritz* values rather than independent measurements. Those splittings therefore live in a curated, cited registry — never scraped from ASD and presented as measurements.

## Tools

| Tool | Returns |
|---|---|
| `get_constant(name)` | `{"status": "ok"\|"ambiguous"\|"not_found", ...}`; the `ok` record names the CODATA release in `source`. |
| `get_levels(species)` | ASD energy levels in **cm⁻¹** (native, no conversion). |
| `get_transitions(species)` | ASD transitions in **nm**, with an **air/vacuum** flag in `unit`. |
| `list_dirac_targets()` | Curated, cited precision splittings for the Dirac-equation exploration. |

### Shared schema

Every quantity-bearing record is `{quantity, value, uncertainty, unit, source, reference}`.

## Install

```bash
uv tool install ./mcp_servers/nist_mcp
```

Verify (deterministic, no stdio hang):

```bash
uv run --project mcp_servers/nist_mcp python -c \
  "from nist_mcp.server import mcp; import asyncio; print(sorted(t.name for t in asyncio.run(mcp.list_tools())))"
# -> ['get_constant', 'get_levels', 'get_transitions', 'list_dirac_targets']
```

## Register with Claude Code

Already wired in this repo's [`.mcp.json`](../../.mcp.json):

```json
{
  "mcpServers": {
    "pyphysics": {"command": "pyphysics-mcp"},
    "nist": {"command": "nist-mcp"}
  }
}
```

After restart, the deferred tools list includes `mcp__nist__*`.

## Caching

ASD lookups are cached on disk under `NIST_MCP_CACHE_DIR` (default `~/.cache/nist_mcp/`), keyed by query + `CACHE_SCHEMA_VERSION`. Pass `refresh=True` to bypass. Clear by deleting the directory.

## Tests

```bash
# Offline (no network).
uv run --project mcp_servers/nist_mcp pytest mcp_servers/nist_mcp/tests -q

# Opt-in live ASD check (hits physics.nist.gov).
NIST_MCP_LIVE=1 uv run --project mcp_servers/nist_mcp pytest mcp_servers/nist_mcp/tests/test_asd.py -q
```

The offline ASD tests run against committed fixtures under `tests/fixtures/`, captured from live H I queries, so a NIST form change surfaces as a test failure.

## Examples

```
Use the nist MCP get_constant tool with name="electron g factor".
Use the nist MCP get_levels tool with species="H I".
Use the nist MCP list_dirac_targets tool.
```
