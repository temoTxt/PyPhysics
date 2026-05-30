# precision-data-mcp

Umbrella MCP server exposing precision-physics experimental data — particle masses (PDG), bound-state QED observables (Lamb shift, hyperfine, g-factor, transition precision), nuclear / particle structure (charge radii incl. proton-radius puzzle, magnetic moments), FLAG lattice-QCD averages, NIST ASD energy levels + CODATA constants, and (deferred) astronomical / cosmological observables — as Claude-callable tools, with every returned quantity carrying **value, uncertainty, unit, and source**.

Companion to [`pyphysics-mcp`](../pyphysics_mcp/). Built for umbrella issue [#92](https://github.com/temoTxt/PyPhysics/issues/92); foundation PR (this PR) closes [#99](https://github.com/temoTxt/PyPhysics/issues/99).

## Architectural posture

Every value returned across every namespace carries `{value, uncertainty, unit, source, retrieved_at, cache_key}`. The `source` field is a `cite_key` from `Roadmapping/History/Bibliography/`. The shared caching + refresh semantics are documented in [`REFRESH_POLICY.md`](REFRESH_POLICY.md); each namespace extends that schema in its own `<namespace>/REFRESH_POLICY.md`.

The umbrella consolidates what would otherwise be multiple separate MCP servers (PDG, bound-state QED, nuclear, FLAG, NIST, ...) into a single runtime with namespaced tool names. The previously-shipped `nist_mcp` is migrated under this umbrella as the `nist.*` namespace per umbrella [§Resolved-decisions #8](https://github.com/temoTxt/PyPhysics/issues/92).

## Tool surface

| Namespace | Sub-issue | Status as of foundation PR |
|---|---|---|
| `nist.*` | [#99](https://github.com/temoTxt/PyPhysics/issues/99) (this PR) | ✅ migrated — 4 tools |
| `pdg.*` | [#96](https://github.com/temoTxt/PyPhysics/issues/96) | planned |
| `qed.*` | [#97](https://github.com/temoTxt/PyPhysics/issues/97) | planned |
| `nuclear.*` + `flag.*` | [#98](https://github.com/temoTxt/PyPhysics/issues/98) | planned |
| `astro.*` | deferred per umbrella #92 | not started |
| `cosmo.*` | deferred per umbrella #92 | not started |

### `nist.*` — migrated from standalone `nist_mcp` ([#90](https://github.com/temoTxt/PyPhysics/issues/90))

| Tool | Returns | Source |
|---|---|---|
| `nist.get_constant(name)` | `{"status": "ok"\|"ambiguous"\|"not_found", ...}` | CODATA via `scipy.constants` (offline) |
| `nist.get_levels(species)` | ASD energy levels in **cm⁻¹** (native, no conversion) | NIST ASD `energy1.pl` |
| `nist.get_transitions(species)` | ASD transitions in **nm**, with air/vacuum flag | NIST ASD `lines1.pl` |
| `nist.list_dirac_targets()` | Curated, cited precision splittings for Dirac-equation exploration | hand-curated registry (will re-route to `qed.*` when [#97](https://github.com/temoTxt/PyPhysics/issues/97) lands) |

The three-way split (CODATA / ASD / curated) is preserved from the original `nist_mcp` design — see the [`nist/REFRESH_POLICY.md`](src/precision_data_mcp/nist/REFRESH_POLICY.md) addendum for source URLs, edition tracking, and refresh routes.

## Shared record schema

Every quantity-bearing record across every namespace is:

```python
{
    "quantity": str,        # e.g. "electron g factor", "1S Lamb shift"
    "value": float,
    "uncertainty": float,
    "unit": str,
    "source": str,          # cite_key from Roadmapping/History/Bibliography/
    "retrieved_at": str,    # ISO-8601 UTC; refreshed on cache miss / refresh=True
    "cache_key": str,       # SHA-256 hexdigest (16-char truncated) of canonical args
}
```

For multi-value disagreements (proton-radius puzzle, Hubble tension, muon g-2 lattice-vs-experiment), the tool returns a *list* of such records — one per measurement / method / theory comparator — never an averaged value. See [`REFRESH_POLICY.md`](REFRESH_POLICY.md) §1, guarantee 4.

## Install

```bash
uv tool install ./mcp_servers/precision_data_mcp
```

## Register

In `.mcp.json` (this repo's root) or `~/.claude/mcp.json`:

```json
{
  "mcpServers": {
    "precision_data": {"command": "precision-data-mcp"}
  }
}
```

The previously-registered `nist` server entry is removed by this PR.

## Run

```bash
precision-data-mcp           # via console_script
python -m precision_data_mcp # equivalent
```

The server speaks MCP over stdio.

## Tests

```bash
cd mcp_servers/precision_data_mcp
uv sync --extra test
uv run pytest
```

Tests for the `nist.*` namespace live under [`tests/nist/`](tests/nist/) and verify the migration preserves the original `nist_mcp` behaviour verbatim.

## Crocco compliance

MCP look-up is **pragmatic AI use** — no prompt-of-record needed for tool invocations. Verification documents that consume tool outputs continue to gate substantive verdicts via `<!-- TODO: human reviews and fills in -->` blocks per the campaign convention. See [`REFRESH_POLICY.md`](REFRESH_POLICY.md) §8.
