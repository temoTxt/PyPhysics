# Umbrella refresh-policy schema for `precision_data_mcp`

This file defines the **shared schema** all per-namespace refresh policies extend. The umbrella server's caching + refresh semantics are uniform across `nist.*`, `pdg.*`, `qed.*`, `nuclear.*`, `flag.*`, `astro.*`, `cosmo.*`. Each namespace's own `<namespace>/REFRESH_POLICY.md` records the namespace-specific source(s), refresh cadence, and audit log; the rules below apply unchanged across all of them.

Established as part of issue [#92](https://github.com/temoTxt/PyPhysics/issues/92) §Resolved-decisions #6 and the [#99](https://github.com/temoTxt/PyPhysics/issues/99) foundation PR.

## 1. The six guarantees every namespace must satisfy

1. **Every returned value carries `{value, uncertainty, unit, source, retrieved_at, cache_key, value_class, safe_for_model_verification}`.** The `source` field is a `cite_key` from `Roadmapping/History/Bibliography/`. Free-form URLs / DOIs / verbal descriptions in `source` are forbidden — they bypass Crocco-rule-3 citation verification.
2. **Cache by default; `refresh=True` opt-in.** No tool may make a live external request unless the caller explicitly passes `refresh=True` or the on-disk cache has no entry for the requested `cache_key`. Cached snapshots remain valid indefinitely; staleness is surfaced at the metadata layer, never silently corrected.
3. **Cached snapshots are flat JSON files under `mcp_servers/precision_data_mcp/<namespace>/cache/`.** No SQLite, no YAML registry, no database backend. The flat-JSON-per-namespace convention matches the existing pattern from `nist_mcp/targets.py` and keeps cache diffs auditable in git for reproducibility.
4. **Disagreement representation is preserved.** When a quantity has multiple authoritative values from competing methods or independent measurements (proton-radius puzzle, Hubble tension, muon g-2 lattice-vs-experiment), the tool returns *all* values with their `method` / `source` labels — never an averaged "best" value that hides the disagreement.
5. **Bibliography stubs accompany every cited primary source.** A `cite_key` returned in a `source` field must correspond to an existing file under `Roadmapping/History/Bibliography/{Primary,Retrospective}/`. New sources require scaffolding via `scaffold_bib_note.py` as part of the namespace's PR — citations are never fabricated.
6. **Theoretical-vs-experimental status is explicit.** Every value carries a `value_class` field placing it on the experimental-vs-theoretical spectrum, plus a derived `safe_for_model_verification` boolean. When unsafe, a `warning` field is auto-attached. Consumers (verification documents, agents) must refuse to use a value with `safe_for_model_verification == False` as the experimental anchor when comparing theoretical predictions — such values are theory comparators only. The taxonomy and the implementation live at [`src/precision_data_mcp/safety.py`](src/precision_data_mcp/safety.py).

### Value-class taxonomy (guarantee 6 detail)

| value_class | safe_for_model_verification | Meaning |
|---|---|---|
| `experimental` | **True** | Single-source direct experimental measurement (e.g. Pohl 2010 muonic-H Lamb shift) |
| `experimental_world_average` | **True** | PDG / equivalent average of multiple independent measurements (e.g. muon mass) |
| `definitional` | **True** | Convention not a measurement target (e.g. electron charge = -1e, spin = 1/2) |
| `codata_adjusted` | **False** | CODATA combines experimental + theoretical inputs; do NOT use as anchor for QED / SM tests |
| `theoretical_lattice_qcd` | **False** | Lattice-QCD calculation (FLAG, BMW); theory comparator only |
| `theoretical_sm_prediction` | **False** | Perturbative or data-driven SM prediction (Theory Initiative a_mu); theory comparator only |
| `theoretical_ritz_value` | **False** | Computed from Rydberg formula + theoretical quantum defects (NIST ASD hydrogen levels); not a measurement |
| `theoretical_calculation` | **False** | Generic theoretical-calculation catch-all; theory comparator only |

## 2. Cache file layout

```
mcp_servers/precision_data_mcp/
├── REFRESH_POLICY.md               (this file — umbrella schema)
├── nist/
│   ├── REFRESH_POLICY.md           (per-namespace addendum: NIST ASD + CODATA)
│   └── cache/
│       └── <cache_key>.json
├── pdg/
│   ├── REFRESH_POLICY.md           (per-namespace addendum: PDG)
│   └── cache/
│       └── <cache_key>.json
├── qed/
│   ├── REFRESH_POLICY.md
│   ├── data.json                   (curated seed JSON for Lamb shift, hyperfine, g-factor, transitions)
│   ├── refresh_log.md              (Crossref/NIST-ASD reconciliation history)
│   └── cache/
│       └── <cache_key>.json
├── nuclear/
│   ├── REFRESH_POLICY.md
│   └── cache/
└── flag/
    ├── REFRESH_POLICY.md
    └── cache/
```

Each cache file is one JSON record. The cache key is a deterministic hash of the tool's argument tuple — see §3.

## 3. Cache-key derivation

A `cache_key` is the SHA-256 hexdigest (truncated to 16 hex chars) of a canonicalised JSON dump of the tool's argument tuple, prefixed with the tool's fully-qualified name:

```python
def cache_key(tool_name: str, args: dict) -> str:
    payload = json.dumps({"tool": tool_name, "args": args}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
```

This makes the cache key reproducible across invocations and across machines, and the truncation keeps filenames tractable. Collisions are vanishingly improbable at 16 hex chars across the entire precision-physics universe; we accept the risk.

## 4. `retrieved_at` semantics

`retrieved_at` is an ISO-8601 UTC timestamp recording when the upstream source was first queried for the cached value. It is **not** refreshed on cache hit — only on cache miss or on `refresh=True`. This makes `retrieved_at` a proper provenance record rather than a last-accessed timestamp.

## 5. `cache_key` semantics

`cache_key` is returned in *every* tool response so that downstream callers (verification documents, agents) can record the exact value they consumed. When two verification documents both cite a quantity, comparing the `cache_key` fields confirms they consumed the same value — no risk of subtle drift between consumers.

## 6. Source-revision tracking

Some upstream sources have versioned releases (PDG yearly *Review of Particle Physics*; CODATA periodic adjustments; NIST ASD updates; FLAG report editions). The cached JSON records the source revision in a `source_revision` sub-field:

- PDG: `pdg_edition: "2024"` — corresponds to a specific PDG MCID database snapshot.
- CODATA: `codata_edition: "2018"` — currently bundled with scipy.
- NIST ASD: `asd_query_timestamp: "<ISO-8601>"` since ASD does not publish edition numbers.
- FLAG: `flag_edition: "2024"`.

When a new edition lands, namespace owners add a corresponding entry to `<namespace>/refresh_log.md` recording the diff against the previous edition before any cached snapshots are refreshed. Verification documents that consumed the old edition retain their old `cache_key` references; they are not silently re-pointed to the new edition.

## 7. The `refresh_log.md` per namespace

Each namespace maintains a chronological log at `<namespace>/refresh_log.md` recording:

- Date of refresh.
- Source-revision diff against the previous cached state.
- Which cache keys were invalidated and which retained.
- Any precision updates that changed a `value` or `uncertainty`.
- Citation diff (any new `cite_key`s scaffolded; any retired).

This is the Crocco-expectation-3 audit artifact: it lets a future reviewer reconstruct *exactly* what data the campaign consumed at any prior point in time.

## 8. Crocco compliance

MCP look-up is **pragmatic AI use** — no prompt-of-record needed for tool invocations. The verification documents that consume tool outputs continue to gate substantive verdicts via `<!-- TODO: human reviews and fills in -->` blocks per the existing campaign convention. This file's rules do not change the campaign's compliance posture; they only sharpen *which* data went into a verdict, which makes substantive-block review easier for the human reviewer.

## 9. Schema evolution

When this schema needs to evolve (e.g., a sixth required field on every returned record, a change to the cache-key derivation), the change goes through a regular PR with a `SCHEMA-CHANGE` label and a migration script under `mcp_servers/precision_data_mcp/migrations/`. Any prior cached snapshots stay valid under the old schema; the umbrella reads both transparently for at least two release cycles before the old schema is retired.
