# `nist.*` namespace — refresh policy addendum

This file extends the umbrella schema at [`mcp_servers/precision_data_mcp/REFRESH_POLICY.md`](../../../REFRESH_POLICY.md) with the namespace-specific source URLs, refresh cadence, and audit-log conventions for the `nist.*` namespace (NIST ASD + CODATA + curated Dirac-target splittings).

Migrated under the umbrella per issue [#99](https://github.com/temoTxt/PyPhysics/issues/99); behavioural contract preserved verbatim from the standalone `nist_mcp` ([issue #90](https://github.com/temoTxt/PyPhysics/issues/90)).

## Sources

| Tool | Upstream source | Edition tracking | Refresh route |
|---|---|---|---|
| `nist.get_constant` | CODATA via `scipy.constants` (offline) | `codata_edition` derived from the installed `scipy` version | bump `scipy>=1.X` in `pyproject.toml` when a new CODATA adjustment lands |
| `nist.get_levels` | NIST ASD `energy1.pl` CGI endpoint | `asd_query_timestamp` (ASD has no edition number) | `refresh=True` on the call |
| `nist.get_transitions` | NIST ASD `lines1.pl` CGI endpoint | `asd_query_timestamp` | `refresh=True` on the call |
| `nist.list_dirac_targets` | Curated, individually-cited primary sources | per-target `cite_key` from `Bibliography/` | hand-edit `tools/targets.py` + the corresponding bib stub |

Note that `nist.get_constant` is **offline** — refresh cadence is the underlying `scipy` package release, not a live network call.

## Bibliography stubs

The two umbrella-level bib stubs covering this namespace's non-paper sources:

- [`nist_asd_2024`](../../../../../../Roadmapping/History/Bibliography/Retrospective/nist_asd_2024.md) — NIST Atomic Spectra Database, current edition. `firstauthor` is the issuing body (NIST).
- [`codata2018_constants`](../../../../../../Roadmapping/History/Bibliography/Retrospective/codata2018_constants.md) — CODATA 2018 fundamental physical constants adjustment.

Additional `cite_key`s appear in `tools/targets.py` for each individually-cited precision splitting (hydrogen fine structure, Lamb shift, 1S hyperfine, helium 3P fine structure, etc.). These were scaffolded during the original [#90](https://github.com/temoTxt/PyPhysics/issues/90) work and are preserved unchanged by this migration.

## Disagreement-representation status

The four `nist.*` tools currently return *single* values per query because their upstream sources are individually authoritative (NIST ASD is the citing body; CODATA values come from a single committee adjustment; the Dirac targets are individually cited).

**Re-routing on issue #97 landing:** the `nist.list_dirac_targets` registry will be retired in favour of `qed.*` tools when [#97](https://github.com/temoTxt/PyPhysics/issues/97) ships its curated JSON database. At that point the *disagreement-representation* rule will apply to the Lamb-shift and 1S-hyperfine values (multiple measurement papers; multiple theory comparators), and the records here will be deprecated with a forwarding pointer rather than deleted. The `cite_key`s themselves stay valid since they live in `Bibliography/`.

## Refresh log

`./refresh_log.md` records the chronological migration + any future precision updates to the Dirac-target registry. The first entry is the [#99](https://github.com/temoTxt/PyPhysics/issues/99) migration itself.
