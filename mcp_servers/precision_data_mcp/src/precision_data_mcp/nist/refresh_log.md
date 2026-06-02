# `nist.*` namespace — refresh log

Chronological record of refreshes, edition bumps, and citation diffs for the `nist.*` namespace. Per the umbrella [`REFRESH_POLICY.md`](../../../REFRESH_POLICY.md) §7, every namespace maintains this log so future reviewers can reconstruct exactly what data the campaign consumed at any prior point in time.

## 2026-05-30 — migration under umbrella `precision_data_mcp` (issue #99)

- **Cause:** structural migration per umbrella issue [#92](https://github.com/temoTxt/PyPhysics/issues/92) §Resolved-decisions #8 (`nist_mcp` becomes a child of the umbrella as `nist.*` namespace).
- **Behavioural diff:** none. The four tools' signatures, cache semantics, return shape, and data values are preserved verbatim from the standalone `nist_mcp` ([#90](https://github.com/temoTxt/PyPhysics/issues/90)). Only the server registration changes.
- **Tool renames:**
  - `get_constant` → `nist.get_constant`
  - `get_levels` → `nist.get_levels`
  - `get_transitions` → `nist.get_transitions`
  - `list_dirac_targets` → `nist.list_dirac_targets`
- **Cache invalidations:** none. Any pre-existing on-disk cache files written by the standalone `nist_mcp` remain valid and are re-used by the umbrella's `nist.*` tools (same cache schema).
- **New bib stubs scaffolded:**
  - [`nist_asd_2024`](../../../../../../Roadmapping/History/Bibliography/Retrospective/nist_asd_2024.md)
  - [`codata2018_constants`](../../../../../../Roadmapping/History/Bibliography/Retrospective/codata2018_constants.md)
- **Caller updates:**
  - [`Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md`](../../../../../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md) — tool references updated to the `nist.*` prefixed form.

Subsequent entries record any precision updates to the Dirac-target registry, CODATA edition bumps via `scipy`, or namespace re-routing (e.g. when [#97](https://github.com/temoTxt/PyPhysics/issues/97) ships and the Dirac-target registry deprecates in favour of `qed.*`).
