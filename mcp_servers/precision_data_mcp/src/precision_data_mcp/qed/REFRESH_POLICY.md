# `qed.*` namespace — refresh policy addendum

Extends the umbrella schema at [`mcp_servers/precision_data_mcp/REFRESH_POLICY.md`](../../../REFRESH_POLICY.md) for bound-state QED precision observables.

Created under issue [#97](https://github.com/temoTxt/csv) (sic — typo intentional; actual link: [#97](https://github.com/temoTxt/PyPhysics/issues/97)). Sub-issue of umbrella [#92](https://github.com/temoTxt/PyPhysics/issues/92). **Highest-impact sub-issue** of the umbrella because the Bethe-Salpeter campaign ([#50](https://github.com/temoTxt/PyPhysics/issues/50)), the Lamb-shift one-loop dual-QED work ([#55](https://github.com/temoTxt/PyPhysics/issues/55)), the hydrogenic Li²⁺ extension ([#78](https://github.com/temoTxt/PyPhysics/issues/78)), the hydrogenic-Z-scan g-factor work ([#82](https://github.com/temoTxt/PyPhysics/issues/82)), and the spectroscopy-collaborator engagement ([#88](https://github.com/temoTxt/PyPhysics/issues/88)) all directly consume this namespace.

## Source — hand-curated JSON

Unlike PDG (which has a JSON API) and NIST ASD (which has tabular HTML), bound-state QED precision data is scattered across PRL / PRA / Nature / Science primary-source papers with no canonical machine-readable database. The seed at [`data.json`](data.json) is therefore hand-curated, with each entry's `source` field a `cite_key` pointing at the primary-source paper.

## Initial seed coverage

| Species | Observables seeded | Primary sources cited |
|---|---|---|
| H | 2S1/2-2P1/2 Lamb shift (multi-value: CODATA + Lundeen direct RF); 1s2S1/2 hyperfine (21-cm); 2s2S1/2 hyperfine (Dirac-paper-cited); 2P3/2-2P1/2 fine structure; 1S-2S transition; 1S-2S H/D isotope shift; 1S-3S; 2S-4P1/2; 2S-8S1/2 | `codata2018_constants`, `lundeen1981_2s2p_microwave`, `essen1957_h_hyperfine`, `heberle1956_h_2s_hyperfine`, `hagley1994_h_fine_structure`, `parthey2011_h_1s2s`, `parthey2010_h_d_isotope_shift`, `fleurbaey2018_1s_3s_spectroscopy`, `beyer2017_2s4p`, `debeauvoir1997_2s_8s_8d` |
| He II (He⁺) | 1s2S1/2 bound-electron g-factor | `kohler2016_he_g_factor` |
| Li III (Li²⁺) | 1s2S1/2 hyperfine | `beckert2007_li2_hyperfine` |
| Si XIV (Si¹³⁺) | 1s2S1/2 bound-electron g-factor | `sturm2013_si13_g_factor` (issue #82 target) |
| muonic_H | 2S1/2-2P3/2 Lamb shift | `pohl2010_muonic_hydrogen` (proton-radius-puzzle trigger) |
| muonium | 1s2S1/2 hyperfine | `liu1999_muonium_hyperfine` |
| positronium | 1s ortho-para splitting | `ishida2014_positronium` |

Adding species or observables is a mechanical follow-up PR: extend `data.json`, scaffold the new `cite_key`'s bib stub, and log the diff in [`refresh_log.md`](refresh_log.md). The seed structure deliberately mirrors how a future automated reconciliation (against Crossref / NIST ASD / individual journal APIs) would shape its data, so the migration is structural rather than contract-breaking.

## Edition tracking

The seed's `$source_revision.data_revision` records the curation date; precision-update events bump it and append to `refresh_log.md`. Existing cached `cache_key` records consumed by prior verification documents remain valid under the old revision per umbrella §6 — they are not silently re-pointed when a higher-precision measurement lands.

## Disagreement-representation cases

The current seed does **not** include multi-value disagreements at the qed.* namespace level (each entry is a single primary-source value). The proton-radius puzzle disagreement lives in the `nuclear.*` namespace under issue [#98](https://github.com/temoTxt/PyPhysics/issues/98), because the puzzle is about *which charge-radius value to extract* from the Lamb-shift measurement, not about the Lamb-shift measurement itself.

Future qed.* multi-value cases (e.g. if a second independent muonic-H Lamb-shift measurement at a competing facility lands and disagrees with Pohl 2010 at >1σ) will follow the umbrella's §Resolved-decisions #5 rule and return `{"values": [...]}` rather than a single averaged value.

## Re-routing of `nist.list_dirac_targets`

Per [`nist/REFRESH_POLICY.md`](../nist/REFRESH_POLICY.md) §Disagreement-representation status, the `nist.list_dirac_targets` registry will be deprecated in favour of `qed.*` tools as the curated data here supersedes it. The deprecation is **not** part of this PR — it lands as a follow-up after the demonstrator refactor of [`Bethe_Salpeter/05_LambShift.md`](../../../../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/05_LambShift.md) validates the qed.* path end-to-end.

## Programmatic refresh

The qed.* values are not refreshed automatically — they live in the hand-curated JSON. The refresh procedure when a higher-precision measurement lands:

1. Identify the new value and its primary-source paper.
2. Scaffold a bib stub for the new paper via `scaffold_bib_note.py`.
3. Update the affected entry in `data.json`; preserve the old cite_key as a comment in the entry's `_note` for audit.
4. Append a chronological entry to `refresh_log.md` recording the diff (old value vs new, old uncertainty vs new, citation change, any campaigns affected).
5. Bump `$source_revision.data_revision` to the current ISO date.

## Bibliography stubs

Seven new primary-source stubs scaffolded as part of [#97](https://github.com/temoTxt/PyPhysics/issues/97):

- `essen1957_h_hyperfine` — Essen et al. 1957 (foundational 21-cm maser measurement)
- `parthey2011_h_1s2s` — Parthey-Hänsch et al. 2011 (H 1S-2S absolute frequency)
- `kohler2016_he_g_factor` — Kohler et al. 2016 (He⁺ g-factor)
- `beckert2007_li2_hyperfine` — Beckert et al. 2007 (Li²⁺ hyperfine)
- `sturm2013_si13_g_factor` — Sturm et al. 2013 (Si¹³⁺ g-factor)
- `pohl2010_muonic_hydrogen` — Pohl et al. 2010 (muonic-H Lamb shift)
- `liu1999_muonium_hyperfine` — Liu et al. 1999 (muonium hyperfine)
- `ishida2014_positronium` — Ishida et al. 2014 (positronium hyperfine)

## Refresh log

See [`refresh_log.md`](refresh_log.md).
