# `pdg.*` namespace — refresh policy addendum

Extends the umbrella schema at [`mcp_servers/precision_data_mcp/REFRESH_POLICY.md`](../../../REFRESH_POLICY.md) for Particle Data Group reference values.

Created under issue [#96](https://github.com/temoTxt/PyPhysics/issues/96). Sub-issue of umbrella [#92](https://github.com/temoTxt/PyPhysics/issues/92).

## Source

The current seed under [`data.json`](data.json) was hand-curated from the **PDG 2024 Review of Particle Physics** (`pdg_edition: "2024"`). Each particle entry's `source` field is `pdg2024_review`; per-measurement entries (electron g-2 from Fan 2023, muon g-2 from Fermilab E989 / BMW lattice QCD / Theory Initiative) cite their primary-source `cite_key` directly.

## Edition tracking

PDG publishes a yearly *Review of Particle Physics*. The seed's `$source_revision.pdg_edition` records which edition the values came from; bumping to a new edition is a deliberate refresh per the umbrella's §6 source-revision rule:

1. Hand-update `data.json` against the new edition's values.
2. Append a chronological entry to [`refresh_log.md`](refresh_log.md) recording the diff (which values changed; which `cite_key`s were retired; which were added).
3. Existing cached `cache_key` records consumed by prior verification documents remain valid under the old edition; they are not silently re-pointed.

For programmatic refreshes against the PDG 2024-vintage JSON API (<https://pdgapi.lbl.gov/>), see §4 below — currently the seed is hand-curated rather than auto-fetched.

## Disagreement-representation cases

The seed surfaces two canonical multi-value disagreements per the umbrella's §Resolved-decisions #5:

1. **Neutron lifetime puzzle.** Ultracold-neutron-bottle measurements (~878 s) disagree with beam measurements (~888 s) at the ~4σ level. The PDG world average masks this; the seed returns both with method labels.
2. **Muon g-2 lattice-vs-experiment tension.** As of 2024 the Fermilab E989 + BNL E821 combined experimental value (a_μ = 1.16592059(22)×10⁻³) disagrees with the Theory Initiative 2020 data-driven Standard Model prediction at ~5σ, *but agrees* with the BMW 2021 lattice-QCD calculation at <1σ. All three are returned with `method` labels via `pdg.get_anomaly(13)`.

Both cases are tracked in `refresh_log.md` for any post-2024 updates (the Theory Initiative 2025 update is expected to reduce the tension; track when it lands).

## Programmatic-access route (deferred)

The current seed is hand-curated. The PDG 2024-vintage JSON API at <https://pdgapi.lbl.gov/> and the `pdg` PyPI package (<https://pypi.org/project/pdg/>) offer automated refresh:

```python
# Sketched, not yet implemented — see issue #96 follow-ups
import pdg
api = pdg.connect()
electron = api.get_particle_by_mcid(11)
# -> {mass, lifetime, ...} from the live PDG database
```

Switching the seed to live-fetched values is mechanical follow-up work; the hand-curated JSON's structure deliberately mirrors what the `pdg` package would return, so the migration is a function rewrite rather than a contract change.

## Bibliography stubs

Non-paper sources cited from the PDG seed:

- [`pdg2024_review`](../../../../../../Roadmapping/History/Bibliography/Retrospective/pdg2024_review.md) — PDG 2024 Review of Particle Physics; issuing body is the PDG collaboration.
- [`fan2023_electron_g_factor`](../../../../../../Roadmapping/History/Bibliography/Retrospective/fan2023_electron_g_factor.md) — Fan et al. 2023 Penning-trap measurement of the electron magnetic moment.
- [`fermilab_e989_2023`](../../../../../../Roadmapping/History/Bibliography/Retrospective/fermilab_e989_2023.md) — Fermilab E989 Muon g-2 collaboration 2023 result combined with BNL E821.
- [`bmw2021_lattice_qcd`](../../../../../../Roadmapping/History/Bibliography/Retrospective/bmw2021_lattice_qcd.md) — Borsanyi-Marquart-Wirth (BMW) collaboration 2021 lattice-QCD calculation of the hadronic vacuum polarisation contribution to a_μ.
- [`theory_initiative_2020`](../../../../../../Roadmapping/History/Bibliography/Retrospective/theory_initiative_2020.md) — Aoyama-Asmussen-Benayoun-… 2020 Theory Initiative consensus paper on the SM prediction for a_μ.

`codata2018_constants` (cited for the electron and proton masses) was scaffolded under issue #99.

## Refresh log

See [`refresh_log.md`](refresh_log.md).
