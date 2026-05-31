# `nuclear.*` namespace — refresh policy addendum

Extends [`mcp_servers/precision_data_mcp/REFRESH_POLICY.md`](../../../REFRESH_POLICY.md) for nuclear / particle structure observables.

Created under issue [#98](https://github.com/temoTxt/PyPhysics/issues/98). Sub-issue of umbrella [#92](https://github.com/temoTxt/PyPhysics/issues/92).

## Source — hand-curated JSON

Like `qed.*`, this namespace's data is scattered across PRL / Nature / Science / PDG / CODATA primary sources with no canonical machine-readable database. The seed at [`data.json`](data.json) is hand-curated; each entry's `source` is a `cite_key`.

## Initial seed coverage

| Nucleus | Observables seeded | Primary sources |
|---|---|---|
| proton | charge radius (5-value puzzle); magnetic moment | `pohl2010_muonic_hydrogen`, `bernauer2010_a1_mainz`, `fleurbaey2018_1s_3s_spectroscopy`, `xiong2019_prad`, `codata2018_constants` |
| neutron | mean-square charge radius `<r^2>_n`; magnetic moment | `pdg2024_review`, `codata2018_constants` |
| deuteron | charge radius (2-value tension); magnetic moment | `pohl2016_muonic_deuterium`, `codata2018_constants` |

## Disagreement-representation cases

**This namespace exists primarily to surface the proton-radius puzzle as the canonical disagreement demonstrator** per the umbrella's §Resolved-decisions #5.

Calling `nuclear.get_charge_radius("proton")` (no `method` argument) returns ALL five values with their `method` labels:

| Method | Value (fm) | Source |
|---|---|---|
| muonic_hydrogen_spectroscopy | 0.84087(39) | `pohl2010_muonic_hydrogen` |
| electron_scattering | 0.879(8) | `bernauer2010_a1_mainz` |
| CODATA_2018_pre_muonic_average | 0.877(7) | `codata2018_constants` |
| 1S-3S_spectroscopy | 0.833(10) | `fleurbaey2018_1s_3s_spectroscopy` |
| PRad_electron_scattering | 0.831(14) | `xiong2019_prad` |

The first three cluster around `~0.88 fm`; the last two and the muonic-H value cluster around `~0.84 fm`. The disagreement is at the ~5σ level depending on which subset is averaged. Calling with `method="muonic_hydrogen_spectroscopy"` (or any explicit method label) returns the single matching value.

The deuteron `nuclear.get_charge_radius("deuteron")` similarly returns 2 values reflecting an analogous ~3σ tension.

## Programmatic refresh

The current seed is hand-curated. Adding a sixth proton-radius measurement (or refreshing a value when a higher-precision result lands):

1. Scaffold a bib stub for the new primary-source paper.
2. Append a new entry to `data.json::nuclei.p.charge_radius.values` with the new measurement's value, uncertainty, method label, and cite_key.
3. Append a chronological entry to [`refresh_log.md`](refresh_log.md) documenting the addition.
4. Do NOT remove the prior values — the disagreement-representation rule preserves the history.

## Bibliography stubs scaffolded

- `pohl2010_muonic_hydrogen` — scaffolded earlier under issue #97 (qed.* namespace uses the same primary source)
- `bernauer2010_a1_mainz` (new this PR)
- `fleurbaey2018_1s_3s_spectroscopy` (new this PR)
- `xiong2019_prad` (new this PR)
- `pohl2016_muonic_deuterium` (new this PR)

## Refresh log

See [`refresh_log.md`](refresh_log.md).
