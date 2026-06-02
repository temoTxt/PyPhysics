# `flag.*` namespace — refresh policy addendum

Extends [`mcp_servers/precision_data_mcp/REFRESH_POLICY.md`](../../../REFRESH_POLICY.md) for FLAG lattice-QCD reference values.

Created under issue [#98](https://github.com/temoTxt/PyPhysics/issues/98). Sub-issue of umbrella [#92](https://github.com/temoTxt/PyPhysics/issues/92).

## Source

FLAG publishes per-edition reports tabulating lattice-QCD averages with systematic-uncertainty breakdowns. The current seed uses FLAG 2024 (`flag_edition: "2024"`, *Eur. Phys. J. C* 84).

## Initial seed coverage

Six quantities from the FLAG 2024 N_f=2+1+1 averages: `m_pi_charged`, `f_pi`, `f_K`, `f_K_over_f_pi`, `g_A`, `sigma_piN`. Adding more is mechanical follow-up via `data.json` + `refresh_log.md`.

## Edition tracking

`$source_revision.flag_edition` records which FLAG edition the values came from. When FLAG 2026 (or later) ships, bumping is a deliberate refresh per umbrella §6:

1. Hand-update `data.json` against the new edition.
2. Append diff to [`refresh_log.md`](refresh_log.md).
3. Existing cached `cache_key` records consumed by prior verification documents remain valid under the old edition.

## Bibliography stub

- `flag2024_averages` (scaffolded as part of this PR; non-paper-source cite-key convention from umbrella §9).

## Refresh log

See [`refresh_log.md`](refresh_log.md).
