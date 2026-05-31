# `nuclear.*` namespace — refresh log

## 2026-05-30 — initial seed (issue #98)

- **Cause:** initial implementation of the `nuclear.*` namespace per umbrella issue [#92](https://github.com/temoTxt/PyPhysics/issues/92) sub-issue [#98](https://github.com/temoTxt/PyPhysics/issues/98).
- **Nuclei seeded:** proton (p), neutron (n), deuteron (d).
- **Proton-radius puzzle:** 5 values seeded with method labels (muonic_hydrogen_spectroscopy 0.84087(39) fm; electron_scattering 0.879(8) fm; CODATA_2018_pre_muonic_average 0.877(7) fm; 1S-3S_spectroscopy 0.833(10) fm; PRad_electron_scattering 0.831(14) fm). All returned by default when method=None per umbrella §5.
- **Deuteron analogue:** 2 values (muonic_deuterium 2.12562(78) fm; CODATA_2018_deuterium_spectroscopy 2.12771(22) fm).
- **Magnetic moments:** proton 2.79284734463(82) μ_N; neutron -1.91304273(45) μ_N; deuteron 0.857438230(24) μ_N. All cite `codata2018_constants`.
- **Bib stubs scaffolded:** `bernauer2010_a1_mainz`, `fleurbaey2018_1s_3s_spectroscopy`, `xiong2019_prad`, `pohl2016_muonic_deuterium`. `pohl2010_muonic_hydrogen` was already scaffolded under issue #97 (qed.* shares the source).
- **Demonstrator status:** the proton-radius puzzle is itself the demonstrator — the disagreement-representation rule is operationally proven by the multi-value return.

Subsequent entries will track:
- Additional nuclei (helion, triton, alpha; eventually heavier nuclei for the campaign's hyperfine-structure work).
- Sixth+ proton-radius measurements as they land (the puzzle is still active as of 2024-2025).
- Any precision-update events touching seeded values.
