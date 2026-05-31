# `qed.*` namespace — refresh log

Chronological record of refreshes, edition bumps, and precision updates for the `qed.*` namespace. Per the umbrella [`REFRESH_POLICY.md`](../../../REFRESH_POLICY.md) §7.

## 2026-05-30 — initial seed (issue #97)

- **Cause:** initial implementation of the `qed.*` namespace per umbrella issue [#92](https://github.com/temoTxt/PyPhysics/issues/92) sub-issue [#97](https://github.com/temoTxt/PyPhysics/issues/97).
- **Species seeded:** H, He II, Li III, Si XIV, muonic_H, muonium, positronium.
- **Observables seeded** (one or more per species):
  - H: 2S1/2-2P1/2 Lamb shift (1057.845(9) MHz, CODATA-2018); 1s2S1/2 hyperfine (1420.405751768(2) MHz, Essen 1957); 1S-2S transition (2466061413187035(10) Hz, Parthey-Hänsch 2011).
  - He II: 1s2S1/2 bound-electron g-factor (1.9999857940(10), Kohler 2016).
  - Li III: 1s2S1/2 hyperfine (11890.018(2) MHz, Beckert 2007).
  - Si XIV: 1s2S1/2 bound-electron g-factor (1.9953395931(6), Sturm 2013) — load-bearing for issue #82.
  - muonic_H: 2S1/2-2P3/2 Lamb shift (49881.88(76) GHz, Pohl 2010) — the proton-radius-puzzle trigger.
  - muonium: 1s2S1/2 hyperfine (4463.302776(51) MHz, Liu 1999).
  - positronium: 1s ortho-para (203389(2) MHz, Ishida 2014).
- **Bib stubs scaffolded** (8 new primary-source stubs):
  - `essen1957_h_hyperfine`
  - `parthey2011_h_1s2s`
  - `kohler2016_he_g_factor`
  - `beckert2007_li2_hyperfine`
  - `sturm2013_si13_g_factor`
  - `pohl2010_muonic_hydrogen`
  - `liu1999_muonium_hyperfine`
  - `ishida2014_positronium`
- **Demonstrator refactor** (sub-issue #97 acceptance criterion #4): [`Bethe_Salpeter/05_LambShift.md`](../../../../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/05_LambShift.md) updated to consume `qed.get_lamb_shift("H", "2S1/2-2P1/2")` rather than the hand-typed `1057.845(9) MHz` value, validating the end-to-end qed.* lookup path.
- **Disagreement-representation cases:** none in this initial seed (each value comes from a single primary source). When a second independent measurement of any seeded quantity lands at a competing facility, it will be added as a `{"values": [...]}` per umbrella §5.

Subsequent entries will track:
- Precision updates to any of the seeded values (e.g., when the next-generation H 1S-2S measurement at MPQ lands a more precise value).
- Additional species added by follow-up PRs (Be IV, B V, C VI, ... for the full hydrogenic Z-scan; multi-electron species beyond hydrogenic).
- Additional observables added by follow-up PRs (more Lamb-shift transitions, hyperfine splittings of higher-n states, etc.).
- Deprecation of `nist.list_dirac_targets` once the qed.* path is fully validated as the canonical access for bound-state QED precision data.
