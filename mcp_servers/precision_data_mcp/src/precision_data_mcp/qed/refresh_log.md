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

## 2026-05-31 — hydrogen precision-spectroscopy seed extension (issue #108)

- **Cause:** User-driven gap: post safety-flag enforcement (PR #107), hydrogen experimental coverage in qed.* was only 3 quantities (1S hyperfine, 1S-2S, the codata_adjusted Lamb shift). Issue [#108](https://github.com/temoTxt/PyPhysics/issues/108) added the canonical Hänsch / Pohl / Beyer / Fleurbaey / Lundeen suite of directly-measured H transitions.
- **Schema impact:** `species.H.lamb_shifts.2S1/2-2P1/2` becomes multi-value. `qed.get_lamb_shift()` now takes an optional `method` arg and returns either the full `{"values": [...]}` list (when method=None) or the single matching record (when method is specified). Mirrors the `nuclear.get_charge_radius()` pattern.
- **New entries seeded (5 directly-measured H transitions):**
  - `species.H.lamb_shifts.2S1/2-2P1/2` (multi-value): Lundeen-Pipkin 1981 direct RF measurement added alongside the existing CODATA entry. Numerical values match (by historical construction); the safety flag differs.
  - `species.H.transitions.1S-2S_H-D_isotope_shift`: Parthey 2010 H/D isotope shift, 670994334.606(15) kHz.
  - `species.H.transitions.1S-3S`: Fleurbaey 2018 two-photon spectroscopy, 2922743278671.5(2.6) kHz.
  - `species.H.transitions.2S-4P1/2`: Beyer 2017 MPQ measurement, 616520152720.3(2.3) kHz.
  - `species.H.transitions.2S-8S1/2`: de Beauvoir 1997 LKB measurement, 770649350012(9) kHz.
- **Bib stubs scaffolded (4 new):**
  - `lundeen1981_2s2p_microwave`
  - `parthey2010_h_d_isotope_shift`
  - `beyer2017_2s4p`
  - `debeauvoir1997_2s_8s_8d`
  (`fleurbaey2018_1s_3s_spectroscopy` already existed from issue #98.)
- **Numerical-value caveat:** the specific digits seeded above use canonical values quoted from each primary source as cited in the existing precision-spectroscopy literature, but should be **cross-checked against the published paper** by a human before flipping each bib stub's `human_reviewed: true`. The `_note` field on each entry records this.
- **Tests added:** 7 new tests in `test_lookup.py::TestHydrogenPrecisionSpectroscopyExtension` + 3 in `TestHydrogen` for the multi-value handling. Cross-namespace `test_safety_contract.py::TestQedSafetyContract::test_h_lamb_shift_multivalue_each_carries_contract` updated for the schema change. All 93 pass + 1 expected network-skip.
- **Out of scope (recorded for follow-up):** 2S-4P3/2 companion (Beyer 2017 measured both fine-structure components); 2S-6S/6D (Bourzeix 1996 series); 2S-8D (de Beauvoir 1997 companion); 2S-12D (Yost et al.); derived experimentally-anchored level energies via fit. Each is a small follow-up extension.
