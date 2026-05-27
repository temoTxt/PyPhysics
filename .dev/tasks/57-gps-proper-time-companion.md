# Implementation plan: GPS relativity — proper-time companion campaign

**Tracks:** issue [#57](https://github.com/temoTxt/PyPhysics/issues/57).
**Status:** plan; per-effect documents to follow in this same PR.
**Sibling:** [#71](https://github.com/temoTxt/PyPhysics/pull/71) `57-gps-standard-derivation` — the standard-method derivation. **This PR depends on #71 merging first** for cross-reference links to resolve; the two PRs are conflict-free at the file level (only `README.md` is touched by #71; this PR adds `README_proper_time.md` and `pt_*.md` files).

This is the **proper-time companion** to the standard-method derivation. The campaign applies the Gill–Zachary proper-time substitution rules ([[two_mathematically_equivalent_versions_of_maxwells_equations]], verified at [`Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md)) to each of the nine effects in the standard derivation. For each effect, the campaign records (a) whether the framework applies as-published, (b) what it predicts, (c) where it would require an extension beyond the verified corpus.

---

## 1. Goal and scope

For each of the nine standard-method effects, work the proper-time formulation explicitly, and produce one of three verdicts:
- **✅ matches standard at GPS precision** — framework applies and gives the same observable; mathematical equivalence at the relevant order.
- **⚠ structural difference (framework extension)** — framework applies but gives an extra term that would be observable in a different regime than GPS.
- **❌ framework doesn't extend natively** — the effect is in curved spacetime (Schwarzschild), and the Gill–Zachary substitution rules are SR + flat-space Maxwell only. The doc carries a *speculative* extension (e.g., `c → b` substitution into the Schwarzschild metric) clearly flagged as substantive AI under [Crocco §1](../../Roadmapping/Tooling/CROCCO_COMPLIANCE.md), and reduces to the standard result at GPS-relevant precision.

### Why this earns its own thread

- **Mathematical equivalence at GPS precision is a real result.** The proper-time framework's claim is that, where it applies, it reproduces standard SR Maxwell at the operational level. GPS is the most stringently tested applied-relativity result available; *demonstrating* the reproduction at the per-effect level is a load-bearing consistency check.
- **The boundary between SR and GR is where the framework needs work.** Five of the nine effects are gravitational (GR-only); two are SR (kinematic and EM-propagation in flat space); two are mixed. Identifying clean boundaries is itself a research outcome.
- **Side-by-side comparability.** The sibling PR #71 produced the standard derivation; this PR mirrors the structure 1:1 so a reader can compare both formulations of each effect.

### Explicit non-goals

- Not a derivation of the proper-time substitution rules — those are verified in [`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) and cited here.
- Not a full extension of the framework to curved spacetime — the speculative GR-extension paragraphs in pt_02, pt_04, pt_07, pt_08 are *minimal hypotheses* (`c → b` in the Schwarzschild metric, with `b → c` at GPS-low-velocity), carried strictly to show what the framework would predict under that extension. Strong substantive-AI disclosure required.
- Not a critique of the framework — agreement with standard at GPS precision is the *expected* outcome; disagreement would be a flag for the framework, not a flag for GPS.

---

## 2. Repository structure

Same campaign folder as the sibling PR #71, with `pt_` prefixed files to avoid file-level collisions:

```
Roadmapping/GPS_Relativity/
├── README.md                                  # standard-method README (PR #71)
├── _template_effect.md                        # standard-method template (PR #71)
├── 01_..._09_*.md                             # standard-method per-effect (PR #71)
│
├── README_proper_time.md                      # ⬅ this PR: companion README
├── _proper_time_GPS_cheatsheet.md             # ⬅ this PR: applicability matrix
├── _template_pt_effect.md                     # ⬅ this PR: PT per-effect template
└── pt_01_..._pt_09_*.md                       # ⬅ this PR: 9 PT per-effect docs
```

Once both PRs merge, a future doc-housekeeping PR can fold `README_proper_time.md` into `README.md` as a second section. Keeping them separate for now avoids merge conflict.

---

## 3. Per-effect document template

The PT per-effect template at [`_template_pt_effect.md`](../../Roadmapping/GPS_Relativity/_template_pt_effect.md) mirrors the standard `_template_effect.md` with three additions:

1. **§0. Framework applicability** — explicit statement of which substitution rules apply, and whether the effect is GR-only.
2. **§5b. GPS-precision-limit equivalence** — explicit reduction `b → c`, `u → v` showing that the PT prediction agrees with the standard prediction at the operational order.
3. **§6b. Where the framework would diverge** — a regime where the framework's prediction is operationally distinguishable from standard (typically: high acceleration, short timescale, or non-circular trajectory).

The Crocco compliance disclosure for this campaign:
- pt_03, pt_05, pt_06 are **pragmatic** — framework applies as-published, derivation is mechanical.
- pt_01, pt_02, pt_04, pt_07, pt_08 are **substantive** — they carry a speculative GR-extension that is not in the verified corpus. Each such doc carries `<!-- TODO: human reviews and fills in -->` on the GR-extension paragraph.

---

## 4. Decision points (deferred to author)

- **D1.** Should pt_02, pt_04, pt_07, pt_08 attempt the `c → b` Schwarzschild extension at all, or just record "framework doesn't extend; standard result stands"? Default: minimal extension with strong substantive-AI disclosure. The reader can then judge whether the extension is the right shape.
- **D2.** When the proper-time and standard formulations make algebraically different expressions that evaluate to the same number at GPS precision, do we record the algebraic difference or only the numerical agreement? Default: algebraic difference, since the algebraic structure is the framework's claim.

---

## 5. Honest scoping

This is a *companion* to the standard derivation, not a critique of either standard relativity or the Gill–Zachary framework. The honest expected outcome at GPS-relevant precision is:
- 4–5 effects: PT framework reproduces standard exactly (mathematical equivalence at this order)
- 4–5 effects: PT framework as-published doesn't extend; speculative extension reduces to standard at GPS precision

A *negative* finding (PT predicts something different at GPS precision that disagrees with measurement) would be a serious problem for the framework. The campaign is structured to be sensitive to such a finding, not to be biased against it.

The author's prior involvement in the framework (co-author of *Dual Relativistic Quantum Mechanics I*, 2021) is acknowledged in [§13 of the Electromagnetism plan](42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix) — the same honest-framing discipline carries here.
