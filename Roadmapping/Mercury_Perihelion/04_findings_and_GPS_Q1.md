# 04 — Findings + cross-reference to GPS author report Q1

**Status:** ✅ flagged.
**Source paper:** [`../Tepper_Gill_Papers/Dual Newtonian Theory.tex`](../Tepper_Gill_Papers/Dual%20Newtonian%20Theory.tex).
**Cross-reference target:** [`../Author_Reports/2026-05_gps_relativity_summary_for_gill.md` §5 Q1](../Author_Reports/2026-05_gps_relativity_summary_for_gill.md).

## 1. What the GPS author report's Q1 asked

The GPS proper-time companion campaign ([PR #73](https://github.com/temoTxt/PyPhysics/pull/73), merged at `ac19d44`) flagged four open questions about how the Gill–Zachary framework should extend to curved spacetime. The first (and most load-bearing):

> **Q1.** What is the correct curved-spacetime extension of the framework? The minimal-extension hypothesis used in pt_01, pt_02, pt_04, pt_07, and pt_08 is `c → b` in the Schwarzschild line element, with `b → c` at low velocity. This is *one* possible extension; we considered three others without resolution. ... For GPS at the `10⁻¹⁰` precision floor, all four hypotheses reduce to the same operational answer; the question matters for next-generation optical-clock GPS (`10⁻¹⁷` stability is at the edge of distinguishing them) and for strong-field regimes where `GM/(rc²)` is no longer `~10⁻¹⁰`.

The GPS campaign noted that **at GPS-relevant precision**, the framework's structural choice did not matter (`b ≈ c` to `10⁻¹⁰`). The question of *which* extension is correct could not be tested by GPS observations.

## 2. What this paper proposes

The *Dual Newtonian Theory* paper (Gill, undated; user-supplied) proposes a **categorically different** kind of answer to Q1 than the four GPS-Q1 candidates. Where those four are all *metric* extensions of Schwarzschild (different ways to put `b` into a curved-spacetime line element), the Mercury paper's approach uses *no* curved-spacetime apparatus at all:

- Keeps the framework's flat-space Hamiltonian `K = π²/(2m) + mc² + V + V²/(2mc²)` (verified ✅ at [`../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` line 41](../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) via the (II.3) "potential-in-the-mass" kernel — *in its EM context*).
- Plugs in `V = −GMm/r` for gravitational interaction (Newtonian gravity as the potential — *the interpretive move* this paper makes; see §0 caveat in [`01_setup_and_verification.md` §3](01_setup_and_verification.md)).
- Derives the modified-Newtonian acceleration `a = −(GM/r²)(1 − GM/(c²r))·ê_r` (Eq. h4) **directly from the `V²/(2mc²)` term in `K`**, without introducing a metric or geodesic equation.

This is a **different *kind* of answer** to Q1: not a refinement of the four metric extensions but a structural alternative ("the framework does not need a curved-spacetime extension at all; the existing flat-space kernel suffices for gravity, under the substitution `V → -GMm/r`"). Putting it on a numbered list as a "5th candidate" alongside the four metric extensions would conflate categorically different approaches; the right reading is that the Mercury paper proposes a *separate path* whose Mercury test is informative independently of the four metric candidates.

## 3. What Mercury's perihelion advance says about this extension

[Document 03](03_numerical_predictions.md) numerically computes the framework's prediction for Mercury under the paper's extension:

| Quantity | Value | Interpretation |
|---|---|---|
| Framework full-dual prediction `Δφ_d` | `+37.79″/century` (positive forward advance) | the prediction of the paper's extension |
| Observed | `≈ +43″/century` | modern residual (Newcomb–Clemence) — *N-body-derived*, see §4 caveat |
| Standard GR | `+42.99″/century` | reproduces observed to `<0.1%` |
| **Framework's `(1 − GM/(c²r))`-induced contribution** (separated from reduced-mass) | **`−6.87″/century`** | *opposite sign from GR's `+42.99″`* — but this is the *contribution*, not the *total* |

**Headline finding for Q1.** This extension **does not reach precision-level agreement** with the observed Mercury perihelion advance. The framework's total prediction is a positive forward advance of `+37.79″/century` (the framework does *not* predict perihelion recession), but it under-predicts the observed `≈43″` by `~12%`. The framework's `(1 − GM/(c²r))`-induced contribution has the opposite sign from GR's relativistic contribution — the modified gravity is *weaker* than Newtonian, *reducing* the reduced-mass-driven advance rather than *adding* to it as GR's effective `1/r³` term does. Whether the `~12%` disagreement constitutes "ruling the extension out" depends on the precision the framework claims; the paper's own text says only that the prediction "well approximates" observation (claiming few-percent agreement, not precision agreement), and a 12% disagreement is consistent with that level of claim without strictly contradicting it.

### 3a. The two-body vs N-body caveat — now resolved (see doc 05)

The 43″ observed residual is computed by subtracting *Newtonian* planetary perturbations from the total observed precession (~5600″/century). The original concern was that the framework's modified gravity would also change the planetary-perturbation calculation, so the two-body comparison might be unfair. [`05_nbody_orbital_mechanics.md`](05_nbody_orbital_mechanics.md) resolves this concern and *sharpens* the finding rather than softening it:

1. **The proper eccentric-orbit perihelion advance is `−7.17″/century` = exactly `−1/6` of GR** (opposite sign, one-sixth magnitude). The paper's `+37.79″` came from a circular-orbit heuristic that conflated orbital-frequency change with ellipse precession; the genuine apsidal-angle calculation gives the `−1/6` structural result.
2. **The N-body planetary corrections cannot close the gap.** The framework's relativistic factor `GM/(c²r)` scales with the *attracting* mass; Jupiter's contribution is `~7×10⁻⁵` of the Sun's, and all planets summed give `~10⁻³ arcsec/century` — three orders of magnitude too small to bridge `−7 → +43`.

So the N-body refinement does *not* rescue the framework's gravitational extension. The discrepancy is structural (sign + factor of 6), not a numerical near-miss. This is the load-bearing finding; the "12% gap" framing below (from the paper's circular-orbit value) is superseded.

## 4. What this resolves and what it does not

**Indicative (not strict rule-out):** the `V²/(2mc²)` extension to gravity proposed in this paper does not reach precision-level agreement with the observed Mercury perihelion advance at the two-body level. Whether this constitutes "ruling the extension out" depends on (a) the precision the framework intends to claim (the paper's text suggests few-percent "well approximates", which 12% is at the edge of), and (b) whether the framework's modified gravity applied N-body-consistently would shift the prediction. Both questions are open.

**Not resolved:** the four candidate extensions from the GPS author report Q1 remain open. None of them is tested by this paper.

**Three load-bearing interpretive questions, all flagged for author review.**

1. **Is the `V²/(2mc²)` kernel intended to apply to gravitational `V`?** The DRQM I (II.3) kernel where this term was verified used `V` as an EM potential. The Mercury paper's substitution `V → -GMm/r` is the author's interpretive move, not a derived consequence. If the author intends the kernel as EM-only, the Mercury application is a separate proposal that this campaign tests but does not endorse.
2. **What precision does the framework claim for the Mercury prediction?** The paper says "well approximates" — at the few-percent level. A 12% disagreement is at the edge of that and depends on which prediction (Corda, full dual, approximate dual) the author considers the framework's intended answer.
3. **N-body consistency.** If the framework's modified Newtonian gravity is applied consistently to all planetary perturbations on Mercury, the residual changes and so does the comparison. A consistent test would settle whether the 12% gap is real or an artifact of the two-body-vs-N-body mismatch.

## 5. Implications for the framework's broader claims

The paper's section §1 (Dual Classical Theory) and §2 introduction explicitly position this work as a "Newtonian-Maxwell relativistic unification". The Mercury computation is the paper's only quantitative test of the gravitational side; the EM side rests on the previously verified Maxwell-paper substitution rules.

If the gravitational extension does not match Mercury at the framework's claimed precision, the unification claim is structurally incomplete: the framework reproduces EM phenomena (as the verified Maxwell paper showed, and as the GPS PT companion confirmed at GPS precision) but does not reproduce gravitational phenomena at the precision Mercury demands.

This is **substantive AI under Crocco rule #5** — flagging it requires the same per-paragraph human-acceptance discipline as other substantive AI findings. The honest reading below is the campaign's framing; whether to convert it into a finding for [`FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md) is the author's call.

<!-- TODO: human reviews and fills in — confirms (a) the "5th candidate extension" framing is honest (the paper proposes something genuinely different from the four GPS-campaign hypotheses, not a refinement of them); (b) the negative finding ("extension does not reproduce Mercury perihelion") is supported by the numerical decomposition in document 03 and is not a misinterpretation of the paper's intent; (c) the GPS author-report cross-reference paragraph (to be added to §5 of that report) accurately conveys this without overclaiming or underclaiming. -->

## 6. Pointer for GPS author report §5 update

The GPS author report at [`../Author_Reports/2026-05_gps_relativity_summary_for_gill.md`](../Author_Reports/2026-05_gps_relativity_summary_for_gill.md) §5 Q1 should be updated with a closing paragraph noting:

- **A categorically different kind of answer to Q1 exists**, proposed in Gill's *Dual Newtonian Theory* paper (the .tex source attached to issue #81): the `V²/(2mc²)` term in the (II.3) "potential-in-the-mass" kernel applied to gravity via `V = −GMm/r`, producing Eq. (h4) `a = −(GM/r²)(1 − GM/(c²r))·ê_r`. Unlike the four Q1 candidates (all metric extensions of Schwarzschild), this is a flat-space modified-Newtonian approach with no curved-spacetime apparatus.
- **The Mercury perihelion test does not reach precision-level agreement** with the observed value: framework full-dual `+37.79″/century` vs observed `≈+43″/century` (12% gap), with the framework's `(1 − GM/(c²r))`-induced contribution opposite-sign from GR's relativistic contribution. The framework's *total* prediction is still positive forward advance, just smaller than GR's. Whether this rules the extension out depends on the precision the framework intends to claim and on N-body consistency (see §3a and §4 above).
- The four metric extensions in Q1 remain open and untested.

The PDF should be rebuilt via `Roadmapping/Author_Reports/build_report.sh` after the §5 update; the page-count check `[3, 7]` should remain satisfied.

## 7. Verdict

❌ The paper's proposed gravitational extension does not reproduce observed Mercury perihelion advance, and the proper orbital-mechanics calculation ([doc 05](05_nbody_orbital_mechanics.md)) shows *why* in a way the original 12% framing did not: the framework's force law gives a perihelion advance of **exactly `−1/6` of GR** (`−7.17″/century`: opposite sign, one-sixth magnitude), and **no N-body refinement can close this** (planetary corrections are `~10⁻³`–`10⁻⁴` too small). This is the classic "1/6 factor" of a force-law-only modification that lacks GR's spatial-metric curvature — connecting to the Van Flandern reference in issue [#81](https://github.com/temoTxt/PyPhysics/issues/81).

For the GPS author report's Q1: the flat-space modified-Newtonian extension is **ruled out for Mercury by a structural argument** (the `−1/6` factor), not merely by a numerical near-miss. The interpretive questions in §4 (is the kernel meant for gravitational `V`? what precision does the framework claim?) remain relevant, but even granting the most charitable reading, the `−1/6` structural factor cannot be reconciled with the observed `+43″` by any orbital-dynamics refinement. Flagged for cross-reference into the GPS author report and (at the author's discretion) into [`../Equation_Verification/FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md) as a potential Finding 4.
