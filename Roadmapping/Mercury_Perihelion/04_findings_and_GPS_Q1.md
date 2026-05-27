# 04 — Findings + cross-reference to GPS author report Q1

**Status:** ✅ flagged.
**Source paper:** [`../Tepper_Gill_Papers/Dual Newtonian Theory.tex`](../Tepper_Gill_Papers/Dual%20Newtonian%20Theory.tex).
**Cross-reference target:** [`../Author_Reports/2026-05_gps_relativity_summary_for_gill.md` §5 Q1](../Author_Reports/2026-05_gps_relativity_summary_for_gill.md).

## 1. What the GPS author report's Q1 asked

The GPS proper-time companion campaign ([PR #73](https://github.com/temoTxt/PyPhysics/pull/73), merged at `ac19d44`) flagged four open questions about how the Gill–Zachary framework should extend to curved spacetime. The first (and most load-bearing):

> **Q1.** What is the correct curved-spacetime extension of the framework? The minimal-extension hypothesis used in pt_01, pt_02, pt_04, pt_07, and pt_08 is `c → b` in the Schwarzschild line element, with `b → c` at low velocity. This is *one* possible extension; we considered three others without resolution. ... For GPS at the `10⁻¹⁰` precision floor, all four hypotheses reduce to the same operational answer; the question matters for next-generation optical-clock GPS (`10⁻¹⁷` stability is at the edge of distinguishing them) and for strong-field regimes where `GM/(rc²)` is no longer `~10⁻¹⁰`.

The GPS campaign noted that **at GPS-relevant precision**, the framework's structural choice did not matter (`b ≈ c` to `10⁻¹⁰`). The question of *which* extension is correct could not be tested by GPS observations.

## 2. What this paper proposes

The *Dual Newtonian Theory* paper (Gill, undated; user-supplied) proposes a **different** curved-spacetime extension than the speculative `c → b` Schwarzschild extension used in the GPS campaign. The paper's extension:

- Keeps the framework's flat-space Hamiltonian `K = π²/(2m) + mc² + V + V²/(2mc²)` (verified ✅ at [`../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` line 41](../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) via the (II.3) "potential-in-the-mass" kernel).
- Plugs in `V = −GMm/r` for gravitational interaction (Newtonian gravity as the potential).
- Derives the modified-Newtonian acceleration `a = −(GM/r²)(1 − GM/(c²r))·ê_r` (Eq. h4) **directly from the `V²/(2mc²)` term in `K`**, without introducing a metric or geodesic equation.

This is a **5th candidate extension** beyond the four considered in the GPS author report Q1 (`c → b` in full metric; `c → b` in `g_{00}` only; null vs timelike geodesics treated differently; `b`-dependent connection). It is structurally simpler — no curved-spacetime apparatus at all; just the framework's existing `V²/(2mc²)` correction applied to gravity.

## 3. What Mercury's perihelion advance says about this extension

[Document 03](03_numerical_predictions.md) numerically computes the framework's prediction for Mercury under the paper's extension:

| Quantity | Value | Interpretation |
|---|---|---|
| Framework full-dual prediction | `+37.79″/century` | the prediction of the paper's extension |
| Observed | `+43″/century` | Le Verrier residual |
| Standard GR | `+42.99″/century` | reproduces observed to `<0.1%` |
| **Framework relativistic part** (separated from reduced-mass) | **`−6.87″/century`** | *opposite sign from GR's `+42.99″`* |

**Headline finding for Q1:** This 5th candidate extension **does not** reproduce observed Mercury perihelion advance. The framework's modified-gravity term has the wrong sign of relativistic correction (the `(1 − GM/(c²r))` factor reduces gravity, producing perihelion *recession*, while GR's effective `1/r³` term enhances gravity and produces perihelion *advance*).

## 4. What this resolves and what it does not

**Resolved (negative):** the `V²/(2mc²)` extension to gravity proposed in this paper is *not* a successful curved-spacetime extension of the framework at the precision Mercury can measure. Whatever the correct framework GR extension is, it is *not* this one.

**Not resolved:** the four candidate extensions from the GPS author report Q1 remain open. None of them is tested by this paper.

**Subtle point.** Both Q1 and this finding are *conditional* on the framework being internally consistent — i.e., on the assumption that the same `K = π²/(2m) + mc² + V + V²/(2mc²)` Hamiltonian governs both EM (where the GPS PT companion verified equivalence to standard at GPS precision) and gravity (where this paper's extension fails to reproduce Mercury). If the framework intends a *different* structural form for gravitational interaction than the simple `V = −GMm/r` substitution into `K`, that's exactly Q1 — and the right next step is to specify which form is intended.

## 5. Implications for the framework's broader claims

The paper's section §1 (Dual Classical Theory) and §2 introduction explicitly position this work as a "Newtonian-Maxwell relativistic unification". The Mercury computation is the paper's only quantitative test of the gravitational side; the EM side rests on the previously verified Maxwell-paper substitution rules.

If the gravitational extension does not match Mercury at the framework's claimed precision, the unification claim is structurally incomplete: the framework reproduces EM phenomena (as the verified Maxwell paper showed, and as the GPS PT companion confirmed at GPS precision) but does not reproduce gravitational phenomena at the precision Mercury demands.

This is **substantive AI under Crocco rule #5** — flagging it requires the same per-paragraph human-acceptance discipline as other substantive AI findings. The honest reading below is the campaign's framing; whether to convert it into a finding for [`FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md) is the author's call.

<!-- TODO: human reviews and fills in — confirms (a) the "5th candidate extension" framing is honest (the paper proposes something genuinely different from the four GPS-campaign hypotheses, not a refinement of them); (b) the negative finding ("extension does not reproduce Mercury perihelion") is supported by the numerical decomposition in document 03 and is not a misinterpretation of the paper's intent; (c) the GPS author-report cross-reference paragraph (to be added to §5 of that report) accurately conveys this without overclaiming or underclaiming. -->

## 6. Pointer for GPS author report §5 update

The GPS author report at [`../Author_Reports/2026-05_gps_relativity_summary_for_gill.md`](../Author_Reports/2026-05_gps_relativity_summary_for_gill.md) §5 Q1 should be updated with a closing paragraph noting:

- **A 5th candidate extension exists**, proposed in Gill's *Dual Newtonian Theory* paper (the .tex source attached to issue #81): the `V²/(2mc²)` term in the (II.3) "potential-in-the-mass" kernel applied to gravity via `V = −GMm/r`, producing Eq. (h4) `a = −(GM/r²)(1 − GM/(c²r))·ê_r`.
- **The Mercury perihelion test rules this extension out** at the `~12%` level — the framework's full-dual prediction under-predicts observed by `5.2″/century`, with the framework's relativistic-correction part having the opposite sign from GR's.
- The four other extensions in Q1 remain open; this finding *narrows the candidate space* by one.

The PDF should be rebuilt via `Roadmapping/Author_Reports/build_report.sh` after the §5 update; the page-count check `[3, 7]` should remain satisfied.

## 7. Verdict

⚠ The paper's proposed gravitational extension does not reproduce observed Mercury perihelion advance at the framework's claimed precision. This is a **substantive finding** for the GPS author report's Q1: the 5th candidate extension is ruled out by Mercury observation. Flagged for cross-reference into the GPS author report and (at the author's discretion) into [`../Equation_Verification/FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md) as a potential Finding 4.
