# Per-problem document template — Quantum Mechanics / Griffiths

This template defines the structure of every per-problem document under [`Griffiths/Ch??_*.md`](Griffiths/). Mirrors the Electromagnetism template ([`Roadmapping/Electromagnetism/_template_problem.md`](../Electromagnetism/_template_problem.md)) with two key differences: **two source citations** (Griffiths 2e + 3e) instead of two unit systems, and the proper-time section uses the **canonical Hamiltonian** `K = H²/(2mc²) + mc²/2` rather than the field equations of [#42](https://github.com/temoTxt/PyPhysics/issues/42).

## Naming convention

Problem headings use `G{N}e-P{C}.{P}`, where `{N}e` is the edition (`2e` or `3e`), `{C}` is the chapter (in the 3e numbering — note Ch. 11/12 renumbered between editions), and `{P}` the 3e problem number. The 2e cross-reference is recorded in the source line. When a problem appears in both editions with matching numbering, the 3e key is canonical and the 2e equivalence is noted.

## Per-problem template (the body of every problem entry)

````markdown
### Problem G3e-P1.4 — short title

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* {one-line reason this problem is in scope for this PR}.
- *Alternatives considered:* {1–3 problems considered and dropped, with reason}.
- *Role in this PR:* {headline-payoff / fluency-builder, per §7.1 of plan}.

**Source:** Griffiths, *Introduction to Quantum Mechanics*, **3e Problem 1.4** (Cambridge 2018; Griffiths & Schroeter). 2e equivalence: Problem 1.4 (Pearson 2005; Griffiths). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** {1–3 sentences. Cite specific equation numbers from Griffiths where relevant.}

**Setup:** {Hilbert-space structure, operators, boundary conditions, given quantities, unknowns.}

#### (a) Griffiths 2e textbook solution

{Standard non-relativistic derivation, step-by-step.}

#### (b) Griffiths 3e textbook solution

{Either the same derivation in 3e form, or — if 3e is unchanged — "Identical to (a)". When 3e reworded the statement, indicate the change.}

#### (c) Proper-time / dual-theory reformulation

Apply the canonical proper-time Hamiltonian from [`_proper_time_K_cheatsheet.md`](_proper_time_K_cheatsheet.md):

- `K = H²/(2mc²) + mc²/2` (DRQM I Eq. I.6)
- `H_0 = √(c²π² + m²c⁴)` (free-particle limit)
- Non-relativistic reduction: `K ≈ mc²/2 + π²/(2m) + mc²/2 = mc² + π²/(2m)` for `π ≪ mc`, so `K - mc² → p²/(2m)` (textbook Schrödinger kinetic energy).

{Proper-time derivation. Mathematica MCP check confirming the reduction.}

**Reduction-to-textbook verdict:** ✅ matches Griffiths in `u \ll c` limit / ⚠ `O((u/c)²)` correction / 🔴 disagreement requiring author review.

**Comparison:**

| Quantity | Griffiths 2e | Griffiths 3e | Proper-time |
|---|---|---|---|
| … | … | … | … |

**Standard-equation comparison:** which textbook reference reproduces the same answer (Sakurai / Shankar / Bethe-Salpeter for hydrogen problems).

**Verdict:** ✅ all three solutions consistent / ⚠ proper-time deviates at `O((v/c)^k)` (see below) / 🔴 inconsistency found.

**Notes for author review (if ⚠ or 🔴):** {description suitable for inclusion in [`FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md).}

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh{NN}_P{C_P}.wl`](../Mathematica_Notebooks/Quantum_Mechanics/).
````

## Branched treatment for `r_e`-finding-touching problems

If a problem's proper-time derivation invokes the anomalous-`g`-factor numerics of DRQM I §III.D (the `r_e ≈ 0.499857...` result that gives `g = -2.0005714` vs measured `-2.00231930...`), the proper-time section carries a **branched treatment** analogous to #42's §5.1 Eq. 24 workflow. The branched subsections are `(c) as-published DRQM I` and `(c′) with corrected r_e`.

The branched treatment is expected to engage primarily in Ch. 7 fine-structure problems and in any future Ch. 4 hydrogen problem that depends on the anomalous `g`-factor. Most non-relativistic Griffiths problems do not touch the `r_e` finding.

## Voice compliance

Prose in (a), (b), (c), and verdict sections follows the voice guide at [`Roadmapping/Tooling/VOICE_MATCH_GILL.md`](../Tooling/VOICE_MATCH_GILL.md). The read-aloud test in §5 of that document is the compliance check.
