# Per-problem document template

This template defines the structure of every per-problem document under [`Jackson/Ch??_*.md`](Jackson/). Each problem is presented in this fixed format so that a reader can locate the same kind of content in the same place across the campaign. The structure was proposed in [§3 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#3-per-problem-template-proposed) and is reproduced verbatim below; refinements identified during PR 0 or PR A will be folded back into this file.

## Naming convention

Problem headings use `J{N}e-P{C}.{P}`, where `{N}e` is the edition (`2e` or `3e`), `{C}` is the chapter, and `{P}` the problem number — e.g., `J3e-P6.1` = Jackson 3rd ed., Chapter 6, Problem 6.1. When a problem appears in both editions with matching numbering, the 3e key is canonical and the 2e equivalence is recorded in the source line.

## Pre-problem header (every chapter file)

Each `Ch??_*.md` opens with a short header recording the chapter's role in the campaign, the PR in which the chapter is being worked, the unit-system regime per [§4](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), and the list of problems contained.

## Per-problem template (the body of every problem entry)

The template below is the canonical structure. Sections marked as conditional (e.g., the SI section in CGS-only chapters, the `(c′)` branch under Eq. 24) appear only when the problem warrants them; the unconditional sections are required for every problem.

````markdown
### Problem J3e-P6.1 — short title

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* {one-line reason this problem is in scope for this PR}.
- *Alternatives considered:* {1–3 problems considered and dropped, with reason}.
- *Role in this PR:* {headline-payoff / fluency-builder, per §7.1 of the plan}.

**Source:** Jackson, *Classical Electrodynamics*, 3e §6.1 Problem 6.1 (and 2e §6.1 Problem 6.1, equivalent). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** {1–3 sentences. Cite specific equation numbers from Jackson where relevant.}

**Setup:** {Geometry, boundary conditions, given quantities, unknowns. Diagram if needed.}

#### (a) Classical solution — Gaussian (CGS)

{Derivation in Gill voice. Mathematica MCP check inline using the same format as the Equation_Verification documents. Each substantive interpretive paragraph carries its own `<!-- TODO: human reviews and fills in -->` block for Ch. 6 (per §13.5 D2 of the plan); granularity is reassessed at PR B.}

#### (b) Classical solution — SI

{Either the full derivation in SI (recommended when the work itself differs), or — if the only difference is constants — a conversion table:

| Quantity | Gaussian → SI |
|---|---|
| `E`-field | `E_SI = E_CGS / (4πε₀)^(1/2)` |
| ... | ... |

This section is omitted in chapters where Jackson uses CGS only (all of 2e; 3e Chs. 11+). See §4 of the plan.}

#### (c) Proper-time reformulation

Apply the substitution rules from [`_proper_time_cheatsheet.md`](_proper_time_cheatsheet.md):

- `c → b = √(c² + u²)`
- `t → τ`, `∂_t → (b/c) ∂_τ`
- `w → u`, `J → (b/c) J`
- `∂_τ(1/b) = −(u·a)/b³`, which re-expresses as `−(u·a)/b⁴ · ∂_τ` once collected as a wave-equation coefficient.

{Derivation in Gill voice. Mathematica MCP check. The derivation works from the underlying Maxwell / force law and re-derives the derived quantity from scratch in the proper-time frame; we do not simply substitute `c → b` into the classical answer.}

#### (c′) Proper-time reformulation — with proposed Eq. 24 correction

*This section appears only when the derivation invokes Maxwell-paper Eq. 24.* Per [§5.1 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#51-branched-treatment-for-eq-24-touching-problems), the derivation is repeated with the proposed correction from [`FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md) (restored factor of `c`, plus the `V²/(2mc²)` term).

#### (c) vs (c′) comparison

*This subsection also appears only when the `(c′)` branch is engaged.* One paragraph summarising how the two branches differ at the answer level. If the branches agree numerically, we record that they do; if they disagree, the disagreement is itself a piece of evidence on which Eq. 24 form is consistent with the rest of the framework.

**Comparison:**

| Quantity | Classical (CGS) | Classical (SI) | Proper-time |
|---|---|---|---|
| ... | ... | ... | ... |

For two-system problems (Chs. 11+ and PR 0 fluency-builders), the SI column is omitted and the table is two-column.

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no / ⚠ yes — {one-line reason}.

**Verdict:** ✅ all sections consistent / ⚠ proper-time deviates (see below) / ❌ inconsistency found.

When the `(c′)` branch is engaged, the verdict is recorded *per branch*.

**Notes for author review (if ⚠ or ❌):** {description suitable for inclusion in [`FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md).}
````

## Companion artifact

Every per-problem document is paired with a companion Wolfram Language notebook at `Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonChNN_PNN.wl`, runnable independent of the Mathematica MCP. The notebook contains the same symbolic checks recorded inline in the document, in a form that another reader can execute end-to-end.

## Voice compliance

Prose in `(a)`, `(b)`, `(c)`, `(c′)`, and the comparison/verdict sections follows the voice guide at [`Roadmapping/Tooling/VOICE_MATCH_GILL.md`](../Tooling/VOICE_MATCH_GILL.md). The read-aloud test in §5 of that document is the compliance check. Selection provenance and verdict-line content are operational and need not match Gill's voice.

## A note on the template's stability

This template is provisional. PR 0 and PR A together will exercise it on ~9 problems across four different chapters and three difficulty levels; refinements identified in either PR are reflected here before PR B starts. Three template choices remain open per [§9 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#9-open-questions--decisions-deferred-to-the-human) and are likely to be resolved by then.
