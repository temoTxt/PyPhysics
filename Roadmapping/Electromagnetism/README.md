# Electromagnetism — Jackson canonical problems × proper-time reformulation

This research thread works the canonical problems of Jackson's *Classical Electrodynamics* (2nd ed., 1975, Gaussian units; 3rd ed., 1998, mixed SI / Gaussian) in three unit systems: CGS, SI, and the proper-time reformulation established by Gill and Zachary in [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]. We follow the substitution rules `w/c = u/b`, `(1/c)∂_t = (1/b)∂_τ`, and `c → b` (with `b² = c² + u²`); we record both forms of the result and, where the proper-time prediction differs from a pure `c → b` redressing of the classical answer, we name the extra term and its physical interpretation.

The two formulations are *mathematically equivalent but not physically equivalent* — the local clock of the source encodes information about the particle's acceleration that the observer's clock does not. This thread records what those rules predict for the problems a graduate student already knows in the classical formulation; it does not propose new physics.

The thread is governed by the plan at [`.dev/tasks/42-electromagnetism-jackson-proper-time.md`](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md) (issue [#42](https://github.com/temoTxt/PyPhysics/issues/42)).

<!-- TODO: human reviews and fills in — confirms the framing above accurately states the thread's scope and conditional posture before per-problem documents start to accumulate -->

## Scope

Selected canonical problems (~5–10 per chapter, ~50–100 total). Not exhaustive coverage. The campaign's value lies in selection, not coverage; see [§7 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list).

## Status

| PR | Scope | Problems | Status |
|---|---|---|---|
| **PR 0** | Chs. 1–2, 5 (fluency warm-up — `u = 0` or steady-current) | 4 | ✅ closed 2026-05-24 (4/4 drafted, retrospective in [§7.4 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#74-pr-0-retrospective-closed-2026-05-24)) |
| **PR A** | Ch. 6 Maxwell + macroscopic media | 5 (6.1, 6.4, 6.5, 6.11, 6.20) | ✅ closed 2026-05-24 (5/5 drafted, retrospective in [§7.5 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#75-pr-a-retrospective-closed-2026-05-24)) |
| **PR B** | Ch. 11 Special Relativity | 5 | ✅ closed 2026-05-24 (5/5 drafted, retrospective in [§7.6 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#76-pr-b-retrospective-closed-2026-05-24)) |
| **PR C** | Ch. 12 Relativistic Dynamics | 5 | ✅ closed 2026-05-24 (5/5 drafted, retrospective in [§7.7 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#77-pr-c-retrospective-closed-2026-05-24)) |
| **PR D** | Ch. 14 Radiation by Moving Charges | 5 | ✅ closed 2026-05-24 (5/5 drafted, retrospective in [§7.8 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#78-pr-d-retrospective-closed-2026-05-24)) |
| **PR E** | Ch. 16 Radiation Damping | 4 | ✅ closed 2026-05-24 (4/4 drafted, retrospective in [§7.9 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#79-pr-e-retrospective-closed-2026-05-24)) |
| **PR F** | Ch. 7 Plane EM Waves | 4 | ✅ drafted 2026-05-24 |
| **PR G** | Ch. 4 Multipoles + Macroscopic Media | 4 | ✅ drafted 2026-05-24 |
| **PR H+** | Remaining backfill: Chs. 1–3 (supplementary), Chs. 5 (supplementary), 8–10, 13, 15 | 4–6 each | planned |

Realistic completion: 9–18 months at part-time pace per [§13.2 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#132-objections-with-no-honest-mitigation).

## Conventions

Three meta-documents govern every per-problem write-up. Authors should read them before drafting:

- [`_template_problem.md`](_template_problem.md) — the per-problem document template. Every `Jackson/Ch??_*.md` mirrors this structure.
- [`_proper_time_cheatsheet.md`](_proper_time_cheatsheet.md) — one-page reference to the substitution rules, extracted from the verified Gill–Zachary Maxwell paper.
- [`Jackson/README.md`](Jackson/README.md) — Jackson-specific status table across all 14 chapters.

Voice and AI-use compliance are governed at the repository level:

- [`Roadmapping/Tooling/VOICE_MATCH_GILL.md`](../Tooling/VOICE_MATCH_GILL.md) — the substantive prose in per-problem documents matches Gill's published-paper voice. The read-aloud test in §5 of that document is the compliance check.
- [`Roadmapping/Tooling/CROCCO_COMPLIANCE.md`](../Tooling/CROCCO_COMPLIANCE.md) — AI-use disclosure. Substantive uses (problem selection, narrative framing, interpretive paragraphs) carry per-paragraph `<!-- TODO -->` blocks; granularity is reassessed at PR B per [§13.5 D2](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Source material

- **Jackson, *Classical Electrodynamics*** — both editions cited. Bibliography entries (to be scaffolded in PR 0): `Roadmapping/History/Bibliography/Retrospective/jackson1975_classical_electrodynamics.md` and `…/jackson1998_classical_electrodynamics.md`. Per-problem documents follow [§4 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling) on which unit system applies to which chapter.
- **Gill and Zachary, *Two Mathematically Equivalent Versions of Maxwell's Equations*** — the source of record for the proper-time substitution rules. Verified at [`Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md). Eqs. (1)–(23) are confirmed ✅; Eq. (24) carries an unresolved finding and triggers branched treatment per [§5.1 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#51-branched-treatment-for-eq-24-touching-problems).

## Related work

- [`Roadmapping/Equation_Verification/`](../Equation_Verification/) — verification thread for the Gill papers' equations. This thread depends on Eqs. (1)–(23) of the Maxwell paper being verified there.
- [`Roadmapping/History/Podcast/`](../History/Podcast/) — five problems from this thread are elevated as podcast questions across episodes 02, 03, 05, and 06. See [§12 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#12-podcast-elevation-candidates).
- Issue [#43](https://github.com/temoTxt/PyPhysics/issues/43) — experimental-comparison sub-investigation against Cole 2018, Poder 2018, and Wistisen 2018 radiation-reaction data.

## Honest limits

This thread is exploratory work in the Gill–Zachary framework, not a validation of it. Per-problem results are conditional predictions of the form "*if* the Gill–Zachary framework is the correct formulation, *then* X." The author has a prior involvement in this framework (co-author of *Dual Relativistic Quantum Mechanics I*, 2021); no amount of careful documentation makes this a neutral evaluation. See [§13 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix) for the full honest framing.
