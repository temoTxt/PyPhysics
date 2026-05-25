# Quantum Mechanics — Griffiths canonical problems × proper-time / dual-theory QM

This research thread works the canonical problems of Griffiths, *Introduction to Quantum Mechanics* (2nd ed. Pearson 2005; 3rd ed. Griffiths & Schroeter, Cambridge 2018) in two editions and three formulations: **Griffiths 2e**, **Griffiths 3e**, and the **proper-time / dual-theory reformulation** using the canonical Hamiltonian

$$K = \frac{H^{2}}{2 m c^{2}} + \frac{m c^{2}}{2}, \qquad H_{0} = \sqrt{c^{2}\boldsymbol\pi^{2} + m^{2}c^{4}}$$

of [`Dual_Relativistic_Quantum_Mechanics_I.md`](../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) (Eqs. I.6, II.1–II.3 ✅) and the analytic Dirac representation of [`Analytic_Representation_of_The_Dirac_Equation.md`](../Equation_Verification/Analytic_Representation_of_The_Dirac_Equation.md).

The two formulations are *mathematically equivalent but not physically equivalent* (Gill–Zachary) at the relativistic level; the campaign demonstrates the reduction `K \to p^{2}/2m` in the non-relativistic limit problem-by-problem, and identifies where the formulations diverge (notably fine structure, Ch. 7).

The thread is governed by the plan at [`.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md`](../../.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md) (issue [#49](https://github.com/temoTxt/PyPhysics/issues/49)). It is a sibling campaign to [#42 Electromagnetism / Jackson](https://github.com/temoTxt/PyPhysics/issues/42) and inherits its [§13 honest-framing discipline](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix).

<!-- TODO: human reviews and fills in — confirms the framing of this thread as a sibling-to-#42 campaign with inherited honest-framing discipline -->

## Scope

Selected canonical problems (~5–10 per chapter, ~60–100 total across 12 Griffiths chapters).

## Status

| PR | Scope | Problems | Status |
|---|---|---|---|
| **PR A** | Ch. 1 Wave Function | 3 | ✅ drafted 2026-05-25 |
| **PR B** | Ch. 2 TI Schrödinger Equation | 4 | ✅ drafted 2026-05-25 |
| **PR C** | Ch. 3 Formalism | 3 | ✅ drafted 2026-05-25 |
| **PR D** ⭐ | Ch. 4 QM in 3D (hydrogen pivot) | 5 | ✅ drafted 2026-05-25 |
| **PR E** | Ch. 5 Identical Particles | 3 | ✅ drafted 2026-05-25 |
| **PR F** | Ch. 6 Symmetries | 3 | ✅ drafted 2026-05-25 |
| **PR G** ⭐ | Ch. 7 TI Perturbation Theory (fine structure) | 5–7 | planned |
| PR H | Ch. 8 Variational Principle (helium) | 3–4 | planned |
| PR I | Ch. 9 WKB | 3–4 | planned |
| PR J | Ch. 10 Scattering | 4–6 | planned |
| PR K | Ch. 11 Quantum Dynamics | 4–6 | planned |
| PR L | Ch. 12 Afterword | 2–3 | planned |

Headline (⭐) chapters carry the load-bearing proper-time content; PRs A–C and the backfill PRs are pedagogical-reduction demonstrations.

## Conventions

- [`_template_problem.md`](_template_problem.md) — per-problem document template (2e + 3e + proper-time three-way).
- [`_proper_time_K_cheatsheet.md`](_proper_time_K_cheatsheet.md) — canonical `K` and dual Dirac substitution rules.
- [`Griffiths/README.md`](Griffiths/README.md) — Griffiths-specific status table across all 12 chapters.

Voice (per [`VOICE_MATCH_GILL.md`](../Tooling/VOICE_MATCH_GILL.md)): substantive QM prose follows Gill's published-paper voice. Griffiths' voice is preserved only in paraphrased quotations of textbook statements.

Crocco compliance (per [`CROCCO_COMPLIANCE.md`](../Tooling/CROCCO_COMPLIANCE.md)): Mathematica + 2e↔3e reconciliation + paraphrasing = pragmatic; physical interpretation of any divergence = substantive.

## Source material

- **Griffiths, *Introduction to Quantum Mechanics*** — both editions cited. Bibliography stubs: [`griffiths2005_qm_2e`](../History/Bibliography/Retrospective/griffiths2005_qm_2e.md) and [`griffiths2018_qm_3e`](../History/Bibliography/Retrospective/griffiths2018_qm_3e.md), both `pdf_status: acquired` (in-copyright, fair-use academic quotation).
- **Verified Gill corpus**: DRQM I (Eqs. I.6, II.1–II.3 ✅) and *Analytic Representation of the Dirac Equation*.

## Honest limits

This thread is exploratory work in the proper-time formulation, not a validation. Per-problem results are conditional predictions of the form "*if* the dual-theory framework is the correct formulation, *then* X." For the majority of Griffiths (non-relativistic problems), the proper-time formulation reduces to textbook QM; the campaign demonstrates the reduction but does not test the framework. The experimental-discrimination work is the subject of the Bethe–Salpeter precision-predictions umbrella, opened concurrently with #49.
