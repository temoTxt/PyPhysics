# Quantum Mechanics — proper-time / dual-theory QM (two campaigns)

The Quantum Mechanics thread runs two coupled campaigns:

1. **[`Griffiths/`](Griffiths/) — pedagogical correspondence** (issue [#49](https://github.com/temoTxt/PyPhysics/issues/49), PR #52): canonical problems of Griffiths, *Introduction to Quantum Mechanics* (2e + 3e), rederived under the proper-time formulation. Demonstrates the reduction `K → p²/2m` at the undergraduate level.
2. **[`Bethe_Salpeter/`](Bethe_Salpeter/) — precision-experiment discrimination** (issue [#50](https://github.com/temoTxt/PyPhysics/issues/50)): canonical results in Bethe & Salpeter, *Quantum Mechanics of One- and Two-Electron Atoms* (Springer 1977), rederived under the proper-time formulation and compared against modern measured values (CODATA-2018 / PDG-2024). Where Griffiths demonstrates the *correspondence*, Bethe–Salpeter is where the proper-time formulation has *measurable* consequences (Lamb shift, hyperfine splitting, fine structure).

Both campaigns use the canonical Hamiltonian

$$K = \frac{H^{2}}{2 m c^{2}} + \frac{m c^{2}}{2}, \qquad H_{0} = \sqrt{c^{2}\boldsymbol\pi^{2} + m^{2}c^{4}}$$

of [`Dual_Relativistic_Quantum_Mechanics_I.md`](../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) (Eqs. I.6, II.1–II.3 ✅) and the analytic Dirac representation of [`Analytic_Representation_of_The_Dirac_Equation.md`](../Equation_Verification/Analytic_Representation_of_The_Dirac_Equation.md).

The two formulations are *mathematically equivalent but not physically equivalent* (Gill–Zachary) at the relativistic level. The Griffiths campaign demonstrates the reduction `K \to p^{2}/2m` in the non-relativistic limit problem-by-problem; the Bethe–Salpeter campaign tests where the two formulations remain experimentally distinguishable.

The threads are governed by the plans at [`.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md`](../../.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md) and [`.dev/tasks/50-bethe-salpeter-precision-predictions.md`](../../.dev/tasks/50-bethe-salpeter-precision-predictions.md). Both inherit the [§13 honest-framing discipline](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix) from [#42 Electromagnetism / Jackson](https://github.com/temoTxt/PyPhysics/issues/42).

<!-- TODO: human reviews and fills in — confirms the framing of the two coupled QM campaigns (Griffiths pedagogical-correspondence + Bethe–Salpeter precision-experiment), and that the latter's Lamb-shift / hyperfine verdicts are the framework's most consequential experimental tests in this repo -->

## Scope

- **Griffiths campaign (#49)**: ~5–10 problems per chapter, ~60–100 total across 12 chapters. Status: ✅ drafted 2026-05-25 (PR #52, PRs A–L).
- **Bethe–Salpeter campaign (#50)**: ~3–6 results per major section, ~30–50 total across 10 PRs (A–J). Status: in progress.

## Griffiths campaign — status

| PR | Scope | Problems | Status |
|---|---|---|---|
| **PR A** | Ch. 1 Wave Function | 3 | ✅ drafted 2026-05-25 |
| **PR B** | Ch. 2 TI Schrödinger Equation | 4 | ✅ drafted 2026-05-25 |
| **PR C** | Ch. 3 Formalism | 3 | ✅ drafted 2026-05-25 |
| **PR D** ⭐ | Ch. 4 QM in 3D (hydrogen pivot) | 5 | ✅ drafted 2026-05-25 |
| **PR E** | Ch. 5 Identical Particles | 3 | ✅ drafted 2026-05-25 |
| **PR F** | Ch. 6 Symmetries | 3 | ✅ drafted 2026-05-25 |
| **PR G** ⭐ | Ch. 7 TI Perturbation Theory (fine structure) | 5 | ✅ drafted 2026-05-25 |
| **PR H** | Ch. 8 Variational Principle (helium) | 3 | ✅ drafted 2026-05-25 |
| **PR I** | Ch. 9 WKB | 3 | ✅ drafted 2026-05-25 |
| **PR J** | Ch. 10 Scattering | 4 | ✅ drafted 2026-05-25 |
| **PR K** | Ch. 11 Quantum Dynamics | 4 | ✅ drafted 2026-05-25 |
| **PR L** | Ch. 12 Afterword | 2 | ✅ drafted 2026-05-25 |

Headline (⭐) chapters carry the load-bearing proper-time content; PRs A–C and the backfill PRs are pedagogical-reduction demonstrations.

## Bethe–Salpeter campaign — status

| PR | Scope | Results | Status |
|---|---|---|---|
| **PR A** | §§1–7 Non-rel hydrogen (scaffold) | 3–4 | in progress |
| **PR B** | §§8–13 Matrix elements + transitions | 3–4 | pending |
| **PR C** ⭐ | §14 Fine structure (pivot) | 3–4 | pending |
| **PR D** | §§15–18 Higher-order rel corrections | 3–4 | pending |
| **PR E** ⭐ | §§19–21 Lamb shift (headline) | 3–4 | pending |
| **PR F** ⭐ | §22 Hyperfine structure (headline) | 2–3 | pending |
| **PR G** | §§23–37 Interaction with radiation | 3–4 | pending |
| **PR H** | §§47–60 Helium ground state | 3–4 | pending |
| **PR I** | §§61–80 Helium excited states | 3–4 | pending |
| **PR J** | Cross-comparison summary | (chapter) | pending |

Headline (⭐) PRs C, E, F carry the campaign's load-bearing experimental discrimination; PR A scaffolds + non-rel reduction, PRs G–I bridge to two-electron precision spectroscopy.

## Conventions

- [`_template_problem.md`](_template_problem.md) — Griffiths per-problem template (2e + 3e + proper-time three-way).
- [`_proper_time_K_cheatsheet.md`](_proper_time_K_cheatsheet.md) — canonical `K` and dual Dirac substitution rules (shared by both campaigns).
- [`Griffiths/README.md`](Griffiths/README.md) — Griffiths-specific status table across all 12 chapters.
- [`Bethe_Salpeter/README.md`](Bethe_Salpeter/README.md) — Bethe–Salpeter campaign status table + per-result template variant (with modern-measurement column).

Voice (per [`VOICE_MATCH_GILL.md`](../Tooling/VOICE_MATCH_GILL.md)): substantive QM prose follows Gill's published-paper voice. Griffiths' voice is preserved only in paraphrased quotations of textbook statements.

Crocco compliance (per [`CROCCO_COMPLIANCE.md`](../Tooling/CROCCO_COMPLIANCE.md)): Mathematica + 2e↔3e reconciliation + paraphrasing = pragmatic; physical interpretation of any divergence = substantive.

## Source material

- **Griffiths, *Introduction to Quantum Mechanics*** — both editions cited. Bibliography stubs: [`griffiths2005_qm_2e`](../History/Bibliography/Retrospective/griffiths2005_qm_2e.md) and [`griffiths2018_qm_3e`](../History/Bibliography/Retrospective/griffiths2018_qm_3e.md), both `pdf_status: acquired` (in-copyright, fair-use academic quotation).
- **Verified Gill corpus**: DRQM I (Eqs. I.6, II.1–II.3 ✅) and *Analytic Representation of the Dirac Equation*.

## Honest limits

These threads are exploratory work in the proper-time formulation, not a validation. Per-problem and per-result claims are conditional predictions of the form "*if* the dual-theory framework is the correct formulation, *then* X."

- **Griffiths (#49)**: for the majority of problems (non-relativistic), the proper-time formulation reduces to textbook QM; the campaign demonstrates the reduction but does not test the framework.
- **Bethe–Salpeter (#50)**: the precision-experiment campaign *does* test the framework, in the regime where measurements distinguish the proper-time predictions from textbook QED. Its honest limits are recorded in [§7 of the #50 plan](../../.dev/tasks/50-bethe-salpeter-precision-predictions.md#7-honest-framing): the Lamb-shift treatment is at Bethe-estimate precision (not full one-loop), and the hyperfine + fine-structure predictions carry a branched treatment on the DRQM I §III.D `r_e` finding.
