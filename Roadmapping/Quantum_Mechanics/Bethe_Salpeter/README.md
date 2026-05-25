# Bethe–Salpeter precision predictions under proper-time / dual-theory formulation

Rederivation of the canonical precision results in **Bethe & Salpeter, *Quantum Mechanics of One- and Two-Electron Atoms*** (Plenum 1957; Springer 1977 reprint) under the proper-time / dual-theory formulation, with explicit comparison against modern measured values (CODATA-2018 / PDG-2024 wherever applicable).

This is the **precision-experiment** companion to the [Griffiths pedagogical campaign](../Griffiths/) (#49). The Griffiths thread demonstrates that the proper-time `K = H₀²/(2mc²) + mc²/2` reduces to non-relativistic QM. This thread tests where the proper-time formulation has *measurable* consequences beyond that reduction.

Governing plan: [`.dev/tasks/50-bethe-salpeter-precision-predictions.md`](../../../.dev/tasks/50-bethe-salpeter-precision-predictions.md) (issue [#50](https://github.com/temoTxt/PyPhysics/issues/50)).

<!-- TODO: human reviews and fills in — confirms the campaign's framing as the precision-experiment counterpart to #49, and the load-bearing role of PRs C / E / F (fine structure, Lamb shift, hyperfine) in determining whether the dual-theory framework is consistent with current measurements -->

## Anchor sources

- **Bethe–Salpeter (1957/1977)** — single canonical edition. PDF held locally per repository policy (`pdf_status: acquired`, in-copyright Springer reprint; only fair-use quoted equations + section refs appear in committed markdown). Bibliography stub: [`bethe1977_one_two_electron_atoms`](../../History/Bibliography/Retrospective/bethe1977_one_two_electron_atoms.md).
- **DRQM I** ([`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md)) — Eqs. I.6 canonical `K`, II.1–II.3 dual Dirac ✅. Eq. III.D `r_e` finding 🔴 propagates into PRs C and F under branched treatment.
- **Analytic Dirac representation** ([`Analytic_Representation_of_The_Dirac_Equation.md`](../../Equation_Verification/Analytic_Representation_of_The_Dirac_Equation.md)) — used for fine-structure derivations.
- **The Classical Electron Problem** ([`The_Classical_Electron_Problem.md`](../../Equation_Verification/The_Classical_Electron_Problem.md)) — supplies the proper-time radiation-reaction structure invoked in the Lamb-shift route (PR E).
- **Proper-time Liénard–Wiechert third term** from [`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) Eq. (7) — used in PR G dipole / multipole work.

## Campaign status

| PR | Bethe–Salpeter sections | Role | Results | Status |
|---|---|---|---|---|
| **PR A** | §§1–7 — Non-rel hydrogen | scaffold + non-rel reduction `K → p²/2m + V₀` | 4 | ✅ drafted 2026-05-25 |
| **PR B** | §§8–13 — Matrix elements + transitions | dipole, multipole, oscillator strengths | 3 | ✅ drafted 2026-05-25 |
| **PR C** ⭐ | **§14 Fine structure (pivot)** | Dirac vs dual Dirac for H 2P₃/₂ – 2P₁/₂ | 3 | ✅ drafted 2026-05-25 (acceptance criterion 2) |
| **PR D** | §§15–18 — Higher-order rel corrections | `(Zα)⁴` terms; FW vs proper-time | 3 | ✅ drafted 2026-05-25 |
| **PR E** ⭐ | **§§19–21 Lamb shift (headline)** | self-energy via `K` + radiation-reaction route | 3 | ✅ drafted 2026-05-25 (acceptance criterion 3) |
| **PR F** ⭐ | **§22 Hyperfine structure (headline)** | H 21-cm line; muonium/positronium deferred to PR I | 2 | ✅ drafted 2026-05-25 |
| **PR G** | §§23–37 — Interaction with radiation | photoionisation, multipole; third-term effect on dipole | 3–4 | pending |
| **PR H** | §§47–60 — Helium ground state | two-electron variational; proper-time energy operator | 3–4 | pending |
| **PR I** | §§61–80 — Helium excited states | two-electron spectroscopy; positronium / muonium where applicable | 3–4 | pending |
| **PR J** | Cross-comparison summary | table of all proper-time vs measured; flagged shifts | (chapter) | pending |

Headline (⭐) PRs are gated by the acceptance criteria in issue #50.

## Acceptance criteria for closing #50

1. PR A merged — `Bethe_Salpeter/` scaffolded + §§1–7 complete (3–4 results), with reduction `K → p²/2m + V₀` verified in Wolfram MCP.
2. PR C merged — §14 fine structure rederived; comparison against measured `2P₃/₂ – 2P₁/₂` splitting recorded.
3. PR E merged — §§19–21 Lamb shift treatment; explicit numerical prediction vs `1057.845(9)` MHz measured.
4. Any ⚠ / 🔴 finding cross-posted to [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md).

Subsequent PRs continue under this campaign but do not require keeping #50 open.

## Conventions

- [`_template_result.md`](_template_result.md) — per-result document template (Bethe–Salpeter textbook + proper-time + modern measurement three-way).
- [`_proper_time_K_cheatsheet.md`](../_proper_time_K_cheatsheet.md) — canonical `K` and dual Dirac substitution rules (shared with Griffiths campaign).

Voice: substantive prose follows [`VOICE_MATCH_GILL.md`](../../Tooling/VOICE_MATCH_GILL.md). The "mathematically equivalent but not physically equivalent" anchor is load-bearing in PRs C, E, F where physical interpretation of any shift is the point of the result.

Crocco compliance (per [`CROCCO_COMPLIANCE.md`](../../Tooling/CROCCO_COMPLIANCE.md)): Mathematica MCP + CODATA / PDG lookups + Bethe–Salpeter transcription = **pragmatic**; physical interpretation of any divergence = **substantive**, with per-paragraph `<!-- TODO -->` blocks throughout headline PRs.

## Honest limits

This campaign tests the dual-theory framework against precision measurements. Three caveats are explicit (per [§7 of the #50 plan](../../../.dev/tasks/50-bethe-salpeter-precision-predictions.md#7-honest-framing)):

1. **The Lamb-shift route is the Bethe (1947) mass-renormalisation estimate, not a one-loop QED calculation.** The proper-time substitution follows Bethe's argument with the dual-Dirac propagator + radiation-reaction structure. Agreement at the Bethe-estimate precision (~few percent of full Lamb) is what we can honestly deliver; sub-MHz precision needs a one-loop calculation that the dual-theory framework has not yet produced.
2. **Hyperfine and fine-structure predictions depend on the anomalous-`g`-factor**, which DRQM I §III.D's flagged `r_e` finding propagates into. PRs C and F carry **branched treatment**: as-published `r_e` (gives `g = -2.0005714`) vs corrected `r_e` (gives measured `g = -2.00231930…`). Both branches are recorded; neither is the campaign's "verdict."
3. **Two-loop corrections, positronium / muonium full precision, heavy-element spectroscopy** are out of scope. The campaign's scope is the one-loop / Bethe-estimate precision floor.

These caveats are the campaign's honest framing. They do not block the work; they set the precision floor against which verdicts are issued.
