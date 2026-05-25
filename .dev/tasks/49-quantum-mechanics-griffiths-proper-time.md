# Implementation plan: Quantum Mechanics — Griffiths × 2e/3e × proper-time / dual-theory QM

**Tracks:** issue [#49](https://github.com/temoTxt/PyPhysics/issues/49).
**Status:** plan only at start. Sibling campaign to [#42 — Electromagnetism / Jackson](42-electromagnetism-jackson-proper-time.md).
**Anchor sources:** verified [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) (Eqs. I.6 canonical `K`, II.1–II.3 dual Dirac ✅) and [`Analytic_Representation_of_The_Dirac_Equation.md`](../../Roadmapping/Equation_Verification/Analytic_Representation_of_The_Dirac_Equation.md).

This plan is **deliberately concise** because it inherits most of its discipline from [§13 of the #42 plan](42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix). Differences from the #42 campaign are flagged below; conventions not flagged are inherited.

---

## 1. Goal and scope

Work the canonical problems of **Griffiths, *Introduction to Quantum Mechanics*** (2nd ed. Pearson 2005; 3rd ed. Griffiths & Schroeter, Cambridge 2018) in two editions and three formulations:

1. **Griffiths 2e** — pre-Schroeter problem set.
2. **Griffiths 3e** — current standard; renumbered Ch. 11 (Quantum Dynamics) and Ch. 12 (Afterword).
3. **Proper-time / dual-theory reformulation** — using the canonical Hamiltonian `K = H²/(2mc²) + mc²/2` with `H_0 = √(c²π² + m²c⁴)` and the dual Dirac equation from DRQM I §II.

Per-problem output: one document showing all three solutions side-by-side, plus an explicit verdict.

### Why this earns its own thread

- **Pedagogical bridge.** Demonstrates that the dual theory passes the undergraduate-QM-correspondence test (proper-time `K` reduces to `p²/2m` at `u \ll c`).
- **Headline payoff at Ch. 7 fine structure.** Proper-time vs Pauli/FW comparison for hydrogen fine structure is the natural falsifiable test.
- **2e/3e dual-edition coverage** mirrors #42's CGS/SI dual coverage.

### Explicit non-goals (per [§1 of #42 plan](42-electromagnetism-jackson-proper-time.md#1-goal-and-scope))

- No exhaustive coverage. Selected canonical problems (~5–10 per chapter, ~60–100 total).
- No reproduction of Griffiths' problem statements verbatim. Paraphrased.
- No new physics claims about the dual theory. Conditional framing throughout.
- Most problems are pedagogical-reduction demonstrations, not experimental discriminators.

---

## 2. Repository structure

```
Roadmapping/Quantum_Mechanics/
├── README.md                              # thread overview + status table
├── _template_problem.md                   # per-problem template (mirrors Electromagnetism's)
├── _proper_time_K_cheatsheet.md           # canonical K and dual Dirac substitution rules
└── Griffiths/
    ├── README.md                          # Griffiths-specific status table across all 12 chapters
    ├── Ch01_Wave_Function.md
    ├── Ch02_TI_Schrodinger.md
    ├── ...
    ├── Ch04_QM_3D.md                      ⭐ pivot (hydrogen via proper-time)
    └── Ch07_TI_Perturbation_Theory.md     ⭐ headline-payoff (fine structure)
```

Sibling to the existing [`Roadmapping/Electromagnetism/`](../../Roadmapping/Electromagnetism/) thread.

---

## 3. Per-problem template

Lives at `Roadmapping/Quantum_Mechanics/_template_problem.md` (to be drafted in PR A). Mirrors the Electromagnetism template; key differences:

- **Two source citations** (Griffiths 2e + 3e) rather than two unit systems.
- **Three solution sections**: (a) Griffiths 2e textbook, (b) Griffiths 3e textbook (or note that it's unchanged), (c) proper-time / dual-theory.
- **Reduction-to-textbook verdict**: ✅ matches in `u \ll c` limit / ⚠ `O((u/c)²)` correction / 🔴 disagreement.

---

## 4. PR sequencing (12 PRs A–L per issue body)

Order = increasing relativistic / proper-time content:

| PR | Griffiths chapter (3e) | Role | Problems |
|---|---|---|---|
| **PR A** | Ch. 1 Wave Function | pedagogical-foundation; reduction `K → p²/2m` | 3 |
| **PR B** | Ch. 2 TI Schrödinger Equation | infinite well, SHO, delta, finite well | 4–6 |
| **PR C** | Ch. 3 Formalism | Hilbert-space recap | 3–4 |
| **PR D** ⭐ | **Ch. 4 QM in 3D (hydrogen)** | **pivot**: hydrogen via canonical `K` | 5–7 |
| PR E | Ch. 5 Identical Particles | Bose/Fermi statistics | 3–4 |
| PR F | Ch. 6 Symmetries | Noether-type derivations | 3–4 |
| **PR G** ⭐ | **Ch. 7 TI Perturbation Theory** | **headline**: fine structure of hydrogen | 5–7 |
| PR H | Ch. 8 Variational Principle | helium ground state | 3–4 |
| PR I | Ch. 9 WKB | tunnelling | 3–4 |
| PR J | Ch. 10 Scattering | Born approximation, partial waves | 4–6 |
| PR K | Ch. 11 Quantum Dynamics | TDPT, emission/absorption | 4–6 |
| PR L | Ch. 12 Afterword | EPR, Bell, measurement | 2–3 |

Total: 50–80 problems across ~12 PRs.

---

## 5. Acceptance criteria (for closing issue #49)

Per the issue body:

1. **PR A merged** — `Quantum_Mechanics/Griffiths/` scaffolded + Ch. 1 complete (3 problems).
2. **PR D merged** — hydrogen via proper-time formulation under canonical `K`; matches DRQM I §III recovery of non-rel limit.
3. **PR G merged** — fine structure of hydrogen reproduced under proper-time; comparison against Pauli/FW recorded.
4. **Root `README.md` updated** to list this campaign alongside Electromagnetism / Verification / History.

---

## 6. Crocco compliance + voice

Inherits the discipline from #42:

- **Pragmatic**: Mathematica MCP, 2e↔3e reconciliation, paraphrasing Griffiths derivations.
- **Substantive**: physical interpretation of any proper-time vs textbook divergence; reduction-to-textbook narrative in Ch. 4 / Ch. 7. Per-paragraph `<!-- TODO -->` for headline content.
- **Bibliography stubs**: Griffiths 2e and 3e via `scaffold_bib_note.py` before any cite. `pdf_status: acquired` (in-copyright).

**Voice**: prose follows the same [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) discipline that #42 established. Gill's voice is the campaign-wide reference for substantive QM prose; Griffiths' voice is the source-text voice and is preserved only in paraphrased quotations.

---

## 7. Honest framing

Per §13 of the #42 plan, all dual-theory predictions are **conditional** on the Gill–Zachary framework being correct. The QM campaign inherits this discipline:

- For non-relativistic Griffiths problems (~80% of the textbook), the proper-time formulation reduces to textbook QM at leading order. The campaign demonstrates the reduction; it does not validate the framework.
- For fine-structure and relativistic-correction problems (Ch. 7, parts of Ch. 4), the proper-time prediction is operationally distinguishable from Pauli/FW by `O((v/c)^4)` terms. Whether the framework gives the correct numerical answer is a separate question.
- DRQM I §III.D records a 🔴 numerical disagreement on the anomalous `g`-factor (`r_e ≈ 0.499857...` gives `g = -2.0005714` vs measured `-2.00231930...`). This is a known flagged finding; Ch. 7 problems that touch it carry a `⚠ blocked on r_e finding` flag analogous to #42's Eq. 24 branched-treatment workflow.

---

## 8. Decision points

Inherited from [§13.5 of #42](42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24) unless overridden:

| Inherited | Override for #49 | Status |
|---|---|---|
| Mechanical deferral of Eq. 24-touching problems | Mechanical deferral of `r_e`-finding-touching problems | confirmed |
| Per-paragraph TODO for headline content | same | confirmed |
| Always full template | same | confirmed |
| Per acceptance criterion: PR A + PR D + PR G | per issue body §"Acceptance criteria" | confirmed |
| No kill switch | same | confirmed |
| Voice = Gill's published-paper voice | same | confirmed |

---

## 9. Out-of-scope (deferred to Bethe–Salpeter umbrella)

The experimental-discrimination side of the dual-theory QM thread — precision tests of hydrogen fine structure, hyperfine splitting, Lamb shift, and `g-2` — is the subject of the **Bethe–Salpeter precision-predictions umbrella** opened concurrently with #49. That umbrella is the experimental-comparison companion to this campaign, mirroring how #43 and #48 relate to #42.
