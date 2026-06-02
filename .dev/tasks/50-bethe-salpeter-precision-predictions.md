# Implementation plan: Bethe–Salpeter precision predictions under proper-time / dual-theory formulation

**Tracks:** issue [#50](https://github.com/temoTxt/PyPhysics/issues/50).
**Status:** plan + scaffold; PRs A–J to follow.
**Sibling campaigns:** [#42 — Electromagnetism / Jackson](42-electromagnetism-jackson-proper-time.md) (PR #51); [#49 — Quantum Mechanics / Griffiths](49-quantum-mechanics-griffiths-proper-time.md) (PR #52, pedagogical prerequisite for this thread).
**Anchor sources:** verified [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) (Eqs. I.6 canonical `K`, II.1–II.3 dual Dirac ✅), [`Analytic_Representation_of_The_Dirac_Equation.md`](../../Roadmapping/Equation_Verification/Analytic_Representation_of_The_Dirac_Equation.md), [`The_Classical_Electron_Problem.md`](../../Roadmapping/Equation_Verification/The_Classical_Electron_Problem.md) (radiative-correction route for Lamb shift), and the proper-time Liénard–Wiechert third term from [`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) Eq. (7).

This plan is **deliberately concise** because it inherits its discipline from [§13 of #42](42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix) and [§7 of #49](49-quantum-mechanics-griffiths-proper-time.md#7-honest-framing). Differences from the sibling campaigns are flagged below.

---

## 1. Goal and scope

Rederive the **canonical precision results in Bethe & Salpeter, *Quantum Mechanics of One- and Two-Electron Atoms*** (Plenum 1957; Springer 1977 reprint) under the proper-time / dual-theory formulation, and compare against current measurements.

Per-result output: one entry showing (a) Bethe–Salpeter's textbook derivation, (b) the proper-time / dual-theory rederivation under `K` and/or the dual Dirac equation, (c) numerical comparison against the modern measured value (CODATA-2018 / PDG-2024 wherever applicable), and (d) an explicit verdict (✅ / ⚠ / 🔴).

### Why this earns its own thread

- **Experimental discrimination.** Where #49 demonstrates pedagogical correspondence (proper-time `K` → `p²/2m` for non-relativistic textbook), Bethe–Salpeter is the textbook where the proper-time formulation has *measurable* consequences: Lamb shift (~MHz, 9-sig-fig measured), hyperfine splitting (~12-sig-fig measured), fine-structure splittings.
- **Headline payoffs at Lamb shift (PR E) and hyperfine (PR F).** Hydrogen 21-cm line is `1420.405751…` MHz; Lamb shift is `1057.845(9)` MHz. The campaign's verdicts on these two numbers determine whether the dual-theory framework is consistent with precision atomic spectroscopy.
- **Most to lose, most to gain.** This is the campaign where the framework can be *falsified* at the current precision floor — and conversely, where agreement to measurement precision is the strongest experimental endorsement available.

### Explicit non-goals

- No exhaustive coverage. Bethe–Salpeter is encyclopedic (~400 sections in the Springer reprint); we triage to the precision-comparable subset (~30–50 results, ~3–6 per major section block).
- No reproduction of Bethe–Salpeter's derivations verbatim. Fair-use quoted equations + section refs only.
- No new physics claims beyond what DRQM I + companion verification docs support. Conditional framing throughout.
- No QED loop calculations from scratch. The Lamb-shift treatment uses the *route* supplied by `The_Classical_Electron_Problem` (radiation reaction → self-energy contribution), not a full one-loop calculation.

---

## 2. Repository structure

```
Roadmapping/Quantum_Mechanics/
├── README.md                              # updated with Bethe–Salpeter campaign row
├── _template_problem.md                   # Griffiths template (inherited)
├── _proper_time_K_cheatsheet.md           # canonical K and dual Dirac (inherited)
├── Griffiths/                             # #49 campaign (pedagogical-correspondence)
└── Bethe_Salpeter/                        # ⭐ this campaign (precision-experiment)
    ├── README.md                          # campaign status table across PRs A–J
    ├── _template_result.md                # per-result template (precision-comparison variant)
    ├── 01_NonRelHydrogen.md               # PR A — §§1–7
    ├── 02_MatrixElements.md               # PR B — §§8–13
    ├── 03_FineStructure.md                # PR C — §14 ⭐ pivot
    ├── 04_HigherOrderRel.md               # PR D — §§15–18
    ├── 05_LambShift.md                    # PR E — §§19–21 ⭐ headline
    ├── 06_Hyperfine.md                    # PR F — §22 ⭐ headline
    ├── 07_RadiationInteraction.md         # PR G — §§23–37
    ├── 08_HeliumGround.md                 # PR H — §§47–60
    ├── 09_HeliumExcited.md                # PR I — §§61–80
    └── 10_CrossComparison.md              # PR J — summary
```

Bethe–Salpeter sits as a **sibling to Griffiths** under `Quantum_Mechanics/`, mirroring how `Jackson/` and `Griffiths/` are siblings in their respective domain folders.

---

## 3. Per-result template

Lives at `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/_template_result.md` (drafted in PR A). Differs from Griffiths' template in three places:

- **One source citation** (single canonical edition — Springer 1977 reprint, identical section/equation numbering to Plenum 1957).
- **Three result blocks**: (a) Bethe–Salpeter textbook (with page + equation ref), (b) proper-time / dual-theory derivation, (c) **modern measurement** (CODATA / PDG / experimental paper).
- **Three-way numerical comparison** column (proper-time prediction vs textbook prediction vs measured value). Residual in MHz / Hz / ppm.

The third column — modern measurement — is the campaign's defining feature; it does not appear in #42 or #49.

---

## 4. PR sequencing (10 PRs A–J per issue body)

Order = Bethe–Salpeter's pedagogical order, weighted toward precision-comparable results:

| PR | Bethe–Salpeter sections | Role | Results |
|---|---|---|---|
| **PR A** | §§1–7 — Non-rel hydrogen | scaffold + non-rel reduction `K → p²/2m + V₀` | 3–4 |
| PR B | §§8–13 — Matrix elements + transitions | dipole, multipole, oscillator strengths | 3–4 |
| **PR C** ⭐ | **§14 Fine structure (pivot)** | Dirac vs dual Dirac for H 2P-3/2 – 2P-1/2 | 3–4 |
| PR D | §§15–18 — Higher-order rel corrections | `(Zα)⁴` terms; FW vs proper-time | 3–4 |
| **PR E** ⭐ | **§§19–21 Lamb shift (headline)** | self-energy via `K` + radiation-reaction route | 3–4 |
| **PR F** ⭐ | **§22 Hyperfine structure (headline)** | H 21-cm line; muonium/positronium deferred to PR I | 2–3 |
| PR G | §§23–37 — Interaction with radiation | photoionisation, multipole; dipole-approx & third term | 3–4 |
| PR H | §§47–60 — Helium ground state | two-electron variational; proper-time energy operator | 3–4 |
| PR I | §§61–80 — Helium excited states | two-electron spectroscopy; positronium / muonium where applicable | 3–4 |
| PR J | Cross-comparison summary | table of all proper-time vs measured; flagged shifts | (chapter) |

Total: ~30–50 results across 10 PRs.

PRs C, E, F are the **headline payoffs**; the campaign's acceptance criteria gate on them.

---

## 5. Acceptance criteria (for closing issue #50)

Per the issue body:

1. **PR A merged** — `Bethe_Salpeter/` scaffolded + §§1–7 non-rel hydrogen complete (3–4 results), with the reduction from `K` to `p²/2m + V₀` verified via Wolfram MCP.
2. **PR C merged** — §14 fine structure of hydrogen rederived under the proper-time formulation; comparison against measured `2P₃/₂ – 2P₁/₂` splitting recorded.
3. **PR E merged** — §§19–21 Lamb shift treatment under the proper-time formulation; explicit numerical prediction vs `1057.845(9)` MHz measured, with verdict.
4. At least one ⚠ or 🔴 finding (if any survives the analysis) cross-posted to [`Roadmapping/Equation_Verification/FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md).

Subsequent PRs (D, F, G, …) continue under this campaign but do not require keeping #50 open.

---

## 6. Crocco compliance + voice

Inherits the discipline from #42 and #49:

- **Pragmatic**: Mathematica MCP verification of non-rel reductions and algebraic identities; transcription of Bethe–Salpeter equations (with page refs); CODATA / PDG lookups.
- **Substantive**: any prose that interprets a proper-time vs measurement divergence (esp. PR C, PR E, PR F). Per-paragraph `<!-- TODO: human reviews and fills in -->` blocks for headline content.
- **Bibliography stub**: `bethe1977_one_two_electron_atoms` via [`scaffold_bib_note.py`](../../Roadmapping/History/Bibliography/scaffold_bib_note.py) before any cite. `pdf_status: acquired` (in-copyright Springer reprint).

**Voice**: prose follows the [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) discipline. The "we plural / connective openers / mathematically equivalent but not physically equivalent" anchor is load-bearing for the Lamb-shift and hyperfine prose where physical interpretation matters most.

---

## 7. Honest framing

The campaign's three honest caveats:

### 7.1 The Lamb-shift route is *not* a one-loop QED calculation

Bethe–Salpeter §§19–21 derive the Lamb shift via Bethe's mass-renormalisation argument (1947), which approximates the QED self-energy. The proper-time treatment we offer follows the **same route**, substituting the dual-Dirac propagator and the proper-time radiation-reaction structure from [`The_Classical_Electron_Problem`](../../Roadmapping/Equation_Verification/The_Classical_Electron_Problem.md) for the standard Bethe estimate. This reproduces the leading log-Bethe structure; a full one-loop dual-QED computation is **out of scope** for this campaign.

What this means for the verdict: agreement at the Bethe-estimate precision (~few percent of the full Lamb shift) is what the campaign can honestly deliver. Agreement to the measured `1057.845(9)` MHz at full precision would require the one-loop calculation, which the dual-theory framework has not yet produced.

### 7.2 Hyperfine splitting depends on the anomalous-`g`-factor

The hyperfine constant `A` depends on the electron's magnetic moment `μ_e = g_e μ_B/2`. DRQM I §III.D's flagged finding on `r_e` propagates directly into the hyperfine prediction. PR F therefore carries a **branched treatment** identical to #49 PR G:

- **Branch (a)** — as-published `r_e ≈ 0.499857`: gives `g = -2.0005714` and a corresponding hyperfine prediction that disagrees with the measured 21-cm line at parts-per-thousand level.
- **Branch (b)** — corrected `r_e ≈ 0.499420`: gives the measured `g = -2.00231930…` and recovers the textbook hyperfine prediction.

Neither branch is the campaign's "verdict"; both are recorded. The flagged finding is unchanged.

### 7.3 Fine structure has the same anomalous-`g` dependence

PR C's fine-structure splitting (especially the magnetic dipole term) inherits the same `r_e` branching. The campaign reports both branches; the experimental verdict is conditional on whichever value of `r_e` the dual-theory framework ultimately adopts.

These three caveats are explicit in PRs C, E, F and surface in the closing PR J cross-comparison table.

---

## 8. Decision points

Inherited from [§13.5 of #42](42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24) and [§8 of #49](49-quantum-mechanics-griffiths-proper-time.md#8-decision-points) unless overridden:

| Inherited | Override for #50 | Status |
|---|---|---|
| Per-paragraph TODO for headline content | same — applies to PRs C, E, F especially | confirmed |
| Always full template | same | confirmed |
| Per-issue acceptance criteria (PR A + PR C + PR E + flagged-finding crosspost) | confirmed | confirmed |
| Mechanical deferral of `r_e`-finding-touching problems | **branched-treatment instead of deferral** for PR C / PR F (precision results cannot be deferred) | confirmed |
| Voice = Gill's published-paper voice | same | confirmed |
| No kill switch | same — full A–J delivery with minimal stoppage | confirmed |

---

## 9. Out-of-scope

- **Two-loop QED corrections** (Lamb shift at sub-MHz precision). The current campaign is a tree-level / one-loop-estimate comparison; sub-MHz precision is deferred to a future thread once the dual-theory framework supplies a renormalised loop calculation.
- **Positronium and muonium spectroscopy** at full precision. PR I touches these where the two-electron Bethe–Salpeter machinery applies; full positronium spectroscopy (analogous to but distinct from H) is its own future thread.
- **Heavy-element precision spectroscopy** (Helium-like ions, Z > 2). Bethe–Salpeter's helium chapters are the upper limit of this campaign's scope.
- **The DRQM I §III.D `r_e` issue itself** is not relitigated. The branched-treatment workflow records both branches; resolving which branch is correct is a separate research question (potentially the subject of a future companion thread).

---

## 10. Linked PRs / branches

- Parent issue: [#50](https://github.com/temoTxt/PyPhysics/issues/50)
- Branch: `50-bethe-salpeter-precision-predictions` (stacks on `49-quantum-mechanics-griffiths-proper-time`)
- PR (to be opened after PR A scaffold): closes #50

---

## 11. Working notes

- Bethe–Salpeter PDF is local-only (in-copyright). All equation references in committed markdown are fair-use quotations + section/equation labels.
- The Springer 1977 reprint preserves the Plenum 1957 numbering; references work for both editions.
- The "load-bearing observation" of the campaign — articulated across all per-result documents — is that the proper-time framework's predictions for precision atomic spectroscopy are (a) computable, (b) numerically distinguishable from textbook QED only at sub-leading orders, and (c) consistent with measurement at the precision the framework can presently deliver. PR J's cross-comparison summary is where this is recorded.
