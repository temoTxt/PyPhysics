# Implementation plan: Electromagnetism research thread — Jackson canonical problems × CGS/SI × proper-time Maxwell

**Tracks:** issue [#42](https://github.com/temoTxt/PyPhysics/issues/42).
**Status:** plan only. No content scaffolded yet. Execution should be staged across multiple focused PRs (see §10).
**Anchor source for the proper-time side:** [`Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) — Eqs. (1)–(23) verified via Wolfram MCP, Eq. (24) flagged for author review.

---

## 1. Goal and scope

Build a **third research campaign** (alongside Equation Verification and the History project) that works the canonical problems from Jackson, *Classical Electrodynamics* (2nd ed., 1975, Gaussian units; 3rd ed., 1998, mixed SI / Gaussian) **in three unit systems**:

1. **CGS / Gaussian** — the system the 2nd ed. uses throughout, and the system the 3rd ed. reverts to in Chs. 11+ (relativity, radiation).
2. **SI / MKS** — the system the 3rd ed. uses for Chs. 1–10.
3. **Proper-time (dual-theory) reformulation** — the same problem rewritten under the substitution rules established in the Gill–Zachary Maxwell paper (Eqs. 1–3′): `c → b`, `t → τ`, `w → u`, with the additional `−(u·a)/b⁴` dissipative term appearing whenever a τ-derivative meets `b`'s τ-dependence.

For each problem the output is one document showing (a) the classical CGS solution, (b) the SI translation, and (c) the proper-time reformulation, with explicit identification of when the dual theory differs from the classical answer by more than a `b/c` redressing of constants.

### Why this campaign

Three reasons it earns its own thread rather than being absorbed into Equation Verification:

- **The Gill Maxwell paper proves only the substitution rules are consistent.** It does not work out *what those rules predict* for the standard problems a graduate student already knows the answer to. This campaign closes that gap.
- **CGS/SI tripping is a known student stumbling block.** A canonical-problem set worked in both systems is its own pedagogical contribution, independent of the dual theory.
- **The dual theory's predictions for radiation-reaction problems (Liénard–Wiechert with the dissipative `(u·a)/b⁴` term) are the most experimentally distinguishable**, and Jackson Chs. 14–16 are the natural venue for that comparison. Putting them on the same template as Ch. 1–10 problems makes the divergence visible.

### Explicit non-goals

- **No exhaustive coverage.** Selected canonical problems (≈5–10 per chapter, ≈50–100 total). Not every numbered problem in either edition.
- **No reproduction of Jackson's problem statements.** Statements are paraphrased; the work product is the *solution*, plus a citation to "Jackson 3e §6.1 Problem 6.1" (etc.) so a reader with the textbook can locate the original. See §6.
- **No new physics claims about the dual theory.** Where the proper-time reformulation predicts something measurably different from classical EM, flag it; do not interpret it. Interpretation is author work (Trey + co-authors of DRQM I).
- **No edition-diffing as a primary deliverable.** The two editions are used as *sources* of canonical problems; we do not maintain a renumbering crosswalk.

---

## 2. Repository structure

New top-level folder, sibling to `Roadmapping/`:

```
Electromagnetism/
├── README.md                          # campaign overview + status table (mirrors Equation_Verification/README.md)
├── _template_problem.md               # per-problem template (see §4)
├── _proper_time_cheatsheet.md         # one-page substitution-rule reference (extracted from Gill Maxwell paper)
├── Jackson/
│   ├── README.md                      # Jackson-specific status table across all 14 chapters
│   ├── Ch01_Introduction_Electrostatics.md
│   ├── Ch02_Boundary_Value_Problems_Electrostatics_I.md
│   ├── …
│   ├── Ch06_Maxwell_Equations_Macroscopic_Media.md   ⭐ primary proper-time pivot
│   ├── …
│   ├── Ch11_Special_Relativity.md     ⭐ velocity-duality direct application
│   ├── Ch12_Relativistic_Dynamics.md  ⭐ proper-time Hamiltonian application
│   ├── Ch14_Radiation_by_Moving_Charges.md   ⭐ dissipative-term test case
│   └── Ch16_Radiation_Damping.md      ⭐ radiation-reaction test case
└── (future) Griffiths/                # placeholder for if/when other textbooks get added
```

Rationale for the `Jackson/` subfolder rather than `Electromagnetism/ChNN_*.md` at the top level:

- Keeps the door open to adding Griffiths or Zangwill problems later without a structural breaking change.
- Makes the Jackson-specific status table the natural place to track per-chapter progress without polluting the top-level README.
- Mirrors the existing `History/{Forward,Bibliography,Podcast}/` subfolder pattern.

Files that will be touched in the campaign but live elsewhere:

- `Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonChNN_PNN.wl` — Wolfram Mathematica notebooks per problem, runnable independent of the MCP (mirrors the existing convention for the verification campaign).
- `Roadmapping/Animations/manim_scenes/jackson_chNN_pNN_*.py` — optional Manim scenes for problems where the proper-time vs classical answer differs visibly (deferred to a later PR).
- Root `README.md` — add the new campaign as a fourth bullet under "Three pillars" → "Four pillars".

---

## 3. Per-problem template (proposed)

Lives at `Electromagnetism/_template_problem.md` and is referenced by every `ChNN_*.md` file. Draft below; the actual template will be refined after PR A (Ch. 6) is reviewed.

````markdown
### Problem JNe-PN.M — short title

**Source:** Jackson, *Classical Electrodynamics*, 3e §6.1 Problem 6.1 (and 2e §6.1 Problem 6.1, equivalent). _Paraphrased; consult the textbook for the precise statement._

**Paraphrased statement:** {1–3 sentences. Cite specific equation numbers from Jackson where relevant.}

**Setup:** {Geometry, boundary conditions, given quantities, unknowns. Diagram if needed.}

#### (a) Classical solution — Gaussian (CGS)

{Derivation, with Mathematica MCP check inline using the same format as the Equation_Verification docs.}

#### (b) Classical solution — SI

{Either the same derivation in SI, or — if the only difference is constants — an explicit conversion table:

| Quantity | Gaussian → SI |
|---|---|
| `E`-field | `E_SI = E_CGS / (4πε₀)^(1/2)` |
| etc. |
}

#### (c) Proper-time reformulation

Apply the substitution rules from [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (3′)]]:

- `c → b = √(c² + u²)`
- `t → τ`, `∂_t → (b/c) ∂_τ` (per Eq. 2)
- `w → u`, `J → (b/c) J` (per Eq. 12)
- Whenever a τ-derivative acts through `1/b`: add `−(u·a)/b³` from `∂_τ(1/b)`.

{Derivation in the proper-time frame. Mathematica MCP check.}

**Comparison:**

| Quantity | Classical (CGS) | Classical (SI) | Proper-time |
|---|---|---|---|
| … | … | … | … |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no / ⚠ yes — {one-line reason}.

**Verdict:** ✅ all three solutions consistent / ⚠ proper-time deviates (see below) / ❌ inconsistency found.

**Notes for author review (if ⚠):** {description suitable for inclusion in `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md`.}
````

Three template fields to debate before locking it (recorded as open questions in §9):

- Should `(a)` and `(b)` be merged when they differ only by constants? (Reduces redundancy but loses pedagogical CGS↔SI training.)
- Where should the Mathematica MCP block live — inline in `(c)` or in a sibling `.wl` file with a one-line "Verified ✅" link? (Inline matches Equation_Verification; sibling files match Mathematica_Notebooks.)
- Should the "Comparison" table be required (forces explicit side-by-side) or optional (only when interesting)?

---

## 4. Unit-system handling

The two editions of Jackson disagree about units in a way that matters for this campaign. Policy:

| Edition | Unit system used | Chapters affected |
|---|---|---|
| 2nd ed. (1975) | Gaussian (CGS) throughout | all (Chs. 1–16) |
| 3rd ed. (1998) | SI for Chs. 1–10, Gaussian for Chs. 11+ | mixed |

For each problem, the **canonical solution language** is determined by which edition's numbering we cite:

- If the problem appears in both editions and the numbering matches → present both CGS and SI explicitly.
- If the problem is in 3rd ed. only and falls in Chs. 1–10 → present in SI, then translate to CGS.
- If the problem is in 3rd ed. only and falls in Chs. 11+ → present in CGS (matches both 2nd ed. and 3rd ed. for that range).
- If the problem is in 2nd ed. only → present in CGS.

The proper-time reformulation always uses **Gaussian units** to match the Gill–Zachary Maxwell paper's conventions. SI translations of proper-time equations are out of scope for this campaign (would require a separate dual-theory SI conversion document; defer).

---

## 5. Proper-time mapping methodology

Extracted from the verified [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] paper. The substitution rules form a **finite list**; every problem reformulation reduces to applying them. The cheat sheet at `Electromagnetism/_proper_time_cheatsheet.md` formalises this.

| Rule | Origin | Application |
|---|---|---|
| `w/c = u/b` | Eq. (1) | Velocity in any kinematic expression |
| `(1/c) ∂_t = (1/b) ∂_τ` | Eq. (2) | Time-derivative in any field equation |
| `∇·E = 4πρ`, `∇·B = 0` | Eq. (3′) | Unchanged from standard Maxwell |
| `∇×E = −(1/b) ∂_τ B` | Eq. (3′) | Faraday with `c → b`, `t → τ` |
| `∇×B = (1/b) ∂_τ E + (4π/b) ρu` | Eq. (3′) | Ampère with `c → b`, `t → τ`, `J → (b/c) J = ρu` |
| `∂_τ(1/b) = −(u·a)/b³` | derived in Eq. (4) section | Whenever a τ-derivative meets `1/b` |
| `b' = γ[b − u·v/c]` | Eq. (11) | Boost of `b` between frames |
| `J' = ` ... | Eq. (12) | Boost of current density |
| Force law: `F = e[E + (u/b)×B] + (V/(mcb)) ...` | Eq. (18) | Lorentz force gets a correction term |

When applying to a Jackson problem, the **decision tree** is:

1. Identify Jackson's classical equation(s) used.
2. Check whether it is a Maxwell equation, a force law, or a derived quantity (potential, energy, momentum, etc.).
3. Apply the substitution rules above to the underlying Maxwell / force law, then re-derive the derived quantity from scratch in the proper-time frame (do NOT just substitute `c → b` into the classical *answer* — derived quantities pick up extra terms from `b(τ)`'s τ-dependence, as Eq. (4) demonstrates).
4. Where the proper-time answer differs from the `c → b` redressing of the classical answer, document the extra term and its physical interpretation (radiation reaction, etc.).

---

## 6. Source-material treatment (copyright + Crocco compliance)

Jackson is in copyright. Treatment:

- **Problem statements are paraphrased**, never reproduced verbatim. Each paraphrase is followed by an unambiguous citation (`Jackson 3e §6.1 Problem 6.1`) so a reader with the textbook can verify.
- **Equations cited from Jackson by number** (e.g. "Jackson Eq. 6.11") are reproduced only when they are short standard expressions a reader would recognise; longer derivations are reproduced from first principles rather than cited verbatim.
- **No PDF of Jackson is committed** to the repo. If a marker-pdf conversion exists locally for authoring convenience, it is gitignored.
- **Per Crocco §2:** every cite of "Jackson 3e" must correspond to an entry in [`Roadmapping/History/Bibliography/Retrospective/jackson1998_classical_electrodynamics.md`](../../Roadmapping/History/Bibliography/Retrospective/) — scaffold the stub via `Roadmapping/History/Bibliography/scaffold_bib_note.py` before the first PR lands.

The campaign's AI-use classification per Crocco §1:

- **Pragmatic:** running Mathematica MCP for symbolic verification, formatting the per-problem template, translating between CGS and SI, paraphrasing Jackson problem statements. _All of this is logged as pragmatic in the per-problem doc; no separate disclosure section required._
- **Substantive:** any prose where Claude proposes a *physical interpretation* of why the proper-time reformulation gives a different answer than classical EM. This is rare (most differences are mechanical — `c → b` plus the dissipative term) but when it happens the per-problem doc must include a `<!-- TODO: human reviews and fills in -->` section per Crocco §5.

---

## 7. Initial chapter selection + canonical-problems list

**First PR (PR A) — Ch. 6 Maxwell equations + macroscopic media.** Why Ch. 6 first:

- It is the natural proper-time pivot — the dual theory's Eq. (3′) is the proper-time version of Jackson's Eqs. (6.6)–(6.9).
- Both editions number Ch. 6 consistently and overlap in problem set, so the CGS/SI dual coverage is most natural here.
- It is the smallest chapter where the three-way comparison (CGS / SI / proper-time) is genuinely interesting on every problem.

**Canonical problems proposed for PR A** (5 problems, all in both editions):

| Jackson 3e # | Jackson 2e # | Short title | Why canonical |
|---|---|---|---|
| 6.1 | 6.1 | Maxwell equations with magnetic monopoles | Tests symmetry under electric ↔ magnetic exchange — does it survive `c → b`? |
| 6.4 | 6.4 | EM momentum of point charge moving uniformly | Direct application of the velocity-duality and `(u/b)×B` force correction |
| 6.5 | 6.5 | Poynting theorem in macroscopic media | Energy-conservation accounting in the proper-time frame |
| 6.11 | 6.11 | Symmetric stress tensor | Lorentz-covariant object — tests the dual theory's Lorentz behaviour |
| 6.20 | 6.20 | Radiation pressure on perfect conductor | Macroscopic boundary-condition problem; tests whether the `(u·a)/b⁴` term contributes |

**Subsequent PRs** prioritise the chapters where the proper-time vs classical difference is most experimentally distinguishable:

- **PR B — Ch. 11 Special Relativity** (5–7 problems). Direct application of the `w/c = u/b` velocity duality. Likely to be the *least* divergent — proper-time SR should agree with classical SR on every kinematic prediction.
- **PR C — Ch. 12 Relativistic Dynamics** (5–7 problems). Direct application of the proper-time Hamiltonian `K = H²/(2mc²) + mc²/2` from DRQM I.
- **PR D — Ch. 14 Radiation by Moving Charges** (5–7 problems). First chapter where the dissipative `(u·a)/b⁴` term *should* contribute non-trivially. Liénard–Wiechert problems are the pivot.
- **PR E — Ch. 16 Radiation Damping** (3–5 problems). The Abraham–Lorentz problem in the dual framework — likely the **headline payoff** of this campaign, since classical EM's runaway/pre-acceleration pathologies are what motivate the dual theory's modification in the first place.
- **PR F+ — Chs. 1–5, 7–10, 13, 15** (4–6 problems each). Backfill of chapters where the proper-time reformulation is structurally identical to classical and exists for pedagogical CGS/SI completeness.

Total scope: ~50–80 problems across ~12 PRs spanning 6–12 months.

---

## 8. Definition of done

Per problem:
- ✅ All three sections (CGS, SI, proper-time) complete with derivations.
- ✅ Mathematica MCP check recorded inline with `Result: …  ✅` annotation (matches Equation_Verification convention).
- ✅ Companion `.wl` notebook committed under `Roadmapping/Mathematica_Notebooks/Electromagnetism/`.
- ✅ Comparison table populated.
- ✅ Verdict assigned (`✅ / ⚠ / ❌`).
- ✅ Citation to Jackson edition + section + problem number.

Per chapter:
- ✅ `ChNN_*.md` table-of-contents updated.
- ✅ `Electromagnetism/Jackson/README.md` status table updated.

Per PR:
- ✅ Top-level `Electromagnetism/README.md` campaign status updated.
- ✅ Root `README.md` "Four pillars" section reflects new campaign progress (after PR A only — not every subsequent PR).
- ✅ Any ⚠ verdicts cross-posted to `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` if they reflect on the Gill papers.

---

## 9. Open questions / decisions deferred to the human

1. **Template debate (§3 bullets).** Three template choices flagged — best resolved after one chapter is fully populated and the trade-offs are concrete.
2. **Mathematica MCP licensing.** Per `Roadmapping/Equation_Verification/README.md`, MCP went online 2026-05-10. Confirmed available for this campaign. No action needed.
3. **Animation companion.** Should each Manim scene for this campaign live under `Roadmapping/Animations/manim_scenes/jackson_*.py`, or get its own subdir `Roadmapping/Animations/manim_scenes/Electromagnetism/`? Defer to first scene's PR.
4. **Cross-referencing to History campaign.** The History chapters on classical EM (Ch. 2: Classical synthesis 1860–1900) cite Maxwell's 1865 paper; should Jackson-problem docs back-reference those historical chapters where they introduce a concept? Recommended: yes, but only for the load-bearing problems (e.g. Ch. 6 displacement current → History Ch. 2).
5. **Future textbooks.** The `Electromagnetism/Jackson/` structure leaves room for `Electromagnetism/Griffiths/` etc. Should the campaign commit to adding Griffiths problems after Jackson is at, say, 50% coverage? Defer to mid-campaign.
6. **Should the canonical-problems list (§7) be locked at plan time, or evolved per-PR?** Recommendation: lock the PR A list (5 problems above); evolve PRs B+ as we learn what's most pedagogically valuable.

---

## 10. PR sequencing

| PR | Scope | Estimated effort | Blocks |
|---|---|---|---|
| **PR A — Scaffold + Ch. 6** | Top-level `Electromagnetism/` dir, README, template, cheat sheet, `Jackson/README.md`, `Ch06_Maxwell_Equations_Macroscopic_Media.md` with 5 problems fully solved. Bibliography stub for Jackson 3e. Root README update. | ~2 weeks of focused work | nothing |
| **PR B — Ch. 11 Special Relativity** | 5–7 problems. Adds animation companion if any problem yields a striking visual. | ~1.5 weeks | PR A |
| **PR C — Ch. 12 Relativistic Dynamics** | 5–7 problems. Cross-refs DRQM I. | ~1.5 weeks | PR A |
| **PR D — Ch. 14 Radiation by Moving Charges** | 5–7 problems. *First chapter where dissipative term contributes.* | ~2 weeks | PR A |
| **PR E — Ch. 16 Radiation Damping** | 3–5 problems. **Headline payoff.** Likely yields one ⚠ for author review. | ~2 weeks | PR D |
| **PR F+** | Backfill Chs. 1–5, 7–10, 13, 15 in priority order. | 1 week each, ≤2 chapters per PR |  |

Total: 12–14 PRs over 6–12 months, depending on cadence.

---

## 11. Acceptance criteria (for closing issue #42)

The issue closes when:
1. PR A is merged (campaign scaffolded, Ch. 6 complete).
2. Either PR D *or* PR E is merged (at least one chapter where the proper-time reformulation contributes a non-trivial new term has been worked).
3. Root `README.md` reflects "Four pillars" with this campaign listed.

Subsequent chapter work continues under PR labels but does not require keeping the issue open.
