# Implementation plan: Electromagnetism research thread — Jackson canonical problems × CGS/SI × proper-time Maxwell

**Tracks:** issue [#42](https://github.com/temoTxt/PyPhysics/issues/42).
**Status:** plan only. No content scaffolded yet. Execution should be staged across multiple focused PRs (see §10).
**Anchor source for the proper-time side:** [`Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) — Eqs. (1)–(23) verified via Wolfram MCP, Eq. (24) flagged for author review.

---

## 1. Goal and scope

Build a **third research campaign** (alongside Equation Verification and the History project) that works the canonical problems from Jackson, *Classical Electrodynamics* (2nd ed., 1975, Gaussian units; 3rd ed., 1998, mixed SI / Gaussian) **in three unit systems**:

1. **CGS / Gaussian** — the system the 2nd ed. uses throughout, and the system the 3rd ed. reverts to in Chs. 11+ (relativity, radiation).
2. **SI / MKS** — the system the 3rd ed. uses for Chs. 1–10.
3. **Proper-time (dual-theory) reformulation** — the same problem rewritten under the substitution rules established in the Gill–Zachary Maxwell paper (Eqs. 1–3′): `c → b`, `t → τ`, `w → u`. Whenever a τ-derivative meets `b`'s τ-dependence, the chain rule introduces an extra `∂_τ(1/b) = −(u·a)/b³` factor; once this is collected as a coefficient of `∂_τ` in the wave equation (after multiplying through by `1/b`), it appears as the dissipative coefficient `−(u·a)/b⁴ · ∂_τ` of Maxwell-paper Eq. (4). The same `(u·a)/b⁴` factor reappears in the Liénard–Wiechert fields of Eq. (7). These are the same physical effect at two stages of one derivation, *not* two independent terms.

For each problem the output is one document showing (a) the classical CGS solution, (b) the SI translation, and (c) the proper-time reformulation, with explicit identification of when the dual theory differs from the classical answer by more than a `b/c` redressing of constants. (Per §4, sections (a) and (b) collapse to a single "Classical solution" section when only one classical unit system applies — true for all of Jackson 2e and for 3e Chs. 11+. Most of the headline PRs (B–E) are therefore two-system "CGS + proper-time" comparisons rather than three-way; the full three-way treatment is reserved for problems that appear in 3e Chs. 1–10.)

### Why this campaign

Three reasons it earns its own thread rather than being absorbed into Equation Verification:

- **The Gill Maxwell paper proves only the substitution rules are consistent.** It does not work out *what those rules predict* for the standard problems a graduate student already knows the answer to. This campaign closes that gap.
- **CGS/SI tripping is a known student stumbling block.** A canonical-problem set worked in both systems is its own pedagogical contribution, independent of the dual theory.
- **The dual theory's predictions for radiation-reaction problems (Liénard–Wiechert with the dissipative `(u·a)/b⁴` term) are the most experimentally distinguishable**, and Jackson Chs. 14–16 are the natural venue for that comparison. Putting them on the same template as Ch. 1–10 problems makes the divergence visible.

### Explicit non-goals

- **No exhaustive coverage.** Selected canonical problems (≈5–10 per chapter, ≈50–100 total). Not every numbered problem in either edition.
- **No reproduction of Jackson's problem statements.** Statements are paraphrased; the work product is the *solution*, plus a citation to "Jackson 3e §6.1 Problem 6.1" (etc.) so a reader with the textbook can locate the original. See §6.
- **No new physics claims about the dual theory.** Where the proper-time reformulation predicts something measurably different from classical EM, flag it as a *conditional* prediction ("*if* the Gill–Zachary framework is the correct formulation, *then* X") and document the size of the predicted deviation; do not interpret it as evidence for the framework. Interpretation is author work (Trey + co-authors of DRQM I). See [§13](#13-devils-advocate-review-and-mitigations) for how this constraint is enforced across the campaign.
- **No edition-diffing as a primary deliverable.** The two editions are used as *sources* of canonical problems; we do not maintain a renumbering crosswalk.

---

## 2. Repository structure

New research thread under `Roadmapping/`, alongside `Equation_Verification/`, `History/`, `Animations/`, and the other existing threads:

```
Roadmapping/
└── Electromagnetism/
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
    └── (future) <other EM textbook>/      # placeholder for a follow-on classical-EM textbook (e.g., Zangwill, Schwinger). NOT Griffiths — see note below.
```

**Note on Griffiths.** Griffiths is canonically associated with quantum mechanics (*Introduction to Quantum Mechanics*), not classical EM. If a Griffiths-textbook campaign is added later, it belongs in a separate `Roadmapping/Quantum/Griffiths/` thread, not here. Griffiths' *Introduction to Electrodynamics* is a legitimate possible follow-on under `Roadmapping/Electromagnetism/`, but the placeholder above is left generic to avoid prejudicing which EM textbook (if any) follows Jackson.

**Why under `Roadmapping/` rather than at the top level.** An earlier draft placed `Electromagnetism/` as a top-level folder sibling to `Roadmapping/`, framed as a "fourth pillar." That framing oversold the campaign — the work is structurally a research thread, not a parallel institutional structure to `Roadmapping/`. Living under `Roadmapping/` puts it alongside the other research threads (Equation Verification, History, Animations) where it naturally belongs, and removes the implicit claim that it deserves the same institutional weight as `Roadmapping/` itself. The "third research campaign" wording in [§1](#1-goal-and-scope) refers to it being the third thread under `Roadmapping/` after Equation Verification and History, not to a top-level pillar.

Rationale for the `Jackson/` subfolder rather than `Roadmapping/Electromagnetism/ChNN_*.md` directly:

- Keeps the door open to adding Griffiths or Zangwill problems later without a structural breaking change.
- Makes the Jackson-specific status table the natural place to track per-chapter progress without polluting the thread's top-level README.
- Mirrors the existing `Roadmapping/History/{Forward,Bibliography,Podcast}/` subfolder pattern.

Files that will be touched in the campaign but live elsewhere:

- `Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonChNN_PNN.wl` — Wolfram Mathematica notebooks per problem, runnable independent of the MCP (mirrors the existing convention for the verification campaign).
- `Roadmapping/Animations/manim_scenes/jackson_chNN_pNN_*.py` — optional Manim scenes for problems where the proper-time vs classical answer differs visibly (deferred to a later PR).
- `Roadmapping/README.md` (if one exists, or root `README.md`) — add a one-line reference to the new thread alongside the existing thread descriptions. **No "Four pillars" or equivalent institutional-weight framing** (see preceding paragraph and [§13.1 O7](#131-objections-with-partial-responses)).

---

## 3. Per-problem template (proposed)

Lives at `Roadmapping/Electromagnetism/_template_problem.md` and is referenced by every `ChNN_*.md` file. Draft below; the actual template will be refined after PR A (Ch. 6) is reviewed.

**Problem-ID convention.** Each problem heading uses `J{N}e-P{C}.{P}`, where `{N}e` is the edition (`2e` or `3e`), `{C}` is the chapter, and `{P}` the problem number — e.g. `J3e-P6.1` = Jackson 3rd ed., Chapter 6, Problem 6.1. When a problem appears in both editions with matching numbering, the 3e key is canonical and the 2e equivalence is noted in the source line.

````markdown
### Problem J3e-P6.1 — short title

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
- Whenever a τ-derivative acts through `1/b`, apply `∂_τ(1/b) = −(u·a)/b³`. When this is collected as a coefficient on `∂_τ` in a second-order wave equation (after multiplying through by `1/b`), the same term re-expresses as `−(u·a)/b⁴ · ∂_τ` — the dissipative coefficient of Eq. (4). The Liénard–Wiechert fields (Eq. 7) carry an `(u·a)/b⁴` factor in their third term for the same reason.

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

- If the problem appears in both editions and the numbering matches → present both CGS and SI explicitly (three-way template).
- If the problem is in 3rd ed. only and falls in Chs. 1–10 → present in SI, then translate to CGS (three-way template).
- If the problem is in 3rd ed. only and falls in Chs. 11+ → present in CGS only (matches both 2nd ed. and 3rd ed. for that range; two-way template).
- If the problem is in 2nd ed. only → present in CGS only (two-way template).

When only one classical unit system applies, the template's (a) and (b) sections merge into a single "Classical solution" section and the comparison table becomes two-column (CGS + proper-time).

The proper-time reformulation always uses **Gaussian units** to match the Gill–Zachary Maxwell paper's conventions. SI translations of proper-time equations are out of scope for this campaign (would require a separate dual-theory SI conversion document; defer).

---

## 5. Proper-time mapping methodology

Extracted from the verified [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] paper. The substitution rules form a **finite list**; every problem reformulation reduces to applying them. The cheat sheet at `Roadmapping/Electromagnetism/_proper_time_cheatsheet.md` formalises this.

| Rule | Origin | Application |
|---|---|---|
| `w/c = u/b` | Eq. (1) | Velocity in any kinematic expression |
| `(1/c) ∂_t = (1/b) ∂_τ` | Eq. (2) | Time-derivative in any field equation |
| `∇·E = 4πρ`, `∇·B = 0` | Eq. (3′) | Unchanged from standard Maxwell |
| `∇×E = −(1/b) ∂_τ B` | Eq. (3′) | Faraday with `c → b`, `t → τ` |
| `∇×B = (1/b) ∂_τ E + (4π/b) ρu` | Eq. (3′) | Ampère with `c → b`, `t → τ`, `J → (b/c) J = ρu` |
| `∂_τ(1/b) = −(u·a)/b³` | derived in Eq. (4) section | Bare chain rule, whenever a τ-derivative meets `1/b` |
| Dissipative coefficient `−(u·a)/b⁴ · ∂_τ` in wave equation | Eq. (4) | Same effect as above, after `× 1/b` to reach standard wave-equation form. Use `ḃ = (u·a)/b` to convert between the two forms freely |
| `b' = γ[b − u·v/c]` | Eq. (11) | Boost of `b` between frames |
| `J' = J + (γ−1)(J·v)v/v² − γ(b/c)ρv` | Eq. (12) | Spatial part of the 4-current `(bρ, J)` boost |
| Force law: `F = e[E + (u/b)×B] + (V/(mcb)) ...` | Eq. (18) | Lorentz force gets a correction term |

When applying to a Jackson problem, the **decision tree** is:

1. Identify Jackson's classical equation(s) used.
2. Check whether it is a Maxwell equation, a force law, or a derived quantity (potential, energy, momentum, etc.).
3. Apply the substitution rules above to the underlying Maxwell / force law, then re-derive the derived quantity from scratch in the proper-time frame (do NOT just substitute `c → b` into the classical *answer* — derived quantities pick up extra terms from `b(τ)`'s τ-dependence, as Eq. (4) demonstrates).
4. Where the proper-time answer differs from the `c → b` redressing of the classical answer, document the extra term and its physical interpretation (radiation reaction, etc.).

### 5.1 Branched treatment for Eq. 24-touching problems

Per [§13.5 D1](#135-decision-points--confirmed-by-author-2026-05-24), any Jackson problem whose proper-time derivation invokes Maxwell-paper Eq. 24 is worked twice: once with Eq. 24 as-published, once with the proposed correction from [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) (restored factor of `c`, plus the `V²/(2mc²)` term).

**Identifying the dependency.** A Jackson problem touches Eq. 24 if its proper-time derivation needs the Pauli-form Hamiltonian `H = (mc² + p²/2m − eϕ + (e/c)A·p + eℏΣ·B/(2m) + V²/(2mc²))` or its non-relativistic reduction. In practice this means:
- Any problem invoking the proper-time Hamiltonian formulation (likely in Ch. 12 dynamics).
- Any spin-orbit or hyperfine-structure problem.
- Any problem deriving an effective potential that includes a `V²/c²` correction.

Most pure-EM Jackson problems (Chs. 1–11, 14, 16) **do not** touch Eq. 24 because they operate at the field-equation level (Eqs. 3′, 4, 7), not the Hamiltonian level. The branched treatment is expected to engage primarily in Ch. 12 problems.

**Per-problem document structure when branched.** Under the proper-time `(c)` section, two parallel subsections appear:

```markdown
#### (c) Proper-time reformulation — as-published Eq. 24

{Derivation using Eq. 24 verbatim from the Gill–Zachary paper. MCP check. Result.}

#### (c′) Proper-time reformulation — with proposed Eq. 24 correction

{Same derivation with the restored c factor and the V²/(2mc²) term added. MCP check. Result.}

#### (c) vs (c′) comparison

{One-paragraph summary of how the two branches differ at the answer level. If they
agree on the answer (within numerical tolerance), say so. If they disagree, the
disagreement is itself a data point on which Eq. 24 form is consistent with the
rest of the framework.}
```

**Verdict semantics under branched treatment.** The per-problem verdict (✅ / ⚠ / ❌) is recorded *per branch*. If both branches yield ✅, the problem is consistent under either Eq. 24 form. If one branch yields ⚠ and the other ✅, the comparison itself is a piece of evidence on which form is the load-bearing one — flagged for the [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) entry on Eq. 24.

**Cost.** A branched problem takes ~1.7× the single-branch time (the second derivation is mostly mechanical given the first). Effort-estimate impact is folded into [§13.1 O4](#131-objections-with-partial-responses).

---

## 6. Source-material treatment (copyright + Crocco compliance)

Jackson is in copyright. Treatment:

- **Problem statements are paraphrased**, never reproduced verbatim. Each paraphrase is followed by an unambiguous citation (`Jackson 3e §6.1 Problem 6.1`) so a reader with the textbook can verify.
- **Equations cited from Jackson by number** (e.g. "Jackson Eq. 6.11") are reproduced only when they are short standard expressions a reader would recognise; longer derivations are reproduced from first principles rather than cited verbatim.
- **No PDF of Jackson is committed** to the repo. If a marker-pdf conversion exists locally for authoring convenience, it is gitignored.
- **Per Crocco §2:** every cite of either edition must correspond to an existing entry under `Roadmapping/History/Bibliography/Retrospective/` — `jackson1998_classical_electrodynamics.md` for the 3rd ed. and `jackson1975_classical_electrodynamics.md` for the 2nd ed. Both stubs must be scaffolded via `Roadmapping/History/Bibliography/scaffold_bib_note.py` before the first PR lands, since the campaign cites both editions from PR A onward.

The campaign's AI-use classification per Crocco §1 (refined in [§13.1 O3](#131-objections-that-materially-shaped-the-plan)):

- **Pragmatic:** running Mathematica MCP for symbolic verification, formatting the per-problem template, translating between CGS and SI, paraphrasing Jackson problem statements. _All of this is logged as pragmatic in the per-problem doc; no separate disclosure section required._
- **Substantive (broader than initially scoped):** any prose where Claude proposes a *physical interpretation*, narrative framing ("headline payoff", "smoking gun"), problem selection rationale, or judgement about which classical answer is the appropriate comparator. Per-problem documents include a `<!-- TODO: human reviews and fills in -->` block for each substantive paragraph and a top-of-document **Selection provenance** note recording which problem was chosen, why, and what alternatives were considered. The §7 / §12.3 "considered and dropped" lists are the campaign-level audit trail.

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

Total scope: ~50–100 problems across ~12–14 PRs spanning 6–12 months (confirmed by author).

### 7.1 Headline-payoff problems vs fluency-builder problems

Each PR's problem list (5–7 items) should be deliberately constructed from two categories, not from a single "best problems we can find" list. **Per [§13.5 D3](#135-decision-points--confirmed-by-author-2026-05-24), all problems in both categories receive the full template treatment** (three-section derivation + MCP verification + Wolfram notebook + Comparison table). The categories below describe *derivation difficulty and narrative value*, not template completeness.

- **Headline-payoff problems** are those where the proper-time-vs-classical contrast is rich, narrative-bearing, and worth elevating into prose or podcast material. The five §12 podcast-elevation picks are the canonical examples. These take the longest per problem (~3–5 days) because the derivation itself is rich and the interpretive write-up is long. They produce the campaign's load-bearing output.
- **Fluency-builder problems** are those where the proper-time-vs-classical contrast is mechanical and the derivation is short. Their purpose is *repetition* — building the author's pattern recognition for how the substitution rules apply across a range of problem geometries. These problems take ~2–3 days each because the derivation is faster, but they still get the full template (Comparison table filled in, MCP verification recorded, Wolfram notebook committed). The shorter timeline reflects shorter derivations, not lighter documentation.

**Why this matters.** The headline-payoff problems are where the framework "does interesting work" (`u·a ≠ 0`, new terms appearing). The fluency-builder problems are where the framework reduces to known classical answers with at most a `c → b` rescaling. Without the fluency builders, the author never develops the muscle memory to recognise when a proper-time derivation has gone off the rails — every problem becomes a fresh adventure, slow and error-prone. With the fluency builders, the author can spot mechanical errors quickly and reserve the long-form interpretive effort for the problems that genuinely warrant it.

**Mix per PR.** Recommended split:

| PR scope (5–7 problems) | Headline-payoff | Fluency-builder | Approx. duration |
|---|---|---|---|
| Minimum | 2 | 2 | ~2 weeks |
| Target | 3 | 3 | ~2.5–3 weeks |
| Maximum (full PR) | 4 | 3 | ~3.5 weeks |

PR A is an exception — it is the campaign's foundation chapter and all 5 of its problems are at least headline-adjacent. From PR B onward the mix applies. The PR A-prequel ([§7.3](#73-pr-a-prequel-adopted)) is a pure fluency-only PR — no headline picks — by design.

### 7.2 Fluency-builder candidates per chapter

These are problem-types, not exact numbered picks; the actual numbered problems are chosen at PR-scope time based on what each Jackson edition includes. Numbers are illustrative.

**Easy proper-time fluency (Chs. 1–2, 5 — for PR F+ but also weavable into earlier PRs as warm-ups):**
- **Ch. 1 electrostatics:** Gauss's law for spherical/cylindrical/planar symmetry (J3e §1.4-style); image-charge problems (J3e §2.x). *Proper-time lesson:* `u = 0` ⇒ `b = c`, ∂_τ-terms vanish, proper-time reduces to classical **exactly**. The work is the CGS↔SI translation. Quickest possible practice.
- **Ch. 2 boundary-value problems:** parallel-plate capacitor energy density; method of images on planar/spherical boundaries. *Same lesson as Ch. 1.* These are the purest fluency builders in the entire campaign.
- **Ch. 5 magnetostatics:** B-field of a current loop (J3e §5.4-style); B-field of a long solenoid; B-field of a rotating charged sphere. *Proper-time lesson:* steady currents ⇒ ∂_τ = 0, only the `J → (b/c)J` rescaling matters. Slightly richer than Chs. 1–2 because of the current-density transformation.

**Intermediate fluency (PR B, Ch. 11 SR — interspersed with Ch. 11 headline picks if any survive §12.3's pruning):**
- Time dilation of a moving clock — `t → τ` substitution in simplest setting.
- Length contraction — same as above.
- Velocity composition `u'_total = (u + v)/(1 + uv/c²)` — direct exercise of `w/c = u/b`.
- Doppler shift for a moving source — `b` vs `c` in frequency transformation.
- 4-momentum boost — pure bookkeeping practice.

**Intermediate fluency (PR C, Ch. 12 dynamics):**
- Lagrangian and equations of motion for a free relativistic particle (J3e §12.1-style).
- E×B drift (J3e §12.3-style) — exercise of the modified force law `F = e[E + (u/b)×B]` from Eq. (18).
- Motion in crossed uniform E and B fields.
- Particle in a plane EM wave (Volkov solution) — `u·a` is computable in closed form.

**Headline-adjacent fluency (PR D, Ch. 14):**
- Liénard–Wiechert *potentials* (J3e §14.1) before the *fields* (the §12 #2 headline) — derive the framework first.
- Synchrotron radiation (relativistic Larmor for circular motion) — `u·a = 0` so the dissipative term vanishes; pure `c → b` rescaling of the classical answer. Pairs naturally with §12 #5 (non-relativistic Larmor) as a "now relativistic but still circular" practice step.
- Bremsstrahlung from a thin foil — `u·a ≠ 0` but the geometry is linear.

**Headline-adjacent fluency (PR E, Ch. 16):**
- Damping of an oscillator by radiation reaction (J3e §16.2-style) — Abraham–Lorentz at the simplest possible kinematic geometry.
- Radiation reaction on a driven oscillator — extends the above with periodic forcing; `u·a` time-varies cleanly.

### 7.3 PR A-prequel (adopted)

Per [§13.5 D7](#135-decision-points--confirmed-by-author-2026-05-24): a **PR A-prequel** is inserted before PR A, containing 3–4 fluency-builder problems where proper-time reduces to classical exactly. Purpose: build muscle memory on low-stakes content and validate the [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) read-aloud test before Ch. 6's richer problems.

**Candidate problems for the PR A-prequel** (4 problems; final selection at PR-scope time, adjustable to whichever appear in both editions):

| # | Jackson 3e | Short title | Fluency lesson |
|---|---|---|---|
| 1 | §2.1 or §2.2 | Method of images, point charge near grounded plane | `u = 0` ⇒ proper-time reduces to classical exactly. The work is pure CGS↔SI translation. Cleanest possible warm-up. |
| 2 | §1.x | Gauss's law for spherical / cylindrical / planar symmetry | Same lesson. Confirms the author's pattern for "static ⇒ no `b`-dependence". |
| 3 | §5.4 or §5.5 | B-field of a steady current loop or solenoid | Steady current ⇒ `∂_τ = 0`, only the `J → (b/c)J` rescaling matters. First exposure to the current-density transformation in a familiar setting. |
| 4 | §5.x | B-field of a rotating charged sphere | Steady (in the rest frame of the rotation) ⇒ same lesson, slightly richer geometry. Optional fourth problem; can be dropped if time pressure exists. |

**Why these four.** Problems 1–2 contain no `b`-dependence whatsoever (electrostatics) and serve as the absolute simplest first exposure. Problems 3–4 introduce the `J → (b/c)J` substitution while keeping `∂_τ = 0`, so the author practises one substitution at a time. None of these problems touches Eq. 24, so the branched-treatment workflow ([§5.1](#51-branched-treatment-for-eq-24-touching-problems)) does not engage in the prequel.

**Deliverable.** A new PR (called "PR 0" or "PR A-prequel" in §10) with:
- 3–4 per-problem documents under `Roadmapping/Electromagnetism/Jackson/Ch01_Introduction_Electrostatics.md` and/or `Ch05_Magnetostatics.md`
- Companion `.wl` notebooks
- A short retrospective at the end of the PR's commit message: "the prequel's purpose was Y; what we learned that informs PR A is Z." This is the author's data point on whether the [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) read-aloud test caught any Gill-voice drift, and whether the full-template treatment (per D3) is sustainable.

**Duration estimate:** ~1 week (4 problems × ~1.5–2 days each — these are the fastest problems in the entire campaign because no proper-time content is exercised).

**Does the prequel gate #42's close?** No. Per [§13.5 D4](#135-decision-points--confirmed-by-author-2026-05-24), only PR A + (PR D OR PR E) gate the issue close. The prequel is preparation, not a campaign deliverable.

### 7.4 PR 0 retrospective (closed 2026-05-24)

PR 0 closed with four per-problem documents drafted and committed:

| Commit | Problem | Verdict | Eq. 24 engaged? |
|---|---|---|---|
| `a346f1b` | J3e-P1.5 — charged-sphere self-energy | ✅ classical = proper-time exactly | no |
| `7e83e9a` | J3e-P2.1 — image method on grounded plane | ✅ classical = proper-time exactly | no |
| `733aeb0` | J3e-P5.4 — B-field of circular current loop | ✅ classical = proper-time exactly (exact `(b/c)`-cancellation) | no |
| `bacc30f` | J3e-P5.13 — rotating-sphere dipole moment | ✅ classical = proper-time at leading order; `O((u/c)^{2})` deviation noted | no |

Empirical observations from PR 0:

- **Effort estimate held.** All four problems completed in one focused session, comfortably within the ~1-week estimate in [§13.1 O4](#131-objections-with-partial-responses). The MCP verification step took 1–2 calls per problem and required no debugging. The per-problem documents averaged ~150 lines of markdown each, with ~30–50 lines of derivation prose and the remainder being scaffolding (Selection provenance, Comparison table, Verdict, companion-notebook link).
- **Voice held under the read-aloud test.** The Gill-voice markers ("we observe that…", "it is easy to show that…", "it follows that…", connective sentence openers, long subordinate-clause sentences with semicolons, the signature "mathematically equivalent but not physically equivalent" anchor) felt natural in derivations of this difficulty. No anti-pattern hyperbole or Claude-voice tells slipped through. The PR A-prequel's "low-stakes content" framing of [§7.3](#73-pr-a-prequel-adopted) appears to have served its purpose.
- **Template held across four problems of varying difficulty.** The Selection provenance + Source + Setup + (a)/(b)/(c) + Comparison + Verdict structure scaled cleanly from the trivial-reduction case (P1.5) to the three-dimensional volume integral (P5.13). The optional `(c')` branched-treatment subsection did not engage for any PR 0 problem — Eq. 24 is a Hamiltonian-level finding, not a field-equation-level one, and pure electrostatic/magnetostatic problems operate entirely at the field-equation level.
- **Per-paragraph `<!-- TODO -->` density was lower than [§13.5 D2](#135-decision-points--confirmed-by-author-2026-05-24) anticipated.** Each PR 0 problem carried exactly one substantive `<!-- TODO -->` block — on the Selection provenance — and no others. Derivations in PR 0 are mechanical (pragmatic AI use); there was no interpretive prose to disclose. For Ch. 6 problems where the dissipative term engages, the density will increase substantially. The decision to reassess granularity at PR B remains the right call.
- **Surprise finding (P5.13):** the rotating-sphere dipole moment carries an `O((u/c)^{2})` proper-time correction that is not a flagged finding but worth recording, because it is the first per-problem result in the campaign where the proper-time formulation produces a non-zero deviation from the classical answer, however small. The correction is far below the observational floor for macroscopic rotating charged spheres. Flagged in the per-problem document's "Notes for author review" as a candidate observation; not posted to `FINDINGS_for_author_review.md` because it is a structural attribute of the velocity-duality bookkeeping rather than an unresolved inconsistency.

**Implications for PR A.** The MCP workflow, the template, and the voice constraint are all validated for the easy-fluency case. PR A's Ch. 6 problems will be the first real test of the dissipative-term machinery and of whether the 3-day-per-headline-problem estimate in [§7.1](#71-headline-payoff-problems-vs-fluency-builder-problems) is accurate. PR A scope is unchanged from [§7](#7-initial-chapter-selection--canonical-problems-list); the 5 Ch. 6 problems begin immediately after this retrospective is committed.

**Implications for the plan.** No changes recommended to §7, §10, §13 based on PR 0 observations alone — the data point is too small. The next decision point is at PR A's mid-point, when the per-paragraph TODO density observation and the headline/fluency cost ratio can be tested against richer content.

### 7.5 PR A retrospective (closed 2026-05-24)

PR A closed with five per-problem documents drafted and committed:

| Commit | Problem | Verdict | Headline observation |
|---|---|---|---|
| `09a550e` | J3e-P6.1 — Maxwell with magnetic monopoles | ✅ duality preserved | **multi-species ambiguity** in proper-time `b` definition (single-source bias of published framework) |
| `b35eb84` | J3e-P6.4 — EM momentum of moving point charge | ✅ classical = proper-time | reproduces the `4/3 factor` exactly; framework consistent with but does not dissolve the Abraham–Lorentz puzzle |
| `e894771` | J3e-P6.5 — Poynting theorem in macroscopic media | ✅ conservation form preserved | strong-equivalence at the conservation-law level (does not generalise to dynamical observables) |
| `2218eed` | J3e-P6.11 — symmetric stress tensor + Lorentz behaviour | ✅ symmetric, trace = −u | **proper-time group ≠ Lorentz group**; stress tensor's covariance properties depend on which group is invoked |
| `0a258d2` | J3e-P6.20 — radiation pressure on perfect conductor | ✅ `P = 2I/c` (classical and proper-time identical) | **null-result canary** for `(u·a)/b⁴` — perfect-conductor idealisation removes the dissipative contribution; first non-null engagement deferred to PR D |

Empirical observations from PR A:

- **Effort estimate held within [§13.1 O4](#131-objections-with-partial-responses)'s ~2.5–3 week budget.** Each Ch. 6 problem took roughly 2× the PR 0 per-problem time (more interpretive content, more substantive paragraphs, richer Mathematica checks). The 3-day-per-headline-problem and 2-day-per-fluency-builder estimates from [§7.1](#71-headline-payoff-problems-vs-fluency-builder-problems) are roughly accurate; PR A's blended cost was closer to 2.5 days/problem averaged across all five. The 5-day upper bound applies only to genuinely interpretive problems with multiple structural observations (e.g., 6.11 stress tensor).
- **Per-paragraph TODO discipline engaged as planned.** Each PR A problem carried 2–4 substantive `<!-- TODO -->` blocks (vs PR 0's typical 1-per-problem). The density is sustainable for Ch. 6 because each problem has multiple genuine interpretive claims that warrant author review. Per [§13.5 D2](#135-decision-points--confirmed-by-author-2026-05-24), reassessment occurs at PR B; the data so far suggests the per-paragraph density is workable and does not yet rise to noise.
- **Three new structural observations surfaced.** None are flagged inconsistencies; all are structural features of the framework's geometric content worth recording: (1) multi-species ambiguity in 6.1, (2) proper-time-group ≠ Lorentz-group in 6.11, (3) null-result for `(u·a)/b⁴` under perfect-conductor idealisation in 6.20. The first two are candidate observations for future Gill–Zachary follow-up papers; the third is a campaign-internal observation about which problems will engage the dissipative term.
- **Two of three headline-payoff Ch. 6 problems ended in null results.** Problems 6.4 (EM momentum) and 6.20 (radiation pressure) reproduce the classical answers exactly. Problem 6.11 surfaces a structural observation (covariance-group distinction) but the *numerical* stress tensor components are unchanged. This was expected per [§7](#7-initial-chapter-selection--canonical-problems-list) and [§13](#13-devils-advocate-review-and-what-we-cannot-honestly-fix): PR A's role is to validate that the framework is consistent with classical EM in cases where the dissipative `(u·a)/b⁴` term cannot engage. The first non-null headline payoff is therefore deferred to PR D (Ch. 14 Liénard–Wiechert), where the source acceleration is essential and the third term of Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] becomes operationally measurable.
- **MCP workflow exercised across many problem types.** PR A's verifications covered symbolic substitution (6.1 duality), volume integrals (6.4 field momentum), conservation-law structure (6.5 Poynting), tensor symmetry and trace identities (6.11), and time-averaged Maxwell stress (6.20). No bugs in the workflow; the per-problem `.wl` notebooks have been written runnable independent of the MCP. Confidence is high that the workflow scales to PR B–E.

**Implications for PR B (Ch. 11 special relativity).** PR A's covariance-group observation in [Problem J3e-P6.11](#problem-j3e-p611--symmetric-stress-tensor-and-lorentz-behaviour) becomes the load-bearing entry point for PR B. Ch. 11 is where the proper-time group of the Maxwell paper §1.3 is exercised in detail; the structural observation surfaced in PR A will sharpen there. If the sharpening produces a genuine inconsistency (vs an extension question), an entry in [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) would be warranted at that time.

**Implications for the plan.** No changes recommended to §7, §10, §13 based on PR A observations. The graceful-degradation provisions in §13.1 O4 were not needed; the campaign is on schedule. The Eq. 24 branched-treatment workflow (per [§5.1](#51-branched-treatment-for-eq-24-touching-problems) and [§13.5 D1](#135-decision-points--confirmed-by-author-2026-05-24)) did not engage in PR A; expected to engage in PR C (Ch. 12 relativistic dynamics, where the Hamiltonian formulation is exercised). Per the [§13.5 D4](#135-decision-points--confirmed-by-author-2026-05-24) criterion, issue #42 remains open pending PR D *or* PR E.

**Headline next problem.** Per [§12.1](#121-per-problem-briefs), the campaign's full dissipative-term test arrives at PR D Problem J3e-Ch14 (Liénard–Wiechert fields of an accelerating point charge). This is where the third term `e(u·a)/(b⁴s³)` of Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] first appears with non-zero numerical content, and where the campaign's experimental-comparison hooks (issue [#43](https://github.com/temoTxt/PyPhysics/issues/43)) become operationally measurable.

### 7.6 PR B retrospective (closed 2026-05-24)

PR B closed with five per-problem documents drafted and committed:

| Commit | Problem | Verdict | Key observation |
|---|---|---|---|
| `ff4baeb` | J3e-P11.1 — relativistic velocity addition | ✅ classical = proper-time | sharpened group distinction at 4-velocity level |
| `6eade8c` | J3e-P11.3 — relativistic Doppler effect | ✅ identical observable | proper-time form `(b+u)/c` is **linear** in `(b, u)`, vs square-root in `\beta` |
| `6b512ce` | J3e-P11.5 — 4-vector position transformation | ✅ Minkowski interval invariant | 4-positions transform identically; group distinction is at 4-velocity level and beyond |
| `675a174` | J3e-P11.6 — Lorentz boost of `\mathbf{E}` and `\mathbf{B}` | ✅ both Lorentz invariants preserved | proper-time form is **linear** in `(b, u)`; multi-source case is where the groups truly diverge |
| `5bd7b21` | J3e-P11.8 — Lorentz invariants of EM field | ✅ identified `F^{\mu\nu}F_{\mu\nu}, F^{\mu\nu}\,{}^{*}F_{\mu\nu}` | closing claim: `F^{\mu\nu}` alone has same invariants in both groups; framework's distinguishing structure lives in invariants involving source 4-velocity |

Empirical observations from PR B:

- **Effort estimate held within [§13.1 O4](#131-objections-with-partial-responses)'s ~2.5–3 week budget.** Five Ch. 11 problems at roughly the same per-problem cost as Ch. 6 problems in PR A — the kinematic content is shorter than Ch. 6's energy-momentum analysis, but the structural-observation content (proper-time group / Lorentz group distinction) is denser. Net cost per problem roughly matches PR A's blended average.
- **Per-paragraph TODO discipline continued sustainably.** Each PR B problem carries 2–4 substantive `<!-- TODO -->` blocks, the same density as PR A. Per [§13.5 D2](#135-decision-points--confirmed-by-author-2026-05-24), the reassessment scheduled for the start of PR B is now overdue — having now seen the density across ten problems (5 in PR A + 5 in PR B), the decision is to **keep the per-paragraph density through PR C** and reassess at PR D, where the dissipative-term content will produce a different interpretive load.
- **One structural pattern emerged across multiple problems.** The proper-time `(b, u)` variables consistently produce **linear** algebraic expressions where the classical `(\gamma, \beta)` variables produce square-roots. Specifically: the Doppler factor `(b+u)/c` vs `\sqrt{(1+\beta)/(1-\beta)}`, the field-boost prefactors `(b E - u B)/c` vs `\gamma(E - \beta B)`. This is a pattern about variable choice, recorded in the per-problem documents but not flagged as a finding because it is a property of the velocity-duality identity itself, not a separate observation.
- **Closing structural claim for the chapter.** `F^{\mu\nu}` alone has the same algebraic Lorentz invariants under both groups; the framework's distinguishing structure lives in invariants involving BOTH `F^{\mu\nu}` and the source 4-velocity `(b, \mathbf{u})`. This claim is the natural bridge to PR C (Ch. 12 relativistic dynamics), where the equations of motion involve `F^{\mu\nu}u_{\nu}` (the modified Lorentz force of Eq. (18) of the Maxwell paper) and where the proper-time group's distinct action will surface at the dynamical level.
- **Wolfram MCP disconnect handled gracefully.** The MCP subprocess dropped mid-PR-B between problems J3e-P11.3 and J3e-P11.5. The campaign halted authoring per the author's directive ("stop until the wolfram MCP is online"), the underlying service was verified to be healthy (HTTP 200 on authenticated `initialize`), the disconnect was diagnosed as Claude Code-side, and after the user reconnected via `/mcp` the campaign resumed cleanly. This is the first operational test of the [§13.2 O4](#132-objections-with-no-honest-mitigation) "no kill switch / accept drift" provision; the campaign successfully paused and resumed without state loss.

**Implications for PR C (Ch. 12 relativistic dynamics).** The closing observation of PR B — that the proper-time framework's distinguishing structure lives in invariants involving both `F^{\mu\nu}` and the source 4-velocity — predicts that PR C will be where the proper-time group's distinct action first becomes operationally visible at the dynamical level. The Hamiltonian formulation in Ch. 12 will exercise Eq. (18) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] (modified Lorentz force), and the Hamiltonian itself depends on `b` through the proper-time kinetic energy. This is also the first place in the campaign where Eq. (24) of the Maxwell paper may engage, triggering the branched-treatment workflow of [§5.1](#51-branched-treatment-for-eq-24-touching-problems) per [§13.5 D1](#135-decision-points--confirmed-by-author-2026-05-24).

**Implications for the plan.** No changes recommended to §7, §10, §13 based on PR B observations. The graceful degradation provisions in §13.2 O4 were exercised (the MCP outage) and worked as designed. Per-paragraph TODO discipline reassessment is deferred from "start of PR B" to "start of PR D", given that the pattern has stabilised over ten problems and the next interpretive load shift is expected at the radiation-reaction chapter.

**Issue #42 acceptance criteria status.** PR A merged ✅, PR D OR PR E pending. PR B is not gating; it is preparatory chapter work for PR C–E.

### 7.7 PR C retrospective (closed 2026-05-24)

PR C closed with five per-problem documents drafted and committed:

| Commit | Problem | Verdict | Key observation |
|---|---|---|---|
| `ab5e661` | J3e-P12.1 — free relativistic Lagrangian | ✅ identical action | proper-time Lagrangian density `L_{\tau} = -mc^{2}` is *constant*; worldline dynamics from constraint `b^{2}-u^{2}=c^{2}` |
| `c1beec2` | J3e-P12.2 — cyclotron motion | ✅ identical observables | **proper-time cyclotron frequency `\omega_{\tau} = qB/(mc)` is γ-INDEPENDENT**; γ-dependence of lab-frame `\omega_{c}` is pure time-dilation. Podcast pick #4 lesson |
| `ea2166f` | J3e-P12.5 — E×B drift | ✅ identical drift velocity | `\mathbf{u}\cdot\mathbf{a} = 0` throughout (drift + cyclotron); dissipative term inactive |
| `13d9bf5` | J3e-P12.10 — Hamiltonian formulation | ✅ identical Hamiltonian | **Eq. 24 anticipated but did NOT engage** (problem is non-spin); branched treatment deferred to PR F+ |
| `801ddff` | J3e-P12.14 — charged particle in plane wave | ⚠ first PR C engagement of `(u·a)/b⁴` | dissipative term first activates with finite time-average; bridges PR C → PR D and issue #43 |

Empirical observations from PR C:

- **Effort estimate held within [§13.1 O4](#131-objections-with-partial-responses)'s ~2.5–3 week budget.** Five Ch. 12 problems at roughly the same per-problem cost as PR A and PR B. Dynamics content (Lagrangian, Hamiltonian, equations of motion) is more derivational than the kinematic content of PR B, but the structural-observation content is comparable. Net cost matches PR A and PR B's blended average.
- **Per-paragraph TODO discipline continued sustainably.** Each PR C problem carries 2–4 substantive `<!-- TODO -->` blocks, the same density as PR A and PR B. Reassessment scheduled for the start of PR D per [§7.6](#76-pr-b-retrospective-closed-2026-05-24) remains the right time; density across 15 problems is now stable and the next interpretive load shift will come at the dissipative-term content of PR D.
- **Eq. 24 branched-treatment workflow did NOT engage in PR C.** [Problem J3e-P12.10](#problem-j3e-p1210--hamiltonian-formulation-for-charged-particle) was anticipated by [§5.1 of the plan](#51-branched-treatment-for-eq-24-touching-problems) as the likely first site of Eq. 24's flagged finding. The actual outcome: Eq. 24 concerns the *Pauli-form* spin-magnetic-moment and `V^{2}/(2mc^{2})` terms, which appear only with explicit spin extension. Jackson Ch. 12 P10 is non-spin classical, so the branched-treatment workflow remains scaffolded but does not engage until a Pauli-reduction problem arrives (likely in PR F+ via a hyperfine or spin-orbit problem). The scaffolding established in [§5.1](#51-branched-treatment-for-eq-24-touching-problems) is correct; it is simply standing in reserve for content that has not yet arrived.
- **Two new headline-level structural observations.**
  - **Cyclotron frequency in proper time is γ-independent** (J3e-P12.2). The cleanest example in the campaign of the proper-time formulation producing an algebraically simpler dynamical statement than the classical formulation. The energy-dependence in the lab-frame cyclotron frequency is pure time-dilation, not a feature of the underlying dynamics.
  - **The dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient first engages in PR C content at J3e-P12.14** (charged particle in plane wave). At low wave intensity (`a_{0} \ll 1`), the contribution recovers the Thomson-scattering rate; at relativistic intensity (`a_{0} \gtrsim 1`), the regime of Cole/Poder/Wistisen 2018 experiments, the comparison is deferred to issue [#43](https://github.com/temoTxt/PyPhysics/issues/43) and to PR D's Liénard–Wiechert content.
- **Bridge to PR D is established.** The closing problem J3e-P12.14 sets up the dissipative-term machinery that PR D (radiation by accelerating charges) will fully exercise. The first non-null engagement of the dissipative coefficient is now operational, and the issue-#43 experimental-comparison work has its first dynamical-side derivation available as a starting point.

**Implications for PR D (Ch. 14 radiation by moving charges).** Per [§12.1](#121-per-problem-briefs), PR D contains podcast picks #2 (Liénard–Wiechert third term) and #5 (non-relativistic Larmor). The dissipative-term engagement first surfaced in J3e-P12.14 is the natural lead-in; PR D will exercise the third term `e(\mathbf{u}\cdot\mathbf{a})/(b^{4}s^{3})` of Eq. (7) directly. The issue #43 comparison work becomes operational with PR D's content.

**Implications for the plan.** No changes recommended to §7, §10, §13 based on PR C observations. The Eq. 24 branched-treatment scaffolding is retained for PR F+ spin-content problems. Per-paragraph TODO reassessment remains scheduled for the start of PR D.

**Issue #42 acceptance criteria status.** PR A merged ✅, PR D OR PR E pending. PR B and PR C are preparatory chapter work; neither is gating per [§13.5 D4](#135-decision-points--confirmed-by-author-2026-05-24).

### 7.8 PR D retrospective (closed 2026-05-24)

**PR D is the campaign's first headline-payoff PR.** All five Ch. 14 problems drafted, all five Wolfram MCP-verified, and the campaign's load-bearing field-level prediction (the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` third term of Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]) is now exhaustively documented.

| Commit | Problem | Verdict | Headline observation |
|---|---|---|---|
| `31c7779` | J3e-P14.1 — Liénard–Wiechert potentials | ✅ identical potentials | foundation: `u/b = v/c` reduces proper-time and classical potentials to the same expressions |
| `437f755` | J3e-P14.2 — LW fields with third term | ⚠ **first qualitatively new field-level prediction** | **proper-time third term `e(u·a)[r×(u×r)]/(b⁴s³)` is absent from classical EM**; predicts longitudinal radiation component when `u·a ≠ 0` (podcast pick #2) |
| `d32a806` | J3e-P14.3 — non-relativistic Larmor | ✅ identical at leading order; `O((u/c)²)` correction | **podcast pick #5's honest moment:** the correction is below floor at velocities where Larmor is valid, and order unity where Larmor must be replaced by Liénard |
| `e6ac1c3` | J3e-P14.5 — synchrotron from circular motion | ✅ identical to classical | third term vanishes by `u·a = 0` (orthogonality); null-result confirmation against precision-tested synchrotron data |
| `8793f90` | J3e-P14.6 — bremsstrahlung from linear deceleration | ⚠ third term engages | non-null counterpoint to synchrotron; **MeV-energy bremsstrahlung identified as candidate new experimental setting** (not previously emphasised in Gill corpus) |

#### Campaign-level observations from PR D

**The dissipative-term machinery is now fully exercised at the field level.** Every dissipative-term scenario in the campaign's plan has been documented:

| Scenario | Engagement | PR D problem |
|---|---|---|
| `\mathbf{u} = 0` (static) | vanishes | not in PR D; documented across PRs 0-A |
| `\mathbf{a} = 0` (uniform motion) | vanishes | not in PR D; documented in PR A (J3e-P6.4) |
| `\mathbf{u} \perp \mathbf{a}` (circular motion) | vanishes | J3e-P14.5 (synchrotron) ✅ null result |
| Steady current | vanishes (averaged) | not in PR D; documented in PR 0 (J3e-P5.4) |
| Perfect-conductor idealisation | vanishes (response is instantaneous) | not in PR D; documented in PR A (J3e-P6.20) |
| `\mathbf{u}\cdot\mathbf{a} \neq 0` (linear acceleration) | **engages** | J3e-P14.6 (bremsstrahlung) ⚠ non-null |
| Plane-wave figure-8 at relativistic intensity | **engages** | not in PR D; first noted in PR C (J3e-P12.14) |
| Radiation-reaction-dominated extreme intensity | **engages** | not in PR D; subject of PR E and issue #43 |

The picture is now complete: the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient engages whenever the source's proper-time velocity has a component along its proper-time acceleration, and vanishes in every kinematic configuration where they are orthogonal or one of them is zero. The campaign has documented every category.

**The campaign's load-bearing finding is articulated.** PR D's J3e-P14.2 is the campaign's first qualitatively new field-level prediction:

> *The proper-time Liénard–Wiechert formulation predicts a **longitudinal radiation component** in any kinematic configuration where `\mathbf{u}\cdot\mathbf{a} \neq 0`. Classical Liénard–Wiechert radiation in the radiation zone is purely transverse; the proper-time third term breaks this purity by introducing components along both `\mathbf{u}` (longitudinal) and `\hat r` (radial).*

Whether this prediction is supported by experiment is the open question of issue [#43](https://github.com/temoTxt/PyPhysics/issues/43). The candidate experimental settings now include:

- **Cole 2018 / Poder 2018 / Wistisen 2018** — GeV electron beam × extreme laser intensity. The original §12 podcast-pick targets.
- **MeV bremsstrahlung at clinical-physics energies** — newly identified in PR D's J3e-P14.6. Precision-experiment-rich setting (medical physics) where the third-term contribution would be order unity and `\sim 1-2\%` measurements achieve distinguishability.

#### Empirical observations from PR D's authoring

- **Effort estimate held within [§13.1 O4](#131-objections-with-partial-responses)'s ~2.5–3 week budget.** Five Ch. 14 problems at roughly the same per-problem cost as PR A–C, with two problems (J3e-P14.2, J3e-P14.6) carrying heavier interpretive load (substantive claims about the third-term contribution) and three lighter (J3e-P14.1, J3e-P14.3, J3e-P14.5).
- **Per-paragraph TODO discipline reassessment is now overdue.** Per [§7.6 of the plan](#76-pr-b-retrospective-closed-2026-05-24), the planned reassessment of per-paragraph density was deferred to "start of PR D"; PR D has now happened. The density remains 2–4 TODO blocks per problem across PR A–D (20 problems), and the noise concern raised in §13.1 O3 has not materialised. **Decision: keep per-paragraph TODO density through PR E** (the final headline PR); reassess at start of PR F+ if backfill content engages.
- **Wolfram MCP performed reliably throughout PR D after the PR B reconnect.** No further disconnects.

#### Implications for PR E (Ch. 16 Radiation Damping)

PR E will engage the radiation-reaction physics of the **Abraham–Lorentz–Dirac** problem, podcast pick #3. The setup is now in place: PR C–D have documented the dissipative term across all kinematic categories; PR E will exercise it in the radiation-back-reaction regime where Cole/Poder/Wistisen 2018 measurements live. The first **non-trivial** radiation-reaction observable — a comparison of proper-time predictions against the classical Lorentz–Abraham–Dirac equation — is the natural PR E content.

Per [§11 acceptance criteria](#11-acceptance-criteria-for-closing-issue-42), PR D being merged satisfies the "PR D OR PR E" requirement; issue #42 can close once PR D is reviewed and merged. PR E becomes optional for #42 closure but remains part of the campaign's scope per §7.

**Implications for the plan.** No changes recommended to §7, §10, §13 based on PR D observations. The per-paragraph TODO discipline decision (keep through PR E, reassess at PR F+) is recorded above. Issue #43's experimental-comparison work has its first complete dynamical-side derivation; the comparison document at `Electromagnetism/Jackson/Experimental_Comparisons/radiation_reaction_2018.md` can now be drafted.

**Issue #42 acceptance criteria status.** PR A merged ✅; **PR D drafted (5/5) and ready for review** — once merged, the second criterion is satisfied. The third criterion (`Roadmapping/README.md` cross-reference) remains pending as a small PR-final commit.

### 7.9 PR E retrospective (closed 2026-05-24)

**PR E is the campaign's second headline-payoff PR and closes the campaign's headline narrative.** All four Ch. 16 problems drafted, all four Wolfram MCP-verified, and the campaign's load-bearing rhetorical claim (the proper-time framework dissolves the classical Abraham–Lorentz pathologies) is now documented with its complete narrative arc and honest caveats.

| Commit | Problem | Verdict | Headline observation |
|---|---|---|---|
| `4fc5aac` | J3e-P16.1 — Abraham–Lorentz + dissolution claim | ⚠ **most rhetorically loaded problem** | **proper-time first-order `\partial_{\tau}` structure → no runaway, no pre-acceleration; classical AL pathologies structurally absent** (podcast pick #3) |
| `0776981` | J3e-P16.2 — RR-damped harmonic oscillator | ✅ identical natural linewidth | concrete-example fluency-builder; `\Gamma = \tau_{0}\omega_{0}^{2}` recovers atomic-transition linewidth `\sim 5 \times 10^{7}` s⁻¹ |
| `eea90f7` | J3e-P16.3 — pre-acceleration IVP analysis | ⚠ pre-acceleration solution form `a(t<0) = (F_{0}/m)e^{t/\tau_{0}}` | detailed-pathology document; structural argument sharpened with three-tier honest framing |
| `5208087` | J3e-P16.5 — proper-time RR for Cole/Poder | ⚠ sketch-level prediction | bridges PR E to issue #43; third-term contribution `~10\%` correction to LL, at boundary of distinguishability |

#### The campaign's headline narrative is now complete

PR E closes the campaign's narrative arc that PR D opened:

1. **PR D (Ch. 14)** established the framework's first qualitatively new field-level prediction: the proper-time third term in the Liénard–Wiechert fields predicts a longitudinal radiation component absent from classical EM.
2. **PR E (Ch. 16)** established the framework's claim to **dissolve the classical Abraham–Lorentz pathologies** (runaway, pre-acceleration) by structural difference (first-order vs third-order in time-derivative).
3. **Issue #43** is the bridge: whether these structural and field-level predictions match the published Cole/Poder/Wistisen 2018 measurements is the campaign's primary experimental test.
4. **Issue #48** (new in PR D) is a secondary bridge: whether the third-term contribution to MeV-energy bremsstrahlung is measurable in clinical-physics precision experiments.

#### Honest caveats — the three-tier framing

The dissolution claim is framed across PR E with three operational caveats that must not be softened:

1. **Structural, not predictive.** The dissolution removes the classical pathologies in the proper-time formulation. It does not prove the proper-time RR force is numerically correct at all energies. Issue #43's comparison is the discriminating test.
2. **Does not address self-energy.** The classical electron's self-energy divergence is a separate pathology not dissolved by the proper-time framework.
3. **Conditional on framework correctness.** Per [§13 of the plan](#13-devils-advocate-review-and-what-we-cannot-honestly-fix), every claim is "*if Gill–Zachary is the correct formulation, then X*". Whether the framework is correct is the open question that the campaign as a whole tests, not pre-judges.

#### Empirical observations from PR E

- **Effort estimate held within [§13.1 O4](#131-objections-with-partial-responses)'s ~2 week budget.** Four Ch. 16 problems at slightly heavier interpretive load than PR D (the dissolution claim required more careful framing), but the per-problem documents averaged roughly the same length as PR A–D.
- **Per-paragraph TODO discipline maintained.** Each PR E problem carries 2–4 substantive `<!-- TODO -->` blocks, consistent with PR A–D's density. The discipline is stable; reassessment scheduled for start of PR F+ remains the right time.
- **Wolfram MCP performed reliably throughout PR E.** No disconnects.
- **Two podcast picks delivered in PR E (#3) and PR D (#2, #5).** Of the five §12 podcast picks, **four have been documented** with their proper-time treatments. The remaining pick — #1 (EM momentum of uniformly-moving charge, J3e-P6.4) — was delivered in PR A.

#### Issue #42 acceptance criteria — NOW FULLY SATISFIED

Per [§11](#11-acceptance-criteria-for-closing-issue-42):

| Criterion | Status |
|---|---|
| PR A merged | ✅ |
| **PR D AND PR E drafted** (either satisfies "D OR E"; both satisfied) | ✅ ✅ |
| Root or `Roadmapping/README.md` cross-reference | pending (small PR-final commit) |

**Issue #42 can close as soon as PR D (or PR E, or both as one large PR) is reviewed and merged**, after the `Roadmapping/README.md` cross-reference is added. The campaign has delivered all the load-bearing content envisioned at planning time.

#### Implications for PR F+ (backfill chapters)

Per [§7 of the plan](#7-initial-chapter-selection--canonical-problems-list), PR F+ covers backfill chapters Chs. 1–5, 7–10, 13, 15 at 4–6 problems each. With PRs A–E complete, the campaign has produced **24 problems** across 5 chapters. The backfill scope is 60–80 additional problems across 10 chapters, totalling 80–100 problems for the full campaign (consistent with [§1's](#1-goal-and-scope) "~50–100 total" estimate).

Per [§13.5 D5](#135-decision-points--confirmed-by-author-2026-05-24) (no kill switch), the campaign can continue into PR F+ indefinitely, or close at the current 24-problem state and declare the headline-payoff content complete. Both are legitimate outcomes; the author's discretion governs the choice.

**Implications for the plan.** No changes recommended to §7, §10, §13 based on PR E observations. The campaign's narrative is complete; PR F+ is pedagogical-completeness work, not headline-payoff work. The Eq. 24 branched-treatment workflow remains scaffolded for any future spin/hyperfine problem in PR F+.

**Recommendation:** consider opening the merge PR for #42 review at this checkpoint. The 24 problems across PR 0 + PR A + PR B + PR C + PR D + PR E constitute a coherent self-contained deliverable, with the headline narrative fully articulated and honestly framed. Issue #42 can close formally on the §11 acceptance criteria being satisfied; PR F+ becomes optional ongoing work without #42's open status as motivation.

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
- ✅ `Roadmapping/Electromagnetism/Jackson/README.md` status table updated.

Per PR:
- ✅ `Roadmapping/Electromagnetism/README.md` thread-status updated.
- ✅ Root or `Roadmapping/README.md` reflects the existence of the new research thread (after PR A only — not every subsequent PR). No "Four pillars" or other institutional-weight framing; a one-line cross-reference alongside the other Roadmapping threads suffices.
- ✅ Any ⚠ verdicts cross-posted to `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` if they reflect on the Gill papers.

---

## 9. Open questions / decisions deferred to the human

1. **Template debate (§3 bullets).** Three template choices flagged — best resolved after one chapter is fully populated and the trade-offs are concrete.
2. **Mathematica MCP licensing.** Per `Roadmapping/Equation_Verification/README.md`, MCP went online 2026-05-10. Confirmed available for this campaign. No action needed.
3. **Animation companion.** Should each Manim scene for this campaign live under `Roadmapping/Animations/manim_scenes/jackson_*.py`, or get its own subdir `Roadmapping/Animations/manim_scenes/Electromagnetism/`? Defer to first scene's PR.
4. **Cross-referencing to History campaign.** The History chapters on classical EM (Ch. 2: Classical synthesis 1860–1900) cite Maxwell's 1865 paper; should Jackson-problem docs back-reference those historical chapters where they introduce a concept? Recommended: yes, but only for the load-bearing problems (e.g. Ch. 6 displacement current → History Ch. 2).
5. **Future textbooks.** The `Roadmapping/Electromagnetism/Jackson/` structure leaves room for a follow-on classical-EM textbook (Zangwill, Schwinger, or Griffiths' *Introduction to Electrodynamics*) under `Roadmapping/Electromagnetism/`. *Griffiths' Introduction to Quantum Mechanics* does **not** belong here — that work, if undertaken, lives under a separate `Roadmapping/Quantum/Griffiths/` thread. Whether and when to add a second EM textbook is deferred to mid-campaign.
6. **Should the canonical-problems list (§7) be locked at plan time, or evolved per-PR?** Recommendation: lock the PR A list (5 problems above); evolve PRs B+ as we learn what's most pedagogically valuable.

---

## 10. PR sequencing

| PR | Scope | Estimated effort | Blocks |
|---|---|---|---|
| **PR 0 — PR A-prequel** (fluency warm-up) | `Roadmapping/Electromagnetism/` scaffold (README, template, cheat sheet, `Jackson/README.md`), bibliography stubs for Jackson 2e + 3e, 3–4 fluency-builder problems from Chs. 1–2 + 5 per [§7.3](#73-pr-a-prequel-adopted). Full template per [§13.5 D3](#135-decision-points--confirmed-by-author-2026-05-24). Closing retrospective on Gill-voice + effort-estimate observations. | ≈1 week (4 problems × ~1.5–2 days; these are the fastest in the campaign because no proper-time content is exercised) | nothing |
| **PR A — Ch. 6** | `Ch06_Maxwell_Equations_Macroscopic_Media.md` with 5 problems fully solved per [§7](#7-initial-chapter-selection--canonical-problems-list) (6.1, 6.4, 6.5, 6.11, 6.20). Full template throughout (per D3). One-line cross-reference in root or `Roadmapping/README.md`. Branched treatment of any Eq. 24-touching problem per [§5.1](#51-branched-treatment-for-eq-24-touching-problems) (none expected in Ch. 6, but the workflow is documented). | ≈2.5–3 weeks (5 problems × ~3 days each, all headline-adjacent so no shorter-fluency-builder mix applies in PR A) | PR 0 |
| **PR B — Ch. 11 Special Relativity** | 5–7 problems. Adds animation companion if any problem yields a striking visual. | ~1.5 weeks | PR A |
| **PR C — Ch. 12 Relativistic Dynamics** | 5–7 problems. Cross-refs DRQM I. | ~1.5 weeks | PR A |
| **PR D — Ch. 14 Radiation by Moving Charges** | 5–7 problems. *First chapter where dissipative term contributes.* | ~2 weeks | PR A |
| **PR E — Ch. 16 Radiation Damping** | 3–5 problems. **Headline payoff.** Likely yields one ⚠ for author review. | ~2 weeks | PR D |
| **PR F+** | Backfill Chs. 1–5, 7–10, 13, 15 in priority order. | 1 week each, ≤2 chapters per PR |  |

Total: 12–14 PRs over 6–12 months, depending on cadence.

---

## 11. Acceptance criteria (for closing issue #42)

The issue closes when:
1. PR A is merged (thread scaffolded under `Roadmapping/Electromagnetism/`, Ch. 6 complete).
2. Either PR D *or* PR E is merged (at least one chapter where the proper-time reformulation contributes a non-trivial new term has been worked).
3. Root `README.md` (or `Roadmapping/README.md` if one exists) cross-references the new thread alongside the existing Roadmapping research threads. No "Four pillars" framing — see [§2](#2-repository-structure) and [§13.1 O7](#131-objections-with-partial-responses).

Subsequent chapter work continues under PR labels but does not require keeping the issue open.

---

## 12. Podcast-elevation candidates

Five Jackson problems are flagged for elevation into the History podcast (issue #7). Selection criterion: the proper-time-vs-classical contrast is rich enough to give the three-voice cast (Historian / Physicist / Experimentalist) genuine material in step 4 of the standard episode structure ([`Roadmapping/History/Podcast/README.md`](../../Roadmapping/History/Podcast/README.md)).

The five form a **pedagogical gradient** — a listener can pick up at any rung without needing earlier episodes. They are numbered by selection order, *not* by where they fall in the gradient; the table is sorted by gradient position.

**All entries in the gradient table are conditional on the Gill–Zachary framework being the correct formulation.** "Effect on the classical answer" describes what the framework *predicts*, not what has been experimentally validated. See [§13](#13-devils-advocate-review-and-what-we-cannot-honestly-fix) for the honest limits of these claims.

| # | Jackson | Problem | `u·a` | Effect on the classical answer | Lesson | Target episode |
|---|---|---|---|---|---|---|
| 4 | 3e §12.2 | Cyclotron motion in uniform B | `= 0` | only `c → b` rescaling | *"the conventions can look identical"* | [`episode_03_old_quantum_theory.md`](../../Roadmapping/History/Podcast/episode_03_old_quantum_theory.md) |
| 1 | 3e §6.4 | EM momentum of uniform-velocity charge | `= 0` | force-law re-bookkeeping via `(u/b)×B` | *"identical observables, different accounting"* | [`episode_02_classical_synthesis.md`](../../Roadmapping/History/Podcast/episode_02_classical_synthesis.md) |
| 5 | 3e §14.2 | Non-relativistic Larmor radiation | `≠ 0`, small | leading correction in `(u/c)²` | *"the dissipative term turns on, and here's how big"* | [`episode_06_synthesis_divergence_map.md`](../../Roadmapping/History/Podcast/episode_06_synthesis_divergence_map.md) |
| 2 | 3e Ch. 14 | Liénard–Wiechert fields of accelerating charge | `≠ 0` | new third term `e(u·a)/(b⁴s³)` in E, B | *"a new measurable structure appears"* | [`episode_06_synthesis_divergence_map.md`](../../Roadmapping/History/Podcast/episode_06_synthesis_divergence_map.md) |
| 3 | 3e Ch. 16 | Abraham–Lorentz radiation reaction | `≠ 0`, large | classical pathology dissolved | *"runaway / pre-acceleration was a coordinate artifact"* | [`episode_05_QED_renormalization_solid_state.md`](../../Roadmapping/History/Podcast/episode_05_QED_renormalization_solid_state.md) |

### 12.1 Per-problem briefs

**#4 — Cyclotron motion (Jackson 3e §12.2).** Cyclotron frequency `ω_c = eB/(γmc)` is something every grad student knows cold. For circular motion `u·a = 0` identically, so the dissipative wave-equation term vanishes; proper-time prediction is `ω_c = eB/(mb)` — same shape, `c → b`. **What's the same** is the lesson. *Historian:* Lorentz force (1895), Lawrence cyclotron (1932), Penning-trap g-2 (Dehmelt 1980s → Fermilab 2021/2023). *Physicist:* derives the dissipative term getting written down and crossed out — primes listener for what `u·a ≠ 0` will look like. *Experimentalist:* Penning-trap g-2 hits ~10⁻¹³ precision; size of `u/c` for a trapped electron puts the proper-time correction below floor. Lowest framing-principle-disclaimer surface area in the five.

**#1 — EM momentum of a uniformly-moving point charge (Jackson 3e §6.4).** Abraham–Minkowski controversy (1908–) is the historical hook: *where does EM momentum live, in the matter or the field?* Uniform motion ⇒ `u·a = 0`, so the dissipative term vanishes again; proper-time effect lives purely in the modified force law `F = e[E + (u/b)×B] + ...` (Eq. 18 of the Maxwell paper). *Historian:* Abraham vs Minkowski, persistence of the debate through Mansuripur (2012). *Physicist:* CGS / SI / proper-time derivations side-by-side. *Experimentalist:* Brevik radiation-pressure experiments, optical-tweezer momentum measurements — *can the `b/c` rescaling distinguish? Probably not at current precision; here's the regime.*

**#5 — Non-relativistic Larmor radiation (Jackson 3e §14.2).** Pedagogical middle rung. `P = (2e²/3c³)|a|²` proper-time-modified to `(2e²/3b³)|a|² × [1 + corrections in u·a/b²]`. Fractional correction scales monotonically in `u/c`: ~10⁻⁴ at chemistry scales, ~10⁻² at 1 keV, order unity at 100 keV. *Historian:* Larmor (1897) + Liénard (1898) as a paired arc. *Physicist:* explicit `u/c` expansion, plot the correction vs velocity. *Experimentalist:* "below floor at the non-relativistic energies where Larmor is valid; at energies where the correction is visible, the right formula is Liénard, not Larmor." Hands off to #2 / #3.

**#2 — Liénard–Wiechert fields (Jackson 3e Ch. 14).** Conditional prediction: *if* the proper-time framework is the correct formulation, the third term in E and B — `e(u·a)[r×(u×r)] / (b⁴s³)` (Maxwell paper Eq. 7) — is a longitudinal radiation component absent from the classical Liénard–Wiechert answer. *Historian:* Liénard (1898) + Wiechert (1900). *Physicist:* explicit proper-time re-derivation, contrast with classical Liénard–Wiechert, computes the predicted-deviation magnitude as a function of `u·a`. `[cue: animation jackson_ch14_lienard_wiechert_third_term]` *Experimentalist:* synchrotron polarization, beamstrahlung at lepton colliders, laser-electron Compton scattering — the Cole / Poder / Wistisen 2018 experiments are the natural place to test whether the conditional prediction is consistent with data. See **issue #43** for the dedicated comparison work.

**#3 — Abraham–Lorentz radiation reaction (Jackson 3e Ch. 16).** Pedagogical pivot, not validation. Classical radiation-reaction theory (Abraham 1903 → Lorentz 1904 → Dirac 1938 → Wheeler–Feynman 1945) carries runaway / pre-acceleration pathology that has resisted clean resolution for ~120 years. *If* the Gill–Zachary framework is the correct formulation, the proper-time wave equation produces radiation-reaction drag from the `(u·a)/b⁴ · ∂_τ` coefficient of Eq. (4) without invoking Lorentz–Dirac, advanced potentials, mass renormalization, or runaway solutions. Whether the framework *is* the correct formulation is a separate, open question — this document derives the conditional prediction only. *Historian:* full Abraham → Dirac → Wheeler-Feynman → Spohn / Rohrlich arc. *Physicist:* derives the classical pathology, shows the proper-time formulation does not exhibit it under the stated substitution rules. *Experimentalist:* Cole et al. PRX **8**, 011020 (2018); Poder et al. PRX **8**, 031004 (2018); Wistisen et al. Nat. Comm. **9**, 795 (2018). Live data at ELI / CALA / Apollon. See **issue #43**. Highest framing-principle-disclaimer surface area in the five; Physicist must deliver the "different mathematical convention, not new physics" disclaimer within the first 60 seconds of the proper-time interlude.

### 12.2 Coordination with the History project

The 9 podcast episode scripts already exist as drafts in [`Roadmapping/History/Podcast/`](../../Roadmapping/History/Podcast/). Elevating these 5 Jackson problems into the corresponding episodes is **edit work on existing scripts**, not new-episode work. It does not depend on the full Electromagnetism campaign being complete — only on the relevant Jackson problem documents (PRs A–E of #42) being far enough along that the Physicist's derivation has a citable source.

Per-problem readiness gates:

| Problem | Episode | Blocked on |
|---|---|---|
| #1 (Momentum) | 02 | PR A (Ch. 6) |
| #4 (Cyclotron) | 03 | PR B or PR C (Ch. 11 or 12) |
| #3 (Radiation reaction) | 05 | PR E (Ch. 16) + issue #43 |
| #5 (Larmor) | 06 | PR D (Ch. 14) |
| #2 (Liénard–Wiechert) | 06 | PR D (Ch. 14) + issue #43 |

Episodes 02 and 03 can therefore be revised after PR A and (B or C) respectively, well ahead of the radiation-reaction chapters. Episodes 05 and 06 are gated on the full radiation-reaction work landing in PRs D–E and on issue #43's comparison document.

### 12.3 Selection rationale (problems considered and dropped)

For traceability, problems considered for elevation but not selected:

- **Jackson 3e §6.1 (magnetic monopoles).** Historically rich (Dirac 1931) but the proper-time reformulation is trivial — monopoles preserve electric–magnetic duality under `c → b`. No experimental hook for the Experimentalist.
- **Jackson 3e §6.11 (symmetric stress tensor).** Too formal; Experimentalist has nothing to say.
- **Jackson 3e Ch. 11 (special relativity).** Too kinematic; proper-time SR agrees with classical SR on every observable. The Physicist's "and here's where they differ" beat doesn't land.

---

## 13. Devil's-advocate review and what we cannot honestly fix

Before implementation begins, this section captures the strongest objections to the campaign and **what the plan does and does not solve about each**. Snake-oil avoidance is the organising principle: where there is no honest mitigation, we say so rather than dress up a partial response as a fix.

Vocabulary used below:
- **Partial response** — something the plan does that helps but does not fully solve the objection. We list what remains unsolved.
- **Honest answer** — no real fix exists; the campaign accepts the constraint.
- **Open risk** — neither solved nor accepted; we don't know yet how bad it is. PR A will tell us.

### 13.1 Objections with partial responses

**O1. The anchor source has an unresolved flagged finding (Maxwell-paper Eq. 24).** [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) flags Eq. 24 as missing a factor of `c` and a `V²/(2mc²)` term.

*Partial response:* the campaign uses Eqs. (1)–(23) only (all Wolfram-verified ✅). Any Jackson problem whose proper-time derivation would invoke Eq. 24 is **deferred** with a `⚠ blocked on Eq. 24` note. The Eq.-dependency mapping is recorded in `_proper_time_cheatsheet.md` so the filter is mechanical, not subjective.

*What this does not solve:* we do not know whether Eq. 24's flaw is a typo or a sign of a deeper issue in the framework. If the latter, the verified equations (1)–(23) may rest on assumptions that quietly break at Eq. 24 — meaning the campaign would be building on a foundation with subtle issues not yet identified. Deferral kicks the can; it does not address the underlying uncertainty.

**O2. "No new physics claims" (§1) is in tension with §12 headline language.**

*Partial response:* all dual-theory results are framed conditionally — "*if* the Gill–Zachary framework is the correct formulation, *then* X." §12.1 briefs for #2 and #3 have been rewritten accordingly.

*What this does not solve:* the act of spending 6–12 months working problems in the framework is itself a signal that the framework is worth this much effort. Even the most careful conditional language cannot cancel that signal. A skeptical reviewer is entitled to read the campaign's existence as soft advocacy regardless of how each individual document is phrased.

**O3. Crocco AI-use classification at scale is genuinely hard.**

*Partial response:* [§6](#6-source-material-treatment-copyright--crocco-compliance) broadened the substantive category to include problem selection, narrative framing, and choice-of-classical-comparator. Per-problem documents require per-paragraph `<!-- TODO -->` blocks for each substantive claim.

*What this does not solve:* once "substantive" is honestly broadened, almost every interpretive sentence in the campaign falls inside it. Per-paragraph TODO blocks will be so frequent that reviewers stop seeing them — disclosure becomes noise. A campaign-level disclosure (one statement per chapter) might be more honest but loses granularity. **We do not yet know which is better.** PR A is the test, and the decision should be revisited at PR B.

**O6. "Three unit systems" oversells what most of the campaign delivers.**

*Partial response:* §1 paragraph 2 and §4 acknowledge that PRs B–E are mostly two-system. Issue title and root README copy should be calibrated when the campaign starts.

*What this does not solve:* casual readers will skim "Electromagnetism / Jackson / CGS / SI / proper-time" and form expectations the campaign cannot meet on every PR. We accept this as residual reader-expectation drift.

**O7. Duplication with Equation Verification.**

*Partial response:* the campaign's value-add over Equation Verification is explicitly (a) CGS↔SI training in problem form; (b) the experimental-comparison ladder (#43+); (c) the podcast-elevation gradient.

*What this does not solve:* if (a) is undercut by O9 (overlap with textbooks), (b) is undercut by O5 (experiments may not discriminate), and (c) is undercut by O8 (small audience), then the value-add reduces to "the author personally finds this useful." That may be a legitimate reason to do it, but it is not the framing the §1 "earns its own thread" justification implies. **If issue #43 produces inconclusive comparisons, the campaign's headline payoff is hollow and the §1 framing must be downgraded.**

### 13.2 Objections with no honest mitigation

**O4. Effort estimate (updated 2026-05-24 with decisions D3 + D1 + D5 + D7 baked in).** Honest per-problem cost under the confirmed decisions:

- Headline-payoff problems with full template (D3): ~3–5 days each.
- Fluency-builder problems with full template (D3): ~2–3 days each (shorter derivation, same documentation burden).
- Eq. 24-branched problems (D1): ~1.7× the single-branch cost — applies to a small number of Ch. 12 / spin-related problems, not the bulk of the campaign.
- PR A-prequel problems (D7): ~1.5–2 days each (no proper-time content exercised, just CGS↔SI + bibliography stubs).

Blended estimate for the full 50–100 problem scope:
- PR A-prequel: 4 problems × ~1.8 days ≈ ~1 week.
- PRs A–F+: 50–100 problems × ~2.7 days blended ≈ 135–270 person-days.
- Branching overhead on the ~5–10 Ch. 12-area problems that touch Eq. 24: ~10–20 extra person-days.

**Total: ~140–290 person-days**, or **~6–14 months at a 30–50%-time hobby budget** for the full 50–100 problem scope.

*Honest answer:* this is a 140–290 person-day campaign against a part-time hobby budget. Realistic completion: **9–18 months for the full scope** under the confirmed decisions (full template throughout per D3, Eq. 24 branching per D1, no graceful degradation). Per [§13.5 D5](#135-decision-points--confirmed-by-author-2026-05-24) the campaign has **no kill switch** — the author has explicitly chosen quality over time-pressure. The trade-off is owned: if the author's available hours evaporate (job change, family events, other research priorities), the campaign may sit at 30–50% complete indefinitely. PRs E and F+ may never land. Issue #42 can still close on its [§11](#11-acceptance-criteria-for-closing-issue-42) criteria (PR A + (PR D OR PR E)) — close at ~20–35% of total scope — but the broader campaign is open-ended by design. No defensive language: the choice was D3 + D5 in full knowledge of the cost, and the cost is documented here.

**O5. The headline radiation-reaction comparison may be experimentally indistinguishable.** Issue #43 may show that proper-time predictions are within the uncertainty bands of one or both classical limits (Landau–Lifshitz, quantum-corrected L-L).

*Honest answer:* the campaign cannot generate new data. Its experimental value depends entirely on what the 2018 papers (and their successors at ELI / CALA / Apollon) can resolve, which is outside the author's control. **We may end with three carefully-derived conditional predictions and no discriminating experiment for any of them.** This is the most likely failure mode of the headline narrative, and the plan has no way to pre-empt it.

**O8. There is no audience clamoring for this work.** Graduate students do not seek out canonical-problem worked solutions from random repos. Equation Verification has DRQM-I co-authors as audience. History has the Crocco methodology and podcast as audience. This campaign has the author + a small dual-theory-curious community + the podcast-elevation pipeline.

*Honest answer:* this is a research workbook by the author for the author, with documentation discipline strong enough that other people can read it if they want to. Pedagogical impact on grad students is hoped-for, not load-bearing. We are not claiming otherwise.

**O9. CGS↔SI training overlaps with Jackson 3e Appendix and Griffiths Ch. 7.**

*Honest answer:* the campaign's CGS↔SI contribution is real but smaller than §1's "its own pedagogical contribution" implies. The distinguishing value-add is that translations are worked alongside the same problem's proper-time derivation. That is a real but narrow contribution; we are not claiming to replace Jackson's appendix.

**O10. The framework could be wrong.**

*Partial mitigation in O2:* conditional framing means predictions are recorded as "the published Gill–Zachary framework predicts X."

*Honest answer:* if the framework proves inconsistent, the campaign documents become *records of what the framework predicted before its inconsistency was found* — useful as historical record, not as validation. **Full mitigation impossible without abandoning the campaign.**

**O11. The author has a conflict of interest.** Trey Morris is a co-author of *Dual Relativistic Quantum Mechanics I* (2021). The framework being explored here is one he has co-authored prior papers about.

*Honest answer:* the campaign cannot be neutral. The author has skin in the game on whether the framework is correct. Claude pairs with the author on derivations and bookkeeping but does not independently validate the framework's correctness; conditional framing is the only honest position Claude can take. **No amount of careful documentation makes this a "neutral evaluation" of the framework**, and we should not market it as one. Readers should be aware of the author's prior involvement when reading any per-problem document.

### 13.3 What this campaign is, plainly

To prevent the framing from drifting toward overselling over the campaign's lifetime:

- **What it is:** a workbook in which the author works Jackson canonical problems in the Gill–Zachary proper-time framework, with CGS / SI translations where applicable. Output is conditional predictions of the form "if the framework is correct, this is what it says about Problem X."
- **What it is not:** a validation of the framework. A discovery of new physics. A pedagogical resource designed for graduate students. A finished textbook. Independent of the author's prior involvement in DRQM I.
- **What success looks like:** the conditional predictions are *correctly derived*, *experimentally-comparable where data exists*, and *clearly framed as conditional* to a reader who is skeptical of the framework. Success **does not** require the predictions to match data, nor the framework to be correct.
- **What failure looks like:** (i) the campaign stalls before PR D/E lands; (ii) the predictions are misframed as validations; (iii) issue #43's experimental comparisons cannot distinguish the framework from the classical limits and the "headline payoff" becomes functionally vacuous. We accept that **all three failure modes are possible** and that the §11 acceptance criteria do not protect against (ii) or (iii).

If a reader walks away from any per-problem document believing the campaign has proved something new about classical EM, the document has failed. The framing-principle disclaimer applies *recursively* — to every problem, every brief, every podcast script.

### 13.4 Open risks (we don't yet know how bad they are)

- **Per-paragraph Crocco disclosure may prove unusable.** See O3. PR A will tell us; reassess at PR B.
- **The author's available hours.** See O4. PR A's completion time is the data point.
- **Whether issue #43 produces discriminating predictions.** See O5. This is the campaign's existential test; without it, the whole headline narrative collapses.
- **Whether deferring Eq. 24-touching problems leaves a coherent campaign.** See O1. If too many Jackson problems route through Eq. 24, the deferral list grows past the point of usability.

### 13.5 Decision points — confirmed by author 2026-05-24

| # | Decision | Confirmed |
|---|---|---|
| 1 | Eq. 24 dependency | **Branched treatment.** Any Jackson problem whose proper-time derivation invokes Maxwell-paper Eq. 24 is worked twice — once with Eq. 24 as-published, once with the proposed correction (missing factor of `c`, plus the `V²/(2mc²)` term). Both derivations appear in the per-problem document under a "Branched on Eq. 24 finding" subsection, with `⚠ branched` flag. See [§5.1](#51-branched-treatment-for-eq-24-touching-problems) for the workflow. |
| 2 | Crocco substantive-disclosure granularity | **Per-paragraph `<!-- TODO -->` blocks for Ch. 6.** Reassess at PR B based on whether the blocks become background noise. |
| 3 | Template treatment per problem | **Always full template.** No "extended set" degradation. Every problem — headline-payoff and fluency-builder alike — gets the full three-section treatment, MCP verification, Wolfram notebook, and Comparison table. The §7.1 headline/fluency distinction now refers to *derivation difficulty*, not template completeness. |
| 4 | Issue-#42 close trigger | **PR A + (PR D OR PR E).** Close cleanly once one headline radiation chapter is done. |
| 5 | Campaign-kill switch | **No kill switch.** The campaign closes only on §11 acceptance criteria; there is no enforced cancellation. Author accepts the trade-off: higher-quality output (combined with D3 above), but the campaign may drift indefinitely if available hours evaporate. See [§13.2 O4](#132-objections-with-no-honest-mitigation) for the honest framing of this risk. |
| 6 | §1 framing downgrade on inconclusive #43 | **Auto-downgrade.** If issue #43 produces inconclusive experimental comparisons (proper-time predictions within uncertainty bands of both classical limits), §1 prose downgrades to "research workbook for the author and DRQM-I co-authors" and the root / `Roadmapping/README.md` cross-reference is amended. No manual review required. |
| 7 | PR A-prequel (fluency warm-up) | **Adopted.** Insert a PR A-prequel of 3–4 easy electrostatic/magnetostatic problems from Chs. 1–2 or 5 before tackling Ch. 6. Purpose: build muscle memory on problems where proper-time reduces to classical exactly, and validate the [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) read-aloud test on low-stakes content. See [§7.3](#73-pr-a-prequel-adopted) for content. |

### 13.6 What this section is not

This is not a list of every conceivable objection — only those that materially shape the plan, warrant explicit acknowledgment as residual risk, or test the campaign's honest framing. It is not a guarantee that the partial responses are sized correctly; PR A's outcome should re-open §13.5. If a reader of this section concludes the campaign should be cancelled rather than launched, that is a legitimate response and the author should consider it on its merits.

---

## 14. Voice and style: matching Gill's published prose

Per-problem documents, the proper-time cheat sheet, and the Physicist voice in elevated podcast scripts are continuous with the Gill corpus and should read as such. Voice continuity is not a cosmetic preference: it signals that the campaign's documents are *work in the framework*, not *commentary about it*. Drift into a third-party voice undermines that signalling and is also one of the easiest ways for AI-generated prose to be spotted as such even when the technical content is correct.

The full style guide lives in the repo at [`Roadmapping/Tooling/VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) — version-controlled, reviewable, and pairing with `CROCCO_COMPLIANCE.md` as a meta-discipline document. This section summarises the load-bearing parts and flags existing plan content that does *not* match Gill's voice and must be revised before elevation into actual documents.

### 14.1 Scope of the style constraint

| Document type | Match Gill's voice? |
|---|---|
| Per-problem documents under `Roadmapping/Electromagnetism/Jackson/Ch*.md` | Yes |
| `_proper_time_cheatsheet.md` | Yes |
| Companion `.wl` notebook header comments | Yes |
| Podcast-script Physicist-voice lines | Yes |
| `Roadmapping/Electromagnetism/README.md` status tables and headers | Short-form, but adopt Gill's diction (no hyperbole) |
| GitHub issue bodies, PR descriptions | No — operational tone |
| `.dev/tasks/` planning docs (including this plan) | No |
| §13 devil's-advocate sections | No — plain operational tone preserves their sharpness |
| TODOs and internal notes | No |

### 14.2 Voice markers to use

Drawn from sampling Gill's prose across [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], [[Dual_Relativistic_Quantum_Mechanics_I]], [[FoundationsII-Classical]], and [[The_Classical_Electron_Problem]]:

- **First person:** "we" as plural of modesty throughout. Never "I". Never address the reader as "you"; use "one can construct", "it is easy to show", "one observes".
- **Connective openers:** "Thus,", "Indeed,", "Furthermore,", "However,", "On the other hand,", "It is clear that…", "It is easy to show that…", "We observe that…", "It follows that…", "In order to see…".
- **Sentence rhythm:** long subordinate-clause sentences linked by semicolons, not fragmentary punchy prose.
- **Attribution:** narrative, by name and year, embedded in prose. "Wheeler and Feynman [13] showed that…", "Liénard, in 1898, computed…". Footnote-style citations supplement but do not replace.
- **Equation→interpretation hand-off:** *separate sentence* after the equation, never inline. Name the mechanism (e.g., "radiation reaction") with attribution to prior work.
- **Reformative framing of classical theory:** "reformulation", "alternative", "extension", "parallel image", "recovers lost physical insight". Acknowledge the success of the classical formulation before contrasting.
- **Hedges:** "to the extent that", "it is not obvious that", "under the assumption that", "if one accepts", "it is not clear how", "one should be able to show".
- **Conceptual anchor phrase:** "*mathematically equivalent but not physically equivalent*" is a Gill signature; reuse it verbatim whenever the contrast between classical and proper-time formulations is being introduced.
- **Parenthetical asides** for clarification: "(in c.g.s. units)", "(what we mean by) physically equivalent". Used sparingly but characteristically.

### 14.3 Anti-patterns

Phrases and rhetorical moves that **do not appear** in Gill's prose and must be removed when content is elevated into per-problem documents:

- Hyperbole: "headline payoff", "smoking gun", "dissolves a 120-year-old pathology", "the framework's superiority", "campaign's emotional arc".
- Casual / idiomatic English: "let's", "the cool thing about", "kind of", "pretty much", "gotcha".
- Markdown-heavy bullet structures for substantive content; Gill writes paragraphs.
- Bold or italic emphasis for rhetorical force; reserve for genuine technical terms.
- "Claude voice" tells: "I'll", "let me", "we should consider", "interestingly", "it's worth noting".

### 14.4 Compliance check

Before any per-problem document is committed, the author (or Claude, when authoring) reads each prose paragraph aloud. If a paragraph could plausibly appear in a Gill paper — same register, same hedging, same connective rhythm — it passes. If it sounds like a third-party summary of a Gill paper, it is revised.

When in doubt, find the closest passage in Gill's corpus that does the same rhetorical work (transition, equation-to-interpretation hand-off, attribution) and pattern-match against it.

### 14.5 Existing plan content that does not match Gill's voice

The §12.1 podcast-elevation briefs use phrases like "headline payoff", "smoking gun", "dissolves a pathology", and "campaign's emotional arc". **These are anti-Gill and were written in this plan deliberately**, since (per §14.1) this is a planning document, not a per-problem document. However, the briefs are scaffolding for content that will be elevated into per-problem documents and podcast scripts. When that elevation happens (PRs C–E), the content must be rewritten to Gill's voice. Specifically:

- Per-problem documents in Chs. 12, 14, 16 must derive the *technical content* of #2, #3, #4, #5 in Gill's voice, drawing on §12.1 only as a scoping note.
- Podcast Physicist-voice lines based on §12 must be rewritten from scratch when they enter `Roadmapping/History/Podcast/episode_*.md`. The Historian and Experimentalist voices follow the podcast README's separate style guidance; only the Physicist voice is bound by §14.

This is logged here so the elevation step does not silently inherit the plan's project-management voice.
