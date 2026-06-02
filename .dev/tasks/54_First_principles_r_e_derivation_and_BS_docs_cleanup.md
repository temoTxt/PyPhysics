# Task 54: First-principles r_e derivation + Bethe-Salpeter docs cleanup at the triangulated value

## Objective

Execute the two remaining strands of issue [#54](https://github.com/temoTxt/PyPhysics/issues/54) after PR [#62](https://github.com/temoTxt/PyPhysics/pull/62)'s empirical resolution of $r_e/r_0 = 0.4994205099128317$ (closes [#61](https://github.com/temoTxt/PyPhysics/issues/61)). Scope 2 (Bethe–Salpeter docs cleanup) is execute-now; Scope 1 (first-principles derivation) is non-urgent and awaits Tepper Gill's input on which framework-internal starting point is natural. Both scopes are tracked in one PR or split as the implementer judges appropriate; this plan documents both so neither is lost.

## Background

**Triangulated value supersedes the bracketing branches.** PR [#62](https://github.com/temoTxt/PyPhysics/pull/62) merged 2026-05-26 with [`Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl`](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl) (notebook), [`Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md`](../../Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md) (follow-up author note), and the Finding 2 update in [`Roadmapping/Equation_Verification/FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md). Tepper's 2026-05-25 author guidance: branches (b) = $0.499857150068631$ and (c) = $0.4994205099128318$ are bracketing guides from a uni-observable search against $g_s$; neither is "intended." The empirical joint fit across six $g_s$-dependent observables resolves to $r_e/r_0 = 0.4994205099128317 \pm 2.5\times10^{-13}$, which matches branch (c) to 16 sig figs.

**Scope 2 — branched-treatment markers are now overdetermined.** Five files under [`Roadmapping/Quantum_Mechanics/Bethe_Salpeter/`](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/) carry branched-treatment language tied to the (b) vs (c) framing:

- [`03_FineStructure.md`](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md) — 13 occurrences. Branched verdict block at [03_FineStructure.md:125-135](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md#L125-L135) for BS-§14.2 (H $2P_{3/2}-2P_{1/2}$).
- [`06_Hyperfine.md`](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/06_Hyperfine.md) — 13 occurrences. Branched verdict at [06_Hyperfine.md:124](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/06_Hyperfine.md#L124) for BS-§22.1 (H 1S hyperfine 21-cm).
- [`07_RadiationInteraction.md`](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/07_RadiationInteraction.md) — 5 occurrences. Verdict at [07_RadiationInteraction.md:113](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/07_RadiationInteraction.md#L113) for BS-§30 (M1 transitions).
- [`09_HeliumExcited.md`](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/09_HeliumExcited.md) — 24 occurrences. Branched verdict at [09_HeliumExcited.md:54](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/09_HeliumExcited.md#L54) for BS-§72 (He ${}^3P_J$) plus positronium ortho-para and muonium hyperfine sections.
- [`10_CrossComparison.md`](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md) — 12 occurrences. Campaign-closing cross-comparison table at [10_CrossComparison.md:25-33](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md#L25-L33) carries the ⚠ / ✅ branched column.

The branched-verdict block structure is established at [03_FineStructure.md:125-135](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md#L125-L135):

```
**Verdict (branched):**
- `(a)` leading Dirac: ✅ ...
- `(b)` as-published `r_e`: ⚠ ...
- `(c)` corrected `r_e`: ✅ ...
The campaign's verdict is conditional on which branch of `r_e` is the intended one. ...
```

The cleanup target structure (un-branched, at the triangulated value):

```
**Verdict:** ✅ at the framework's Bethe-estimate / sub-leading-QED precision floor, at the triangulated $r_e/r_0 = 0.4994205099128317$ (per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62)). The honest reading remains: the prediction is `(g_s/-2)^n × textbook leading-g_s term`, which matches measurement by construction once `g_s` is the back-fit value — a self-consistency check at leading-`g_s` precision, not an independent discrimination from standard QED (per [`10_CrossComparison.md` §2](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md#2-the-r_e-back-fit-self-consistency-across-six-g_s-dependent-observables)).
```

**Scope 1 — derivation companion docs.** The triangulation notebook ([`r_e_triangulation.wl`](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl)) and the `BetheSalpeter_S3.wl` companion ([`BetheSalpeter_S3.wl`](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S3.wl)) establish the `.wl` style: top docstring block citing companion verification doc, `ClearAll`/`Print`/`FullSimplify` per-cell, honest-scope discipline inline. The three candidate starting points (Section 1 below) are documented in the follow-up author note at [`2026-05_re_triangulation_followup_for_gill.md`](../../Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md). A fourth candidate originally considered — retrieving an unpublished working-notebook derivation step — is empty per Tepper's 2026-05-25 author guidance (the original cutoff was a numerical search, not a derivation). The derivation's companion verification doc would extend the Eq. (III.21–23) section of [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) or live as a sibling under `Roadmapping/Equation_Verification/`.

**FINDINGS Finding 2 current state.** Verdict marker is ⚠ (was 🔴) per the PR #62 update at [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md). The marker moves to ✅ on Scope 1 completion if the derivation agrees with the triangulated value; stays ⚠ (with the derivation result recorded) if a third refinement is exposed.

## Implementation Plan

The plan is split by scope. Scope 2 has no external dependencies and can ship first as a self-contained PR. Scope 1 is dependent on Tepper's input and may sit indefinitely.

### Scope 2 — BS-docs cleanup at the triangulated value (execute-now)

1. **For each of the four per-PR files** (`03_FineStructure.md`, `06_Hyperfine.md`, `07_RadiationInteraction.md`, `09_HeliumExcited.md`):
   - Replace the per-result `**Verdict (branched):**` block with an un-branched `**Verdict:**` block per the cleanup target structure above.
   - Keep the prediction tables intact (the branch (c) column values are the un-branched-verdict numerical predictions; do not re-do the numerics).
   - Drop branch-(b) prediction rows from the per-result prediction tables, leaving the framework prediction (= branch (c) values) and the measurement column. Where the table explicitly labels rows "Proper-time `(b)` as-published `r_e`" / "Proper-time `(c)` corrected `r_e`", collapse to a single "Proper-time prediction (at triangulated `r_e`)" row.
   - Rewrite prose that says "branch (b) / branch (c)" to remove the bracketing-guide framing — replace with references to "the triangulated $r_e/r_0$ value" and cross-references to [PR #62](https://github.com/temoTxt/PyPhysics/pull/62) and the [follow-up author note](../../Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md).
   - Preserve the honest-scope language ("back-fit self-consistency, not independent corroboration") — that framing is unchanged by the cleanup; it applies to the un-branched verdict just as it applied to the branched one.
   - Preserve existing per-result `<!-- TODO: human reviews and fills in -->` blocks and add a new short TODO to each modified verdict block confirming the un-branching is faithful.

2. **`10_CrossComparison.md` campaign-closing chapter** (most extensive change):
   - The §2 prediction table at [10_CrossComparison.md:25-33](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md#L25-L33) replaces the ⚠ / ✅ "branched" column with a single ✅ "verdict" column at the triangulated value. Numerical predictions in the table are the branch (c) values, unchanged.
   - The §2 prose ("Reframe note", "Honest reading", "flagged finding's resolution path") is rewritten: drop "branch (b) vs branch (c)" framing; record Tepper's bracketing-guide guidance and PR #62's triangulation as the empirical resolution; preserve the honest reading ("one back-fit applied six times, not six independent corroborations").
   - §1's PR-by-PR result inventory at [10_CrossComparison.md:25-33](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md#L25-L33) updates the six "⚠ / ✅ branched" verdict markers to ✅ with a `(triangulated)` annotation.
   - §4 closing observation: rewrite "whether the resolution of the `r_e` finding falls on branch (b) or branch (c)" to record that the resolution is at the triangulated value; the conditionality ("if branch (c) is intended ...") becomes unconditional at the triangulated value, modulo the Scope 1 first-principles thread.

3. **Crocco compliance for the cleanup pass.** The relabeling itself is *pragmatic* AI (mechanical reframing of verdict markers; no numerical predictions change). The interpretive prose updates in §2 and §4 of `10_CrossComparison.md` are *substantive* AI (the honest-reading framing is rewritten); each substantive update carries a `<!-- TODO: human reviews and fills in -->` block.

4. **FINDINGS doc consistency check.** No Finding 2 changes needed — the PR #62 update already records the triangulated value and the verdict shift 🔴 → ⚠. The Scope 2 cleanup pass adds a one-line reference from each modified BS-docs verdict block back to the FINDINGS Finding 2 entry, so the trail is bidirectional.

### Scope 1 — Candidates summary document for the first-principles derivation (immediate; derivation work itself is deferred and conditional)

**Conceptual reframe (2026-05-26):** the published $r_e$ value was a cutoff fixed by numerical search against $g_s$, not a derived quantity. PR #62's triangulation widened the search from one observable to a joint fit across six observables that share $(g_s/-2)^n \times \text{textbook}$ structure, which confirms the framework's structural prescription is self-consistent under a single cutoff but does not promote that cutoff to a derived quantity. There are no pre-existing derivation notebooks to wire up; a first-principles derivation would be genuinely new theory work, dependent on Tepper's input on which framework-internal route is natural. Scope 1's immediate deliverable is therefore a **candidates summary document** that lays out the three remaining candidate routes (with 1–3 paragraph explanations each) for Tepper to react to; the derivation work itself is deferred to a downstream task contingent on his selection. (A fourth candidate originally considered — retrieving an unpublished working-notebook derivation step — was ruled out by Tepper's 2026-05-25 author guidance that the original cutoff value was a numerical search alone.)

5. **Create the candidates summary document.** Draft `Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md`, following the Author_Reports naming convention. The document covers:
   - **Context** — the honest reading that the value remains a cutoff, not a derived quantity; the joint-fit-self-consistency reading per `10_CrossComparison.md §2` (one structural fact applied to six manifestations, not six independent constraints); the open derivation question.
   - **Candidate 1 — Proper-time self-energy integral.** 1–3 paragraphs covering the framework's analog of the one-loop electron self-energy diagram in the proper-time formulation; what the route would compute and what the output would look like.
   - **Candidate 2 — Variational determination via the renormalised dual-Dirac equation.** 1–3 paragraphs on fixing $r_e$ by a mass-condition at the cutoff, analogous to standard-QED renormalisation conditions on $\alpha(\mu)$; the additional consistency conditions that may be needed for closure.
   - **Candidate 3 — Closed-form Schwinger identification of the triangulated value.** 1–3 paragraphs on $r_e/r_0 = (2 - \alpha/(2\pi))/(4 + \alpha/\pi)$ reproducing the Schwinger one-loop anomalous moment $g_e = -2 - \alpha/\pi$ via the DRQM I cutoff formula. Recasts the derivation question from "what value should the cutoff have?" to "why does the framework specify this particular cutoff prescription?"
   - **Closing** — explicit honest framing: this thread is non-urgent; the triangulated value already serves at the framework's precision floor; we defer to Tepper on whether any of the three candidates is natural to pursue; the campaign continues with the triangulated value as the empirical cutoff regardless.

6. **Build the PDF deliverable.** Run `cd Roadmapping/Author_Reports && REPORT_DATE=2026-05-26 ./build_report.sh 2026-05_re_derivation_candidates_for_gill` to produce the `.tex` and `.pdf`. Target page count: 3–5 within the build script's [3, 7] guard.

7. **Cross-reference from the triangulation note and the issue.** Add a short pointer from [`Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md`](../../Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md)'s "Where this leaves the open question" section to the new candidates summary doc, so the trail from triangulation note → candidates summary is documented. (Optional: skip if the user prefers the two docs stand alone.)

8. **Deferred work — derivation notebook + companion verification doc + FINDINGS final update.** Conditional on Tepper selecting a candidate and confirming it is pursueable, the downstream task would:
   - Draft `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation_<candidate>.wl` following the `.wl` style precedent at [`BetheSalpeter_S3.wl`](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S3.wl) (top docstring block citing the companion verification doc + the chosen candidate; per-cell symbolic derivation + numerical evaluation; closing cell comparing derived $r_e/r_0$ against the triangulated $0.4994205099128317$).
   - Append a "First-principles derivation of $r_e$" subsection to [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md)'s Eq. (III.21–23) treatment.
   - Append a "Finding 2 final" subsection to [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) recording verdict ⚠ → ✅ (agreement with triangulation) or ⚠ stays (third refinement exposed).
   - Crocco compliance: the derivation work is *substantive* AI end-to-end; per-section `<!-- TODO: human reviews and fills in -->` blocks throughout.
   - This step is **out of scope** for the immediate Scope 1 deliverable and is recorded here only so the downstream task is documented; pick up as a new task when Tepper's selection lands.

## Files to Modify

| File | Change |
|---|---|
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md` | Scope 2 — replace branched verdict block for BS-§14.2 with un-branched ✅ at the triangulated value; rewrite prose to drop branch-(b)/(c) framing; cross-reference PR #62 + follow-up note |
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/06_Hyperfine.md` | Scope 2 — same treatment for BS-§22.1 (H 1S hyperfine 21-cm) |
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/07_RadiationInteraction.md` | Scope 2 — same treatment for BS-§30 (M1 transitions) |
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/09_HeliumExcited.md` | Scope 2 — same treatment for BS-§72 (He ${}^3P_J$), positronium ortho-para, muonium hyperfine |
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md` | Scope 2 — most extensive: §1 PR-by-PR inventory verdicts (six ⚠ / ✅ → ✅); §2 prediction table column collapse; §2 + §4 honest-reading prose rewrite; campaign-closing observation updated to record triangulated resolution |
| `Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md` | Scope 1 — create; candidates summary doc (1–3 paragraphs per candidate; markdown source-of-truth) |
| `Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.tex` | Scope 1 — create (derived); pandoc-built LaTeX from the markdown |
| `Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.pdf` | Scope 1 — create (derived); built PDF deliverable for Tepper, target 3–5 pages |
| `Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md` | Scope 1 — optional one-line cross-reference from "Where this leaves the open question" to the new candidates summary doc |

## Dependencies

- **Scope 2** has no external dependencies — executes against `main` post-PR #62 merge using files already present.
- **Scope 1** is dependent on Tepper Gill's input on which of the three candidate starting points is natural for the framework. No external libraries.
- **Wolfram MCP** for running the Scope 1 derivation notebook (already documented in [`CLAUDE.md`](../../CLAUDE.md)). The three MCP gotchas (single-line code; avoid `V` reserved for Vanadium; avoid `e` resolving to Euler's number; non-commutative `Dot`) apply.

## Acceptance Criteria

### Scope 2 (BS-docs cleanup)

- [ ] All five Bethe-Salpeter files (`03_FineStructure.md`, `06_Hyperfine.md`, `07_RadiationInteraction.md`, `09_HeliumExcited.md`, `10_CrossComparison.md`) carry un-branched ✅ verdicts at the triangulated value $r_e/r_0 = 0.4994205099128317$ in place of the prior ⚠ / ✅ branched markers.
- [ ] No numerical prediction values are changed in any prediction table — the branch (c) values are the un-branched values, unchanged.
- [ ] Each modified verdict block in the four per-PR files (`03`, `06`, `07`, `09`) carries a cross-reference to [PR #62](https://github.com/temoTxt/PyPhysics/pull/62) and to [`Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md`](../../Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md).
- [ ] `10_CrossComparison.md` §1 inventory's six "⚠ / ✅ branched" verdict markers are all replaced with ✅ + `(triangulated)` annotation; §2 prediction-table branch-(b)/(c) columns collapsed to a single verdict column; §2 + §4 prose rewritten to drop branch-(b)/(c) framing while preserving the honest-reading framing ("one back-fit applied six times").
- [ ] No `branched` or `branch (b)` or `branch (c)` strings remain in any of the five files after the cleanup (verify via `grep -c branched ... | grep -v ':0$'` returning empty).
- [ ] Each substantive prose update in `10_CrossComparison.md` carries a `<!-- TODO: human reviews and fills in -->` block.

### Scope 1 (candidates summary document; immediate)

- [ ] `Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md` exists, carries YAML frontmatter (title, author, date, subject), and covers all three remaining candidates with 1–3 paragraph explanations each (Candidate 4 ruled out per Tepper's 2026-05-25 author guidance).
- [ ] The document opens with explicit honest framing: the $r_e$ value remains a cutoff (not a derived quantity); the triangulation increased its credibility by widening the constraint set; no derivation notebooks exist at present.
- [ ] The document closes with deferral to Tepper on which candidate (if any) is natural for the framework; thread explicitly non-urgent; the triangulated value stands as the campaign's $r_e$ disposition regardless.
- [ ] `.tex` and `.pdf` built via `build_report.sh`; PDF page count within `[3, 7]`; PDF metadata reads correctly (Title from YAML, Author + Date pinned).
- [ ] Crocco compliance: the document is *substantive* AI (it makes claims about what each candidate would compute and what its plausibility appears to be); a `<!-- TODO: human reviews and fills in -->` block is present near the top or inline per the per-paragraph TODO discipline.
- [ ] Deferred-work step (step 8 above — derivation notebook + Eq. III.21–23 update + FINDINGS final) is **out of scope** for this PR and recorded only as the downstream task once Tepper's selection lands.

## Testing

Commands a reviewer should run:

```bash
# Scope 2 — cleanup pass invariants

# 1. No branched markers remain in the five BS files after cleanup
for f in 03_FineStructure 06_Hyperfine 07_RadiationInteraction 09_HeliumExcited 10_CrossComparison; do
  count=$(grep -cE 'branched|branch \(b\)|branch \(c\)' "Roadmapping/Quantum_Mechanics/Bethe_Salpeter/$f.md" || true)
  echo "$f: $count branched references (target: 0)"
done

# 2. Each modified BS file references PR #62
for f in 03_FineStructure 06_Hyperfine 07_RadiationInteraction 09_HeliumExcited 10_CrossComparison; do
  grep -l "PR #62\|pull/62" "Roadmapping/Quantum_Mechanics/Bethe_Salpeter/$f.md" || echo "MISSING: $f"
done

# 3. Triangulated value appears in each modified file
grep -l "0.4994205099128317" Roadmapping/Quantum_Mechanics/Bethe_Salpeter/*.md

# Scope 1 — candidates summary doc (immediate)

# 4. Candidates doc exists with YAML frontmatter
test -f Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md
head -6 Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md | grep -q '^title:' && echo 'YAML title present'

# 5. All three remaining candidates are present as section headings (Candidate 4 ruled out per Tepper's 2026-05-25 guidance)
grep -cE '^## Candidate [1-3]' Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md
# expected: 3

# 6. PDF builds and lands in the [3, 7] page guard
REPORT_DATE=2026-05-26 ./Roadmapping/Author_Reports/build_report.sh 2026-05_re_derivation_candidates_for_gill 2>&1 | tail -5
pdfinfo Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.pdf | grep -E '^(Title|Author|Pages):'

# 7. Honest-framing keywords present in the document (cutoff vs derived; triangulation credibility; non-urgent)
grep -cE 'cutoff|credibility|non-urgent|triangulat' Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md

# Scope 1 — deferred-work check (not run as part of this PR; documents what the future task would test)
# - r_e_derivation_<candidate>.wl static + Wolfram MCP run
# - Dual_Relativistic_Quantum_Mechanics_I.md "First-principles derivation" subsection present
# - FINDINGS Finding 2 "final" subsection present (verdict ⚠ → ✅ or ⚠ stays)
```

Tests added: none in a test suite — this repository has no formal test framework (per [`CLAUDE.md`](../../CLAUDE.md) "There are no tests"). Validation is via the grep checks (Scope 2) and the Wolfram MCP run (Scope 1).
