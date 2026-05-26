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

**Scope 1 — derivation companion docs.** The triangulation notebook ([`r_e_triangulation.wl`](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl)) and the `BetheSalpeter_S3.wl` companion ([`BetheSalpeter_S3.wl`](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S3.wl)) establish the `.wl` style: top docstring block citing companion verification doc, `ClearAll`/`Print`/`FullSimplify` per-cell, honest-scope discipline inline. The four candidate starting points (Section 1 below) are documented in the follow-up author note at [`2026-05_re_triangulation_followup_for_gill.md`](../../Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md). The derivation's companion verification doc would extend the Eq. (III.21–23) section of [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) or live as a sibling under `Roadmapping/Equation_Verification/`.

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

### Scope 1 — First-principles r_e derivation (non-urgent; awaits Tepper)

5. **Wait for Tepper's input** on which starting point is natural for the framework. The four candidates are documented in [`2026-05_re_triangulation_followup_for_gill.md`](../../Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md):
   - proper-time self-energy integral as the regulator scale at which mass renormalisation closes,
   - variational determination via the renormalised dual-Dirac equation at the cutoff,
   - structural constant of the dual representation (e.g., $b$-factor projection structure),
   - working-notebook derivation step not in the published prose.

6. **Once a starting point is identified**, draft `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation.wl` following the `.wl` style precedent at [`BetheSalpeter_S3.wl`](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S3.wl):
   - top docstring block citing the companion verification doc, the framework's renormalisation prescription, and the chosen starting point;
   - per-cell `ClearAll` / symbolic derivation / `Print` numerical evaluation with explicit honest-scope labels (algebraic identity vs derived numerical result vs conditional prediction);
   - closing cell that compares the derived numerical $r_e/r_0$ against the triangulated $0.4994205099128317$ and reports whether the derivation agrees within the triangulation's $\sigma_r = 2.5\times10^{-13}$, refines further, or exposes a third refinement.

7. **Companion verification doc.** Append a new subsection "First-principles derivation of $r_e$ — 2026-MM-DD" to [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md)'s Eq. (III.21–23) treatment, recording: the chosen starting point, the symbolic-derivation summary, the derived numerical $r_e/r_0$, the comparison against the triangulated value, and the verdict-marker disposition (⚠ → ✅ if agreement; ⚠ stays with the derivation result recorded if a third refinement is exposed).

8. **FINDINGS Finding 2 final update.** Append a "Finding 2 final — 2026-MM-DD" subsection to [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) recording the verdict-marker disposition per step 7's outcome.

9. **Crocco compliance for Scope 1.** The derivation is *substantive* AI end-to-end (it makes claims about what the framework's internal logic specifies for $r_e$). Per-section `<!-- TODO: human reviews and fills in -->` blocks throughout the notebook + verification doc + FINDINGS update.

## Files to Modify

| File | Change |
|---|---|
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md` | Scope 2 — replace branched verdict block for BS-§14.2 with un-branched ✅ at the triangulated value; rewrite prose to drop branch-(b)/(c) framing; cross-reference PR #62 + follow-up note |
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/06_Hyperfine.md` | Scope 2 — same treatment for BS-§22.1 (H 1S hyperfine 21-cm) |
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/07_RadiationInteraction.md` | Scope 2 — same treatment for BS-§30 (M1 transitions) |
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/09_HeliumExcited.md` | Scope 2 — same treatment for BS-§72 (He ${}^3P_J$), positronium ortho-para, muonium hyperfine |
| `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md` | Scope 2 — most extensive: §1 PR-by-PR inventory verdicts (six ⚠ / ✅ → ✅); §2 prediction table column collapse; §2 + §4 honest-reading prose rewrite; campaign-closing observation updated to record triangulated resolution |
| `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation.wl` | Scope 1 — create; first-principles derivation notebook (filename placeholder; rename to reflect chosen starting point once Tepper confirms) |
| `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` | Scope 1 — append "First-principles derivation of $r_e$" subsection to the Eq. (III.21–23) treatment |
| `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` | Scope 1 — append "Finding 2 final" subsection recording the verdict-marker disposition (⚠ → ✅ or stay ⚠) |

## Dependencies

- **Scope 2** has no external dependencies — executes against `main` post-PR #62 merge using files already present.
- **Scope 1** is dependent on Tepper Gill's input on which of the four candidate starting points is natural for the framework. No external libraries.
- **Wolfram MCP** for running the Scope 1 derivation notebook (already documented in [`CLAUDE.md`](../../CLAUDE.md)). The three MCP gotchas (single-line code; avoid `V` reserved for Vanadium; avoid `e` resolving to Euler's number; non-commutative `Dot`) apply.

## Acceptance Criteria

### Scope 2 (BS-docs cleanup)

- [ ] All five Bethe-Salpeter files (`03_FineStructure.md`, `06_Hyperfine.md`, `07_RadiationInteraction.md`, `09_HeliumExcited.md`, `10_CrossComparison.md`) carry un-branched ✅ verdicts at the triangulated value $r_e/r_0 = 0.4994205099128317$ in place of the prior ⚠ / ✅ branched markers.
- [ ] No numerical prediction values are changed in any prediction table — the branch (c) values are the un-branched values, unchanged.
- [ ] Each modified verdict block in the four per-PR files (`03`, `06`, `07`, `09`) carries a cross-reference to [PR #62](https://github.com/temoTxt/PyPhysics/pull/62) and to [`Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md`](../../Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md).
- [ ] `10_CrossComparison.md` §1 inventory's six "⚠ / ✅ branched" verdict markers are all replaced with ✅ + `(triangulated)` annotation; §2 prediction-table branch-(b)/(c) columns collapsed to a single verdict column; §2 + §4 prose rewritten to drop branch-(b)/(c) framing while preserving the honest-reading framing ("one back-fit applied six times").
- [ ] No `branched` or `branch (b)` or `branch (c)` strings remain in any of the five files after the cleanup (verify via `grep -c branched ... | grep -v ':0$'` returning empty).
- [ ] Each substantive prose update in `10_CrossComparison.md` carries a `<!-- TODO: human reviews and fills in -->` block.

### Scope 1 (first-principles derivation; non-urgent)

- [ ] `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation.wl` (or per-chosen-starting-point filename) exists, runs end-to-end via Wolfram MCP, derives a numerical $r_e/r_0$ from the chosen framework-internal starting point, and reports the comparison against the triangulated $0.4994205099128317$.
- [ ] [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md)'s Eq. (III.21–23) treatment carries the new "First-principles derivation" subsection with: the chosen starting point, the symbolic-derivation summary, the derived numerical value, and the verdict disposition.
- [ ] [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) Finding 2 carries the "Finding 2 final" subsection recording verdict ⚠ → ✅ (agreement with triangulation) or ⚠ stays (third refinement exposed).
- [ ] Crocco compliance: the derivation is *substantive* AI; per-section `<!-- TODO: human reviews and fills in -->` blocks throughout the notebook + verification doc + FINDINGS update.
- [ ] All Wolfram cells obey the three MCP gotchas in [`CLAUDE.md`](../../CLAUDE.md) (single-line; no `V` or `e` as symbols; non-commutative `Dot` handled).

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

# Scope 1 — derivation notebook (after Tepper input)

# 4. Notebook static check
uv run python -c "
import pathlib
p = pathlib.Path('Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation.wl')
if not p.exists(): raise SystemExit('not yet authored (Scope 1 pending Tepper input)')
s = p.read_text()
assert s.startswith('(* ::Package:: *)'), 'missing Package marker'
assert 'Companion to' in s, 'missing companion-to citation in header'
assert '0.4994205099128317' in s, 'missing comparison against triangulated value'
print('notebook static checks: OK')
"

# 5. Run the notebook via Wolfram MCP (standard verification workflow in CLAUDE.md).

# 6. FINDINGS doc shows the "Finding 2 final" subsection
grep -n 'Finding 2 final' Roadmapping/Equation_Verification/FINDINGS_for_author_review.md
```

Tests added: none in a test suite — this repository has no formal test framework (per [`CLAUDE.md`](../../CLAUDE.md) "There are no tests"). Validation is via the grep checks (Scope 2) and the Wolfram MCP run (Scope 1).
