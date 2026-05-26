# Task 61: Triangulate optimal r_e/r_0 across 6 precision observables

## Objective

Find the single value of $r_e/r_0$ that best fits all six precision atomic-physics observables tabulated in [Roadmapping/Author_Reports/2026-05_interim_for_gill.md:71-78](../../Roadmapping/Author_Reports/2026-05_interim_for_gill.md#L71-L78) jointly, build a runnable Wolfram-Language notebook that performs the search, and update [Roadmapping/Equation_Verification/FINDINGS_for_author_review.md:51-87](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md#L51-L87) Finding 2 with the triangulated value plus per-observable residuals. The work follows Tepper Gill's 2026-05-25 author guidance that branches (b) and (c) are bracketing guides rather than theoretical predictions — the cutoff is a numerically-searched parameter.

## Background

**The guiding formula and the discrepancy that motivates the search.** [Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md:454](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md#L454) records DRQM I Eq. (III.21–22):

$$g_r = 2\!\left[1 - \frac{4 r_0}{2 r + r_0}\right]$$

evaluated at a cutoff $r_e$. [Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md:480](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md#L480) reports the sensitivity $dg_r/d(r_e/r_0) \approx 4.0046$ at the relevant cutoff, so an $r_e$ change at the $\sim 10^{-14}$ level moves $g_e$ by $\sim 10^{-13}$ (the experimental precision on $g_e$).

**Finding 2's current state.** [Roadmapping/Equation_Verification/FINDINGS_for_author_review.md:51-87](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md#L51-L87) records two values:
- branch (b) = $0.499857150068631$ (as-published in DRQM I §III.D; reproduces neither $g_e$ nor the other five observables to measurement precision)
- branch (c) $\approx 0.4994205099128318$ (inverse-solved to reproduce $g_e$; reproduces the other five observables only at the Bethe-estimate precision floor)

Per Tepper's 2026-05-25 guidance, neither is "the" answer — they are bracketing values from a uni-observable search that should be generalised to a joint fit over all six observables.

**The six observables and where their measurements live.** The table at [Roadmapping/Author_Reports/2026-05_interim_for_gill.md:71-78](../../Roadmapping/Author_Reports/2026-05_interim_for_gill.md#L71-L78) is the input:

| Observable | Measurement | Stated uncertainty |
|---|---|---|
| Electron $g_s$ | $-2.00231930\ldots$ | $\sim 10^{-12}$ |
| H $2P_{3/2}-2P_{1/2}$ | $10{,}969.13(10)$ MHz | $\sim 10^{-8}$ |
| H 1S hyperfine (21-cm) | $1{,}420.405\,751\,768(2)$ MHz | $\sim 10^{-12}$ |
| He ${}^3P_0-{}^3P_1$ | $29{,}616.952$ MHz | $\sim 10^{-9}$ |
| Positronium ortho-para | $203{,}389(2)$ MHz | $\sim 10^{-5}$ |
| Muonium hyperfine | $4{,}463.302\,776(51)$ MHz | $\sim 10^{-8}$ |

**Source for the per-observable formulas (sequencing constraint).** The per-observable derivations linking each measurement to $r_e/r_0$ live in the Bethe-Salpeter campaign content under `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/`. That directory **does not yet exist on `main`** — it lives on the open PR branch `origin/50-bethe-salpeter-precision-predictions` (PR [#53](https://github.com/temoTxt/PyPhysics/pull/53), state: OPEN at planning time). This issue's implementation should be sequenced after PR #53 merges, OR the implementer can rebase this branch onto `origin/50-bethe-salpeter-precision-predictions` to inherit those files. The notebook itself can be drafted from DRQM I Eqs. (III.21–22) plus standard QED structure for the five non-$g_e$ observables, but the in-repo cross-references in the notebook's docstring will not resolve until PR #53 is on main.

**Mathematica-notebook style precedent.** [Roadmapping/Mathematica_Notebooks/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.wl:1-44](../../Roadmapping/Mathematica_Notebooks/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.wl#L1-L44) establishes the style: a `(* ::Package:: *)` header, a docstring comment block citing the companion verification doc, conventions (units, vector handling), and per-cell `(* PENDING *)` markers as the MCP wires in. The new notebook follows the same shape with cells that run end-to-end under Wolfram MCP.

**Crocco compliance posture.** The choice of objective function (precision-weighted $\chi^2$ vs unweighted least-squares vs minimax), the choice of weighting (pure measurement variance vs measurement variance plus framework-precision-floor noise term), and the interpretation of the resulting residual structure are *substantive* AI decisions per [CLAUDE.md](../../CLAUDE.md) "Crocco compliance" §1. The notebook header and the Finding 2 update therefore carry per-section `<!-- TODO: human reviews and fills in -->` blocks; the mechanical scan/refine arithmetic is *pragmatic*.

## Implementation Plan

1. **Create the subdirectory and README stub.** Add `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/` with a one-paragraph `README.md` documenting the subdirectory's scope (Quantum-Mechanics-side `.wl` notebooks, companion to `Roadmapping/Quantum_Mechanics/`).

2. **Draft `r_e_triangulation.wl` skeleton.** Header block following the precedent at [Roadmapping/Mathematica_Notebooks/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.wl:1-22](../../Roadmapping/Mathematica_Notebooks/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.wl#L1-L22). Convention block declares: precision-weighted $\chi^2$ as default objective; measurement variances pulled from the table above; framework precision floor (Bethe-estimate) handled as an optional noise term for fine-structure / Lamb-derived inputs.

3. **Encode the six observables as functions of $r_e/r_0$.**
   - $g_s$: direct from $g_r = 2[1 - 4 r_0/(2 r + r_0)]$ at $r = r_e$.
   - The five spectroscopic observables: each is a sum of a leading dual-theory term linear in $g_r$ (or quadratic via $g_s^2$) plus standard QED contributions. The implementer transcribes the per-observable formulas from `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/` (available after PR [#53](https://github.com/temoTxt/PyPhysics/pull/53) merges) and reduces each to $f_i(r_e/r_0)$ with all other QED inputs frozen at CODATA values. Each cell is `(* PENDING *)` until tested under Wolfram MCP per the precedent.

4. **Define the joint objective.** Default form, expressed in single-line Wolfram per CLAUDE.md MCP gotcha rule 1:

    ```
    chi2[x_] := Sum[((fi[x] - mi)/sigmai)^2, {i, 1, 6}]
    ```

    The notebook records the alternative (unweighted, minimax) for comparison and documents *why* the chosen default is the default. This is substantive AI — the choice is recorded in a `<!-- TODO: human reviews and fills in -->` block.

5. **Run the search.** Two-pass: (a) `Plot[chi2[x], {x, 0.499, 0.5}]` for an over-the-bracket sweep (branch (b) and branch (c) bracket the minimum); (b) `FindMinimum[chi2[x], {x, x0}]` with $x_0$ at the visual minimum, refining to machine precision. Record best-fit $r_e/r_0$ and the value of $\chi^2$ at the optimum.

6. **Propagate uncertainty.** Compute the one-sigma uncertainty on the optimum via the Hessian: $\sigma_x = \sqrt{2 / (d^2\chi^2/dx^2)}_{x=x_{\rm opt}}$. Compare against the precision of the published $r_e$ (15 sig figs) and report whether the search-defined precision is consistent.

7. **Diagnostic table.** For each of the six observables, evaluate $f_i(x_{\rm opt})$ and report the residual in measurement-σ. Flag any observable with $|f_i - m_i|/\sigma_i > 3$ as in tension with the joint fit (a "stretched-fit" signal).

8. **Update Finding 2.** Append a new subsection "Finding 2 update — 2026-MM-DD: empirical triangulation per Tepper guidance (2026-05-25)" to [Roadmapping/Equation_Verification/FINDINGS_for_author_review.md](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md). The subsection records: Tepper's bracketing-guide guidance, the triangulated value, the one-sigma uncertainty, the per-observable residual table, and the implication for the campaign's verdicts on the six observables (now recorded at the triangulated value, not the bracketing values). Per Crocco §5, the subsection carries a `<!-- TODO: human reviews and fills in -->` block.

9. **Update DRQM I §III.D verdict.** Append an "Update 2026-MM-DD" block to the Eq. (III.21–23) section of [Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md:451-505](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md#L451-L505) pointing to the FINDINGS update. Verdict marker shifts from 🔴 to ⚠ (the discrepancy is now characterised, not unresolved); ✅ vs 🔴 is deferred to the first-principles rederivation in [#54](https://github.com/temoTxt/PyPhysics/issues/54).

10. **Draft the follow-up author note.** Create `Roadmapping/Author_Reports/2026-MM_followup_for_gill_re_triangulation.md` — a short (≤ 1 page) note summarising the triangulated value, the per-observable residuals, and one or two specific questions for Tepper if the residuals expose a stretched-fit signal. The note follows the [Roadmapping/Author_Reports/README.md](../../Roadmapping/Author_Reports/README.md) naming convention and the [Roadmapping/Author_Reports/build_report.sh](../../Roadmapping/Author_Reports/build_report.sh) pipeline (markdown → LaTeX → PDF) only if/when the next interim report is built.

11. **Crocco-pragmatic markers** in the notebook (the scan + arithmetic) and **Crocco-substantive markers** in the objective-choice rationale, the residual-interpretation cells, and the follow-up note's interpretive paragraphs. Each substantive block carries `<!-- TODO: human reviews and fills in -->` per [CLAUDE.md](../../CLAUDE.md) Crocco §5.

## Files to Modify

| File | Change |
|---|---|
| `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/README.md` | Create — one-paragraph subdirectory scope note |
| `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl` | Create — runnable triangulation notebook per implementation steps 2–7 |
| `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` | Append "Finding 2 update" subsection per implementation step 8 |
| `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` | Append "Update" block to the Eq. (III.21–23) section per implementation step 9 |
| `Roadmapping/Author_Reports/2026-MM_followup_for_gill_re_triangulation.md` | Create — follow-up author note per implementation step 10 (replace `MM` with actual month at write time) |

## Dependencies

- **Sequencing dependency on PR [#53](https://github.com/temoTxt/PyPhysics/pull/53).** The per-observable formulas for the five non-$g_e$ observables live on `origin/50-bethe-salpeter-precision-predictions`. Implementation either waits for #53 to merge or rebases onto that branch. No new external libraries.
- **Wolfram MCP** for running the notebook (already a documented system dependency per [CLAUDE.md](../../CLAUDE.md) "How equation verification works"). The three MCP gotchas (single-line code, `V` reserved as Vanadium, `e` resolving to Euler's number) apply per [CLAUDE.md](../../CLAUDE.md).
- **No new Python or system packages.**

## Acceptance Criteria

- [ ] `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl` exists, runs end-to-end under Wolfram MCP, and reports a best-fit $r_e/r_0$ with a one-sigma uncertainty derived from the objective's Hessian.
- [ ] The six predictions at the optimum are recorded in a diagnostic cell alongside their measurements with residuals expressed in measurement-σ.
- [ ] The diagnostic cell flags any observable whose residual exceeds $3\sigma$ as a stretched-fit signal.
- [ ] [Roadmapping/Equation_Verification/FINDINGS_for_author_review.md](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) Finding 2 carries an "Update" subsection recording Tepper's bracketing-guide guidance (2026-05-25), the triangulated value, the per-observable residual table, and the campaign-verdict implications.
- [ ] [Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) Eq. (III.21–23) section carries an "Update" block pointing to the FINDINGS update; the verdict marker shifts from 🔴 to ⚠ pending the [#54](https://github.com/temoTxt/PyPhysics/issues/54) first-principles rederivation.
- [ ] A follow-up author note exists under `Roadmapping/Author_Reports/` per the [Roadmapping/Author_Reports/README.md](../../Roadmapping/Author_Reports/README.md) naming convention, summarising the triangulation result for Tepper.
- [ ] The notebook header, the Finding 2 update, and the follow-up note each carry the appropriate Crocco-pragmatic-vs-substantive marker per [CLAUDE.md](../../CLAUDE.md) Crocco §5; substantive blocks carry `<!-- TODO: human reviews and fills in -->`.
- [ ] All Wolfram cells obey the three MCP gotchas in [CLAUDE.md](../../CLAUDE.md) (single-line code; no `V` or `e` as symbols; non-commutative `Dot` handled).

## Testing

Commands a reviewer should run:

```bash
# 1. Notebook integrity — header docstring, package marker, no syntax errors at the wolframscript level.
uv run python -c "
import re, pathlib
p = pathlib.Path('Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl')
s = p.read_text()
assert s.startswith('(* ::Package:: *)'), 'missing Package marker'
assert 'Companion to' in s, 'missing companion-to citation in header'
assert '0.499' in s, 'optimum bracket window not present'
print('notebook static checks: OK')
"

# 2. Run the notebook via Wolfram MCP (per the standard verification workflow in CLAUDE.md).
# The reviewer reads through cells, confirms each prints the expected diagnostic, and
# records the best-fit r_e/r_0 + one-sigma uncertainty.

# 3. FINDINGS doc grep — the update subsection is present and links to issue #61.
grep -n "Tepper" Roadmapping/Equation_Verification/FINDINGS_for_author_review.md | head
grep -n "#61\|issue/61\|issues/61" Roadmapping/Equation_Verification/FINDINGS_for_author_review.md | head

# 4. DRQM I doc grep — verdict marker shift from 🔴 to ⚠ at Eq. (III.21–23).
grep -n "III.21\|III.22\|III.23" Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md | head

# 5. Author-report follow-up note — file exists and follows naming convention.
ls Roadmapping/Author_Reports/2026-*_followup_for_gill_re_triangulation.md
```

Tests added: none in a test suite — this repository has no formal test framework (per [CLAUDE.md](../../CLAUDE.md) "What this repo is: ... There are no tests"). Validation is via the Wolfram MCP run in step 2 above and the doc-grep checks in steps 3–5.
