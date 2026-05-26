# Overnight research brief — Candidate 3 (issue #66)

**Branch:** `66-theory-candidate-3-closed-form-schwinger-identification-of-the-triangulated-r_e`
**Issue:** https://github.com/temoTxt/PyPhysics/issues/66
**Master:** #67 — https://github.com/temoTxt/PyPhysics/issues/67
**Triangulated target:** $r_e/r_0 = 0.499\,420\,509\,912\,831\,7 \pm 2.5\times10^{-13}$ (PR #62)
**Closed-form Schwinger:** $(2 - \alpha/(2\pi))/(4 + \alpha/\pi) = 0.499\,419\,632\,156\ldots$ ($\sim 10^{-6}$ match)

## What this branch is doing

Determine whether the framework's cutoff prescription reproduces the Schwinger closed-form $r_e/r_0 = (2 - \alpha/(2\pi))/(4 + \alpha/\pi)$ as a derived identity (encoding the Schwinger one-loop QED anomalous moment $g_e = -2 - \alpha/\pi$), or whether the $\sim 10^{-6}$ match to the triangulated value is coincidental at the framework's precision floor.

Per issue #66, this candidate has **two paths**:

1. **Tepper-confirmation path (lowest cost)** — async; not actionable overnight. Recorded as a follow-up question, but the loop does not block on it.
2. **Empirical-test path (no author input required)** — extend the joint $\chi^2$ fit (the PR #62 triangulation, notebook `r_e_triangulation.wl`) with higher-precision observables (e.g., Penning-trap $g_e$ all-orders, deuterium hyperfine, improved positronium ortho-para), and test whether the residual to the closed-form $(2-\alpha/(2\pi))/(4+\alpha/\pi)$ systematically tracks the Karplus–Kroll two-loop QED corrections. **This is the overnight focus.**

See issue #66 for full method and acceptance criteria — this brief does not duplicate them.

## Source-of-record (read in order on first iteration)

1. `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl` — the PR #62 triangulation notebook; baseline to extend
2. `Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md` — companion summary doc, Candidate 3 section
3. `Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md` — follow-up author note
4. `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` §III.D (Eqs. III.18–III.23) — cutoff prescription
5. `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` Finding 2 — verdict to be updated
6. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md` — the cross-comparison observables already in the joint fit
7. `CLAUDE.md` — repo conventions (Gaussian units, Wolfram MCP gotchas, `uv run` only)

## Iteration protocol

Each `/loop` iteration:

1. Read `STATE.md` to see where the last iteration stopped.
2. Identify the next pending step. On iteration 0, that is "read source-of-record §1, record the six observables in the existing triangulation + their precisions".
3. Execute **one** substantive step. Typical step shapes:
   - Read a primary source not yet read; record key formulas and observable definitions in STATE.md.
   - Set up or extend a Wolfram MCP cell in `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_schwinger_residual_test.wl` (`.wl` style).
   - Identify a candidate higher-precision observable to add to the joint fit; document its measured value + uncertainty + framework-prediction formula.
   - Add the observable to the extended joint $\chi^2$; refit $r_e/r_0$.
   - Compute the residual from the extended fit to the closed-form $(2-\alpha/(2\pi))/(4+\alpha/\pi)$; test whether it tracks the Karplus–Kroll two-loop QED corrections at the framework's precision floor.
   - Draft a section of `Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` §III.D append ("Schwinger identification — empirical residual test").
   - Maintain a "Question for Tepper" stub in STATE.md so the human orchestrator can lift it into a #66 comment in the morning.
4. Append a dated entry to `STATE.md` (ISO timestamp) summarising what advanced + what is queued + outcome-matrix branch per issue #66 / master #67.
5. `git add` the changed files explicitly (no `git add -A`). New commit only (never `--amend`). Push to this branch.
6. **Self-pace next iteration via ScheduleWakeup** (no fixed interval):
   - 1200–1800s for substantive math/Wolfram/data-curation work,
   - 60–270s only for quick continuations that don't cross a 5-min cache boundary.
7. Stop the loop (omit ScheduleWakeup) only when:
   - All acceptance criteria for the empirical-test path in issue #66 can be checked, **or**
   - A `BLOCKED: <reason>` state is recorded — typically: an observable's framework-prediction formula not yet derived, or a data-curation block on measured-value uncertainties.

## Guardrails

- **Wolfram MCP gotchas** per CLAUDE.md: single-line code; avoid `V` or `e` as symbols (collide with Wolfram built-ins); non-commutative `Dot` handled explicitly.
- **Crocco compliance**: substantive AI use throughout; per-section `<!-- TODO: human reviews and fills in -->` blocks in every interpretation paragraph + the Wolfram notebook header + the FINDINGS update.
- **Python via `uv run`** only (per CLAUDE.md). Do not invoke bare `python`/`pip`.
- **Commit hygiene**: never `--amend`; one commit per iteration; clear one-line message starting with `candidate-3:`.
- **Push to this branch only**; never to `main` or another candidate's branch.
- **Never run destructive git ops** (`reset --hard`, `push --force`, `branch -D`, `clean -f`). If a step would require one, halt and record `BLOCKED: <reason>` in STATE.md.
- **Do not open PRs or post issue comments** during overnight iterations; the human orchestrator (in `PyPyshics-Main`) handles GitHub-side updates after morning review.
- **Measured-value provenance**: every observable added to the extended fit must record its measurement source (paper + year + uncertainty) in STATE.md. No bare numbers without provenance.
- **Distinguish "Karplus–Kroll consistency"** from "match to within precision": the empirical-test path requires the *residual* to track Karplus–Kroll, not just to be small. Record the residual prediction formula and the data residual side-by-side.

## Done criteria (for halting the loop)

- Acceptance criteria for the empirical-test path in issue #66 can be checked, **or**
- A `BLOCKED: <reason>` state is recorded that cannot be advanced without author input, **or**
- The empirical residual test has produced a definite verdict: Karplus–Kroll consistent (Schwinger encoding intentional) → outcome B, or inconsistent (encoding coincidental) → outcome C.

## Outcome-matrix (per master #67)

The iteration's terminal record in STATE.md should classify the result as one of:

- **A** — Empirical-test path refits to $r_e \approx 0.4994205099 \cdot r_0$ (no significant shift from PR #62 baseline) → triangulated value confirmed; Finding 2 candidate ✅.
- **B** — Empirical-test path shows residual to closed-form Schwinger tracks Karplus–Kroll two-loop corrections → Schwinger encoding intentional; Finding 2 candidate ✅.
- **C** — Empirical-test path shows residual *does not* track Karplus–Kroll → Schwinger match coincidental; Finding 2 → ⚠ with the residual structure recorded.
- **D** — Empirical-test path intractable (data unavailable / formula-blocker) → Tepper-confirmation path becomes the only route; record as the Question-for-Tepper.
