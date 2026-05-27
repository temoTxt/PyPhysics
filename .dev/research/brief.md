# Overnight research brief — Candidate 2 (issue #65)

**Branch:** `65-theory-candidate-2-variational-determination-of-r_e-via-renormalised-dual-dirac-equation`
**Issue:** https://github.com/temoTxt/PyPhysics/issues/65
**Master:** #67 — https://github.com/temoTxt/PyPhysics/issues/67
**Triangulated target:** $r_e/r_0 = 0.499\,420\,509\,912\,831\,7 \pm 2.5\times10^{-13}$ (PR #62)
**Closed-form Schwinger reference:** $(2 - \alpha/(2\pi))/(4 + \alpha/\pi) = 0.499\,419\,632\,156\ldots$ ($\sim 10^{-6}$ match)

## What this branch is doing

Fix $r_e/r_0$ by demanding the renormalised dual-Dirac equation reproduces the physical electron mass eigenvalue at the cutoff: $E(r_e) = m_e c^2$ plus framework-internal binding/self-energy contributions. Solve the coupled eigenvalue / renormalisation-condition system; identify additional closure conditions framework-internally; cross-check the resulting $r_e/r_0$ against the triangulated value and the Schwinger closed-form. See issue #65 for full method and acceptance criteria — this brief does not duplicate them.

## Source-of-record (read in order on first iteration)

1. `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` — especially §II (dual Dirac equation setup) and §III (Dirac eigenvalue problem, Eqs. III.1–III.23); §III.D introduces the cutoff $r_e$ as a free parameter
2. `Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md` — companion summary doc, Candidate 2 section
3. `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` Finding 2 — verdict to be updated
4. `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S3.wl` — `.wl` style template
5. `CLAUDE.md` — repo conventions (Gaussian units, Wolfram MCP gotchas, `uv run` only)

## Iteration protocol

Each `/loop` iteration:

1. Read `STATE.md` to see where the last iteration stopped.
2. Identify the next pending step. On iteration 0, that is "read source-of-record §1, record the eigenvalue-problem structure with $r_e$ as a parameter".
3. Execute **one** substantive step. Typical step shapes:
   - Read a primary source not yet read; record key equations and renormalisation-relevant identities in STATE.md.
   - Set up or extend a Wolfram MCP cell in `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation_variational.wl` (`.wl` style).
   - Execute / debug a Wolfram cell.
   - Enumerate candidate closure conditions (gauge invariance at cutoff, current conservation at radial boundary, magnetic-moment relation) and assess each for framework-internal motivation vs ad-hoc.
   - Solve the coupled eigenvalue / renormalisation-condition system symbolically.
   - Draft a section of `Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` §III.D append ("First-principles derivation of $r_e$ — variational route").
   - Cross-check resulting $r_e/r_0$ against triangulated $0.4994205099128317$ and Schwinger $(2-\alpha/(2\pi))/(4+\alpha/\pi)$.
4. Append a dated entry to `STATE.md` (ISO timestamp) summarising what advanced + what is queued + outcome-matrix branch per issue #65 / master #67.
5. `git add` the changed files explicitly (no `git add -A`). New commit only (never `--amend`). Push to this branch.
6. **Self-pace next iteration via ScheduleWakeup** (no fixed interval):
   - 1200–1800s for substantive math/Wolfram work,
   - 60–270s only for quick continuations that don't cross a 5-min cache boundary.
7. Stop the loop (omit ScheduleWakeup) only when:
   - All acceptance criteria in issue #65 can be checked, **or**
   - A `BLOCKED: <Tepper input needed>` state is recorded — typically: which closure conditions are framework-internal vs ad-hoc, or the form of framework-specified binding/self-energy contributions at the cutoff (both flagged in issue #65 dependencies).

## Guardrails

- **Wolfram MCP gotchas** per CLAUDE.md: single-line code; avoid `V` or `e` as symbols (collide with Wolfram built-ins); non-commutative `Dot` handled explicitly.
- **Crocco compliance**: substantive AI use throughout; per-section `<!-- TODO: human reviews and fills in -->` blocks in every derivational paragraph + the Wolfram notebook header + the FINDINGS update.
- **Python via `uv run`** only (per CLAUDE.md). Do not invoke bare `python`/`pip`.
- **Commit hygiene**: never `--amend`; one commit per iteration; clear one-line message starting with `candidate-2:`.
- **Push to this branch only**; never to `main` or another candidate's branch.
- **Never run destructive git ops** (`reset --hard`, `push --force`, `branch -D`, `clean -f`). If a step would require one, halt and record `BLOCKED: <reason>` in STATE.md.
- **Do not open PRs or post issue comments** during overnight iterations; the human orchestrator (in `PyPyshics-Main`) handles GitHub-side updates after morning review.
- **Closure-condition honesty**: when proposing a closure condition, the STATE log must explicitly mark it as `framework-internal` (load-bearing motivation traced to a specific equation/identity in DRQM I or a companion verification doc) or `ad-hoc` (no such trace). Ad-hoc closures may be explored as scaffolding but cannot be the basis of a terminal outcome-A/B result.

## Done criteria (for halting the loop)

- All acceptance criteria checkboxes in issue #65 can be checked, **or**
- A `BLOCKED: <reason>` state is recorded that cannot be advanced without author input, **or**
- A definite $r_e/r_0$ value has been computed with framework-internal closure and the outcome-matrix branch (per master #67) determined.

## Outcome-matrix (per master #67)

The iteration's terminal record in STATE.md should classify the result as one of:

- **A** — Derivation reproduces $r_e \approx 0.4994205099 \cdot r_0$ at framework precision → Finding 2 candidate ✅.
- **B** — Derivation reproduces Schwinger closed-form $(2-\alpha/(2\pi))/(4+\alpha/\pi)$ at one-loop precision → confirms Schwinger encoding (cross-references Candidate 3); Finding 2 candidate ✅.
- **C** — Derivation reproduces a different definite value → new finding; Finding 2 → ⚠ with the new value recorded.
- **D** — Derivation intractable / Tepper-blocker → other candidates may still resolve.
