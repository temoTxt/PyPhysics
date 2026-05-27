# Research brief — Li²⁺ Bethe–Salpeter Z-extension (issue #78)

**Branch:** `78-bethe-salpeter-z-extension-li2plus`
**Issue:** https://github.com/temoTxt/PyPhysics/issues/78
**Parent:** #50 Bethe–Salpeter campaign (closed; this extends it)
**Sibling:** #75 framework-specs author engagement (parallel, not gating)

## What this branch is doing

**Tepper's 2026-05-27 suggestion** — apply the dual-theory / proper-time framework to the single-electron Li²⁺ hydrogenic ion at $Z=3$, predicting four precision observables and comparing to measurement. The thread tests whether the framework's $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ cutoff (Z=1 triangulation, PR #62) is **universal in Z** or **Z-scaled** — the structural complement to PR #70's lepton-axis muon test (same physics question, Z-axis instead of $m_\ell$-axis).

See issue #78 for full goal, method, acceptance criteria — this brief does not duplicate them.

## Four observables to predict and verify (in order of headline-precision)

1. **Bound-electron g-factor in Li²⁺.** Sturm 2014 *Nature* **506**, 467: $g_e^{\rm bound}({}^7\text{Li}^{2+}) = 2.000\,025\,170\,7(10)$, fractional $\sim 5\times 10^{-10}$. Apply (III.22) with Z-scaled (or Z-universal) cutoff; compare. **Headline precision test.**
2. **2S₁/₂–2P₁/₂ Lamb shift in Li²⁺.** Schiffer 1995 *PRL* **74**, 2188: $\Delta E_{2S-2P}({}^7\text{Li}^{2+}) = 62\,765(21)$ MHz. Apply BS-§19 Bethe-estimate apparatus at Z=3. **Z⁴-enhanced precision-floor test.**
3. **Fine structure 2P₃/₂–2P₁/₂ in Li²⁺.** ~$7\,367$ MHz from Bayfield/Riis era. Apply BS-§14 at Z=3.
4. **Li-7 1s hyperfine splitting.** ~$12.7$ GHz from Beckmann 1974. Apply BS-§22 with Li-7 nuclear spin $I=3/2$.

## Source-of-record (read in order on first iteration)

1. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/05_LambShift.md` — BS-§19 Lamb-shift apparatus (Z=1 baseline; Z-scaling to be applied)
2. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md` — BS-§14 fine-structure apparatus (Z=1 baseline)
3. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/06_Hyperfine.md` — BS-§22 hyperfine apparatus (Z=1 baseline)
4. `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` — especially §III.D Eqs. (III.22)/(III.23) — note (III.23) gives muon/proton analogues at the *same* dimensionless $r_e/r_0$ but with $r_0^{\mu,p} = e^2/(m_{\mu,p}c^2)$; this is the closest precedent for Z-extension
5. `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl` — Pass A / Pass B joint-fit pattern to mirror for Z=3
6. `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` Finding 2 — verdict updates from #69/#70/#72 (just landed); this branch may add a Z-extension update

## Per-observable deliverable shape (modeled on BS-§N template)

For each of the four observables:

- **As measured:** value, uncertainty, source (paper + year + DOI when available).
- **QED Z-expansion:** standard textbook treatment, with the leading $(Z\alpha)^n$ terms identified.
- **Framework prediction:** apply the relevant BS-§N formula with Z=3 scaling. Two readings to compute and compare:
  - **(Z-i) Universal cutoff:** $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ (Z=1 triangulated value used directly).
  - **(Z-ii) Z-scaled cutoff:** if and only if a clear framework-internal Z-scaling emerges from the algebra; otherwise reported as not derivable.
- **Wolfram MCP check:** symbolic + numerical confirmation per CLAUDE.md gotchas (single-line cells, `ee` / `potV` / etc.).
- **Verdict:** ✅ at framework precision / ⚠ at framework precision floor / 🔴 disagreement.

## Iteration protocol

Each `/loop` iteration does **one** substantive step:

1. Read `STATE.md` to find the last iteration's pending-next.
2. Execute one step. Step shapes (each iteration picks one):
   - Read a source-of-record document not yet read; record key formulas + Z-scaling in STATE.md.
   - Set up or extend `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_joint_fit.wl` (new notebook this branch is creating).
   - Execute / debug a Wolfram cell.
   - Draft a per-observable result section in `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/11_Li2plus_HydrogenicIon.md` (new doc this branch is creating).
   - Compute the joint $\chi^2$ at Z=3 once all four predictions are in; report Li²⁺-optimal $r_e/r_0$ and compare to Z=1 triangulated.
3. Append a dated entry to `STATE.md` (ISO timestamp) recording: what advanced, what's queued next, current observable focus (#1/#2/#3/#4), and any BLOCKED state.
4. `git add` changed files explicitly. New commit only (no `--amend`). Push to this branch.
5. Self-pace next iteration via ScheduleWakeup (1200–1800s for substantive math; 60–270s for cache-warm continuations).
6. Stop the loop (omit ScheduleWakeup) only when:
   - All four observables have results + a Z-axis verdict is recorded in STATE.md, **or**
   - A hard BLOCKED state is recorded (e.g., a framework formula's Z-extension is genuinely ambiguous and needs Tepper input).

## Guardrails

- **Wolfram MCP gotchas** per CLAUDE.md: single-line cells; avoid `V` / `e` as symbols; non-commutative `Dot` handled.
- **Crocco compliance:** substantive AI; `<!-- TODO: human reviews and fills in -->` blocks in every interpretation paragraph + the Wolfram notebook header + the FINDINGS / per-observable doc updates. **Substantive AI tag on every Z-extension decision** (e.g., "we apply the BS-§14 formula at Z=3 using ... — substantive choice").
- **Python via `uv run`** only.
- **Commit hygiene:** new commits only; one commit per substantive iteration; clear one-line message starting with `li2plus:`.
- **Push to this branch only**; never to `main` or another branch.
- **Never run destructive git ops.** If a step would require one, halt and record `BLOCKED` in STATE.md.
- **Measurement provenance:** every observable value added must record its measurement source (paper + year + DOI/journal) in STATE.md and in the per-observable doc.
- **Z-i vs Z-ii honesty:** report both readings where they differ; do not silently pick one. The "which is intended" question is for Tepper (#75 sibling).
- **Cross-particle cross-check (echo of PR #70):** for the g-factor observable specifically, the framework's universal-vs-particle-specific cutoff was already constrained by the muon at $\sim 57{,}000\,\sigma_{a_\mu}$. The Li²⁺ test extends this along the Z axis; report whether Li²⁺ adds a new constraint.

## Done criteria

- All four per-observable result documents exist + Wolfram-MCP-verified.
- `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/11_Li2plus_HydrogenicIon.md` (or equivalent) integrates the four with a cross-comparison table.
- `r_e_Li2plus_joint_fit.wl` notebook reports the Z=3 joint optimum + residuals.
- A clear Z-axis verdict is recorded: cutoff Z-universal / Z-scaled / per-Z-back-fit.
- FINDINGS Finding 2 update appended (Z-extension result), with `<!-- TODO: human reviews and fills in -->` block.

## Outcome-matrix (for STATE.md classification per iteration)

- **A** — All four predictions match measurement at framework precision under the **same** $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ cutoff used for Z=1. Z-universal cutoff confirmed; framework's apparatus passes the Z-axis cross-consistency check at Z=3.
- **B** — Predictions match under a **Z-scaled** cutoff with a framework-internally-derivable form; Z-dependence characterised.
- **C** — Predictions match only under a **per-Z back-fit** (each Z requires its own cutoff with no derivable scaling); analogous to PR #70's lepton-axis verdict.
- **D** — Predictions cannot be made tractably under either reading; framework apparatus is structurally inadequate for Z>1. BLOCKED on Tepper input.
