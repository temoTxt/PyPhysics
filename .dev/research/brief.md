# Research brief — Li²⁺ spectroscopy (Lamb shift + fine structure) — issue #78

**Branch:** `78-li2plus-spectroscopy`
**Issue:** https://github.com/temoTxt/PyPhysics/issues/78
**Parent campaign:** [#50](https://github.com/temoTxt/PyPhysics/issues/50) Bethe–Salpeter (closed; #78 extends it)
**Parallel branches (do not overlap with these):**
- `78-bethe-salpeter-z-extension-li2plus` (running; owns observable #1 g-factor + integration / joint χ² + `r_e_Li2plus_joint_fit.wl`)
- `78-li2plus-hyperfine` (running; owns observable #4 hyperfine)
**This branch owns:** observables #2 (2S–2P Lamb shift) and #3 (2P₃/₂–2P₁/₂ fine structure)
**Output file (distinct):** `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/12_Li2plus_Spectroscopy.md`

## What this branch is doing

Predict the Li²⁺ Lamb shift and fine structure under the dual-theory framework, compare to measurement, and contribute Z=3 prediction values to the joint χ² fit owned by the parallel Self-Energy branch.

**Priority signal from Self-Energy iter-1** (recorded 2026-05-27): the **2S–2P Lamb shift is a weak discriminator** of the cutoff-universality question — $r_e$ enters only at sub-leading order via the anomalous-moment piece (per `Bethe_Salpeter/05_LambShift.md` BS-§20 line 114), so (Z-i) and (Z-ii) readings give identical predictions at the Bethe-estimate precision floor. **Therefore prioritize observable #3 fine structure** (which engages $r_e$ through the anomalous-g mechanism and is a real Z-axis discriminator) over observable #2 Lamb shift (which can be reported as ✅-by-inheritance trivially, no Z-axis-test value).

## Observable scope

### #2 — 2S₁/₂–2P₁/₂ Lamb shift in Li²⁺ *(weak discriminator; document briefly, do not over-invest)*

- **Measurement:** Schiffer 1995 *PRL* **74**, 2188: $\Delta E_{2S-2P}({}^7\text{Li}^{2+}) = 62\,765(21)$ MHz.
- **Z-scaling:** leading $\propto (Z\alpha)^4 m_e c^2 / n^3$; Bethe-log $\log(K/|E_m-E_n|)$ shrinks with Z (≈9.83 at Z=1 → ≈7.64 at Z=3). Net Z=3 prediction ≈ 59.3 × 1057.845 ≈ 62,729 MHz at framework's Bethe-estimate precision floor.
- **Framework apparatus:** apply BS-§19 / BS-§20 formula at Z=3 with the Z-scaled Bethe log; verdict will be ✅-by-inheritance at the framework's precision floor (which is ~250 MHz at Z=3, vs the 21 MHz measurement uncertainty — so the framework's precision is worse than measurement, meaning "agreement at framework precision" is the honest verdict).

### #3 — Fine structure 2P₃/₂–2P₁/₂ in Li²⁺ *(real Z-axis discriminator; primary focus)*

- **Measurement:** Riis et al. (1994) and earlier Bayfield era; current value $\approx 7\,367$ MHz with kHz–MHz precision depending on source.
- **Z-scaling:** Sommerfeld–Dirac leading term $\propto (Z\alpha)^4 m_e c^2$; Z=3 gives $81 \times$ the H value of 10,969 MHz divided by Bohr-scaling factors → ~7.4 GHz. Apply BS-§14 formula.
- **Anomalous-g coupling to $r_e$:** the fine structure splitting depends on $g_s$ (Thomas-precession reduction); under the framework's $(g_s/-2)^n \cdot \text{textbook}$ scaling with $n=1$ for fine structure (per `Bethe_Salpeter/10_CrossComparison.md §2`), $r_e$ enters the prediction through $g_s = g_r(r_e/r_0)$. **This is the Z-axis discriminator we want.**
- **Framework prediction under both readings:**
  - **(Z-i) Universal cutoff:** $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ (Z=1 triangulated). Predicts Z=3 FS using the universal cutoff with Z=3 anchor formula.
  - **(Z-ii) Z-scaled cutoff:** if a framework-internal Z-scaling emerges from §III.D-extension reasoning, apply it; otherwise report as not derivable.

## Source-of-record (read in order)

1. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md` — BS-§14 (Z=1 baseline + Z⁴ scaling identification) **← PRIMARY focus**
2. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/05_LambShift.md` — BS-§19/§20 (Z=1 baseline). Already partly read by Self-Energy iter-1; pull key Z-scaling identities only.
3. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md` §2 — the $(g_s/-2)^n \cdot \text{textbook}$ scaling pattern; confirm $n_{\rm FS} = 1$ for fine structure.
4. `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` §III.D Eqs. (III.22)/(III.23) — anomalous-g formula + muon/proton analogue (the framework's only published Z-precedent).
5. `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl` — Z=1 joint-fit notebook for pattern reference (Self-Energy owns the Z=3 analogue, do not duplicate the notebook here).

## Iteration protocol

Each `/loop` iteration:

1. Read `STATE.md` to find the last iteration's pending-next.
2. Execute **one** substantive step. Step shapes:
   - Read a source-of-record document not yet read; record key Z=1 formulas + Z-scaling identities in STATE.md.
   - Set up or extend a Wolfram MCP cell in `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_spectroscopy.wl` (new notebook this branch creates, distinct from Self-Energy's `r_e_Li2plus_joint_fit.wl`).
   - Execute / debug a Wolfram cell to compute Z=3 framework prediction for #3 or #2.
   - Draft a per-observable section in `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/12_Li2plus_Spectroscopy.md` (new doc this branch creates).
   - Report (Z-i)/(Z-ii) reading differences explicitly when they exist.
3. Append a dated entry to `STATE.md` (ISO timestamp) recording: what advanced, what's queued next, current observable focus (#2 or #3), the current outcome-matrix branch (A/B/C/D), and any BLOCKED state with full measurement provenance for any value added.
4. `git add` changed files explicitly. New commit only (no `--amend`). Commit message starting with `li2plus-spec:`. Push to this branch.
5. Self-pace next iteration via ScheduleWakeup — **aggressive cadence: 60–180s** (cache-warm; favor 60s the runtime floor unless the next step is explicitly idle-pending; the user has requested fastest viable iteration).
6. Stop the loop (omit ScheduleWakeup) only when:
   - Both observables (#2 + #3) have framework predictions documented in `12_Li2plus_Spectroscopy.md` with verdicts, **or**
   - A hard BLOCKED state requiring Tepper input is recorded.

## Guardrails

- **Wolfram MCP gotchas** per CLAUDE.md: single-line cells; avoid `V` / `e` as symbols; non-commutative `Dot` handled.
- **Crocco compliance:** substantive AI; `<!-- TODO: human reviews and fills in -->` blocks in every interpretation paragraph + the Wolfram notebook header + the per-observable doc updates.
- **Python via `uv run`** only.
- **Commit hygiene:** new commits only; clear one-line messages starting with `li2plus-spec:`.
- **Push to this branch only** (`78-li2plus-spectroscopy`); never to `main` or another branch.
- **Never run destructive git ops.** Halt and record `BLOCKED` if a step would require one.
- **Do not edit files owned by other branches.** Specifically: do **not** create or edit `Bethe_Salpeter/11_Li2plus_HydrogenicIon.md`, `Bethe_Salpeter/13_Li2plus_Hyperfine.md`, or `r_e_Li2plus_joint_fit.wl` — those are owned by Self-Energy / Hyperfine branches respectively. If you need to modify shared files (`Dual_Relativistic_Quantum_Mechanics_I.md`, `FINDINGS_for_author_review.md`), make scoped additions only; the orchestrator will reconcile at merge time.
- **Measurement provenance:** every measured value added must record its source (paper + year + DOI/journal) in STATE.md and in `12_Li2plus_Spectroscopy.md`.
- **(Z-i)/(Z-ii) honesty:** report both readings where they differ; do not silently pick one.

## Done criteria

- `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/12_Li2plus_Spectroscopy.md` contains per-observable result sections for #2 + #3, each with measurement + framework prediction + verdict.
- `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_spectroscopy.wl` exists with Wolfram-MCP-verified cells for both observables.
- STATE.md records the outcome-matrix branch for each observable separately + an overall verdict for this branch's two observables.

## Outcome-matrix (per observable, recorded in STATE.md)

- **A** — Framework prediction matches measurement under the Z=1 triangulated universal cutoff.
- **B** — Match under a Z-scaled cutoff with framework-internal derivation.
- **C** — Match only under per-Z back-fit.
- **D** — Intractable; BLOCKED on Tepper input.
