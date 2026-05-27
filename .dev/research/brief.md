# Research brief — Li²⁺ hyperfine — issue #78

**Branch:** `78-li2plus-hyperfine`
**Issue:** https://github.com/temoTxt/PyPhysics/issues/78
**Parent campaign:** [#50](https://github.com/temoTxt/PyPhysics/issues/50) Bethe–Salpeter (closed; #78 extends it)
**Parallel branches (do not overlap with these):**
- `78-bethe-salpeter-z-extension-li2plus` (running; owns observable #1 g-factor + integration / joint χ² + `r_e_Li2plus_joint_fit.wl`)
- `78-li2plus-spectroscopy` (running; owns observables #2 Lamb shift + #3 fine structure)
**This branch owns:** observable #4 — Li-7 1s hyperfine splitting (nuclear $I = 3/2$)
**Output file (distinct):** `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/13_Li2plus_Hyperfine.md`

## What this branch is doing

Predict the Li-7 1s hyperfine splitting in Li²⁺ under the dual-theory framework, comparing to measurement and contributing the Z=3 hyperfine prediction value to the joint χ² fit owned by the parallel Self-Energy branch.

Hyperfine is a **real Z-axis discriminator** of the cutoff-universality question — the Fermi-contact term involves $\boldsymbol\sigma \cdot \mathbf{I}$ matrix elements that scale with $g_s = g_r(r_e/r_0)$ to power $n=1$ (per `Bethe_Salpeter/10_CrossComparison.md §2`), so $r_e$ enters the prediction directly. This is the analog of the H 21-cm BS-§22 result.

## Observable scope

### #4 — Li-7 1s hyperfine splitting

- **Measurement:** Beckmann 1974 era: $\Delta\nu_{\rm HFS}({}^7\text{Li}^{2+}, 1s) \approx 12\,732$ MHz. Refined by Riis 1994 and modern groups; precision is at the kHz level for the most recent measurements. (Look up exact CODATA-2018 or current best value in iter 1; record DOI + year.)
- **QED breakdown:** Fermi contact $\propto (Z\alpha)^4 (m_e/m_p) \mu_N \mu_I$ with nuclear spin $I = 3/2$ for Li-7 (vs $I = 1/2$ for the proton in H). Nuclear-spin angular structure is therefore *different* from H 21-cm, not just a Z-rescaling.
- **Z-scaling:** leading Fermi-contact $\propto Z^3 g_s (m_e/m_p)$; for Li²⁺ with $Z=3$ that's a $27\times$ enhancement over H's $1420$ MHz baseline before nuclear-magnetic-moment ratio corrections. Final ~12.7 GHz emerges after the Li-7 nuclear $g_I$ and reduced-mass factors.
- **$r_e$ coupling:** through $g_s = g_r(r_e/r_0)$ at power $n=1$. **(Z-i)** universal cutoff vs **(Z-ii)** Z-scaled cutoff give different predictions.
- **Framework apparatus:** apply BS-§22 Fermi-contact formula with $Z=3$ scaling and Li-7 nuclear-spin structure $I=3/2$. The $I=3/2$ structure changes the hyperfine multiplet from a single splitting (H proton's two-level $F=0,1$) to a multi-level coupling — the $1s_{1/2}$ ground state couples to $I=3/2$ to give $F=1,2$ (two levels, with the splitting being the headline observable).

## Source-of-record (read in order)

1. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/06_Hyperfine.md` — BS-§22 (Z=1 H 21-cm baseline + Fermi-contact apparatus) **← PRIMARY focus**
2. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md` §2 — confirm $n_{\rm HFS} = 1$ in the $(g_s/-2)^n$ scaling pattern.
3. `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` §III.D Eqs. (III.22)/(III.23) — anomalous-g formula + the muon/proton analogue (the framework's only published Z/particle-mass precedent).
4. Look up: CODATA / NIST atomic-spectroscopy reference value for Li²⁺ 1s hyperfine. Beckmann 1974 *Z. Phys.* **270**, 173 is the standard early reference; modern values refined by Schweikhard 1991 and others.

## Iteration protocol

Each `/loop` iteration:

1. Read `STATE.md` to find the last iteration's pending-next.
2. Execute **one** substantive step. Step shapes:
   - Read a source-of-record document not yet read; record key Z=1 Fermi-contact formula + Z/nuclear-spin scaling identities in STATE.md.
   - Set up or extend a Wolfram MCP cell in `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_hyperfine.wl` (new notebook this branch creates, distinct from other branches' notebooks).
   - Execute / debug a Wolfram cell to compute Z=3 + $I=3/2$ framework prediction.
   - Draft a per-observable section in `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/13_Li2plus_Hyperfine.md` (new doc this branch creates).
   - Report (Z-i)/(Z-ii) reading differences explicitly when they exist; record measurement provenance.
3. Append a dated entry to `STATE.md` (ISO timestamp) recording: what advanced, what's queued next, current observable focus (always #4 for this branch), the current outcome-matrix branch (A/B/C/D), and any BLOCKED state with measurement provenance.
4. `git add` changed files explicitly. New commit only (no `--amend`). Commit message starting with `li2plus-hfs:`. Push to this branch.
5. Self-pace next iteration via ScheduleWakeup — **aggressive cadence: 60–180s** (cache-warm; favor 60s the runtime floor unless the next step is explicitly idle-pending; the user has requested fastest viable iteration).
6. Stop the loop (omit ScheduleWakeup) only when:
   - Observable #4 has a framework prediction documented in `13_Li2plus_Hyperfine.md` with a verdict, **or**
   - A hard BLOCKED state requiring Tepper input is recorded.

## Guardrails

- **Wolfram MCP gotchas** per CLAUDE.md: single-line cells; avoid `V` / `e` as symbols; non-commutative `Dot` handled.
- **Crocco compliance:** substantive AI; `<!-- TODO: human reviews and fills in -->` blocks in every interpretation paragraph + the Wolfram notebook header + the per-observable doc updates. **Substantive AI tag on the Li-7 $I=3/2$ nuclear-spin coupling decision** specifically — extending H 21-cm's $I=1/2$ machinery to $I=3/2$ is a non-trivial substantive choice.
- **Python via `uv run`** only.
- **Commit hygiene:** new commits only; clear one-line messages starting with `li2plus-hfs:`.
- **Push to this branch only** (`78-li2plus-hyperfine`); never to `main` or another branch.
- **Never run destructive git ops.** Halt and record `BLOCKED` if a step would require one.
- **Do not edit files owned by other branches.** Specifically: do **not** create or edit `Bethe_Salpeter/11_Li2plus_HydrogenicIon.md`, `Bethe_Salpeter/12_Li2plus_Spectroscopy.md`, `r_e_Li2plus_joint_fit.wl`, or `r_e_Li2plus_spectroscopy.wl` — those belong to Self-Energy / Spectroscopy branches. If you need scoped additions to shared files (`Dual_Relativistic_Quantum_Mechanics_I.md`, `FINDINGS_for_author_review.md`), make them additive only; orchestrator reconciles at merge time.
- **Measurement provenance:** every measured value must record its source (paper + year + DOI/journal) in STATE.md and in `13_Li2plus_Hyperfine.md`. The Li-7 nuclear magnetic moment $\mu_I$ and nuclear $g_I$ enter the prediction — record their CODATA / NUBASE provenance too.
- **(Z-i)/(Z-ii) honesty:** report both readings where they differ; do not silently pick one.
- **Nuclear-structure caveat:** Li-7 has hyperfine anomaly (Bohr–Weisskopf effect) and quadrupole structure. These are *out of scope* per the issue #78 acceptance criteria, but record them as framework-floor caveats in the verdict prose.

## Done criteria

- `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/13_Li2plus_Hyperfine.md` contains a per-observable result section for #4 with measurement + framework prediction (both Z-i and Z-ii readings if they differ) + verdict.
- `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_hyperfine.wl` exists with Wolfram-MCP-verified cells for the Fermi-contact prediction at Z=3 with $I=3/2$ structure.
- STATE.md records the outcome-matrix branch for observable #4.

## Outcome-matrix (recorded in STATE.md)

- **A** — Framework prediction matches measurement under the Z=1 triangulated universal cutoff at framework precision.
- **B** — Match under a Z-scaled cutoff with framework-internal derivation.
- **C** — Match only under per-Z back-fit (each Z requires its own cutoff).
- **D** — Intractable; BLOCKED on Tepper input (specifically: the $I = 1/2 \to I = 3/2$ extension of BS-§22's Fermi-contact apparatus may need framework-internal guidance).
