# Research brief — Hydrogenic Z-scan g-factor (issue #82)

**Branch:** `82-hydrogenic-z-scan-g-factor`
**Issue:** https://github.com/temoTxt/PyPhysics/issues/82
**Sibling:** #78 (Li²⁺ four-observable spot-check; this ticket extends to multi-Z g-factor only)
**Parallel branches under #78 (do not overlap):**
- `78-bethe-salpeter-z-extension-li2plus` (owns Li²⁺ g-factor + joint χ²)
- `78-li2plus-spectroscopy` (owns Li²⁺ Lamb shift + fine structure)
- `78-li2plus-hyperfine` (owns Li²⁺ hyperfine)
**This branch owns:** bound-electron g-factor at multiple Z's (He⁺, Be³⁺, C⁵⁺, Si¹³⁺, Ca¹⁹⁺ — Li²⁺ comes from #78's parallel work, not re-computed here).
**Output files (distinct):** `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/14_HydrogenicIon_Zscan.md` + `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Zscan_fit.wl`

## What this branch is doing

Turn the Z=3 spot-check (running on #78's three branches) into a **Z-scan curve fit**: collect bound-electron g-factor measurements at 4–6 Z values, compute framework predictions under (Z-i) universal-cutoff and (Z-ii) Z-scaled-cutoff readings, and fit $r_e^{(Z)}/r_0$ to the data. Verdict A/B/C/D classifies the Z-dependence.

See issue #82 for full goal, method, acceptance criteria — this brief does not duplicate them.

## Z-scan ion list (target — refine per-iteration with DOI/year provenance)

| Ion | $Z$ | Best $g_e^{\rm bound}$ value | Primary reference |
|---|---|---|---|
| ³He⁺ | 2 | $2.000\,008\,021(15)$ | Hoffmann 1989 / Köhler 2015 update |
| ⁹Be³⁺ | 4 | (verify in iter 1) | various Penning-trap groups |
| ¹²C⁵⁺ | 6 | $2.001\,041\,590\,18(3)$ | Sturm 2011 *PRL* **107**, 023002 |
| ²⁸Si¹³⁺ | 14 | $1.995\,348\,958\,7(5)$ | Sturm 2013 *PRL* **110**, 263002; Köhler 2016 *Nat. Comm.* **7**, 10246 |
| ⁴⁰Ca¹⁹⁺ | 20 | sub-ppb | Glazov 2019 / Köhler-Langes 2018 |

Li²⁺ (Z=3) value $g_e^{\rm bound}({}^7\text{Li}^{2+}) = 2.000\,025\,170\,7(10)$ (Sturm 2014 *Nature* **506**, 467) — already being worked on `78-bethe-salpeter-z-extension-li2plus`; this branch *imports* that value for the joint fit, does not re-derive.

## Source-of-record (read in order)

1. `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` §III.D Eqs. (III.22)/(III.23) — anomalous-g formula + muon/proton analogue (framework's published precedent for particle/scale variation)
2. `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl` — Z=1 joint-fit pattern to mirror at multi-Z
3. `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` Finding 2 — recently-updated verdict structure (this branch will append a multi-Z update)
4. `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md` §2 — the $(g_s/-2)^n \cdot \text{textbook}$ scaling; for g-factor itself the prediction is $g_r(r_e/r_0)$ directly

## Per-ion deliverable shape

For each ion in the Z-scan:

- **As measured:** value, uncertainty, source (paper + year + DOI).
- **QED prediction:** standard textbook bound-state QED expansion $g_e^{\rm bound}(Z) = g_e^{\rm free} \cdot [1 - (Z\alpha)^2/3 - (Z\alpha)^4/12 + \ldots] + \text{QED corrections}$.
- **Framework prediction under both readings:**
  - **(Z-i) Universal cutoff:** $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ (Z=1 triangulated value, universal). Predict $g_r(r_e/r_0) = -2.00231930\ldots$ — *the same value at every Z*. Residual = measured $g_e^{\rm bound}(Z) - (-2.00231930)$.
  - **(Z-ii) Z-scaled cutoff:** invert the g-formula at each measured $g_e^{\rm bound}(Z)$ to get $r_e^{(Z)}/r_0 = (2 - a_e^{\rm bound}(Z))/(2(2 + a_e^{\rm bound}(Z)))$. Tabulate the resulting per-Z back-fit values; check for systematic Z-dependence (linear in $Z\alpha$? quadratic? log? Or random scatter consistent with measurement noise + framework precision?).
- **Per-ion verdict:** ✅ at framework precision / ⚠ at framework precision floor / 🔴.

## Joint χ² + Z-scaling fit (the headline deliverable)

After per-ion values are recorded, the Wolfram MCP notebook `r_e_Zscan_fit.wl` performs:

1. **(Z-i) test:** with $r_e/r_0$ fixed at the Z=1 triangulated value, compute $\chi^2$ across all Z's. If $\chi^2 \ll N_{\rm ions}$ → Z-universal cutoff confirmed; **outcome A.**
2. **(Z-ii) per-Z back-fit:** report the per-Z back-fit values $r_e^{(Z)}/r_0$ for each ion. Plot or tabulate vs $Z\alpha$. Two sub-cases:
   - Systematic Z-dependence with a clean functional form (e.g., $r_e^{(Z)}/r_0 = a + b(Z\alpha)^2 + \ldots$ with $\chi^2_{\rm fit} \ll N_{\rm ions}$): **outcome B** with the derived form.
   - No clean form; per-Z back-fit values scatter or follow QED's bound-state $a_e(Z\alpha)$ inheritance: **outcome C** (framework's cutoff inherits QED bound-state structure per-Z, not derivable from framework's own apparatus).
3. **Cross-check with PR #70's lepton-axis verdict:** PR #70 showed the cutoff is particle-specific via $a_\ell$. The Z-axis result should be consistent with that: $r_e^{(Z)}/r_0 = (2 - a_e^{\rm bound}(Z\alpha))/(2(2 + a_e^{\rm bound}(Z\alpha)))$ would inherit QED's bound-state $a_e(Z\alpha)$ identically — extending PR #70's "particle-specific through $a_\ell$" verdict to "particle-and-Z-specific through $a_\ell^{\rm bound}(Z\alpha)$".

## Iteration protocol

Each `/loop` iteration:

1. Read `STATE.md` to find the last iteration's pending-next.
2. Execute **one** substantive step. Step shapes:
   - Read a source-of-record document not yet read; record key formulas + Z-scaling in STATE.md.
   - Look up + record the precise measurement value + DOI/year for one ion in the Z-scan table (one ion per iter to keep steps small).
   - Set up or extend `r_e_Zscan_fit.wl` (per-Z prediction cell or joint-fit cell).
   - Draft a per-ion section in `14_HydrogenicIon_Zscan.md`.
   - Once all ions are catalogued, execute the joint χ² + Z-scaling fit.
3. Append a dated entry to `STATE.md`: what advanced, what's queued next, current ion focus (He⁺/Be³⁺/C⁵⁺/Si¹³⁺/Ca¹⁹⁺ or "joint fit"), outcome-matrix branch tentative, any BLOCKED state.
4. `git add` changed files explicitly. New commit only. Commit message starting with `zscan:`. Push to this branch.
5. Self-pace next iteration via ScheduleWakeup — **aggressive cadence: 60–180s** (cache-warm; favor 60s the runtime floor unless the next step is explicitly idle-pending; the user has requested fastest viable iteration).
6. Stop the loop only when:
   - All ions in the Z-scan have measured values + framework predictions documented + joint χ² + Z-scaling fit reported, **or**
   - A hard BLOCKED state is recorded.

## Guardrails

- **Wolfram MCP gotchas** per CLAUDE.md: single-line cells; avoid `V` / `e` as symbols (use `potV` / `ee`); non-commutative `Dot` handled.
- **Crocco compliance:** substantive AI; `<!-- TODO -->` blocks throughout. Substantive-AI tag on the Z-scaling functional-form choice (if Outcome B emerges).
- **Python via `uv run`** only.
- **Commit hygiene:** new commits only; `zscan:` prefix.
- **Push to this branch only** (`82-hydrogenic-z-scan-g-factor`); never to `main` or another branch.
- **Never run destructive git ops.** Halt and record `BLOCKED` if needed.
- **Do not edit files owned by #78 branches:** `11_Li2plus_HydrogenicIon.md`, `12_Li2plus_Spectroscopy.md`, `13_Li2plus_Hyperfine.md`, `r_e_Li2plus_joint_fit.wl`, `r_e_Li2plus_spectroscopy.wl`, `r_e_Li2plus_hyperfine.wl` all belong to parallel branches.
- **Measurement provenance is non-negotiable** for every ion: paper + year + DOI/journal in both STATE.md and `14_HydrogenicIon_Zscan.md`. No bare values.
- **Li²⁺ value is imported, not re-derived:** use the value Self-Energy branch's work produces; do not duplicate that work.

## Done criteria

- 5+ ions catalogued with measurement provenance (target: He⁺ Z=2, Li²⁺ Z=3 imported, Be³⁺ Z=4, C⁵⁺ Z=6, Si¹³⁺ Z=14; bonus Ca¹⁹⁺ Z=20 if precision data lookup is fast).
- Per-ion framework prediction under (Z-i) and (Z-ii) readings in `14_HydrogenicIon_Zscan.md`.
- Joint χ² + Z-scaling fit reported in `r_e_Zscan_fit.wl` with Wolfram-MCP-verified numerics.
- Z-axis verdict (A/B/C/D) recorded in both STATE.md and `14_HydrogenicIon_Zscan.md`.
- FINDINGS Finding 2 multi-Z update appended.

## Outcome-matrix

- **A** — Z-universal cutoff fits all 4–6 ions at framework precision (single $r_e/r_0$).
- **B** — Z-scaled cutoff $r_e^{(Z)}/r_0 = f(Z\alpha)$ with a framework-internally derivable form fits the data; Z-dependence characterised.
- **C** — Per-Z back-fit values follow QED's bound-state $a_e(Z\alpha)$ inheritance (no framework-internal Z-derivation; analogous to PR #70's lepton-axis verdict on the Z-axis).
- **D** — Predictions cannot be made tractably under either reading. BLOCKED on Tepper input.
