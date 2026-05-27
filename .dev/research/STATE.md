# Li²⁺ spectroscopy (Lamb shift + fine structure) — state log

## Iteration 0 — 2026-05-27 — initialized

- Branch `78-li2plus-spectroscopy` checked out from main (post-#67 r_e thread closure).
- `.dev/research/brief.md` written (Lamb shift + fine structure scope; fine structure prioritized per Self-Energy iter-1 finding that Lamb shift is a weak Z-axis discriminator).
- `.dev/research/loop_prompt.md` written.
- No prediction work yet.
- **Current observable focus:** none yet; **#3 fine structure prioritized over #2 Lamb shift** per the brief's priority signal from the parallel Self-Energy branch (iter-1 finding: $r_e$ enters Lamb shift only at sub-leading order, below the Bethe-estimate precision floor).
- **Next:** read source-of-record §1 (`Bethe_Salpeter/03_FineStructure.md`, BS-§14) — record the Z=1 fine-structure formula, the $(g_s/-2)^{n=1}$ scaling, and the $(Z\alpha)^4$ leading Z-scaling. Then proceed to §4 (DRQM I §III.D) to confirm the anomalous-g coupling chain.
- **Outcome-matrix:** not yet determinable.
- **Status:** READY.

## Iteration 1 — 2026-05-27 — read source-of-record §1 (BS-§14 fine structure)

**Step taken:** read `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md` (BS-§14, the Z=1 fine-structure source-of-record). One step (read + record); no notebook/doc work this iteration.

**Key Z=1 formulas + Z-scaling identities recorded:**

- **Sommerfeld–Dirac level formula** (BS-§14, line 5–7):
  $E_{nj} = -\frac{mc^2(Z\alpha)^2}{2n^2}\left[1 + \frac{(Z\alpha)^2}{n^2}\left(\frac{n}{j+1/2}-\frac{3}{4}\right) + \mathcal{O}((Z\alpha)^4)\right]$.
- **Leading Dirac 2P₃/₂–2P₁/₂ splitting at Z=1:** $\Delta E_{FS} = m_e c^2\,\alpha^4/32 \approx 10\,949$ MHz (BS-§14.2, line 84, 90). Leading term $\propto (Z\alpha)^4 m_e c^2$ ⇒ **splitting scales as $Z^4$** (this is the leading Z-scaling identity for #3).
- **Anomalous-g correction** (BS-§14.2, line 94–96): $\Delta E_{anom} = \frac{g_s-2}{2}\cdot\Delta E_{leading}$ — the $n=1$ case of the $(g_s/-2)^n\cdot\text{textbook}$ pattern. At $g_s=-2.00231930$ (triangulated $r_e/r_0=0.4994205099128317$) ⇒ $\approx 12.7$ MHz at Z=1. **This is the $r_e$ entry point — the Z-axis discriminator for #3.**
- **Z=1 total:** $10\,949 + 12.7 = 10\,962$ MHz vs CODATA 10,969(1) MHz; residual $-7$ MHz (recoil + two-loop, out of Bethe-estimate scope).
- **Z=1 measurement provenance:** CODATA-2018 $\Delta E_{FS}(2P_{3/2}-2P_{1/2}) = 10\,969(1)$ MHz; experimental $10\,969.13(10)$ MHz (Hagley & Pipkin 1994 + refinements; BS-§14.2 line 86).
- **Honesty flag inherited from source (BS-§14.2 line 100, 121):** the "✅ at triangulated $r_e$" is **back-fit self-consistency, NOT independent corroboration** — the triangulated $r_e$ is by construction the value giving the measured $g_s$, so $\frac{g_{s,\text{meas}}-2}{2}\times E_{\text{leading}}$ is identical to textbook QED's leading-$g_s$ formula. The 7 MHz residual is ~70σ from the 0.1 MHz experimental precision; "✅" is relative to the Bethe-estimate floor only. Carry this framing into 12_Li2plus_Spectroscopy.md #3.

**⚠ NUMERICAL FLAG (queued for reconciliation before #3 verdict):** the brief's #3 figures are inconsistent with the $Z^4$ scaling derived from the source. Brief says Z=3 FS ≈ 7.4 GHz and "measurement ≈ 7,367 MHz" (brief.md lines 28–29). But $Z^4=81$ at Z=3 ⇒ $81\times 10\,949 \approx 887$ GHz leading Dirac — **two orders of magnitude larger** than the brief's 7,367 MHz. The 7,367 MHz figure is likely a mis-transcription or refers to a different splitting/species (possibly neutral-Li / Li⁺, not hydrogenic Li²⁺). **Do not adopt 7,367 MHz as the measurement** until provenance is verified. Reconcile Li²⁺ 2P FS measurement provenance (literature lookup) at the step before drafting #3's verdict; do NOT edit brief.md (research input).

**Queued next:** read source-of-record §3 (`10_CrossComparison.md §2`) to confirm $n_{FS}=1$ in the $(g_s/-2)^n\cdot\text{textbook}$ scaling pattern; then set up `r_e_Li2plus_spectroscopy.wl` cell 1 computing Z=3 leading Dirac FS $= m_e c^2 (Z\alpha)^4/32$ at $Z=3$ with ⁷Li reduced-mass, to fix the correct ~887 GHz scale and supersede the brief's suspect 7.4 GHz.

- **Current observable focus:** #3 fine structure (primary).
- **Outcome-matrix:** not yet determinable (Z=3 prediction not yet computed).
- **(Z-i)/(Z-ii):** not yet diverged — both readings share the leading Dirac $Z^4$ term; any divergence enters only through the anomalous-g piece's cutoff treatment, not yet computed.
- **Status:** READY.

## Iteration 2 — 2026-05-27 — read source-of-record §3 (10_CrossComparison §2/§3)

**Step taken:** read `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md` (PR J cross-comparison; §2 scaling pattern + §3 Lamb-shift reframe). One step (read + record).

**Confirmed:**

- **$n_{FS}=1$ CONFIRMED** for fine structure. §2 (line 57): $f_i(g_s)$ is "linear for spin–orbit and Fermi contact contributions; quadratic for two-fermion spin–spin contributions like positronium ortho-para." The 2P₃/₂–2P₁/₂ splitting is a **spin–orbit** observable ⇒ linear ⇒ $n=1$. ($n=2$ applies only to two-fermion spin-spin, e.g. positronium ortho-para.)
- **Framework prediction structure for #3:** $\Delta E_{FS}^{\rm fw}(Z) = (g_s/-2)^1 \times \Delta E_{FS}^{\rm leading\,Dirac}(Z)$. Numerically at Z=1: $(g_s/-2)=(-2.00231930/-2)=1.00115965$, $\times 10\,949 = 10\,962$ MHz — reproduces BS-§14.2's 10,962. The $(g_s/-2)$ factor = the anomalous-moment multiplier; $r_e$ enters via $g_s = g_r(r_e/r_0)$ (brief line 30).
- **Triangulated cutoff (Z=1) re-confirmed across sources:** $r_e/r_0 = 0.4994205099128317$ (PR #62 joint-best-fit; gives $g_s=-2.00231930$ by construction).
- **#3 honesty disposition (inherited, carry into doc):** §2 lines 70–72 + §1 line 48 — the six $g_s$-dependent "✅ at triangulated $r_e$" results are **one back-fit applied six times, not six independent corroborations**. #3's framework "✅" is back-fit self-consistency at the textbook leading-$g_s$ precision floor.
- **#2 Lamb-shift weak-discriminator CONFIRMED from source (§3 line 90):** "The Lamb shift result is $r_e$-independent because the leading log-Bethe contribution is $g_s$-symmetric." ⇒ **(Z-i) and (Z-ii) give identical #2 predictions** — validates the brief's priority signal. #2 is ✅-by-inheritance / reproduction-not-endorsement (§3 lines 86–89); document briefly, no Z-axis-test value.

**Queued next:** scaffold `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_spectroscopy.wl` (new notebook, Crocco header + `<!-- TODO -->`), cell 1 = Z=3 leading Dirac FS $= m_e c^2 (Z\alpha)^4/32$ at $Z=3$ with ⁷Li reduced-mass correction, to fix the correct ~887 GHz scale (supersedes brief's suspect 7.4 GHz / 7,367 MHz). Use single-line Wolfram cells; symbols `ee`/`potV` per CLAUDE.md gotchas.

- **Current observable focus:** #3 fine structure (primary).
- **Outcome-matrix:** #3 not yet determinable (Z=3 prediction not yet computed). #2 trending **A** (✅ by inheritance under universal cutoff, $r_e$-independent ⇒ (Z-i)=(Z-ii)).
- **(Z-i)/(Z-ii):** #2 — **no divergence** (Lamb shift $r_e$-independent, source-confirmed). #3 — leading Dirac $Z^4$ shared; potential divergence only in the anomalous-g cutoff treatment (not yet computed).
- **Status:** READY.

## Iteration 3 — 2026-05-27 — scaffold r_e_Li2plus_spectroscopy.wl + cell 1 (Z=3 leading Dirac FS)

**Step taken:** created new notebook `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_spectroscopy.wl` (this branch's notebook; distinct from Self-Energy's `r_e_Li2plus_joint_fit.wl`). Crocco-compliant header (substantive AI, author line, `<!-- TODO human review -->`), records the $n_{FS}=1$ scaling + back-fit honesty scope + (Z-i)/(Z-ii) framing. **Cell 1 written** (single-line, `ee`/`potV` gotchas respected): $\Delta E_{FS}(Z) = m_e c^2 (Z\alpha)^4/32$ with sanity Z=1→~10949 MHz and Z=3→~887 GHz prediction, plus $Z^4$ ratio check.

**Hand-computed expected values (to be MCP-verified next iter):** constants $m_ec^2=510998.95$ eV, $h=4.135667696\times10^{-15}$ eV·s, $\alpha=1/137.035999084$ ⇒ Z=1: $\approx 10\,949$ MHz (matches BS-§14.2 ✓); Z=3: $\approx 886\,800$ MHz $= 886.8$ GHz; ratio $=81$ ✓. **This confirms the iter-1 flag: the brief's 7.4 GHz / 7,367 MHz Z=3 figures are wrong by ~2 orders of magnitude; the leading-Dirac Li²⁺ 2P FS scale is ~887 GHz.**

**Queued next:** execute cell 1 via Wolfram MCP (`mcp__wolfram__WolframLanguageEvaluator`) to verify the three printed values; record "Result ✅" inline in the notebook. Then add cell 2 = framework prediction $\Delta E_{FS}^{\rm fw}(Z{=}3)=(g_s/-2)\times\Delta E_{\rm leading}(Z{=}3)$ at $g_s=-2.00231930$ (Z-i universal cutoff) and compute the anomalous offset at Z=3.

- **Current observable focus:** #3 fine structure (primary).
- **Outcome-matrix:** #3 not yet determinable (need measurement provenance + framework value). #2 trending **A**.
- **(Z-i)/(Z-ii):** unchanged — #2 no divergence; #3 divergence (if any) only via anomalous-g cutoff, not yet computed. NOTE: a measurement value for Li²⁺ 2P FS is still UNRESOLVED — brief's 7,367 MHz refuted by $Z^4$ scaling; literature provenance lookup still required before #3 verdict (do not adopt brief figure).
- **Status:** READY.

## Iteration 4 — 2026-05-27 — execute cell 1 via Wolfram MCP (Z=3 leading Dirac FS verified)

**Step taken:** ran cell 1 of `r_e_Li2plus_spectroscopy.wl` through `mcp__wolfram__WolframLanguageEvaluator`. Recorded the result inline in the notebook.

**Wolfram MCP result (2026-05-27):** `{10949.283786302933, 886891.9866905379, 81.00000000000003}` ✅
- Z=1 leading Dirac FS = **10,949.28 MHz** — matches BS-§14.2's ~10,949 MHz sanity value ✓.
- Z=3 leading Dirac FS = **886,892 MHz = 886.9 GHz** (infinite-nuclear-mass).
- $Z^4$ ratio = **81.00** ✓.
- (Harmless `Symbol::undefined2` warning on the `ZZ` pattern var in `SetDelayed`; does not affect the numeric output.)

**Conclusion:** the Li²⁺ 2P₃/₂–2P₁/₂ leading-Dirac scale is **~887 GHz**, definitively refuting the brief's draft 7.4 GHz / 7,367 MHz figures (wrong by ~2 orders of magnitude). The ~887 GHz figure is now the verified scale of record for #3.

**Queued next:** add cell 2 to the notebook — framework prediction $\Delta E_{FS}^{\rm fw}(Z{=}3) = (g_s/-2)^1 \times \Delta E_{\rm leading}(Z{=}3)$ at $g_s = -2.00231930$ (the Z=1-triangulated universal cutoff, reading **(Z-i)**); compute the Z=3 anomalous offset = $\frac{g_s-2}{2}\times 886\,892$ MHz $\approx +1\,029$ MHz, framework total $\approx 887\,921$ MHz. Then execute it. (Reduced-mass + recoil refinement for ⁷Li is a sub-% correction, ~−200 MHz; defer as a refinement note.) After that, the literature-provenance lookup for the Li²⁺ 2P FS measurement is the gating item before drafting #3's verdict in the doc.

- **Current observable focus:** #3 fine structure (primary).
- **Outcome-matrix:** #3 not yet determinable (framework value pending cell 2; measurement provenance still unresolved). #2 trending **A**.
- **(Z-i)/(Z-ii):** unchanged — #2 no divergence; #3 leading $Z^4$ verified, anomalous-g cutoff divergence not yet computed.
- **Status:** READY.

## Iteration 5 — 2026-05-27 — cell 2: framework #3 prediction at Z=3 (Z-i), Wolfram-verified

**Step taken:** added cell 2 to `r_e_Li2plus_spectroscopy.wl` and executed via Wolfram MCP. Framework fine-structure prediction for Li²⁺ under reading **(Z-i)** (universal Z=1-triangulated cutoff), $n_{FS}=1$: $\Delta E_{FS}^{\rm fw}(Z{=}3) = (g_s/-2)\times\Delta E_{\rm leading}(Z{=}3)$.

**Wolfram MCP result (2026-05-27):** `{1.00115965, 886891.9866905379, 1028.4842923656106, 887920.4709829035}` ✅
- $(g_s/-2)$ factor = **1.00115965** (at $g_s=-2.00231930$, triangulated $r_e/r_0=0.4994205099128317$).
- leading Dirac Z=3 = 886,892 MHz; **anomalous offset Z=3 = +1,028.5 MHz**; **framework total (Z-i) = 887,920 MHz = 887.92 GHz**.
- Infinite-nuclear-mass; ⁷Li reduced-mass + recoil refinement ~−0.02% (~−200 MHz), within Bethe-estimate floor.

**(Z-i)/(Z-ii) status:** (Z-i) framework value now in hand = **887,920 MHz**. **(Z-ii) Z-scaled cutoff NOT YET DERIVABLE** — pending the DRQM I §III.D-extension reading (source-of-record §4); if no framework-internal Z-scaling emerges there, **(Z-ii) ≡ (Z-i)** (the anomalous factor $a_e=0.00115965$ is tiny, so even a modest Z-scaling of the cutoff shifts only the ~1,028 MHz anomalous piece, i.e. sub-0.1% of the total — #3 is a *weak* discriminator at the leading scale but the anomalous piece is where any Z-axis signal lives).

**Queued next:** read source-of-record §4 — `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` §III.D Eqs. (III.22)/(III.23) — the anomalous-g / $r_e$ formula + muon/proton analogue (the framework's only published Z-precedent). Determine whether §III.D's reasoning implies a Z-scaled cutoff [→ (Z-ii) derivable, outcome **B**] or only a Z-universal one [→ (Z-ii)≡(Z-i), outcome **A**]. THEN literature-provenance lookup for the Li²⁺ 2P FS measurement (gating item before #3 verdict). THEN draft #3 section in `12_Li2plus_Spectroscopy.md`.

- **Current observable focus:** #3 fine structure (primary).
- **Outcome-matrix:** #3 framework value computed (887,920 MHz, Z-i); branch determination pending §III.D (A vs B) + measurement provenance. #2 trending **A**.
- **Status:** READY.

## Iteration 6 — 2026-05-27 — read source-of-record §4 (DRQM I §III.D) — (Z-ii) RESOLVED → collapses to (Z-i)

**Step taken:** read `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` Eqs. (III.21)–(III.23) + §III.D-extension (lines 451–588). One step (read + determination).

**DECISIVE (Z-i)/(Z-ii) DETERMINATION:**

- **§III.D cutoff formula (III.21–22):** $g_r = 2\left[1 - \dfrac{4r_0}{2r+r_0}\right]$, with $r_0 = e^2/(m_e c^2)$ — the **electron's own classical radius**. **There is NO nuclear charge $Z$ anywhere in the formula.** The cutoff $r_e/r_0$ is fixed by matching the **free-electron** anomalous moment $g_s=-2.00231930$ (back-fit; published $0.499857150068631$ is wrong, triangulated $0.4994205099128317$ is the value that reproduces $g_s$).
- **Cross-particle (III.23):** muon/proton get their OWN $r_0^\ell = e^2/(m_\ell c^2)$ — the cutoff is **particle-mass-dependent** (per-lepton), and the muon back-fit ($r_\mu/r_0^\mu=0.499417379$) differs from the electron's at $>57\,000\sigma$. But this dependence is on the **bound particle's identity/mass, NOT on the external nuclear charge Z.**
- **⇒ The framework provides NO Z-scaling of the cutoff.** The electron's $g_s$ is intrinsic and Z-independent (same in H, Li²⁺, free). **Reading (Z-ii) "Z-scaled cutoff" is NOT DERIVABLE from §III.D — it collapses to (Z-i).** The SAME $r_e/r_0 = 0.4994205099128317$ (hence same $g_s$, same $(g_s/-2)=1.00115965$ factor) applies at Z=3. **Outcome A for #3.**

**Honest consequence (important reframe):** because the framework cutoff is a free-electron property with no Z-handle, **#3 is a WEAKER Z-axis discriminator than the brief hoped.** The Li²⁺ prediction is just $Z^4$-scaled leading Dirac × the same $(g_s/-2)$ factor — no NEW Z-axis signal beyond the trivial $Z^4$ scaling. The brief's premise (#3 "engages $r_e$ through the anomalous-g mechanism and is a real Z-axis discriminator") holds only in that $r_e$ enters the prediction, but it does NOT vary with Z, so it cannot discriminate a Z-universal vs Z-scaled cutoff — there is no Z-scaled cutoff on offer. Both this branch's observables (#2, #3) end up **outcome A / (Z-i)=(Z-ii)**.

**Caveat to document:** a real bound-state-QED effect exists (bound-electron g-factor gains $(Z\alpha)^2$ corrections in the nuclear field), but §III.D does NOT model it (it models the free anomalous moment). At Z=3 this correction to the ~1,028 MHz anomalous piece is $\sim(Z\alpha)^2 a_e \times$(that piece) ≈ sub-MHz — within the Bethe-estimate floor, so it doesn't change the verdict.

**Queued next:** literature-provenance lookup for the actual measured Li²⁺ 2P₃/₂–2P₁/₂ fine-structure interval (gating item before #3 verdict; brief's 7,367 MHz already refuted, need the real ~887 GHz-scale measurement + source). Use WebSearch / known refs (e.g. Riis et al. 1994; precision Li²⁺ FS literature). THEN draft #3 section in `12_Li2plus_Spectroscopy.md` with measurement + framework 887,920 MHz + verdict (outcome A, back-fit-self-consistency framing, weak-discriminator note).

- **Current observable focus:** #3 fine structure (primary).
- **Outcome-matrix:** **#3 → A** (Z-universal cutoff; (Z-ii)≡(Z-i); framework value 887,920 MHz). **#2 → A** (Lamb shift $r_e$-independent). Verdict-drafting pending measurement provenance.
- **(Z-i)/(Z-ii):** **RESOLVED — no divergence for either observable.** #2: Lamb shift $r_e$-independent. #3: §III.D cutoff has no Z-dependence ⇒ (Z-ii) not derivable, ≡ (Z-i).
- **Status:** READY.

## Iteration 7 — 2026-05-27 — measurement-provenance lookup #3 (WebSearch): scale confirmed ~888 GHz, brief figure refuted

**Step taken:** WebSearch for the hydrogen-like Li²⁺ 2P₃/₂–2P₁/₂ fine-structure interval.

**Findings:**
- **Scale CONFIRMED:** hydrogenic Li²⁺ 2P FS ≈ H value × $Z^4$(=81) ≈ 0.366 cm⁻¹ × 81 ≈ **29.6 cm⁻¹ ≈ 888 GHz** (a generic source rounded to "0.4 cm⁻¹ × 81 = 32.4 cm⁻¹"). Consistent with our Wolfram-verified leading-Dirac 886,892 MHz / framework 887,920 MHz. **The brief's 7,367 MHz / 7.4 GHz figures are definitively refuted** (~2 orders of magnitude low).
- **Brief's confusion source identified:** **neutral lithium** (3-electron, "Li I") 2p₃/₂–2p₁/₂ fine structure is ≈ 10 GHz (arXiv:1405.4111 "QED corrections to the 2P fine splitting in Li" treats neutral Li, ~10 GHz). The brief conflated hydrogenic **Li²⁺** (1-electron, Z=3, FS ~888 GHz) with neutral/many-electron Li (FS ~10 GHz). These are different systems; #3's target is the hydrogenic ~888 GHz interval.
- **Direct high-precision EXPERIMENTAL value of the hydrogenic Li²⁺ 2P FS interval was NOT surfaced** in this search. The canonical Li²⁺ 2s–2p precision experiment is **Riis et al., Phys. Rev. A 49, 207 (1994)** (measured the 2S–2P Lamb shift, i.e. #2; the FS interval is less directly/precisely measured for hydrogenic Li²⁺). #2's measurement is already provenanced in the brief: **Schiffer et al., PRL 74, 2188 (1995), $\Delta E_{2S-2P}=62\,765(21)$ MHz** (verify Riis vs Schiffer attribution when drafting).

**Provenance status:** #3 — framework prediction in hand (887,920 MHz, outcome A); precise experimental comparison value still to be pinned. NOT a hard BLOCKED state: framework prediction is computed and outcome A is determined; the open item is the experimental reference value, which is a literature-completeness matter, not a Tepper-input dependency.

**Queued next:** one targeted lookup of Riis et al. 1994 (PRA 49, 207) for the precise hydrogenic Li²⁺ 2S–2P / 2P-FS values to pin #3's experimental comparison and cross-check #2's Schiffer-1995 value. THEN draft the #3 section in `12_Li2plus_Spectroscopy.md` (measurement + framework 887,920 MHz + verdict A, back-fit-self-consistency + weak-discriminator framing).

- **Current observable focus:** #3 fine structure (primary).
- **Outcome-matrix:** **#3 → A** (framework 887,920 MHz; scale confirmed ~888 GHz). **#2 → A**. Both verdicts pending only doc-drafting + final measurement-value citation.
- **(Z-i)/(Z-ii):** RESOLVED — no divergence for either observable (unchanged from iter-6).
- **Status:** READY.

**Sources (iter-7):** arXiv:1405.4111 (Li 2P fine splitting, neutral Li); generic Z⁴-scaling fine-structure references. Riis et al., Phys. Rev. A 49, 207 (1994) [to verify next iter]. Schiffer et al., PRL 74, 2188 (1995) [#2, from brief].

## Iteration 8 — 2026-05-27 — Riis et al. 1994 lookup: it's Li⁺ (2-electron), NOT hydrogenic Li²⁺ — second brief mis-attribution

**Step taken:** targeted WebSearch for Riis et al., PRA 49, 207 (1994).

**Finding (correction):** Riis et al. 1994 measures **Li⁺** (singly-ionized, **two-electron** 1s2s ³S₁–1s2p ³P_J manifold), NOT hydrogenic **Li²⁺** (one-electron, Z=3). It is "the best two-electron Lamb shift measurement in the literature" (±0.5 MHz); Li⁺ ⁷Li spin-averaged effective shift = 37,429.40(39) MHz. **This is a He-like (Z=3, 2-electron) system — a different observable from this branch's hydrogenic Li²⁺ targets.** So the brief's "Riis et al. (1994)" attribution for #3 is a **second mis-attribution** (after the 7,367 MHz / neutral-Li confusion in iter-7).

**Consequence for #3 measurement provenance — characterized:**
- **No clean high-precision direct EXPERIMENTAL value of the hydrogenic Li²⁺ 2²P₃/₂–2²P₁/₂ fine-structure interval was located** in two targeted searches (iter-7 generic + iter-8 Riis). The literature the brief pointed to is for *different species* (neutral Li ~10 GHz; He-like Li⁺ 1s2s–1s2p).
- **The comparison value for #3 is therefore the THEORETICAL hydrogenic-ion value** (Sommerfeld–Dirac + QED, $Z^4$-scaled): leading Dirac 886,892 MHz; full-QED ≈ 887,9xx MHz (leading + anomalous + recoil + two-loop), well-known for hydrogenic ions. Framework's 887,920 MHz reproduces this at the Bethe-estimate floor — **outcome A**, comparison is framework-vs-textbook-QED (both $Z^4$ Dirac × same $(g_s/-2)$), NOT framework-vs-direct-experiment.
- **NOT a hard BLOCKED state:** the framework prediction is computed and outcome A determined; the limitation is experimental-data sparseness for this specific interval (a literature fact, recorded honestly in the verdict), not a Tepper-input dependency.
- **#2 provenance note:** keep the brief's Schiffer et al., PRL 74, 2188 (1995), $\Delta E_{2S-2P}(^7\text{Li}^{2+})=62\,765(21)$ MHz for the Lamb shift; flag "verify Schiffer is hydrogenic Li²⁺ 2S₁/₂–2P₁/₂" as a drafting check (it is consistent with the brief's Z-scaled ≈62,729 MHz estimate).

**Queued next:** **DRAFT the #3 fine-structure section** in `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/12_Li2plus_Spectroscopy.md` (create the file with header + #3 section): framework prediction 887,920 MHz (leading 886,892 + anomalous 1,028), outcome A, (Z-i)≡(Z-ii) with the §III.D no-Z-dependence reasoning, back-fit-self-consistency caveat, weak-Z-discriminator reframe, comparison to theoretical ~888 GHz with the honest "no direct precision experiment for this interval located" note + Crocco `<!-- TODO -->` blocks. THEN draft #2 (brief) and overall branch verdict.

- **Current observable focus:** #3 fine structure (primary) — moving to documentation.
- **Outcome-matrix:** **#3 → A**, **#2 → A**. Both framework predictions computed; ready to document.
- **(Z-i)/(Z-ii):** RESOLVED — no divergence (unchanged).
- **Status:** READY.

**Sources (iter-8):** Riis et al., *Phys. Rev. A* **49**, 207 (1994) — Li⁺ (2-electron) Lamb shift/HFS [journals.aps.org/pra/abstract/10.1103/PhysRevA.49.207].
