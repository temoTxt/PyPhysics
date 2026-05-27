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
