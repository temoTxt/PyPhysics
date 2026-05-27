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
