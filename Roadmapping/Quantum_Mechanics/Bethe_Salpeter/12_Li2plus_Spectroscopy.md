# §12 — Li²⁺ spectroscopy: Lamb shift + fine structure (Z=3 extension)

**Issue [#78](https://github.com/temoTxt/PyPhysics/issues/78), branch `78-li2plus-spectroscopy`.** This document extends the Bethe–Salpeter precision-predictions campaign (parent [#50](https://github.com/temoTxt/PyPhysics/issues/50)) from hydrogen (Z=1) to **hydrogen-like Li²⁺** (Z=3, **one electron**). It covers two observables:

- **#2 — 2S₁/₂–2P₁/₂ Lamb shift** in Li²⁺ *(weak Z-axis discriminator; documented briefly).*
- **#3 — 2P₃/₂–2P₁/₂ fine-structure splitting** in Li²⁺ *(primary focus).*

(Observable #1, the bound-electron g-factor / joint χ² integration, is owned by the parallel Self-Energy branch in [`11_Li2plus_HydrogenicIon.md`](11_Li2plus_HydrogenicIon.md); observable #4, hyperfine, by [`13_Li2plus_Hyperfine.md`](13_Li2plus_Hyperfine.md). This document does not duplicate those.)

> **Species scope (load-bearing).** "Li²⁺" here means the **hydrogenic** (one-electron, nuclear charge Z=3) lithium ion. It is **not** neutral lithium (Li I, 3-electron, 2p fine structure ≈ 10 GHz) nor He-like Li⁺ (Li II, 2-electron, 1s2s–1s2p; the system measured by Riis et al. 1994). The hydrogenic Li²⁺ 2P fine-structure interval scales as $Z^4$ from hydrogen and sits at **≈ 888 GHz**, not the ~10 GHz of the many-electron lithium species.

Companion notebook: [`r_e_Li2plus_spectroscopy.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_spectroscopy.wl) (Wolfram-MCP-verified cells).

## Results

| Result | Status | Outcome | Role |
|---|---|---|---|
| [§12.3 — 2P₃/₂–2P₁/₂ fine structure (Z=3)](#1232--2p2p-fine-structure-z3-primary) | drafted | **A** | **primary — Z-axis test, ✅ at framework precision** |
| §12.2 — 2S₁/₂–2P₁/₂ Lamb shift (Z=3) | pending | (A) | weak discriminator — ✅-by-inheritance |

---

### §12.3 — 2P₃/₂–2P₁/₂ fine structure (Z=3) *(primary)* <a id="1232--2p2p-fine-structure-z3-primary"></a>

**Source-of-record:** [`03_FineStructure.md`](03_FineStructure.md) BS-§14 (Z=1 baseline), [`10_CrossComparison.md` §2](10_CrossComparison.md#2-the-r_e-self-consistency-across-six-g_s-dependent-observables-at-the-triangulated-value) (the $(g_s/-2)^n$ scaling), [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) §III.D (anomalous-g cutoff). *Substantive AI.*

**Sommerfeld–Dirac leading term.** The fine-structure splitting is the $(Z\alpha)^4$ leading Dirac interval

```math
\Delta E_{FS}(2P_{3/2}-2P_{1/2}) = \frac{m_e c^2\,(Z\alpha)^4}{32},
```

which scales as $Z^4$. Wolfram-MCP evaluation (notebook cell 1) gives, with $m_ec^2 = 510998.95$ eV, $\alpha = 1/137.035999084$:

| Z | leading Dirac $\Delta E_{FS}$ | check |
|---|---|---|
| 1 (H) | 10,949.28 MHz | matches BS-§14.2 leading value ✓ |
| 3 (Li²⁺) | **886,891.99 MHz ≈ 886.9 GHz** | $Z^4$ ratio = 81.00 ✓ |

(Infinite-nuclear-mass; the ⁷Li reduced-mass + recoil refinement is ~−0.02%, ~−200 MHz, within the framework's Bethe-estimate precision floor.)

**Framework prediction (anomalous-g correction).** Fine structure is a **spin–orbit** observable, hence **linear** in the electron gyromagnetic ratio: $n_{FS}=1$ in the campaign's $(g_s/-2)^n\times\text{textbook}$ structure (confirmed against [`10_CrossComparison.md` §2](10_CrossComparison.md#2-the-r_e-self-consistency-across-six-g_s-dependent-observables-at-the-triangulated-value): linear for spin–orbit, quadratic only for two-fermion spin–spin). The framework prediction is therefore

```math
\Delta E_{FS}^{\text{framework}}(Z) = \left(\frac{g_s}{-2}\right)\,\Delta E_{FS}^{\text{leading Dirac}}(Z),
```

evaluated at the Z=1-triangulated cutoff $r_e/r_0 = 0.4994205099128317$ ([PR #62](https://github.com/temoTxt/PyPhysics/pull/62)), which gives $g_s = -2.00231930$ (so $g_s/-2 = 1.00115965$). Wolfram-MCP (notebook cell 2):

| quantity | value |
|---|---|
| $(g_s/-2)$ factor | 1.00115965 |
| leading Dirac (Z=3) | 886,892 MHz |
| anomalous offset (Z=3) | +1,028.5 MHz |
| **framework total (Z=3)** | **887,920 MHz ≈ 887.92 GHz** |

**(Z-i) vs (Z-ii) — RESOLVED: no divergence.** This branch's central question was whether the framework cutoff is Z-universal (reading **(Z-i)**) or Z-scaled (reading **(Z-ii)**). Inspection of DRQM I §III.D (Eqs. III.21–III.23) settles it: the cutoff formula

```math
g_r = 2\!\left[1 - \frac{4 r_0}{2 r + r_0}\right], \qquad r_0 = \frac{e^2}{m_e c^2},
```

contains **no nuclear charge $Z$**. The cutoff $r_e/r_0$ is keyed to the **free electron's own** classical radius $r_0=e^2/(m_ec^2)$ and is fixed by matching the **free-electron** anomalous moment $g_s$, which is intrinsic to the electron and identical in H, Li²⁺, and free space. The §III.D cross-particle extension (Eq. III.23) makes the cutoff *particle-mass*-dependent (muon, proton get their own $r_0^\ell$), but never *nuclear-charge*-dependent. **Therefore reading (Z-ii) is not derivable from the framework — it collapses to (Z-i):** the same $g_s=-2.00231930$ and the same prediction above apply at Z=3.

> **Consequence — #3 is a *weaker* Z-axis discriminator than anticipated.** Because the framework offers no mechanism for a Z-dependent cutoff, the Li²⁺ prediction is just the $Z^4$-scaled leading Dirac interval times the same $(g_s/-2)$ factor used at Z=1. It carries no new Z-axis signal that could distinguish a Z-universal from a Z-scaled cutoff — there is no Z-scaled cutoff on offer.

**Caveat (bound-state QED).** A genuine bound-electron-$g$-factor effect exists (the in-medium $g$ gains $(Z\alpha)^2$ corrections in the nuclear field), but §III.D does not model it — it models the *free* anomalous moment. At Z=3 this correction to the ~1,028 MHz anomalous piece is $\sim(Z\alpha)^2 a_e\times$(that piece) ≈ sub-MHz, well within the Bethe-estimate floor, so it does not alter the verdict.

**Measurement provenance.** No clean high-precision *direct experimental* value of the hydrogenic Li²⁺ 2²P₃/₂–2²P₁/₂ fine-structure interval was located (two targeted literature searches, 2026-05-27). The references commonly cited in this context are for *different species* — neutral Li (≈10 GHz) and He-like Li⁺ (Riis et al., *Phys. Rev. A* **49**, 207 (1994), 1s2s–1s2p). The comparison value for #3 is therefore the **theoretical** hydrogenic-ion result: the $Z^4$-scaled Sommerfeld–Dirac + QED interval, ≈ 888 GHz (leading Dirac 886,892 MHz + anomalous + recoil + two-loop), which is itself extremely well-established for hydrogenic ions.

**Verdict: ✅ outcome A — at the framework's Bethe-estimate precision, read as back-fit self-consistency, not independent corroboration.**

- The framework reproduces the leading Dirac $Z^4$ interval exactly (886,892 MHz, identical to textbook).
- The anomalous-$g$ correction uses the textbook leading-$g_s$ formula $(g_s/-2)\times E_{\text{leading}}$ with the *measured/triangulated* $g_s$ — so the agreement is the same self-consistency flagged in [BS-§14.2](03_FineStructure.md#result-bs-142--2p--2p-fine-structure-splitting-branched) and [`10_CrossComparison.md` §2](10_CrossComparison.md#2-the-r_e-self-consistency-across-six-g_s-dependent-observables-at-the-triangulated-value), now extrapolated to Z=3 under the (Z-universal) cutoff. It is **not** an independent discrimination of the dual theory from standard QED.
- Two-loop + recoil corrections (out of Bethe-estimate scope) set the residual; the comparison is framework-vs-textbook-QED (both $Z^4$ Dirac × the same $(g_s/-2)$), limited by the sparseness of direct experimental data for this specific interval.

<!-- TODO: human reviews and fills in — confirms (a) the n_FS=1 spin-orbit scaling and the 887,920 MHz framework value, (b) the (Z-ii)≡(Z-i) determination from §III.D's Z-independent cutoff and the resulting "weak Z-discriminator" reframe, (c) the measurement-provenance disposition (no direct hydrogenic-Li²⁺ 2P-FS precision experiment located; comparison vs theoretical ~888 GHz), and (d) the back-fit-self-consistency verdict framing carried over from BS-§14.2. -->

---
