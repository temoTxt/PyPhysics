# Hydrogenic Z-scan g-factor — state log

## Iteration 0 — 2026-05-27 — initialized

- Branch `82-hydrogenic-z-scan-g-factor` checked out from main (post-#67 closure).
- Folder `PyPyshics-Thread_4` cloned fresh from origin for this thread.
- `.dev/research/brief.md` written (Z-scan across 4–6 hydrogenic ions, g-factor only; Li²⁺ imported from #78 parallel work).
- `.dev/research/loop_prompt.md` written.
- No prediction work yet.
- **Current ion focus:** none yet.
- **Next:** read source-of-record §1 (`Dual_Relativistic_Quantum_Mechanics_I.md` §III.D Eqs. III.22/III.23) — record the anomalous-g formula and the framework's published precedent for particle/scale variation. Then look up + record the precise He⁺ ($Z=2$) g-factor measurement value with DOI/year provenance.
- **Outcome-matrix:** not yet determinable.
- **Status:** READY.

## Iteration 1 — 2026-05-27 — read source-of-record §1 (DRQM I §III.D g-formula + particle-variation precedent)

**Step taken:** Read `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` Eqs. (III.21)–(III.23) + §III.D-extension (lines 451–674). Recorded the g-formula, the cutoff–anomaly inversion, and the framework's published precedent for particle/scale variation.

**Key formulas recorded:**
- **g-formula (III.21–22):** $g_r = 2\left[1 - \dfrac{4r_0}{2r+r_0}\right] = 2\left[1 - \dfrac{4}{2x+1}\right]$, $\;x \equiv r/r_0$, with $r_0 = e^2/(mc^2)$ (classical-radius scale). Limit checks: $g_r(x{=}1/2) = -2$ (tree-Dirac), $g_r(x{\to}0) = -6$. Both verified in-doc.
- **Cutoff–anomaly inversion (§III.D-extension):** identifying $g_r(r_e/r_0) = -2(1+a_e)$ with $a_e \equiv (|g_e|-2)/2$ gives the closed form $\boxed{r_e/r_0 = (2-a_e)/(2(2+a_e))}$. With CODATA-full $a_e^{\rm expt}$ → $r_e/r_0 = 0.4994205099128317$ ($\sigma_r = 2.5\times10^{-13}$). This is the **Z=1 universal-cutoff value** to use for the (Z-i) test.
- **Schwinger one-loop closed form (Branch B):** $r_e/r_0|_{\rm Schwinger} = (2-\alpha/(2\pi))/(4+\alpha/\pi) = 0.49941963215699$, off triangulated by $\Delta r = +8.78\times10^{-7}$.
- **Sensitivity:** $dg_r/dx \approx 4.0046$ at the cutoff — so $\sigma_x \approx \sigma_g/4$. Carries to the per-Z back-fit error propagation.

**Particle/scale-variation precedent (the load-bearing read for the Z-axis):**
- **Eq. (III.23):** muon and proton get their **own free cutoffs** $r_\mu, r_p$ with $r_0^\mu = e^2/(m_\mu c^2)$, $r_0^p = e^2/(m_p c^2)$; same dimensionless g-formula, but paper specifies **no numerical $r_\mu, r_p$** — they are left as separate free parameters per particle.
- **PR #70 cross-lepton test (lines 548–552):** back-fitting $r_\mu/r_0^\mu$ from $a_\mu^{\rm exp}$ (FNAL 2023, $a_\mu = 116592059(22)\times10^{-11}$, PRD 108 092009) gives $0.499417379350$ vs electron $0.499420509913$ — differ by $3.13\times10^{-6}$ (~57 kσ in muon units). **Conclusion: a universal dimensionless cutoff $r/r_0$ across particles is ruled out at >57 kσ; the cutoff is particle-specific through $a_\ell$.** Not a falsification of the paper (which leaves them free), but it constrains any future closed-form to be mass-dependent.

**Z-axis mapping (how the lepton-axis precedent transposes to Z):** for a bound electron in hydrogenic ion of charge $Z$, the analogue is $r_e^{(Z)}/r_0 = (2 - a_e^{\rm bound}(Z\alpha))/(2(2 + a_e^{\rm bound}(Z\alpha)))$. The (Z-i) test asks whether the *same* $0.4994205099$ fits all Z (→ A); (Z-ii) inverts per-Z. By analogy to PR #70's lepton-axis verdict, the prior expectation is that the back-fit inherits QED's bound-state $a_e(Z\alpha)$ structure per-Z (→ C), but this is **not yet determinable for the Z-axis** until the measured values are fitted.

- **Current ion focus:** none (source-read iteration); He⁺ (Z=2) is next.
- **Next:** look up + record the precise ³He⁺ (Z=2) bound-electron g-factor measurement value with full DOI/year provenance (brief table lists $2.000\,008\,021(15)$, Hoffmann 1989 / Köhler 2015 update — verify which is the current best and its DOI). Enter it into the Z-scan table in STATE.md; per-ion section in `14_HydrogenicIon_Zscan.md` to follow once 2+ values are in hand.
- **Outcome-matrix tentative:** leaning **C** by analogy to PR #70 lepton-axis (particle-specific through $a_\ell$ → Z-specific through $a_e(Z\alpha)$), but **A** is the clean falsifiable test and remains open until data are fitted. Not yet determinable.
- **Status:** READY (1 source-of-record doc read; 0/5 ions catalogued).
