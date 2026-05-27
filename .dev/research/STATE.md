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

## Iteration 2 — 2026-05-27 — ³He⁺ (Z=2) g-factor measurement + provenance (ion 1/5)

**Step taken:** Looked up and recorded the precise ³He⁺ bound-electron g-factor with full provenance (one-ion measurement step). Verified via web search against the primary literature.

**⚠ Brief-table correction:** `brief.md` listed `2.000 008 021(15)` for ³He⁺ — this is a **transcription error** (unphysical: implies a binding shift of $\sim2\times10^{-3}$ off the free-electron value, but the He⁺ Breit/Dirac binding shift is only $-\tfrac13(Z\alpha)^2 \approx -7.1\times10^{-5}$, giving $|g|\approx2.00218$, not $2.000008$). Corrected value below. The brief's "Hoffmann 1989 / Köhler 2015" pointer is also superseded — the definitive direct measurement is Schneider 2022.

**³He⁺ (Z=2) — recorded value:**
- **Experimental:** $g_e^{\rm exp}(^3{\rm He}^+) = -2.002\,177\,415\,79(34)_{\rm stat}(30)_{\rm sys}$ → combined $\sigma = 45\times10^{-11}$, i.e. $-2.00217741579(45)$.
- **Theory (bound-state QED):** $g_e^{\rm theo} = -2.002\,177\,416\,252\,23(39)$ (exp–theory agree at $\sim5\times10^{-10}$).
- **Source:** A. Schneider, B. Sikora, S. Dickopf, M. Müller, N. S. Oreshkina, A. Rischka, I. A. Valuev, S. Ulmer, J. Walz, Z. Harman, C. H. Keitel, A. Mooser, K. Blaum, *"Direct measurement of the ³He⁺ magnetic moments,"* **Nature 606, 878–883 (2022)**. DOI: **10.1038/s41586-022-04761-7**. (First direct high-precision Penning-trap measurement; ~10× better than prior indirect results.)

**Framework-relevant derived quantities (to be Wolfram-verified at the per-ion / joint-fit step):**
- $a_e^{\rm bound}(Z{=}2) = (|g_e|-2)/2 = 0.001\,088\,707\,895$ — note this is **smaller** than free-electron $a_e = 0.001\,159\,652\,18$, by the binding correction $-\tfrac13(Z\alpha)^2$ ($(Z\alpha)^2 = 2.130\times10^{-4}$ at Z=2).
- (Z-ii) back-fit estimate: $r_e^{(Z=2)}/r_0 = (2-a)/(2(2+a)) \approx 0.499\,456$ (hand estimate; **must verify in `r_e_Zscan_fit.wl`**), vs Z=1 triangulated $0.499\,420\,510$ — an **upward drift of $\sim3.5\times10^{-5}$**.
- **Direction of evidence:** the back-fit $r_e^{(Z)}/r_0$ rises with Z because bound-state QED *reduces* the effective anomaly. This is the QED-bound-state $a_e(Z\alpha)$ inheritance signature — consistent with Outcome **C** (per-Z inheritance, no framework-internal Z-derivation), or **B** if the drift fits a clean $a + b(Z\alpha)^2$ form (it should, since the leading binding term is exactly $-\tfrac13(Z\alpha)^2$).

- **Current ion focus:** ³He⁺ (Z=2) — done.
- **Next:** look up + record ⁹Be³⁺ (Z=4) bound-electron g-factor with DOI/year provenance (likely a Mainz/MPIK Penning-trap result; verify the current best value + source). Per the brief, "(verify in iter 1)" — no value pre-listed, so this needs a genuine lookup.
- **Outcome-matrix tentative:** **C-leaning** (with B as a live possibility given the clean $-\tfrac13(Z\alpha)^2$ leading term); A appears disfavoured already given the Z=2 drift, but not falsified until the joint χ² is run.
- **Status:** READY (1/5+ ions catalogued: ³He⁺ Z=2 ✅; Li²⁺ Z=3 to be imported from #78).

## Iteration 3 — 2026-05-27 — ⁹Be³⁺ (Z=4) resolved: no precision data → drop + substitute ²⁰Ne⁹⁺ (Z=10)

**Step taken:** Attempted the ⁹Be³⁺ (Z=4) bound-electron g-factor lookup (the brief's "(verify in iter 1)" placeholder). Web-searched for a Be³⁺ Penning-trap measurement; none located. Cross-checked against my knowledge of the measured-ion set.

**Finding — Be³⁺ has no precision bound-electron g-factor measurement:**
- Hydrogen-like ions with high-precision Penning-trap bound-electron g-factors (the established set): ¹²C⁵⁺ (Häffner 2000), ¹⁶O⁷⁺ (2004), ²⁸Si¹³⁺ (Sturm 2011/2013), ⁴⁰Ca¹⁹⁺ (2019), ²⁰Ne⁹⁺ (Sailer 2022), ³He⁺ (Schneider 2022), ¹¹⁸Sn⁴⁹⁺ (Morgner 2023), ²⁰⁸Pb⁸¹⁺ (2023/24).
- **⁹Be³⁺ is NOT in this set.** Beryllium appears in trap physics as the singly-charged ⁹Be⁺ laser-cooling / quantum-logic coolant ion, not as hydrogenic Be³⁺. The brief's "various Penning-trap groups" pointer for Be³⁺ does not resolve to an actual measurement.

**Decision (substantive — ion-list refinement):** **Drop ⁹Be³⁺ (Z=4)** from the Z-scan for lack of measurement provenance (the brief explicitly authorises per-iteration refinement of the target list). **Substitute ²⁰Ne⁹⁺ (Z=10).** Rationale: (a) it has the **most precise** bound-electron g-factor measured to date (coupled-ion differential method, $\sim10^{-11}$ relative), giving the fit a high-weight mid-Z anchor; (b) it improves the Z lever-arm distribution — revised scan is **Z = {2, 3, 6, 10, 14, 20}**, a near-even spread for resolving any $(Z\alpha)^2$ curvature, vs the original {2,3,4,6,14,20} which clustered low.
- Provisional source for the substitute (to verify + record next iter): T. Sailer, V. Debierre, Z. Harman, F. Heiße, C. König, J. Morgner, B. Tu, A. V. Volotka, C. H. Keitel, K. Blaum, S. Sturm, *"Measurement of the bound-electron g-factor difference in coupled ions,"* **Nature 606, 479–483 (2022)**, DOI 10.1038/s41586-022-04807-w.

- **Current ion focus:** Be³⁺ (resolved → dropped); ²⁰Ne⁹⁺ queued.
- **Next:** look up + record the precise ²⁰Ne⁹⁺ (Z=10) bound-electron g-factor absolute value + uncertainty + DOI (Sailer 2022). Note: that paper's headline is a *difference/isotope-shift* measurement — confirm whether it reports an absolute ²⁰Ne⁹⁺ g-factor or whether the absolute value must be sourced from the companion theory/earlier work; record provenance accordingly.
- **Outcome-matrix tentative:** **C-leaning** (B live), unchanged from iter 2.
- **Status:** READY (1/5+ ions catalogued; revised target Z = {2, 3, 6, 10, 14, 20}; Be³⁺ dropped, Ne⁹⁺ substituted).

## Iteration 4 — 2026-05-27 — ²⁰Ne⁹⁺ resolved: differential-only → finalize 6-ion absolute set Z={2,3,6,8,14,20}

**Step taken:** Looked up the ²⁰Ne⁹⁺ (Z=10) absolute bound-electron g-factor (queued from iter 3). Two web searches. Finding forces an ion-set decision.

**Finding — Sailer 2022 is a *difference* measurement; no usable absolute ²⁰Ne⁹⁺ g-factor:**
- Sailer et al., *Nature* **606**, 479 (2022), DOI 10.1038/s41586-022-04807-w, measures the **isotope-shift difference** $g(^{20}\mathrm{Ne}^{9+}) - g(^{22}\mathrm{Ne}^{9+})$ to $5.6\times10^{-13}$ via the coupled-ion / common-magnetron-orbit method. That precision is *differential* — it cancels the common magnetic-field systematic. The **absolute** individual g-factors are limited by field calibration and are **not published at a precision useful as an absolute Z-scan anchor**.
- Consequence: Ne⁹⁺ cannot serve the (Z-ii) per-Z back-fit, which requires an absolute measured $g_e^{\rm bound}(Z)$ to invert $r_e^{(Z)}/r_0 = (2-a)/(2(2+a))$. **Drop ²⁰Ne⁹⁺** as a back-fit anchor.

**Decision (substantive — ion set now FINAL):** Use the six hydrogenic ions with **published absolute** Penning-trap bound-electron g-factors:

| Ion | Z | Absolute g status | Source (to transcribe with full provenance) |
|---|---|---|---|
| ³He⁺ | 2 | ✅ recorded (iter 2) | Schneider 2022, Nature 606, 878 |
| ⁷Li²⁺ | 3 | import from #78 | (Self-Energy branch supplies) |
| ¹²C⁵⁺ | 6 | pending transcribe | Sturm 2011, PRL 107, 023002 (brief table) |
| ¹⁶O⁷⁺ | 8 | **new — substitute for Ne⁹⁺** | Verdú et al. 2004, PRL 92, 093002 |
| ²⁸Si¹³⁺ | 14 | pending transcribe | Sturm 2013, PRL 110, 263002 (brief table) |
| ⁴⁰Ca¹⁹⁺ | 20 | pending lookup | Köhler-Langes/Sturm 2016 |

Final Z-scan: **Z = {2, 3, 6, 8, 14, 20}** — all absolute, well-distributed. ¹⁶O⁷⁺ (Z=8) replaces the dropped Ne⁹⁺ for mid-Z coverage; it has a clean absolute measurement (Verdú 2004) and bridges the C⁵⁺(6)–Si¹³⁺(14) gap.

**Rationale note for the writeup:** the differential-vs-absolute distinction is itself relevant to the framework test — the (Z-ii) back-fit is a per-ion *absolute* inversion, so only absolute measurements qualify. (The Sailer isotope-shift result tests nuclear-recoil QED, orthogonal to the cutoff-radius question here.)

- **Current ion focus:** Ne⁹⁺ (resolved → dropped, differential-only); O⁷⁺ queued.
- **Next:** look up + record the ¹⁶O⁷⁺ (Z=8) absolute bound-electron g-factor + uncertainty + DOI (Verdú, Djekić, Stahl, Valenzuela, Vogel, Werth, Beier, Kluge, Quint, PRL 92, 093002 (2004), "Electronic g Factor of Hydrogenlike Oxygen ¹⁶O⁷⁺").
- **Outcome-matrix tentative:** **C-leaning** (B live), unchanged.
- **Status:** READY (1/6 ions catalogued with absolute values; ion set FINAL at Z={2,3,6,8,14,20}).

## Iteration 5 — 2026-05-27 — ¹⁶O⁷⁺ (Z=8) g-factor + provenance (ion 2/6)

**Step taken:** Looked up and recorded the ¹⁶O⁷⁺ absolute bound-electron g-factor (queued from iter 4). Verified via web search against the primary source.

**¹⁶O⁷⁺ (Z=8) — recorded value:**
- **Experimental:** $g_e^{\rm exp}(^{16}{\rm O}^{7+}) = 2.000\,047\,025\,4(15)_{\rm stat}(44)_{\rm sys}$ → combined $\sigma \approx 46\times10^{-10}$, i.e. $2.0000470254(46)$.
- **Theory (bound-state QED):** $g_e^{\rm theo} = 2.000\,047\,020\,2(6)$ (exp–theo agree at 1.1σ; 0.25% BS-QED test).
- **Source:** J. L. Verdú, S. Djekić, S. Stahl, T. Valenzuela, M. Vogel, G. Werth, T. Beier, H.-J. Kluge, W. Quint, *"Electronic g Factor of Hydrogenlike Oxygen ¹⁶O⁷⁺,"* **Phys. Rev. Lett. 92, 093002 (2004)**. DOI: **10.1103/PhysRevLett.92.093002**. (Single ¹⁶O⁷⁺ ion in a Penning trap; first calculated resonance line shape.)

**Framework-relevant derived quantities (Wolfram-verify at joint-fit step):**
- $a_e^{\rm bound}(Z{=}8) = (g-2)/2 = +0.000\,023\,512\,7$. The Z=8 binding term $-\tfrac13(Z\alpha)^2 \approx -1.136\times10^{-3}$ nearly cancels the free-electron anomaly $a_e=+1.1597\times10^{-3}$ → net $g \approx 2.0000$.
- (Z-ii) back-fit estimate: $r_e^{(Z=8)}/r_0 = (2-a)/(2(2+a)) \approx 0.499\,994$ (hand estimate; verify in `.wl`).

**Z-trend across the catalogued + projected points (the headline emerging pattern):**
| Z | source | $g_e^{\rm bound}$ | $a_e^{\rm bound}=(g-2)/2$ | back-fit $r_e^{(Z)}/r_0$ (est.) |
|---|---|---|---|---|
| 1 (free) | CODATA | 2.00231930 | +0.00115965 | 0.499420510 (triangulated) |
| 2 | Schneider22 | 2.00217742 | +0.00108871 | ≈0.499456 |
| 8 | Verdú04 | 2.00004703 | +0.00002351 | ≈0.499994 |
| 14 (proj.) | Sturm13 | 1.99534896 | −0.00232552 | ≈0.501164 |

- The back-fit $r_e^{(Z)}/r_0$ **monotonically sweeps through 0.5** as Z increases (crossing 0.5 near Z≈8–9, where the binding correction exactly cancels the free anomaly so $g=2$ and $r_e/r_0=1/2$ — the tree-Dirac value). This is a **large, clean, monotonic Z-dependence**, not scatter.
- **Outcome A (Z-universal single cutoff) is now effectively ruled out** by inspection: a fixed $r_e/r_0=0.499421$ predicts $g=-2.00231930$ at every Z, but measured $|g|$ ranges 2.0000–2.0022 across the scan — residuals of order $10^{-3}$–$10^{-4}$, i.e. $10^6$–$10^7\sigma$. The joint χ² will confirm quantitatively.
- The back-fit curve is exactly the inversion of QED's bound-state $g(Z\alpha)=2[1-\tfrac13(Z\alpha)^2-\ldots]+\tfrac{\alpha}{\pi}+\ldots$ → **Outcome C** (per-Z inheritance of QED's $a_e(Z\alpha)$; the framework supplies the g-formula but not the $-\tfrac13(Z\alpha)^2$ binding term). Distinguishing **C vs B** hinges on whether the framework can *derive* the $(Z\alpha)^2$ coefficient internally — it cannot (g-formula leaves each state's cutoff free), so C is strongly favoured; the joint fit's $a+b(Z\alpha)^2$ form-fit will characterise the curve either way.

- **Current ion focus:** ¹⁶O⁷⁺ (Z=8) — done.
- **Next:** transcribe ¹²C⁵⁺ (Z=6) absolute g-factor with full provenance — brief table lists $2.001\,041\,590\,18(3)$, Sturm 2011 *PRL* **107**, 023002; verify value + DOI (DOI 10.1103/PhysRevLett.107.023002).
- **Outcome-matrix tentative:** **C** (firming up; A effectively excluded by the Z=8 point; B requires an internal $(Z\alpha)^2$ derivation the framework lacks).
- **Status:** READY (2/6 ions catalogued: He⁺ Z=2 ✅, O⁷⁺ Z=8 ✅; remaining C⁵⁺ Z=6, Si¹³⁺ Z=14, Ca¹⁹⁺ Z=20 + Li²⁺ Z=3 import).

## Iteration 6 — 2026-05-27 — ¹²C⁵⁺ (Z=6) g-factor + provenance fix (ion 3/6)

**Step taken:** Transcribed/verified the ¹²C⁵⁺ absolute bound-electron g-factor (queued from iter 5). Two web searches resolved a value-vs-source ambiguity.

**⚠ Brief-table source correction (provenance):** The brief assigns C⁵⁺ to "Sturm 2011 *PRL* **107**, 023002" — **that paper is the ²⁸Si¹³⁺ measurement, not C⁵⁺.** The brief's *value* ($2.00104159018(3)$) is correct and comes from the 2014 *Nature* electron-mass paper, not PRL 107. (Flag for the Si¹³⁺ iteration: re-verify that source too, since the brief lists Si¹³⁺ as "Sturm 2013 *PRL* **110**, 263002" — likely also needs checking.)

**¹²C⁵⁺ (Z=6) — recorded value:**
- **Experimental:** $g_e^{\rm bound}(^{12}{\rm C}^{5+}) = 2.001\,041\,590\,18(3)$ (relative $\sigma \approx 1.5\times10^{-11}$). Supersedes Häffner et al. 2000 ($2.001\,041\,596(5)$, *PRL* **85**, 5308) — the 2014 value is ~1.2σ lower and ~170× more precise.
- **Source:** S. Sturm, F. Köhler, J. Zatorski, A. Wagner, Z. Harman, G. Werth, W. Quint, C. H. Keitel, K. Blaum, *"High-precision measurement of the atomic mass of the electron,"* **Nature 506, 467–470 (2014)**, DOI **10.1038/nature13026**. Refined analysis: F. Köhler et al., *"The electron mass from g-factor measurements on hydrogen-like carbon ¹²C⁵⁺,"* **J. Phys. B 48, 144032 (2015)**, DOI 10.1088/0953-4075/48/14/144032.

**Framework-relevant derived quantities (Wolfram-verify at joint-fit step):**
- $a_e^{\rm bound}(Z{=}6) = (g-2)/2 = +0.000\,520\,795\,09$. Binding term $-\tfrac13(Z\alpha)^2 \approx -6.39\times10^{-4}$ partially cancels free $a_e$.
- (Z-ii) back-fit estimate: $r_e^{(Z=6)}/r_0 = (2-a)/(2(2+a)) \approx 0.499\,739$ — sits between Z=2 (0.499456) and Z=8 (0.499994), as required by the monotonic sweep.

**Updated Z-trend table (3 absolute points + free + Si projection):**
| Z | source | $g_e^{\rm bound}$ | $a_e^{\rm bound}$ | back-fit $r_e^{(Z)}/r_0$ (est.) |
|---|---|---|---|---|
| 1 (free) | CODATA | 2.00231930 | +0.00115965 | 0.499420510 |
| 2 | Schneider22 | 2.00217742 | +0.00108871 | ≈0.499456 |
| 6 | Sturm14 | 2.00104159 | +0.00052080 | ≈0.499739 |
| 8 | Verdú04 | 2.00004703 | +0.00002351 | ≈0.499994 |
| 14 (proj.) | (verify) | 1.99534896 | −0.00232552 | ≈0.501164 |

- **Current ion focus:** ¹²C⁵⁺ (Z=6) — done.
- **Next:** look up + record ²⁸Si¹³⁺ (Z=14) absolute g-factor with full provenance. Brief value $1.995\,348\,958\,7(5)$; **verify the source** (brief cites "Sturm 2013 *PRL* **110**, 263002" and also "Köhler 2016 *Nat. Comm.* **7**, 10246" — determine the correct primary; note Sturm 2011 *PRL* **107**, 023002 is the original Si¹³⁺ paper).
- **Outcome-matrix tentative:** **C** (A excluded; B needs an internal $(Z\alpha)^2$ derivation the framework lacks).
- **Status:** READY (3/6 ions catalogued: He⁺ Z=2 ✅, C⁵⁺ Z=6 ✅, O⁷⁺ Z=8 ✅; remaining Si¹³⁺ Z=14, Ca¹⁹⁺ Z=20 + Li²⁺ Z=3 import).

## Iteration 7 — 2026-05-27 — ²⁸Si¹³⁺ (Z=14) g-factor + provenance fix (ion 4/6)

**Step taken:** Looked up and recorded the ²⁸Si¹³⁺ absolute bound-electron g-factor (queued from iter 6). Web search confirmed value + resolved the source attribution.

**⚠ Brief-table source correction (provenance):** The brief assigns Si¹³⁺ to "Sturm 2013 *PRL* **110**, 263002" — the **correct primary is Sturm et al., *PRL* 107, 023002 (2011)** (the same paper the brief had wrongly attributed to C⁵⁺ in iter 6). Net: the brief swapped/mis-cited both C⁵⁺ and Si¹³⁺ sources; both now corrected (C⁵⁺ → Nature 506 2014; Si¹³⁺ → PRL 107 2011).

**²⁸Si¹³⁺ (Z=14) — recorded value:**
- **Experimental:** $g_e^{\rm bound}(^{28}{\rm Si}^{13+}) = 1.995\,348\,958\,7(5)_{\rm stat}(3)(8)_{\rm sys}$ → combined $\sigma \approx 10\times10^{-10}$, i.e. $1.9953489587(10)$ (relative $\sim5\times10^{-10}$; 10 significant digits).
- **Theory (BS-QED, 2-loop):** $g_e^{\rm theo} = 1.995\,348\,958\,0(17)$ (excellent agreement; "most stringent test of BS-QED" at the time).
- **Source:** S. Sturm, A. Wagner, B. Schabinger, J. Zatorski, Z. Harman, W. Quint, G. Werth, C. H. Keitel, K. Blaum, *"g Factor of Hydrogenlike ²⁸Si¹³⁺,"* **Phys. Rev. Lett. 107, 023002 (2011)**. DOI: **10.1103/PhysRevLett.107.023002**.

**Framework-relevant derived quantities (Wolfram-verify at joint-fit step):**
- $a_e^{\rm bound}(Z{=}14) = (g-2)/2 = \mathbf{-0.002\,325\,520\,65}$ — **NEGATIVE**: at Z=14 the binding term $-\tfrac13(Z\alpha)^2 \approx -7.43\times10^{-3}$ exceeds the free anomaly $+1.16\times10^{-3}$, so $g<2$.
- (Z-ii) back-fit estimate: $r_e^{(Z=14)}/r_0 = (2-a)/(2(2+a)) \approx 0.501\,164$ — **>0.5**, matching the iter-5/6 projection exactly. **The measured-data back-fit confirms the curve crosses 0.5** (at the Z where binding cancels the free anomaly, $g=2 \Leftrightarrow r_e/r_0=1/2$, near Z≈8–9).

**Updated Z-trend table (4 absolute points + free):**
| Z | source | $g_e^{\rm bound}$ | $a_e^{\rm bound}=(g-2)/2$ | back-fit $r_e^{(Z)}/r_0$ (est.) |
|---|---|---|---|---|
| 1 (free) | CODATA | 2.00231930 | +0.00115965 | 0.499420510 |
| 2 | Schneider22 | 2.00217742 | +0.00108871 | ≈0.499456 |
| 6 | Sturm14 | 2.00104159 | +0.00052080 | ≈0.499739 |
| 8 | Verdú04 | 2.00004703 | +0.00002351 | ≈0.499994 |
| 14 | Sturm11 | 1.99534896 | −0.00232552 | ≈0.501164 |
| 20 (pending) | (lookup) | ≈1.989 | (neg, larger) | (>0.5018) |

- **Current ion focus:** ²⁸Si¹³⁺ (Z=14) — done.
- **Next:** look up + record ⁴⁰Ca¹⁹⁺ (Z=20) absolute g-factor + uncertainty + DOI (brief cites "Glazov 2019 / Köhler-Langes 2018"; verify the primary measurement — likely Köhler-Langes/Sturm/Vogel; note theory by Glazov et al.).
- **Outcome-matrix tentative:** **C** (A excluded — measured data now span $g$=1.9953→2.0022, a $4.7\times10^{-3}$ spread vs the (Z-i) single-value prediction; B needs an internal $(Z\alpha)^2$ derivation the framework lacks).
- **Status:** READY (4/6 ions catalogued: He⁺ Z=2 ✅, C⁵⁺ Z=6 ✅, O⁷⁺ Z=8 ✅, Si¹³⁺ Z=14 ✅; remaining Ca¹⁹⁺ Z=20 + Li²⁺ Z=3 import).
