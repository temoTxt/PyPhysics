# Bethe–Salpeter §11–§14 extension: Hydrogenic-ion bound-electron g-factor Z-scan

**Issue:** [#82](https://github.com/temoTxt/PyPhysics/issues/82) (sibling to [#78](https://github.com/temoTxt/PyPhysics/issues/78); extends [#50](https://github.com/temoTxt/PyPhysics/issues/50) Bethe–Salpeter campaign to multi-Z g-factor data).
**Companion notebook:** [`../../Mathematica_Notebooks/Quantum_Mechanics/r_e_Zscan_fit.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/r_e_Zscan_fit.wl).
**Branch:** `82-hydrogenic-z-scan-g-factor`.
**State:** updated 2026-05-27 (loop iter 11/12) — **5 hydrogenic ions catalogued with verified absolute g-factors + DOIs (Z = 2, 6, 8, 14, 50); joint χ² + Z-scaling fit executed in Wolfram MCP.** Earlier draft was a steel-man revision after the 2026-05-27 devil's-advocate review on PR [#87](https://github.com/temoTxt/PyPhysics/pull/87). Final ion set substituted ¹¹⁸Sn⁴⁹⁺ (Z=50) for the originally-targeted ⁴⁰Ca¹⁹⁺ (which has **no published hydrogenic g-factor** — the literature "Ca¹⁹⁺ g-factor" results are lithium-like, 3-electron); Li²⁺ (Z=3) import from #78 remains pending (the brief's placeholder value is unphysical). The 5-ion set spanning Z=2–50 is verdict-complete.

**Substantive AI** (Crocco rule #1). `<!-- TODO: human reviews and fills in -->` blocks throughout.

## 1. Why this doc exists — the Z-axis test the campaign needs

The original $r_e$ triangulation ([PR #62](https://github.com/temoTxt/PyPhysics/pull/62)) fitted a single $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ to six $g_s$-dependent observables, all at Z ≤ 2. The lepton-axis cross-check ([PR #70](https://github.com/temoTxt/PyPhysics/pull/70)) showed the framework's cutoff is not universal across leptons — it inherits QED's $a_\ell$ per particle. **This doc is the Z-axis analog**: does the same $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ work across the hydrogenic-ion sequence, or does it inherit QED's bound-state $a_e(Z\alpha)$ per Z?

Per Eq. (III.22), the framework's cutoff formula has **no $Z$**: $r_e/r_0 = (2 - a_e)/(2(2 + a_e))$ depends only on the free-electron anomaly $a_e$ and the free-electron $r_0 = e^2/(m_e c^2)$. The framework's (Z-i) reading is therefore that $g_e^{\rm bound}$ is the *same* at every Z — $g_e^{\rm bound}(Z) = g_r(r_e/r_0) = -2.00231930\ldots$ independent of Z.

**The empirical test:** are measured bound-electron g-factors Z-invariant?

## 2. Ion catalog — final (5 verified, Z = 2–50; Li²⁺ pending import)

| Ion | $Z$ | Measured $g_e^{\rm bound}$ | Uncertainty | Source | DOI |
|---|---|---|---|---|---|
| ³He⁺ | 2 | $-2.002\,177\,415\,79$ | $\pm 4.5\times 10^{-10}$ | Schneider et al., *Nature* **606**, 878 (2022) | [10.1038/s41586-022-04761-7](https://doi.org/10.1038/s41586-022-04761-7) |
| ⁷Li²⁺ | 3 | _[import from #78 Self-Energy branch — pending; brief placeholder $2.0000251707$ is unphysical, see below]_ | _[pending]_ | _[from #78]_ | _[from #78]_ |
| ¹²C⁵⁺ | 6 | $-2.001\,041\,590\,18$ | $\pm 1.5\times 10^{-11}$ | Sturm et al., *Nature* **506**, 467 (2014); refined Köhler et al., *J. Phys. B* **48**, 144032 (2015) | [10.1038/nature13026](https://doi.org/10.1038/nature13026) |
| ¹⁶O⁷⁺ | 8 | $-2.000\,047\,025\,4$ | $\pm 4.6\times 10^{-9}$ | Verdú et al., *Phys. Rev. Lett.* **92**, 093002 (2004) | [10.1103/PhysRevLett.92.093002](https://doi.org/10.1103/PhysRevLett.92.093002) |
| ²⁸Si¹³⁺ | 14 | $-1.995\,348\,958\,7$ | $\pm 1.0\times 10^{-9}$ | Sturm et al., *Phys. Rev. Lett.* **107**, 023002 (2011) | [10.1103/PhysRevLett.107.023002](https://doi.org/10.1103/PhysRevLett.107.023002) |
| ¹¹⁸Sn⁴⁹⁺ | 50 | $-1.910\,562\,059$ | $\pm 1.0\times 10^{-9}$ | Morgner et al., *Nature* **622**, 53 (2023) | [10.1038/s41586-023-06453-2](https://doi.org/10.1038/s41586-023-06453-2) |

**Provenance discipline (iter-by-iter web-verified against primary literature):**
- **Source corrections vs the brief's table:** the brief swapped two attributions — ¹²C⁵⁺'s precise value is from the 2014 *Nature* electron-mass paper (not "PRL 107, 023002", which is the ²⁸Si¹³⁺ paper), and ²⁸Si¹³⁺ is "PRL 107, 023002 (2011)" (not "PRL 110, 263002 (2013)"). Both corrected here.
- **Dropped ions:** ⁹Be³⁺ (no hydrogenic g-factor measurement exists — Be is a coolant ion, used as ⁹Be⁺); ²⁰Ne⁹⁺ (Sailer 2022 is a *differential* isotope-shift measurement — no absolute g-factor at useful precision); ⁴⁰Ca¹⁹⁺ (the brief's "Glazov 2019 / Köhler-Langes 2018" sources are both **lithium-like**, not hydrogenic — no published H-like Ca¹⁹⁺ g-factor).
- **Substitutions:** ¹⁶O⁷⁺ (Z=8) for mid-Z coverage; ¹¹⁸Sn⁴⁹⁺ (Z=50) as the high-Z anchor (extends the lever arm, sharpens the curvature test).
- **High-Z caveat:** at Z=50, $(Z\alpha)^2 = 0.133$, so the perturbative $g = 2[1-\tfrac13(Z\alpha)^2-\ldots]$ series is invalid; the fit uses the *measured* $g$ directly (the back-fit $r_e^{(Z)}/r_0 = (2-a)/(2(2+a))$ is exact in $a = (|g|-2)/2$ regardless of expansion order).

<!-- TODO: human reviews and fills in — confirms the ion-catalog provenance (each row's DOI verified against primary source), the three drop decisions (Be³⁺/Ne⁹⁺/Ca¹⁹⁺), and the Li²⁺ value to import from #78 once that branch's PR merges. -->

## 3. Per-ion comparison vs framework (Z-i) prediction

Framework (Z-i) prediction is the *same* value $|g_r(r_e/r_0)| = 2.002\,319\,304\,362\,56$ at every Z (no $Z$ in the cutoff formula). MCP-computed residuals and σ-counts (`r_e_Zscan_fit.wl` §S3):

| $Z$ | Measured $|g_e^{\rm bound}|$ | Residual (meas − pred) | Residual / $\sigma_{\rm meas}$ |
|---|---|---|---|
| 1 (free) | $2.002\,319\,30$ | 0 | 0 (triangulated by construction) |
| 2 | $2.002\,177\,42$ | $-1.42\times 10^{-4}$ | $-3.2\times 10^5\,\sigma$ |
| 6 | $2.001\,041\,59$ | $-1.28\times 10^{-3}$ | $-4.3\times 10^7\,\sigma$ |
| 8 | $2.000\,047\,03$ | $-2.27\times 10^{-3}$ | $-4.9\times 10^5\,\sigma$ |
| 14 | $1.995\,348\,96$ | $-6.97\times 10^{-3}$ | $-7.0\times 10^6\,\sigma$ |
| 50 | $1.910\,562\,06$ | $-9.18\times 10^{-2}$ | $-9.2\times 10^7\,\sigma$ |

**Joint test statistic (5 ions, 0 free parameters):** $\boxed{\chi^2_{\rm (Z\text{-}i)} = 1.03\times 10^{16}}$ (Wolfram MCP, `r_e_Zscan_fit.wl` §S3). **Branch A (Z-universal cutoff) is empirically refuted** — the framework's (Z-i) prediction is *the same value at every Z* and matches measurement only at Z=1 (where it was triangulated). The residual grows monotonically with Z, reaching $9\%$ at Sn⁴⁹⁺.

## 4. Back-fit $r_e^{(Z)}/r_0$ inherits QED's bound-state $a_e(Z\alpha)$ — Branch C verdict

Inverting the framework's g-formula per Z: $r_e^{(Z)}/r_0 = (2 - a_e^{\rm bound}(Z))/(2(2 + a_e^{\rm bound}(Z)))$ where $a_e^{\rm bound}(Z) = (|g_e^{\rm bound}(Z)| - 2)/2$.

MCP-computed (`r_e_Zscan_fit.wl` §S2):

| $Z$ | $a_e^{\rm bound}(Z) = (|g| - 2)/2$ | back-fit $r_e^{(Z)}/r_0$ |
|---|---|---|
| 1 (free) | $+0.001\,159\,652$ | $0.499\,420\,510$ (triangulated) |
| 2 | $+0.001\,088\,708$ | $0.499\,455\,942$ |
| 6 | $+0.000\,520\,795$ | $0.499\,739\,670$ |
| 8 | $+0.000\,023\,513$ | $0.499\,988\,244$ (sweeps through 0.5 here) |
| 14 | $-0.002\,325\,521$ | $0.501\,164\,114$ |
| 50 | $-0.044\,718\,970$ | $0.522\,870\,866$ |

The back-fit $r_e^{(Z)}/r_0$ **monotonically sweeps through 0.5** as Z increases — exactly the QED $-\frac{1}{3}(Z\alpha)^2 + \mathcal{O}((Z\alpha)^4)$ Breit/Dirac binding correction inheritance. At $Z \approx 8$, $-(Z\alpha)^2/3 \approx -1.1\times 10^{-3}$ roughly cancels the free-electron $a_e = +1.16\times 10^{-3}$, giving $a_e^{\rm bound} \approx 0$ and $r_e/r_0 \approx 1/2$ (the tree-Dirac value); by Sn⁴⁹⁺ (Z=50) the cutoff has moved $+2.9\%$ to $0.5229$.

## 5. Operational Branch B (data fits clean form) vs derivational Branch B (framework supplies the coefficient)

The back-fit data fits a clean $c_0 + c_2(Z\alpha)^2 + c_4(Z\alpha)^4$ form — this is **operational Branch B**: the Z-dependence is functionally clean and follows the textbook bound-state QED expansion. But the framework's published apparatus has **no Z-dependence at all** in the cutoff formula — the framework cannot supply the $c_2$ coefficient from its own structure. This is **derivational Branch B failing**: the data exhibits a Z-dependence the framework's apparatus is silent about.

MCP form-fit (`r_e_Zscan_fit.wl` §S4), back-fit $r_e^{(Z)}/r_0$ vs $(Z\alpha)^2$:
- **linear:** $r_e^{(Z)}/r_0 = 0.499\,383\,590 + 0.176\,393\,13\,(Z\alpha)^2$ (residuals up to $6\times10^{-5}$ — structured, needs the quartic);
- **quadratic:** $r_e^{(Z)}/r_0 = 0.499\,420\,608 + 0.166\,275\,68\,(Z\alpha)^2 + 0.074\,154\,07\,(Z\alpha)^4$.

Two diagnostics from the quadratic fit pin the inheritance precisely:
- **Intercept $c_0 = 0.499\,420\,608$** recovers the independently-determined Z=1 triangulated / free-electron cutoff $0.499\,420\,510$ **to $9.8\times10^{-8}$** — i.e. the Z-scan extrapolates to the free-electron value as $Z\to0$, as inheritance requires.
- **Slope $c_2 = 0.166\,275\,68 \approx 1/6 = 0.166\,667$** (to $2.3\times10^{-3}$) is exactly QED's leading prediction: $\dfrac{d(r_e/r_0)}{d(Z\alpha)^2} = \dfrac{d(r_e/r_0)}{da}\cdot\dfrac{da}{d(Z\alpha)^2} = \left(-\tfrac12\right)\left(-\tfrac13\right) = +\tfrac16$ at $a\approx0$ (the $-\tfrac13$ is QED's Breit/Dirac binding coefficient; the $-\tfrac12$ is $d(r_e/r_0)/da$ of the framework's inversion at $a\approx0$).

The distinction matters for the synthesis verdict:

- **If "Branch B" means "data fits a clean form":** ✅ confirmed — the $c_0 + c_2(Z\alpha)^2 + c_4(Z\alpha)^4$ fit is excellent, $c_0$ = free cutoff, $c_2 \approx 1/6$.
- **If "Branch B" means "framework derives the coefficient":** ✗ refuted — the $-\tfrac13$ in $c_2$ is **QED's** bound-state binding coefficient, inherited through the inversion; the framework supplies only $g_r[x]$ and has no internal mechanism to produce $-(Z\alpha)^2/3$.

The campaign's master outcome-matrix [(per #67)](https://github.com/temoTxt/PyPhysics/issues/67) uses Branch B to mean *derivational* (the framework produces the form). On that definition, **Branch C is the correct verdict** — per-Z inheritance from QED's bound-state $a_e(Z\alpha)$, with no framework-internal Z-derivation.

## 6. Verdict — Branch C (per-Z inheritance), Branch A structurally refuted

**Headline verdict: Branch C.** The hydrogenic-ion Z-scan data (5 ions, Z=2–50) forces the framework into per-Z back-fits that inherit QED's bound-state $a_e(Z\alpha)$ structure. The framework's (Z-i) universal cutoff is empirically refuted at **$\chi^2 = 1.03\times10^{16}$** (5 ions, 0 free params; MCP `r_e_Zscan_fit.wl` §S3), with per-Z residuals from $3\times10^5\,\sigma$ (He⁺) to $9\times10^7\,\sigma$ (Sn⁴⁹⁺). The back-fit values follow the QED $-(Z\alpha)^2/3$ binding correction (form-fit slope $c_2\approx1/6$, intercept = free cutoff to $10^{-7}$; §S4) with **no framework-internal derivation** of the coefficient.

**Same verdict shape as [PR #70's](https://github.com/temoTxt/PyPhysics/pull/70) lepton-axis verdict.** The framework's apparatus says the cutoff is universal-in-X (across leptons / across Z); the data shows the cutoff is X-specific through QED's anomalous-moment inheritance ($a_\ell$ / $a_e^{\rm bound}(Z\alpha)$). The framework's apparatus is structurally inadequate to capture this Z-dependence.

## 7. Cross-PR reconciliation with the three #78 branches

Three parallel branches ([PR #84](https://github.com/temoTxt/PyPhysics/pull/84) Self-Energy, [PR #85](https://github.com/temoTxt/PyPhysics/pull/85) Spectroscopy, [PR #86](https://github.com/temoTxt/PyPhysics/pull/86) Hyperfine) concluded **"Branch A by construction"** for the Li²⁺ four-observable set, on the structural observation that no $Z$ appears in the framework's cutoff formula.

This branch's Z-scan data appears to refute that finding ($10^5$–$10^7\,\sigma$). The two findings are reconcilable under proper framing:

- **The three #78 branches are correct about the framework's apparatus** — the cutoff formula has no $Z$; the framework's stated position is Branch A. They are *correct as statements about what the framework's published structure asserts*.
- **This branch is correct about the empirical content** — the multi-Z data refutes the universal cutoff at $> 10^5\,\sigma$; **Branch C is the operationally correct verdict** for any genuine empirical test.
- **The reconciliation:** the framework's apparatus is Z-trivial (the three #78 branches' "by construction" finding); the data demonstrates that this Z-triviality forces Branch C (this branch's empirical finding). Saying "Branch A" *about the apparatus* and "Branch C" *about the empirical comparison* is the honest joint verdict.

**The three #78 branches' verdicts should be read as:** "the framework's apparatus has no Z-content to test, so the Z=3 Li²⁺ measurements don't independently discriminate the framework's $r_e$." This is true and useful. **It is *not* a finding that the framework agrees with measurement at Z=3** — the cutoff-engaging Li²⁺ observables in the three #78 branches had no valid measurements, so the comparison was BLOCKED, not ✅.

## 8. Open items

**Done (loop iters 1–12):**
- ✅ Joint $\chi^2$ computed explicitly in `r_e_Zscan_fit.wl` §S3: $\chi^2_{\rm (Z\text{-}i)} = 1.03\times10^{16}$ (5 ions, replaces the prior "$10^5$–$10^7\,\sigma$" estimates with a single citable statistic).
- ✅ Back-fit fitted to $c_0 + c_2(Z\alpha)^2 + c_4(Z\alpha)^4$ (§S4): $c_0 = 0.499\,420\,608$ (= free cutoff to $10^{-7}$), $c_2 = 0.166\,28 \approx 1/6$, $c_4 = 0.074\,15$ — confirms QED $-\tfrac13$ Breit inheritance.

**Remaining (not blocking PR merge):**
- Li²⁺ Z=3 import from [PR #84](https://github.com/temoTxt/PyPhysics/pull/84) once that branch merges (the brief's placeholder $2.0000251707$ is unphysical — the Z=3 binding trend requires $g\approx2.00200$; the correct value must come from #78, not be re-derived here). Slotting it in is a precision refinement; it does not change the Z=2–50 verdict.
- Append the Z-scan verdict to [`FINDINGS_for_author_review.md` Finding 2](../../Equation_Verification/FINDINGS_for_author_review.md) (parallel to the lepton-axis update from [PR #70](https://github.com/temoTxt/PyPhysics/pull/70)). _(Next loop step.)_
- Cross-link to [#75](https://github.com/temoTxt/PyPhysics/issues/75) (framework-specs, Tepper engagement) — the derivational-Branch-B question is the open theoretical thread that would resolve C → B.

<!-- TODO: human reviews and fills in — confirms (a) the Branch C verdict and the cross-PR reconciliation framing, (b) the operational vs derivational Branch B distinction, (c) the parallel to PR #70's lepton-axis verdict, and (d) the link to #75 as the open derivational question. -->
