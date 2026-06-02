# §13 (Z-extension) — Li²⁺ 1s hyperfine splitting (observable #4)

**Issue [#78](https://github.com/temoTxt/PyPhysics/issues/78), observable #4.** This chapter extends the hydrogen 21-cm Fermi-contact apparatus of [`06_Hyperfine.md` BS-§22.1](06_Hyperfine.md#result-bs-221--fermi-contact-term-and-21-cm-line-branched) from Z=1 (proton, nuclear spin $I=1/2$) to hydrogenic **⁷Li²⁺** (Z=3, nuclear spin $I=3/2$), evaluating the dual-theory framework prediction for the 1s ground-state hyperfine splitting.

It is the hyperfine companion to the Bethe–Salpeter Z-extension campaign (parallel branches own observable #1 g-factor / joint χ² fit, #2 Lamb shift, #3 fine structure). The Z=3 hyperfine prediction value derived here is available to the joint-χ² fit owned by the Self-Energy branch; this chapter does **not** edit that fit.

Companion Wolfram notebook: [`r_e_Li2plus_hyperfine.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_hyperfine.wl) (Cells 1–5, all Wolfram-MCP-verified).

---

## Result #4 — ⁷Li²⁺ 1s hyperfine splitting

**Selection provenance.** Hyperfine is the campaign's most precision-sensitive `r_e`-dependent observable (the leading Fermi-contact term is linear in `g_s`, [`10_CrossComparison.md §2`](10_CrossComparison.md#2-the-r_e-self-consistency-across-six-g_s-dependent-observables-at-the-triangulated-value)). Li²⁺ is the campaign's **first out-of-sample** `g_s`-dependent observable: Z=3 was *not* in the six-observable joint fit that defined the triangulated cutoff. *Substantive AI; the `I=1/2 → I=3/2` angular-structure extension is flagged substantive below.*

**Source.** [Bethe–Salpeter §22](06_Hyperfine.md) Fermi-contact Hamiltonian, applied at Z=3 with Li-7 nuclear spin $I=3/2$:

```math
H_{HF} = \frac{8\pi}{3}\,g_{I}\mu_{N}\,g_{s}\mu_{B}\,\delta^{3}(\mathbf{r})\,\mathbf{I}\cdot\mathbf{S},
```

with $g_I$ the **Li-7 nuclear g-factor** (replacing the proton's $g_p$), and $\langle 1s|\delta^3(\mathbf r)|1s\rangle = Z^3/(\pi a_0^3)$ carrying the $Z^3$ hydrogenic enhancement.

### The $I=3/2$ angular structure — *substantive-AI extension of BS-§22's $I=1/2$ apparatus*

H (proton, $I=1/2$) couples $S=1/2$ to give a single two-level splitting $F=0,1$; the 21-cm line is that one interval. Li-7 ($I=3/2$) coupled to $S=1/2$ gives **$F=1,2$**, and the headline observable is the $F=2\leftrightarrow1$ interval. With $H_{HF}=a\,\mathbf I\cdot\mathbf S$ and $\langle\mathbf I\cdot\mathbf S\rangle = \tfrac12[F(F+1)-I(I+1)-S(S+1)]$ (Wolfram-MCP, Cell 4):

| $F$ | $\langle\mathbf I\cdot\mathbf S\rangle$ |
|---|---|
| 2 | $+3/4$ |
| 1 | $-5/4$ |

so the $F=2\leftrightarrow1$ splitting $= 2a = (I+\tfrac12)\,a$ (interval rule). This $(I+\tfrac12)$ angular factor, together with the nuclear-moment factor $\mu_I/I$, is the substantive change from H — *not* a pure $Z$-rescaling.

<!-- TODO: human reviews and fills in — confirms the I=1/2→3/2 angular-structure extension (F=1,2; headline interval = 2a) is the correct generalization of BS-§22's proton I=1/2 two-level apparatus, and that the nuclear-moment factor μ_I(2I+1)/(2I) is the right scaling carrier. -->

### Measurement provenance — *no direct Li²⁺ experiment exists*

There is **no direct experimental ⁷Li²⁺ (hydrogenic) 1s HFS measurement**. The authoritative comparator is the QED-theory value of **Pachucki, Patkóš & Yerokhin, *Hyperfine splitting in ⁶,⁷Li⁺*, [arXiv:2309.00436](https://arxiv.org/abs/2309.00436) (2023), Table VII**, which predicts the Li²⁺ hfs from the *experimental Li⁺* hfs plus QED theory (its uncertainty derives from the Li⁺ measurement). The early atomic-beam reference Beckmann, Boeklen & Elke, *Z. Phys.* **270**, 173 (1974) concerns neutral-Li / Li⁺, not hydrogenic Li²⁺.

Li-7 nuclear data: $I=3/2$, $\mu_I = +3.256427\,\mu_N$ (NUBASE/CODATA); $g_I = \mu_I/I = 2.1709$.

<!-- TODO: human reviews and fills in — extract the exact ⁷Li²⁺ 1s HFS value (and ⁶Li²⁺) from Pachucki et al. 2023 Table VII (PDF) and insert into the residual table below; confirm the "no direct measurement; comparator is QED-theory-from-Li⁺" provenance framing. -->

### Framework prediction (dual-theory, at the triangulated `r_e`)

The framework prediction is $(g_s/{-2})^{1}\times$ (textbook leading Fermi-contact term), evaluated at the triangulated $r_e/r_0 = 0.4994205099128317$ ([PR #62](https://github.com/temoTxt/PyPhysics/pull/62)), giving $g_s = -2.00231930$. The Li/H scaling of the textbook leading term (Wolfram-MCP, Cell 2):

```math
\frac{\Delta\nu_{\rm Li}}{\Delta\nu_{\rm H}} = \frac{\mu_I^{\rm Li}(2I_{\rm Li}+1)/(2I_{\rm Li})}{\mu_p(2I_p+1)/(2I_p)}\times Z^3 = \frac{4.3419}{5.5857}\times 27 = 20.988 .
```

Scaling the $g_s{=}{-}2$ leading-Fermi H baseline ($1418.4$ MHz, BS-§22.1) and applying the single $(g_s/{-}2)=1.001160$ factor (Cell 3) gives the **minimal** framework prediction. Adding the standard $r_e$-independent corrections the full Dirac solution carries (relativistic factor $1/[\gamma(2\gamma-1)]$, $\gamma=\sqrt{1-(Z\alpha)^2}$, $+0.072\%$; reduced-mass ratio³, $+0.140\%$; Cell 5) gives the prediction at standard-QED point-nucleus completeness.

**Wolfram MCP check** (notebook Cells 2,3,5):

```text
Cell2: scaleRatio (Li/H) = 20.9878 ; textbook Li2+ (g_s=-2 leading) = 29769.1 MHz
Cell3: framework Li2+ = 29803.6 MHz  (cross-check scaling framework-H 1420.04: 29803.5 MHz — agree to 0.1 MHz) ✅
Cell5: standard-QED point-nucleus comparator = 29866.8 MHz ; framework(minimal) − comparator = −63.2 MHz
```

### Numerical comparison

| Quantity | $\Delta\nu_{\rm HFS}$(⁷Li²⁺, 1s) | Note |
|---|---|---|
| Textbook leading Fermi ($g_s=-2$) | `29\,769` MHz | $Z^3\times$ H baseline $\times$ nuclear factor |
| **Framework, minimal ($\times g_s/{-2}$, Z-i)** | **`29\,804` MHz** | matches campaign H methodology (leading only) |
| Framework $+$ relativistic Dirac $+$ reduced-mass | `29\,867` MHz | framework inherits Dirac (BS-§14.1) + kinematics |
| **Standard-QED point-nucleus comparator** | **`29\,867` MHz** | $=$ framework with same corrections, *by construction* |
| Pachucki et al. 2023 (incl. Bohr–Weisskopf) | _[human extract from Table VII]_ | external; nuclear structure reduces by $\sim0.5$–$1\%$ |

Residual (framework minimal − standard-QED point-nucleus) $= -63$ MHz ($-0.21\%$): entirely the relativistic + reduced-mass corrections the minimal leading route omits. **Both are $r_e$-independent**, so the residual does not bear on the verdict (exact analog of the H 21-cm leading-vs-full gap).

### Verdict — **Outcome-matrix branch A** (structural self-consistency)

✅ **at the campaign's leading-`g_s` / Bethe-estimate precision floor**, under the **(Z-i) universal cutoff**: the framework reproduces the standard-QED Li²⁺ 1s HFS at Z=3.

**Structural-self-consistency caveat (load-bearing).** The framework's *only* $r_e$-dependent content is the $(g_s/{-2})=+0.116\%$ factor. That anomalous-moment enhancement is the **free-electron $g_s$ — it is Z-independent and identical to what textbook QED applies at every Z.** Under the universal cutoff (which was *fit* to reproduce the measured $g_s$), the framework applies $g_s=-2.00231930$ at Z=3 exactly as textbook QED does. So this Li²⁺ result, like the H 21-cm result, is **agreement by construction at the leading-`g_s` floor — not an independent Z=3 discrimination** of dual-theory content from standard QED. See the same caveat at [`06_Hyperfine.md §22.1`](06_Hyperfine.md#result-bs-221--fermi-contact-term-and-21-cm-line-branched) and [`10_CrossComparison.md §2`](10_CrossComparison.md#2-the-r_e-self-consistency-across-six-g_s-dependent-observables-at-the-triangulated-value).

**(Z-i)/(Z-ii) honesty.** The out-of-sample-at-Z=3 framing would be a genuine independent test *only if* the framework's $r_e$ ran with Z (**(Z-ii)**), making $g_s(Z{=}3)\neq-2.00231930$. **The framework asserts a universal cutoff (Z-i), not a Z-running one** — so (Z-ii) is a hypothetical the framework does not claim; constructing it would require framework-internal guidance (would be outcome-branch D). Branches B/C (Z-scaled or per-Z back-fit cutoff) are likewise not posited by the framework.

**Nuclear-structure floor caveat.** Li-7 has a hyperfine anomaly (Bohr–Weisskopf, finite nuclear magnetization) and electric-quadrupole structure that *reduce* the point-nucleus HFS by $\sim0.5$–$1\%$. These are **out of scope** per issue #78 acceptance criteria; the framework's point-nucleus leading-Fermi prediction does not capture them, and they sit above the framework's precision floor — recorded here as a floor caveat, not folded into the verdict.

<!-- TODO: human reviews and fills in — confirms (a) the branch-A assignment with the structural-self-consistency caveat (g_s enhancement is Z-independent, identical to textbook QED) is the correct honest disposition, mirroring the H 21-cm verdict; (b) the (Z-ii) Z-running reading is correctly identified as NOT the framework's stated claim; (c) the Bohr–Weisskopf/quadrupole floor caveat is faithfully scoped out per issue #78; and (d) the framework prediction 29 804 MHz (minimal) / 29 867 MHz (with standard corrections) is the value to contribute to the joint-χ² fit. -->

**Notes for author review.** No direct ⁷Li²⁺ 1s HFS measurement exists; the comparison is against QED theory (Pachucki et al. 2023). The dual-theory framework predicts $\approx 29.8$ GHz, reproducing the standard-QED point-nucleus value by construction (both apply the same universal $g_s$ and the same Dirac/Fermi structure). This is the hyperfine consequence of the same `r_e` finding recorded in [`FINDINGS_for_author_review.md` Finding 2](../../Equation_Verification/FINDINGS_for_author_review.md), now exercised out-of-sample at Z=3 — where it remains a self-consistency statement, since the $r_e$-dependent $g_s$ enhancement does not vary with Z under the framework's universal-cutoff stance.

---

### Cross-PR reconciliation with [PR #87](https://github.com/temoTxt/PyPhysics/pull/87) (hydrogenic-ion Z-scan, issue [#82](https://github.com/temoTxt/PyPhysics/issues/82)) — *steel-man revision 2026-05-27*

The "Branch A by structural self-consistency" verdict above assumes the framework's universal cutoff means $g_s = -2.00231930$ is the right multiplier *at every Z*. The parallel hydrogenic-ion g-factor scan ([PR #87](https://github.com/temoTxt/PyPhysics/pull/87)) shows this assumption fails empirically.

**The Z-scan data the structural argument needs to confront:**

| $Z$ | Source | Measured $g_e^{\rm bound}$ |
|---|---|---|
| 1 (free) | CODATA | $-2.00231930$ |
| 2 (³He⁺) | Schneider 2022 | $-2.00217742(45)$ |
| 6 (¹²C⁵⁺) | Sturm 2011 | $-2.00104159$ |
| 8 (¹⁶O⁷⁺) | Verdú 2004 | $-2.00004703(46)$ |

The bound-electron $g_s$ at Z=2 is already $-2.00218$, not $-2.00232$. **For Li²⁺ at Z=3, the bound-electron $g_s$ is $\sim -2.00210$** (interpolating the $-(Z\alpha)^2/3$ binding correction). The framework's HFS prediction uses $-2.00232$ — i.e., the framework's $0.116\%$ enhancement is ~$0.011$ percentage-points too large in $|g_s|$ at Z=3.

**Concrete effect on #4:** the framework's 29,804 MHz minimal prediction uses $(g_s/-2) = 1.00116$. The Z=3 bound-electron value gives $(g_s/-2) \approx 1.00105$. The framework's prediction is therefore ~$3$ MHz too high compared to what bound-state QED would predict for Li²⁺'s actual $g_e^{\rm bound}$. This is far below the Bohr–Weisskopf nuclear-structure floor ($\sim 0.5$–$1\%$ ≈ 150–300 MHz), so it doesn't shift the verdict — but it does mean the "by construction" framing is more honest as "by construction at the Bethe-estimate floor; below that floor, the framework's $g_s$ is empirically wrong at every Z > 1."

**Reconciliation.** The verdicts are mutually consistent under proper framing:

- **Branch A is the framework's stated position** — universal cutoff, identical $g_s$ at every Z; reproduces standard-QED point-nucleus HFS at the leading-$g_s$ floor.
- **Branch C is what high-precision Z=3 data would force** — the bound-state $g_s$ at Z=3 differs from the free-electron value by $\sim 10^{-4}$, which is below the BW floor for #4 but above the Bethe-estimate floor for the underlying $r_e$-engagement.
- **Same verdict shape as the lepton axis** ([PR #70](https://github.com/temoTxt/PyPhysics/pull/70)): apparatus is X-trivial (A); data requires X-specific inheritance from QED (C).

**Sharpened verdict.** "Branch A by structural self-consistency" is **correct as a statement about the framework's apparatus**; it is **not** an independent Z=3 discrimination, and the $g_s$-Z-invariance assumption it relies on is empirically falsified at Z=2 (Schneider 2022) and would fail at Z=3 if a high-precision Li²⁺ HFS measurement were available. The "agreement by construction at the Bethe-estimate floor" framing remains honest because the framework's apparatus precision doesn't extend below where the Z-mismatch shows up — but the *reason* it works is that the framework's apparatus is structurally too coarse to expose its own Z-blindness on this observable.

<!-- TODO: human reviews and fills in — confirms (a) the cross-PR reconciliation framing (Branch A structurally + Branch C operationally is the honest joint verdict), (b) the concrete Z=3 g_s-mismatch effect on #4 (~3 MHz, below BW floor but indicative), (c) the parallel to PR #70's lepton-axis verdict, and (d) the sharpened verdict prose that distinguishes "framework's apparatus says A" from "data requires C." -->
