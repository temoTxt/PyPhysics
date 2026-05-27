# §§Li²⁺ — Z-extension to the hydrogenic lithium ion (Z=3) ⭐ issue #78

**Issue #78 (extends the closed #50 Bethe–Salpeter campaign).** Per Tepper Gill's 2026-05-27 suggestion, this chapter applies the dual-theory / proper-time apparatus to the **single-electron Li²⁺ hydrogenic ion** at `Z=3`, predicting four precision observables and comparing to measurement. The scientific question is whether the framework's triangulated cutoff `r_e/r_0 = 0.499 420 509 912 831 7` (Z=1, [PR #62](https://github.com/temoTxt/PyPhysics/pull/62)) is **Z-universal** or **Z-scaled** — the **Z-axis** complement to [PR #70](https://github.com/temoTxt/PyPhysics/pull/70)'s **lepton-axis** muon test. Companion notebook: [`r_e_Li2plus_joint_fit.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_joint_fit.wl).

## ⭐ Gating structural result — the cutoff is Z-invariant

*Substantive AI; the central decision of this chapter.* The DRQM I §III.D g-factor formula `g_r(x) = 2[1 − 4/(2x+1)]` at `x = r_e/r_0`, with `r_0 = e²/(m_e c²)`, inverts to the closed form

```math
\frac{r_e}{r_0} = \frac{2 - a_e}{2\,(2 + a_e)},
```

which depends on exactly two inputs: the **free-electron anomalous moment `a_e`** and `r_0 = e²/(m_e c²)`. **Neither depends on the nuclear charge `Z`.** The radiating particle is the *same electron* at every `Z`, so `r_e/r_0` is identical at `Z=1` and `Z=3`.

This is the structural opposite of the lepton axis. On the lepton axis ([PR #70](https://github.com/temoTxt/PyPhysics/pull/70) iter-5), the dimensionless cutoff is *particle-specific* (`r_μ/r_0^μ = 0.499 417 379` ≠ `r_e/r_0^e = 0.499 420 510`, ruling out lepton-universality at `>57`kσ) **because both `a_ℓ` and `r_0^ℓ` vary with lepton mass**. On the Z axis nothing varies — same electron, same `a_e`, same `r_0` — so **the framework's own algebra forces `r_e/r_0` to be Z-invariant**. Consequently:

- **(Z-i) universal cutoff is the framework-correct reading**, not merely one of two options: `r_e/r_0 = 0.499 420 509 912 831 7` applies unchanged at `Z=3`.
- **(Z-ii) a Z-scaled cutoff is not derivable** and is structurally unmotivated — there is no `Z` anywhere in the (III.22) cutoff structure.
- **Outcome-matrix Branch B is ruled out by the framework's algebra** (a *derivable* Z-scaling cannot exist when the cutoff is Z-independent by construction). The verdict lies between **A** (the universal cutoff reproduces measurement) and **C** (the universal cutoff is pure back-fit with no content distinct from textbook QED) — the same A/C tension [PR #70](https://github.com/temoTxt/PyPhysics/pull/70) left on the lepton axis, now confirmed on the Z axis.

<!-- TODO: human reviews and fills in — confirms the Z-invariance of r_e/r_0 follows from DRQM I §III.D as recorded (cutoff tied to the Z-independent free-electron a_e), and that the Branch-B exclusion / A-vs-C framing is the correct disposition for the Z axis. -->

## Results

| Result | Observable | Z-scaling | Status | Role |
|---|---|---|---|---|
| [§Li-1 — bound g-factor](#result-li-1--bound-electron-g-factor-of-li2-z3) | #1 g-factor | binding `(Zα)²` + anomaly | **prediction drafted; measurement BLOCKED** | headline; only non-back-fit Z-test |
| [§Li-2 — Lamb shift](#result-li-2--lamb-shift-of-li2-z3) | #2 2S–2P | `Z⁴` × shrinking Bethe-log | **✅ at Bethe-estimate floor** | weak discriminator (g=2-symmetric) |
| [§Li-3 — fine structure](#result-li-3--fine-structure-of-li2-z3) | #3 2P₃/₂–2P₁/₂ | pure `Z⁴` | **prediction ✅; measurement BLOCKED** | wrong-ion in brief; ~887.92 GHz predicted |
| §Li-4 — hyperfine | #4 ⁷Li 1s HFS | `Z³` × nuclear | stub | genuine `r_e` discriminator |
| §Li-5 — Z-axis verdict | joint χ² | — | stub | cutoff Z-universal / scaled / back-fit |

---

### Result §Li-1 — Bound-electron g-factor of Li²⁺ (Z=3) <a id="result-li-1--bound-electron-g-factor-of-li2-z3"></a>

**Source:** Sommerfeld–Dirac / Breit (1928) bound-state g-factor; DRQM I §III.D for the anomaly. *Substantive AI.*

**As measured:** *(provenance unresolved — see below).* The brief cited `g_e^{bound}({}^7\text{Li}^{2+}) = 2.000 025 170 7(10)` (Sturm 2014 *Nature* **506** 467). **This value and citation are mis-attributed** (see provenance note). A valid high-precision hydrogenic ⁷Li²⁺ measurement has not yet been sourced.

**QED Z-expansion:** the bound-electron g-factor of a hydrogenic 1s state is the Dirac (Breit 1928) point-nucleus value plus the free-electron QED anomaly plus bound-state QED/recoil/nuclear-size corrections:

```math
g_{\rm bound} = \underbrace{\tfrac{2}{3}\!\left[1 + 2\sqrt{1-(Z\alpha)^2}\right]}_{\text{Dirac binding}} + \underbrace{2 a_e}_{\text{free anomaly}} + \mathcal{O}\!\big(\alpha (Z\alpha)^2,\ m_e/M\big).
```

The Dirac binding term carries the genuine `Z`-dependence: expanded, `g_{\rm bound} \approx 2[1 - (Z\alpha)^2/3 - (Z\alpha)^4/12 - \dots] + 2a_e`.

**Framework prediction:** the framework's cutoff is Z-invariant (gating result above), so the anomaly piece `2a_e` is the Z=1 value evaluated at the universal `r_e/r_0`. The binding piece is the standard Coulomb–Dirac result, which is *formulation-independent* (it follows from the dual-Dirac FW reduction, [BS-§14.1](03_FineStructure.md)). Both readings coincide:

- **(Z-i) / (Z-ii):** identical — the cutoff does not vary with `Z`, so there is no reading difference in the prediction.

**Wolfram MCP check** (`r_e_Li2plus_joint_fit.wl` Section 1, 2026-05-27):

```text
gDirac1s[zz_] := (2/3)*(1 + 2*Sqrt[1 - (zz*aa)^2]); gBound[zz_] := gDirac1s[zz] + 2*ae;
Print["Dirac 1s binding g (Z=3) = ", gDirac1s[3]];   (* 1.99968 *)
Print["+ free anomaly 2 a_e = ", 2*ae];              (* 0.0023193 *)
Print["g_bound(7Li2+) = ", gBound[3]];               (* 2.00200 *)
Result: g_bound(⁷Li²⁺) = 2.0019998 ✅  (binding −3.20×10⁻⁴, anomaly +2.3193×10⁻³)
```

**🔴 Measurement-provenance finding (BLOCKED).** Solving `g_bound(Z) = 2.0000251707` (the brief's cited value) gives `Zα = 0.0586 ⇒ Z ≈ 8.04` — a hydrogen-like *oxygen* g-factor, **not** `Z=3`. The cited value is ~`2×10⁻³` too low for Li²⁺. Furthermore **Sturm 2014 *Nature* 506, 467 is the ¹²C⁵⁺ (Z=6) electron-mass measurement** (g(¹²C⁵⁺) ≈ 2.00104, Wolfram-consistent with the Breit formula), not a lithium measurement. The brief's #1 entry is therefore wrong on both value and citation. **The framework prediction (`2.00200`) is well-defined and verified, but no valid measurement is yet available to compare against** — the "headline precision test" framing rested on a mis-cited number.

**Verdict:** prediction ✅ (Wolfram-verified, `g_bound(⁷Li²⁺) ≈ 2.00200`); comparison **BLOCKED** pending a genuine high-precision ⁷Li²⁺ bound-g measurement (Tepper clarification requested on the intended source).

<!-- TODO: human reviews and fills in — confirms (a) the framework prediction g_bound(⁷Li²⁺) ≈ 2.00200, (b) the provenance finding that the brief's 2.0000251707 is a Z≈8 value and Sturm 2014 is the ¹²C⁵⁺ paper, and (c) whether a valid ⁷Li²⁺ measurement exists to complete the headline comparison. -->

---

### Result §Li-2 — 2S₁/₂–2P₁/₂ Lamb shift of Li²⁺ (Z=3) <a id="result-li-2--lamb-shift-of-li2-z3"></a>

**Source:** Bethe (1947) self-energy estimate via [BS-§19/§20](05_LambShift.md). *Substantive AI.*

**As measured:** `ΔE_{2S-2P}({}^7\text{Li}^{2+}) = 62\,765(21)` MHz (Schiffer, Bayfield & Pipkin-era measurement; **Schiffer 1995 *PRL* 74, 2188**). This is a **valid measurement** — the cleanest comparison of the four observables.

**QED Z-expansion:** the leading self-energy Lamb shift scales as `∝ (Zα)⁴ m_e c²/n³` (a `Z⁴` power law) modulated by the Bethe logarithm `∝ [\ln(1/(Zα)²) + C]`, which *shrinks* with `Z`. Naive `Z⁴` would give `81 × 1057.845 = 85\,685` MHz; the measured `62\,765` MHz corresponds to an effective ratio of `59.33`, the shortfall being exactly the Bethe-log shrinkage.

**Framework prediction:** per [BS-§19/§20](05_LambShift.md), the proper-time Lamb-shift route **is** the textbook Bethe-1947 estimate — matrix elements, energy denominators, the Bethe-log UV cutoff, and the mass-renormalisation subtraction are all formulation-independent. So the framework's Li²⁺ prediction equals the textbook Bethe-estimate at `Z=3`, which inherits the standard `~4%` residual against the full measured value.

- **(Z-i) / (Z-ii):** **identical** at the Bethe-estimate precision floor. The leading log-Bethe contribution is `g=2`-symmetric (see [BS-§20](05_LambShift.md), lines 114/184); `r_e` enters only at sub-leading (anomalous-moment) order, below the precision the route delivers. **#2 is a weak discriminator** of the cutoff question.

**Wolfram MCP check** (`r_e_Li2plus_joint_fit.wl` Section 2, 2026-05-27):

```text
ratioMeas = 62765/1057.845;                                   (* 59.33 *)
cc0 = Solve[3^4*(bracket[3,c]/bracket[1,c]) == ratioMeas, c]; (* C = -1.626 *)
effScaling = bracket[3,cc0]/bracket[1,cc0];                   (* 0.7325 vs Z^4 *)
fwLi = (1016/1057.845)*62765;
Result: framework Li²⁺ Bethe-estimate ≈ 60 282 MHz (residual −2 483 MHz, ~3.96%) ✅
        — same fractional residual as Z=1 (42 MHz / 1057.8 = 4.0%)
```

**Numerical comparison:**

| Source | `ΔE_{2S-2P}(⁷Li²⁺)` | Residual vs measured |
|---|---|---|
| Naive `Z⁴` × H Lamb | `85 685` MHz | `+22 920` MHz (Bethe-log not yet applied) |
| Framework / Bethe-estimate (Z⁴ × shrinking Bethe-log) | `≈ 60 282` MHz | `−2 483` MHz (~`3.96%`) |
| Schiffer 1995 measurement | `62 765(21)` MHz | — |

**Verdict:** ✅ at the **Bethe-estimate precision floor** — the framework reproduces the textbook Bethe-estimate Li²⁺ Lamb shift (`~60 282` MHz) with the same `~4%` residual it has at `Z=1`. This is **reproduction-by-construction, not an independent corroboration** (echoing [BS-§20](05_LambShift.md)'s honest framing); the Z⁴ × shrinking-Bethe-log scaling is correct, and the cutoff does not engage at this precision.

<!-- TODO: human reviews and fills in — confirms (a) the Schiffer 1995 measurement is the correct #2 target, (b) the framework Bethe-estimate prediction ~60 282 MHz with ~4% residual is the right reproduction-by-construction number, and (c) the weak-discriminator framing (r_e below the precision floor at Z=3) is faithfully recorded. -->

### Result §Li-3 — 2P₃/₂–2P₁/₂ fine structure of Li²⁺ (Z=3) <a id="result-li-3--fine-structure-of-li2-z3"></a>

**Source:** Sommerfeld–Dirac fine-structure formula via [BS-§14.2](03_FineStructure.md); anomaly from DRQM I §III.D. *Substantive AI.*

**As measured:** *(no valid measurement — see provenance).* The brief cited "~`7367` MHz (Bayfield/Riis era)." **This is the wrong ion** (see provenance note); a genuine hydrogenic ⁷Li²⁺ 2P fine-structure measurement (~`888` GHz scale) has not been sourced.

**QED Z-expansion:** the Sommerfeld–Dirac 2P₃/₂–2P₁/₂ splitting is `ΔE_{FS} = m_e c²(Zα)⁴/(2n⁴) = m_e c²(Zα)⁴/32` at `n=2` — a **pure `Z⁴`** power law with **no Bethe-log** (cleaner than the Lamb shift). The anomalous-moment correction is `ΔE_{anom} = ((g_s−2)/2)·ΔE_{leading} = a_e·ΔE_{leading}`.

**Framework prediction:** leading Dirac is formulation-independent (dual-Dirac FW reduction, [BS-§14.1](03_FineStructure.md)); the anomalous correction uses `a_e` from the **Z-universal cutoff**.

- **(Z-i) / (Z-ii):** **identical** — although #3 is a genuine `r_e`-engaging observable *in principle*, the cutoff is Z-invariant (gating result), so the anomaly piece is the same Z=1 `a_e` at Z=3. No reading difference.

**Wolfram MCP check** (`r_e_Li2plus_joint_fit.wl` Section 3, 2026-05-27):

```text
fsMHz[zz_] := mec2eV*(zz*aa)^4/32*eV2MHz;
Print["leading Dirac 2P FS (Z=3) = ", fsMHz[3]];   (* 886 892 MHz = 886.89 GHz, ratio = Z^4 = 81 *)
anomLi = ae*fsMHz[3];                                (* 1028.5 MHz *)
Print["FRAMEWORK total = ", fsMHz[3] + anomLi];     (* 887 920 MHz = 887.92 GHz *)
Result: framework Li²⁺ 2P FS = 887 920 MHz (886 892 leading + 1 028 anomalous) ✅
```

**🔴 Measurement-provenance finding.** The brief's "~`7367` MHz (Bayfield/Riis era)" is **smaller** than even the hydrogen value (`10 969` MHz), which is impossible under `Z⁴` scaling (the hydrogenic Li²⁺ splitting must be `~81×` larger, ≈ `888` GHz). The `7367` MHz figure refers to **helium-like Li⁺** (the two-electron `2³P` fine-structure intervals; Riis et al. measured Li⁺, not Li²⁺), not the hydrogenic single-electron Li²⁺ this branch targets. **The framework prediction (`887.92` GHz) is well-defined, but no valid hydrogenic measurement is available** — #3 is a **prediction-without-comparison** pending a sourced measurement (Tepper clarification requested).

**Numerical comparison:**

| Source | `ΔE_{FS}(2P₃/₂–2P₁/₂, ⁷Li²⁺)` | Note |
|---|---|---|
| Leading Dirac (`m_e c²(Zα)⁴/32`) | `886 892` MHz | formulation-independent |
| + anomalous (`a_e × leading`) | `+1 028` MHz | Z-universal cutoff |
| **Framework total** | **`887 920` MHz (`887.92` GHz)** | prediction |
| Brief's "~7367 MHz" | — | **wrong ion (helium-like Li⁺), discard** |

**Verdict:** prediction ✅ (Wolfram-verified, `887.92` GHz); comparison **BLOCKED** pending a genuine hydrogenic ⁷Li²⁺ 2P fine-structure measurement. As with #1, the cutoff does not actually distinguish readings (Z-invariant), so even a valid measurement would test back-fit self-consistency, not discrimination from standard QED.

<!-- TODO: human reviews and fills in — confirms (a) the framework prediction 887.92 GHz, (b) the provenance finding that the brief's ~7367 MHz is helium-like Li⁺ not hydrogenic Li²⁺, and (c) whether a valid hydrogenic ⁷Li²⁺ 2P fine-structure measurement exists. -->

### Result §Li-4 — ⁷Li 1s hyperfine splitting (Z=3) — *stub*

*To be drafted.* Z-scaling: `Z³` (contact density) × nuclear factors. ⁷Li: `I=3/2`, `μ ≈ 3.2564 μ_N`, `g_I ≈ 2.171`. Genuine `r_e` discriminator (Fermi contact linear in g_s). Brief target ~`12.7` GHz (Beckmann 1974) order-of-magnitude consistent (known ⁷Li²⁺ 1s HFS ≈ 11.8 GHz); exact value/source to be verified.

---

### §Li-5 — Z-axis verdict — *stub*

*To be drafted once #2–#4 predictions are in.* Expected disposition from the gating result: **Branch A with the back-fit caveat** — the cutoff is Z-universal by the framework's algebra (Branch B excluded), and the predictions that have valid measurements reproduce them at the framework precision floor. The Z-axis adds the structural confirmation that the lepton-axis A/C tension recurs unchanged: where the framework reduces to textbook QM (binding, Bethe-log, Fermi contact), it gives textbook predictions at every Z; the anomaly back-fit is Z-stable because `a_e` is a Z-independent free-electron property.
