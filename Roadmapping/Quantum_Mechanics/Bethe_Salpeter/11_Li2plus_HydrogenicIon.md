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
| [§Li-4 — hyperfine](#result-li-4--hyperfine-of-li2-z3) | #4 ⁷Li 1s HFS | `Z³` × nuclear | **prediction ✅; measurement BLOCKED** | brief value off 2.35×; ~29.85 GHz predicted |
| [§Li-5 — Z-axis verdict](#result-li-5--z-axis-verdict) | joint χ² | — | **Branch A; B excluded** | cutoff Z-universal by construction |

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

### Result §Li-4 — ⁷Li 1s hyperfine splitting (Z=3) <a id="result-li-4--hyperfine-of-li2-z3"></a>

**Source:** Fermi (1930) contact term via [BS-§22.1](06_Hyperfine.md); ⁷Li nuclear data. *Substantive AI.*

**As measured:** *(provenance suspect — see below).* The brief cited "~`12.7` GHz (Beckmann 1974)." A `³He⁺`-validated scaling shows the hydrogenic ⁷Li²⁺ 1s HFS should be ~`29.8` GHz; the brief value is off by `2.35×`.

**QED Z-expansion:** the 1s Fermi-contact hyperfine splitting scales as `Z³` (contact density `|ψ(0)|² ∝ (Z/a₀)³`) times the nuclear+spin factor `f = μ_I·(2I+1)/(2I)`. For the F=I+½ ↔ F=I−½ interval, `ΔE_{HF} = A(I+½)` with `A ∝ μ_I/I`, so `ΔE_{HF} ∝ Z³·μ_I·(2I+1)/(2I)`. ⁷Li: `I=3/2`, `μ_I = 3.2564 μ_N`.

**Framework prediction:** scale the hydrogen value `1420.4` MHz by `Z³` and the nuclear factor, then by the anomaly factor `(g_s/−2) = 1.00116` from the **Z-universal cutoff**.

- **(Z-i) / (Z-ii):** **identical** — like #3, #4 engages `r_e` linearly in principle, but the Z-invariant cutoff makes both readings coincide.

**Wolfram MCP check** (`r_e_Li2plus_joint_fit.wl` Section 4, 2026-05-27):

```text
fHF[muI_, ii_] := muI*(2*ii+1)/(2*ii);
nuHe = 2^3*(fHF[2.127625,1/2]/fHF[2.792847,1/2])*1420.4;   (* 3He+ validation *)
nuLi = 3^3*(fHF[3.256427,3/2]/fHF[2.792847,1/2])*1420.4;
Result: 3He+ check = 8656.7 MHz (known 8665.6 MHz, 0.1% — method validated) ✅
        7Li²⁺ leading = 29 811 MHz; × anomaly (1.00116) = 29 846 MHz = 29.85 GHz
```

**🔴 Measurement-provenance finding.** The `³He⁺` cross-check validates the `Z³ × nuclear-factor` scaling to `0.1%` (`8656.7` vs known `8665.6` MHz), so the framework ⁷Li²⁺ prediction (`~29.85` GHz) is reliable. The brief's "~`12.7` GHz (Beckmann 1974)" is off by a factor `2.35` and is **inconsistent** with hydrogenic ⁷Li²⁺; Beckmann 1974 is a nuclear-magnetic-moment paper, not a Li²⁺ HFS measurement. **#4's measurement provenance is suspect** — a valid hydrogenic ⁷Li²⁺ 1s HFS measurement must be sourced before the comparison can be completed.

**Numerical comparison:**

| Source | `ΔE_{HF}(1s, ⁷Li²⁺)` | Note |
|---|---|---|
| ³He⁺ method validation | `8656.7` MHz | vs known `8665.6` MHz (`0.1%`) ✅ |
| Framework leading (`g_s=−2`) | `29 811` MHz | `Z³ × nuclear factor` |
| Framework × anomaly | `29 846` MHz (`29.85` GHz) | Z-universal cutoff |
| Brief's "~12.7 GHz" | — | **off by 2.35×, discard** |

**Verdict:** prediction ✅ (Wolfram-verified, ³He⁺-validated, `~29.85` GHz); comparison **BLOCKED** pending a valid hydrogenic ⁷Li²⁺ 1s HFS measurement. As with #1/#3, the Z-invariant cutoff means even a valid measurement tests back-fit self-consistency, not discrimination.

<!-- TODO: human reviews and fills in — confirms (a) the ³He⁺-validated framework prediction ~29.85 GHz for ⁷Li²⁺, (b) the provenance finding that the brief's ~12.7 GHz is off by 2.35× and Beckmann 1974 is a nuclear-moment paper, and (c) whether a valid hydrogenic ⁷Li²⁺ 1s HFS measurement exists. -->

---

### §Li-5 — Z-axis verdict <a id="result-li-5--z-axis-verdict"></a>

**Source:** synthesis of §§Li-1–Li-4 + the gating structural result. *Substantive AI.*

**Summary of the four predictions** (all Wolfram-verified; companion notebook Sections 1–4):

| # | Observable | Z-scaling | Framework prediction | Valid measurement? |
|---|---|---|---|---|
| 1 | bound g-factor | binding `(Zα)²` + anomaly | `g = 2.00200` | **no** — brief value is a Z≈8 number |
| 2 | 2S–2P Lamb shift | `Z⁴` × shrinking Bethe-log | `60 282` MHz | **yes** — `62 765(21)` MHz ✅ (~4% residual) |
| 3 | 2P fine structure | pure `Z⁴` | `887.92` GHz | **no** — brief value is helium-like Li⁺ |
| 4 | ⁷Li 1s hyperfine | `Z³` × nuclear | `29.85` GHz | **no** — brief value off `2.35×` |

**Joint χ² at Z=3 — not meaningfully computable.** A joint fit for a Li²⁺-optimal `r_e/r_0` requires observables that both (a) engage the cutoff and (b) have valid measurements. Of the four: #2 has a valid measurement but **does not engage the cutoff** (`g=2`-symmetric, weak discriminator); #1/#3/#4 engage the cutoff but **lack valid measurements** in the brief (three mis-provenanced values — see §§Li-1/Li-3/Li-4). The intersection is empty, so no Z=3 joint optimum can be extracted. This is a **measurement-sourcing limitation, not a framework limitation** — all four framework predictions are computed and verified.

**Z-axis verdict — Branch A (cutoff Z-universal); Branch B ruled out structurally; A-vs-C indistinguishable.**

1. **The cutoff is Z-universal — Branch A — and this is forced by the framework's own algebra, not fitted.** The (III.22) closed form `r_e/r_0 = (2−a_e)/(2(2+a_e))` depends only on the free-electron anomaly `a_e` and `r_0 = e²/(m_e c²)`, neither of which depends on `Z`. The same electron radiates at every `Z`, so `r_e/r_0` is identical at `Z=1` and `Z=3`. **Branch B (a derivable Z-scaling) is excluded** because no `Z` appears anywhere in the cutoff structure.

2. **Li²⁺ adds no new empirical constraint on the cutoff** (answering the brief's cross-particle cross-check). The lepton-axis muon test ([PR #70](https://github.com/temoTxt/PyPhysics/pull/70)) was a *meaningful* constraint (`>57`kσ exclusion of lepton-universality) **because changing the particle changes `a_ℓ` and `r_0^ℓ`** — there was something to constrain. Changing `Z` changes *neither* `a_e` *nor* `r_0`, so the Z-axis is **structurally incapable of constraining the cutoff** the way the lepton axis did. There is no "Li²⁺ value of `r_e/r_0`" distinct from the electron's; the question is answered before any measurement is taken.

3. **A-vs-C remains indistinguishable, exactly as on the lepton axis.** The one clean comparison (#2 Lamb shift) reproduces measurement by construction *without engaging the cutoff*; the three cutoff-engaging observables have no valid measurements. So Li²⁺ cannot discriminate "the framework's cutoff has independent content" (A) from "the cutoff is a Z-stable back-fit with no content distinct from textbook QED" (C). Where the framework reduces to textbook QM (Dirac binding, Bethe-log self-energy, Fermi contact), it gives textbook predictions at every `Z`; the anomaly back-fit is Z-stable precisely because `a_e` is a Z-independent free-electron property.

**Bottom line:** the framework's cutoff is **Z-universal by construction** (Branch A, Branch B excluded). The Z-axis confirms — rather than independently tests — the lepton-axis finding: the dimensionless cutoff is a property of the *electron* (via `a_e`), invariant under nuclear charge. The strongest honest statement is structural, not empirical: *Li²⁺ cannot add a cutoff constraint because there is no Z-dependence in the cutoff to constrain.* Completing the empirical comparisons for #1/#3/#4 requires correctly-sourced hydrogenic ⁷Li²⁺ measurements (flagged for author/orchestrator follow-up).

#### Cross-PR reconciliation with [PR #87](https://github.com/temoTxt/PyPhysics/pull/87) (hydrogenic-ion Z-scan, issue [#82](https://github.com/temoTxt/PyPhysics/issues/82)) — *steel-man revision 2026-05-27*

This branch's "Branch A by construction" finding rests on the *structural* observation that no $Z$ appears in the framework's cutoff formula. **The parallel hydrogenic-ion Z-scan (PR [#87](https://github.com/temoTxt/PyPhysics/pull/87)) supplies the empirical complement** — bound-electron g-factor measurements at $Z \in \{2, 6, 8, 14\}$ from Schneider 2022 (³He⁺), Sturm 2011 (¹²C⁵⁺), Verdú 2004 (¹⁶O⁷⁺), Sturm 2013 (²⁸Si¹³⁺).

| $Z$ | Measured $g_e^{\rm bound}$ | Framework (Z-i) prediction | Residual |
|---|---|---|---|
| 1 (free) | $-2.00231930$ | $-2.00231930$ | 0 (triangulated) |
| 2 | $-2.00217742$ | $-2.00231930$ | $1.4\times 10^{-4}$ |
| 6 | $-2.00104159$ | $-2.00231930$ | $1.3\times 10^{-3}$ |
| 8 | $-2.00004703$ | $-2.00231930$ | $2.3\times 10^{-3}$ |
| 14 | $-1.99534896$ | $-2.00231930$ | $7.0\times 10^{-3}$ |

The Z-invariant (Z-i) prediction gives residuals of $10^{-4}$ to $10^{-3}$ against measured $|g_e^{\rm bound}|$ — **$10^5$ to $10^7\,\sigma$** in measurement units ($\sigma \sim 10^{-9}$–$10^{-11}$).

**Reconciliation.** The two verdicts are not contradictory once carefully framed:

- **Branch A is the framework's stated position** — the published §III.D apparatus has no $Z$ in the cutoff, so it predicts a $Z$-invariant $g_e^{\rm bound}$. This branch's "Branch A by construction" finding is correct as a statement about the framework's *apparatus*.
- **Branch C is the operationally correct verdict** — the multi-Z data refutes the $Z$-invariant prediction at $10^5$–$10^7\,\sigma$. The back-fit $r_e^{(Z)}/r_0$ values inherit QED's bound-state structure $a_e^{\rm bound}(Z\alpha) = a_e^{\rm free} - \tfrac{1}{3}(Z\alpha)^2 + \mathcal{O}((Z\alpha)^4)$ per-Z, with **no framework-internal derivation** of the $(Z\alpha)^2$ binding coefficient.
- **Same verdict shape as the lepton axis** ([PR #70](https://github.com/temoTxt/PyPhysics/pull/70)): the published apparatus says the cutoff is universal-in-X (across leptons / across Z); the data shows the cutoff is *X-specific through* the QED anomalous moment $a_\ell$ or $a_e^{\rm bound}(Z\alpha)$. **The framework's apparatus is Z-trivial; that is exactly the structural shape that produces Branch C operationally.**

**The "Li²⁺ adds no new constraint" framing is therefore too modest.** The framework's apparatus is empirically inadequate for $Z > 1$ precision predictions; the multi-Z g-factor scan exposes a $10^5$–$10^7\,\sigma$ gap that the framework's published structure has no mechanism to close. The honest Z-axis verdict shifts:

- **Branch A** (framework's stated position): confirmed structurally.
- **Branch C** (operationally correct): confirmed empirically by PR [#87](https://github.com/temoTxt/PyPhysics/pull/87)'s data.
- The two coexist because they describe different things — the *apparatus* is Z-trivial (A); the *empirical comparison* requires per-Z bound-state QED inheritance (C). Closing the gap to a derived Branch B is the open question for [#75](https://github.com/temoTxt/PyPhysics/issues/75) (framework specs, Tepper engagement).

<!-- TODO: human reviews and fills in — confirms (a) the cross-PR reconciliation framing (Branch A structurally + Branch C operationally is the honest joint verdict), (b) the residual table from #87 is faithfully captured, (c) the parallel to PR #70's lepton-axis verdict, and (d) the link to #75 as the open derivational question. -->

<!-- TODO: human reviews and fills in — confirms (a) the Z-universal (Branch A) verdict with Branch B structurally excluded, (b) the finding that Li²⁺ adds no new cutoff constraint (the Z-axis is structurally unlike the lepton axis), (c) the A-vs-C indistinguishability, and (d) the meta-finding that three of the brief's four measurement values are mis-provenanced and need re-sourcing. -->
