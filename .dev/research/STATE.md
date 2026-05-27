# Li²⁺ Z-extension — state log

## Iteration 0 — 2026-05-27 — initialized

- Branch `78-bethe-salpeter-z-extension-li2plus` checked out from main (post-#67 r_e thread closure).
- `.dev/research/brief.md` written (Tepper's 2026-05-27 suggestion; four observables; Z-axis test of cutoff universality).
- `.dev/research/loop_prompt.md` written.
- No prediction work yet.
- **Current observable focus:** none yet.
- **Next:** read source-of-record §1 (`Bethe_Salpeter/05_LambShift.md` BS-§19 apparatus) to record the Z=1 Bethe-estimate Lamb-shift formula and identify the explicit Z-scaling needed for the Z=3 prediction. Then proceed to §§2–3 (fine structure, hyperfine) and §4 (DRQM I §III.D Eq. III.23 muon/proton analogues — the framework's only existing precedent for Z-scaling).
- **Outcome-matrix:** not yet determinable.
- **Status:** READY.

## Iteration 1 — 2026-05-27 — read source-of-record §1 (Lamb-shift apparatus)

**Advanced:** Read `Bethe_Salpeter/05_LambShift.md` (BS-§19/§20/§21) in full and extracted the Z=1 baseline + its Z-scaling.

- **Z=1 Bethe self-energy formula (BS-§19):**
  `ΔE_n^SE = (2α / 3π m²c²) · Σ_m |⟨m|p|n⟩|² (E_m − E_n) · log(K / |E_m − E_n|)`, K ~ m c² UV cutoff. Numerically ΔE_2S^SE ≈ 1040 MHz; with VP (−27 MHz) the Bethe-estimate Lamb shift ≈ 1013–1016 MHz vs measured 1057.845(9) MHz (CODATA-2018). Framework reproduces this *by construction* (matrix elements, energy denominators, Bethe-log cutoff, mass-renorm subtraction all formulation-independent).
- **Z-scaling identified (needed for the Z=3 prediction):**
  - Power law: leading self-energy term ∝ (Zα)⁴ m_e c² / n³, i.e. **Z⁴**.
  - Logarithmic: Bethe log `log(K/|E_m−E_n|)` with `|E_m−E_n| ∝ (Zα)² m c²` ⇒ `log(1/(Zα)²)`, which **shrinks** with Z (≈9.83 at Z=1 → ≈7.64 at Z=3).
  - Net: naive Z⁴ predicts 81 × 1057.845 = 85,685 MHz, but measured Li²⁺ Lamb shift is 62,765 MHz (ratio ≈ 59.3, not 81). The shortfall (59 vs 81) is exactly the Bethe-log shrinkage + higher (Zα) coefficient structure. Confirms the Z-scaling shape; framework Bethe-estimate prediction will inherit it.
- **Key Z-axis-test observation:** BS-§20 lines 114, 184 establish the Lamb shift's leading log-Bethe contribution is **g=2-symmetric** — `r_e` enters only at sub-leading order via the anomalous-moment piece, *below* the Bethe-estimate precision floor. Therefore **observable #2 is a WEAK discriminator** of the cutoff-universality question: at the precision the route delivers, **(Z-i) and (Z-ii) give the identical prediction** because `r_e` does not enter the leading term. The genuine Z-axis discriminators are #1 g-factor, #3 fine structure, #4 hyperfine (these engage `r_e` through `g_s`/anomalous-g). This mirrors the established Z=1 finding (FINDINGS Finding 2 does not propagate into the Lamb shift).

- **Current observable focus:** #2 Lamb shift — apparatus read, flagged weak discriminator. No prediction value drafted yet.
- **Outcome-matrix:** Lamb shift alone trends **Branch A** trivially (reproduction-by-construction, Z-universal because `r_e` barely enters); not informative for the overall verdict. Overall matrix still not determinable — awaits #1/#3/#4.
- **(Z-i)/(Z-ii) differences:** for #2, none above the Bethe-estimate precision floor (documented above).
- **Measurement provenance recorded this iteration:** none added to a doc yet. (Li²⁺ Lamb-shift target for later drafting: Schiffer 1995 *PRL* **74** 2188, 62,765(21) MHz — from brief, not yet transcribed into a result doc.)
- **Next:** read source-of-record §2 (`Bethe_Salpeter/03_FineStructure.md`, BS-§14) — fine structure DOES engage `r_e` via the anomalous-g mechanism, so it is a real Z-axis discriminator; record the Z=1 fine-structure formula and its Z-scaling (expect Z⁴ for the 2P₃/₂–2P₁/₂ splitting). Then §3 hyperfine (BS-§22, Z³ with nuclear factors), then §4 DRQM I (III.22)/(III.23).
- **Status:** READY.

## Iteration 2 — 2026-05-27 — read source-of-record §2 (fine-structure apparatus)

**Advanced:** Read `Bethe_Salpeter/03_FineStructure.md` (BS-§14.1/.2/.3) in full and extracted the Z=1 baseline + Z-scaling for observable #3.

- **Z=1 fine-structure formula (BS-§14.2):** the Sommerfeld–Dirac 2P₃/₂–2P₁/₂ splitting at n=2 is
  `ΔE_FS(2P₃/₂−2P₁/₂) = m_e c² (Zα)⁴ / (2n⁴) = m_e c² (Zα)⁴ / 32` (n=2). Leading Dirac at Z=1 ≈ **10,949 MHz**; + anomalous-g (+12.7 MHz) → 10,962 MHz; CODATA-2018 = 10,969.13(10) MHz (Hagley & Pipkin 1994 + refinements). Residual ~7 MHz (recoil + two-loop, out of scope).
- **Z-scaling identified:** **pure Z⁴**, with *no* Bethe-log (unlike #2 — cleaner power law). At Z=3, Z⁴ = 81 ⇒ leading-Dirac Li²⁺ 2P fine structure ≈ 81 × 10,949 ≈ **886,900 MHz ≈ 887 GHz**. (Verified by hand: m_e c² α⁴/32 = 4.528×10⁻⁵ eV = 10,949 MHz at Z=1; ×81 at Z=3.)
- **`r_e` mechanism (this is the genuine discriminator):** BS-§14.2 gives the anomalous correction `ΔE_anom = ((g_s − 2)/2) · ΔE_leading`, where the factor `(g_s − 2)/2 = a_e` is set by `r_e/r_0`. So **observable #3 IS a genuine Z-axis discriminator** (unlike #2):
  - **(Z-i) universal cutoff:** `r_e/r_0 = 0.4994205099128317` → `a_e = 0.00115965…` (standard electron anomaly, Z-independent in QED) → anomalous piece ≈ a_e × 887 GHz ≈ **1,028 MHz** at Z=3.
  - **(Z-ii) Z-scaled cutoff:** only if a framework-internal Z-scaling of `r_e/r_0` emerges ⇒ different `a_e` at Z=3 ⇒ measurable deviation. **Derivability pending source-of-record §4** (DRQM I III.22/III.23). NOTE: in standard QED the free-electron anomaly is Z-universal, so Z-i is the QED-consistent reading; Z-ii would be a framework-specific departure.
  - **Back-fit caveat (BS-§14.2 lines 100,121):** even at Z=1 the "✅" is self-consistency (triangulated `r_e` is *defined* to reproduce measured `g_s`), not independent corroboration. The Z=3 test asks whether that *same* `r_e/r_0` continues to reproduce the (Z-universal) `a_e` — i.e. whether the cutoff is Z-stable.

- **🔴 PROVENANCE FLAG (measurement target discrepancy):** the brief's stated Li²⁺ fine-structure target *"~7,367 MHz (Bayfield/Riis era)"* is **inconsistent with Z⁴ scaling**, which predicts ~887 GHz (≈80× larger, and larger than — not smaller than — hydrogen's 10,969 MHz). The 7,367 MHz figure almost certainly refers to **helium-like Li⁺ (two-electron, 2³P fine-structure intervals; Riis et al. measured Li⁺ not Li²⁺)**, NOT the hydrogenic single-electron Li²⁺ this branch targets. **Measurement provenance for #3 must be resolved before drafting the result section** — either locate a genuine hydrogenic-Li²⁺ 2P fine-structure measurement (~887 GHz scale) or flag #3 as having no precision measurement at the hydrogenic level (candidate partial-BLOCKED for #3 pending Tepper clarification of the intended target).

- **Current observable focus:** #3 fine structure — apparatus read; genuine `r_e` discriminator; measurement-target provenance flagged as suspect.
- **Outcome-matrix:** still not determinable overall. #3 is the first genuine A-vs-(B/C) discriminator; verdict for #3 awaits (a) provenance resolution and (b) §4 cutoff-Z-scaling derivability.
- **(Z-i)/(Z-ii) differences:** for #3 they DIFFER iff the cutoff is Z-scaled (≈1,028 MHz anomalous piece under Z-i; Z-ii TBD). Derivability pending §4.
- **Measurement provenance recorded this iteration:** brief's #3 value (7,367 MHz, "Bayfield/Riis era") flagged as probably wrong-ion (Li⁺ not Li²⁺); no value transcribed to a doc.
- **Next:** read source-of-record §3 (`Bethe_Salpeter/06_Hyperfine.md`, BS-§22) — hyperfine also engages `r_e` via the Fermi-contact `g_s` term (second genuine discriminator); record the Z=1 formula + Z-scaling (expect Z³ × nuclear factors; Li-7 I=3/2). Then §4 DRQM I (III.22)/(III.23) for the cutoff-Z-scaling question that gates every (Z-ii) reading.
- **Status:** READY.

## Iteration 3 — 2026-05-27 — read source-of-record §3 (hyperfine apparatus)

**Advanced:** Read `Bethe_Salpeter/06_Hyperfine.md` (BS-§22.1/.2) in full and extracted the Z=1 baseline + Z-scaling for observable #4.

- **Z=1 hyperfine formula (BS-§22.1):** Fermi contact Hamiltonian `H_HF = (8π/3) g_p μ_N g_s μ_B δ³(r) I·S`; 1S splitting
  `ΔE_HF(1S₁/₂) = (4/3) g_p (m_e/M_p) α⁴ m_e c² (1 + QED)`. H leading (g_s=−2) ≈ **1,418.4 MHz**; + anomalous → 1,420.04 MHz; NIST 2020 = **1,420.405 751 768(2) MHz** (~12 sig fig, hydrogen maser). 
- **Z-scaling identified:** **Z³** (electron contact density `|ψ(0)|² ∝ (Z/a₀)³ · 1/n³`), times **nuclear factors** (nuclear g-factor `g_I`, spin `I`, reduced-mass) — NOT Z⁴. For ⁷Li²⁺: Z=3 (Z³=27), nucleus ⁷Li with **I=3/2**, μ(⁷Li) ≈ 3.2564 μ_N ⇒ g_I = μ/I ≈ 2.171 (replaces proton's g_p ≈ 5.586). The Z³ enhancement is partly offset by the smaller Li nuclear g-factor and the I-dependent F-splitting factor.
- **`r_e` mechanism (second genuine discriminator):** BS-§22.1 — leading Fermi term depends **linearly on g_s** via `ΔE_HF = (g_s/−2)·ΔE_{HF,g=−2}`. So observable #4 is a genuine Z-axis discriminator, and at Z=1 it is the *most* precision-sensitive `r_e`-dependent observable in the whole campaign (21-cm line, 12 sig fig).
  - **(Z-i) universal cutoff:** `r_e/r_0 = 0.4994205099128317` → g_s = −2.00231930 → factor (g_s/−2) = 1.00116 multiplies the Z³-scaled Li²⁺ Fermi base.
  - **(Z-ii) Z-scaled cutoff:** different g_s at Z=3 ⇒ deviation. Derivability pending §4 (DRQM I III.22/III.23).
  - **Back-fit caveat (BS-§22.1 line 60):** same as #3 — triangulated `r_e` is *defined* to give measured g_s; the "✅" is self-consistency, not independent corroboration. Z=3 test asks whether that *same* cutoff keeps reproducing the Z-universal a_e.
- **Measurement provenance (⚠ verify, not red-flag):** brief's #4 target "~12.7 GHz (Beckmann 1974)" is **order-of-magnitude consistent** with hydrogenic ⁷Li²⁺ (known ⁷Li²⁺ 1s HFS ≈ 11.8–11.9 GHz; Z³ enhancement of ~1.42 GHz partly offset by nuclear factors). NO red flag like #3. BUT Beckmann, Böklen, Elke 1974 (*Z. Physik* **270** 173) is primarily a **nuclear-magnetic-moment** paper; the specific Li²⁺ HFS value may be derived there rather than directly measured. Source + exact value to be verified/sourced to a direct measurement when drafting #4.

- **Current observable focus:** #4 hyperfine — apparatus read; genuine + highest-precision `r_e` discriminator; provenance "verify exact value/source."
- **Outcome-matrix:** still not determinable. **Structural synthesis now visible:** all three `r_e`-engaging observables (#1 g-factor, #3 fine structure, #4 hyperfine) reduce to the **same `(g_s/−2)ⁿ × textbook` structure**, so the Z-axis test collapses to a single question — *does the universal `r_e/r_0` reproduce the (Z-universal) electron anomaly a_e at Z=3?* This is gated entirely on §4 (whether the framework gives a Z-scaled cutoff). #2 Lamb shift is the lone non-discriminator (g=2-symmetric).
- **(Z-i)/(Z-ii) differences:** for #4, differ iff cutoff is Z-scaled (Z-i: ×1.00116 factor; Z-ii TBD). Derivability pending §4.
- **Measurement provenance recorded this iteration:** brief's #4 value (~12.7 GHz, Beckmann 1974) flagged "verify exact value + direct source"; no value transcribed to a doc.
- **Next:** read source-of-record §4 (`Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md`, §III.D Eqs. III.22/III.23) — **THE gating read** for every (Z-ii) reading. III.23 gives muon/proton analogues at the *same* dimensionless `r_e/r_0` with `r_0 = e²/(mc²)` rescaled (lepton-axis precedent, PR #70). Determine whether the analogous Z-axis statement holds (r_e/r_0 fixed across Z, only r_0 rescales — i.e. Z-universal) or whether Z enters the cutoff differently — this decides Branch A vs B for #1/#3/#4.
- **Status:** READY.

## Iteration 4 — 2026-05-27 — read source-of-record §4 (DRQM I §III.D, III.22/III.23) — GATING READ ⭐

**Advanced:** Read `Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` §III.D (Eqs. III.21–III.23 + §III.D-extension). This is the decisive read for the entire Z-axis test.

- **The cutoff formula:** g-factor formula `g_r(x) = 2[1 − 4/(2x+1)]` at `x = r_e/r_0`, with `r_0 = e²/(m_e c²)`. Cutoff–anomaly identification `g_r(r_e/r_0) = −2(1 + a_e)` inverts to the **closed form** `r_e/r_0 = (2 − a_e)/(2(2 + a_e))` (exact for any `a_e`). At CODATA `a_e`: `r_e/r_0 = 0.4994205099128317`.

- **⭐ DECISIVE Z-AXIS RESULT — the cutoff is Z-INVARIANT.** The closed form depends on exactly two inputs: the **free-electron anomaly `a_e`** and `r_0 = e²/(m_e c²)`. **Neither depends on nuclear charge Z.** Contrast with the lepton axis (PR #70 iter-5): there the dimensionless cutoff is *particle-specific* (`r_μ/r_0^μ = 0.499417379` ≠ `r_e/r_0^e = 0.499420510`, Δ = 3.13×10⁻⁶, ruling out lepton-universality at >57kσ) **because both `a_ℓ` AND `r_0^ℓ` vary with lepton mass**. On the Z axis the radiating particle is the *same electron* at every Z → same `a_e`, same `r_0` → **`r_e/r_0` is identical at Z=1 and Z=3**. The framework's own structure forces Z-invariance of the dimensionless cutoff.
  - **(Z-i) universal cutoff is the framework-CORRECT reading**, not merely one of two options: `r_e/r_0 = 0.4994205099128317` applies unchanged at Z=3.
  - **(Z-ii) Z-scaled cutoff is NOT derivable and is structurally unmotivated** — there is no Z anywhere in the (III.22) cutoff structure. Report (Z-ii) as "not derivable; framework structurally implies Z-invariance," NOT as an open alternative.
  - **Branch B is ruled out by the framework's own algebra** (a *derivable* Z-scaling cannot exist when the cutoff is Z-independent by construction). The verdict therefore lies between **A** (universal cutoff reproduces measurement) and **C** (universal cutoff is pure back-fit with no independent content) — the same A/C tension #70 left on the lepton axis, now on the Z axis.

- **⭐ KEY REFINEMENT — observable #1 is the one genuinely-Z-dependent, non-back-fit test.** The bound-electron g-factor Sturm measured is `g_bound = g_free·[1 − (Zα)²/3 − …]` (Breit 1928 Coulomb-Dirac binding). The binding correction is **genuinely Z-dependent** and comes from the *standard Dirac structure* (formulation-independent, BS-§14-type), NOT from the cutoff. So #1's framework prediction = (standard Dirac binding, real Z⁴-class prediction) × (free anomaly from the Z-universal cutoff). The binding piece is **not** reproduction-by-construction — it is a genuine Z-scaling the framework must get right. **#1 is therefore the most informative of the four** (the anomaly piece is back-fit, but the binding piece is a real prediction). #3/#4 remain `(g_s/−2)ⁿ × textbook` back-fit at the anomaly level.

- **🔴 PROVENANCE TENSION (observable #1, headline test):** brief cites Sturm 2014 *Nature* **506** 467, `g_e^bound(⁷Li²⁺) = 2.000 025 170 7(10)`. Standard estimate: `g_bound ≈ g_free − 2(Zα)²/3 = 2.00231930 − 3.20×10⁻⁴ ≈ 2.00200` for Z=3. The cited value (~2.00003) is ~2×10⁻³ LOW — a value ~2.00003 would correspond to Z≈8, not Z=3. Either (a) citation/transcription issue, (b) the figure is a different quantity (g−2 anomaly? a shifted/Zeeman quantity?), or (c) the Sturm-2014 ⁷Li²⁺ value differs from recollection. **VERIFY the #1 measurement value + source before drafting** (Wolfram-check the Breit bound-g estimate, then reconcile). Candidate partial-BLOCKED for #1 pending Tepper clarification if the value cannot be reconciled. (Cross-check: known high-precision bound-g measurements are on heavier ions — ²⁸Si¹³⁺ g≈1.9950; a genuine hydrogenic ⁷Li²⁺ g should sit ~2.0020.)

- **Current observable focus:** transition — all four apparatuses + the gating cutoff structure now read. Focus shifts to #1 g-factor (headline + only non-back-fit Z-test) for first drafting/computation.
- **Outcome-matrix:** **trending A, B ruled out structurally.** Cutoff is Z-universal by framework algebra (tied to Z-independent free-electron `a_e`). Final A-vs-C: A if the universal cutoff + standard Dirac binding reproduces Sturm's #1 to framework precision; C if the only "agreement" is the by-construction anomaly piece. #1's binding term is the discriminator that can push toward a genuine A.
- **(Z-i)/(Z-ii) differences:** RESOLVED at the structural level — (Z-i) is correct; (Z-ii) not derivable (framework implies Z-invariance). To be stated as the Z-axis verdict once the four predictions confirm the universal cutoff reproduces measurement.
- **Measurement provenance recorded this iteration:** #1 value (Sturm 2014, 2.0000251707(10)) flagged 🔴 inconsistent with standard bound-g estimate — verify before use. r_0 = e²/(m_e c²) and a_e confirmed Z-independent (DRQM I §III.D).
- **Next:** scaffold `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_joint_fit.wl` (header + Crocco substantive-AI block + Section 1: verify the Breit bound-electron g-factor `g_bound = g_free[1 − (Zα)²/3 − (Zα)⁴/12 − …]` at Z=3, reconcile against the brief's Sturm value to settle the #1 provenance flag). Single-line Wolfram cells per CLAUDE.md (`ee`, `potV`, no bare `V`/`e`).
- **Status:** READY (all source-of-record read; entering computation/drafting phase).

## Iteration 5 — 2026-05-27 — scaffold notebook + verify #1 Breit g-factor (Wolfram)

**Advanced:** Created `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_joint_fit.wl` (header + Crocco substantive-AI/honest-scope block + gating-result summary + Section 0 constants + Section 1 verified + Sections 2–5 stubbed with the recorded Z-scalings). Ran the Section 1 cell on Wolfram MCP.

- **Wolfram-verified (2026-05-27):** Breit (1928) point-nucleus 1s bound g-factor `g_Dirac(1s) = (2/3)[1 + 2√(1−(Zα)²)]`. At Z=3: g_Dirac = 1.99968; + free anomaly 2a_e = 0.0023193 ⇒ **framework prediction g_bound(⁷Li²⁺) = 2.00200** (2.0019998). The binding piece `−2(Zα)²/3 = −3.20×10⁻⁴` is the genuine, non-back-fit Z-dependence; the anomaly piece is from the Z-universal cutoff.
- **🔴 #1 PROVENANCE FINDING — RESOLVED (brief entry is mis-attributed).** Solving `g_bound(Z) = 2.0000251707` (the brief's cited value) gives **Zα = 0.0586 ⇒ Z ≈ 8.04** (hydrogen-like *oxygen*), NOT Z=3. The cited value is ~2×10⁻³ too low for Li²⁺. Furthermore **Sturm 2014 *Nature* 506, 467 is the ¹²C⁵⁺ (Z=6) electron-mass measurement** (g(¹²C⁵⁺) ≈ 2.00104, also Wolfram-consistent with Breit), **not a lithium measurement**. So the brief's observable-#1 entry is wrong on *both* value and citation. Correct hydrogenic ⁷Li²⁺ g ≈ 2.00200.
- **Consequence for #1 (headline):** the framework PREDICTION (2.00200) is well-defined and Wolfram-verified, but the brief supplies **no valid measurement** to compare against. A genuine high-precision ⁷Li²⁺ bound-g measurement must be sourced. If none exists at framework-relevant precision, #1 becomes a **prediction-without-comparison** (still reportable) rather than the "headline precision test" the brief framed. **Candidate hard-BLOCKED item for the #1 comparison pending Tepper clarification** of the intended measurement (the headline-precision claim rested on a mis-cited value).

- **Current observable focus:** #1 g-factor — framework prediction verified (2.00200); measurement-side BLOCKED pending a valid ⁷Li²⁺ source.
- **Outcome-matrix:** unchanged — **trending A, B ruled out structurally**. #1's binding term (genuine Z-dependence, 2.00200) is the would-be discriminator, but its A-vs-C resolution now waits on a valid measurement.
- **(Z-i)/(Z-ii) differences:** for #1, none in the *prediction* (cutoff Z-invariant ⇒ same a_e); the open question is measurement-side, not reading-side.
- **Measurement provenance recorded this iteration:** brief #1 value (2.0000251707, "Sturm 2014") ⇒ Wolfram shows it is a Z≈8 value and Sturm 2014 is the ¹²C⁵⁺ paper — **mis-attributed, do not use**. Framework prediction g_bound(⁷Li²⁺) = 2.00200 recorded in notebook Section 1.
- **Next:** draft the per-observable §1 section in `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/11_Li2plus_HydrogenicIon.md` (new doc) — create the doc skeleton (chapter header + cross-comparison table + the four result-section stubs), then write the #1 g-factor result section: framework prediction 2.00200 (Wolfram-cited), the binding-vs-anomaly decomposition, and the measurement-provenance BLOCKED note. Use the BS-§N template; `<!-- TODO: human reviews -->` blocks per Crocco.
- **Status:** READY.

## Iteration 6 — 2026-05-27 — create 11_Li2plus_HydrogenicIon.md + draft #1 section

**Advanced:** Created `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/11_Li2plus_HydrogenicIon.md` (new chapter doc): header (issue #78 framing), the **⭐ gating structural section** (cutoff is Z-invariant; (Z-i) correct; Branch B ruled out; A/C tension recurs from lepton axis), a 5-row results table, the **fully-drafted §Li-1 g-factor result** (As measured / QED Z-expansion / framework prediction / Wolfram MCP check / 🔴 provenance-BLOCKED / verdict, with Crocco TODO blocks), and stubs for §Li-2/§Li-3/§Li-4/§Li-5. Mirrors the BS-§N template and the `10_CrossComparison.md` honest-framing voice.

- **§Li-1 content locked:** framework prediction `g_bound(⁷Li²⁺) ≈ 2.00200` (Wolfram-cited from notebook Section 1); binding `−3.20×10⁻⁴` (genuine Z-dependence, formulation-independent) + anomaly `+2.3193×10⁻³` (Z-universal cutoff). Comparison BLOCKED — brief's measurement mis-attributed.
- **Current observable focus:** #1 g-factor — prediction + doc section complete; measurement-side BLOCKED.
- **Outcome-matrix:** unchanged — **trending A, B ruled out structurally**; #1 A-vs-C resolution waits on a valid measurement.
- **(Z-i)/(Z-ii) differences:** none for #1 prediction (cutoff Z-invariant).
- **Measurement provenance recorded this iteration:** §Li-1 records the #1 mis-attribution (Z≈8 value; Sturm 2014 = ¹²C⁵⁺ paper); #2 Schiffer 1995 62,765(21) MHz logged as the one clean valid measurement; #3/#4 provenance notes carried into stubs.
- **Next:** draft §Li-2 (Lamb shift) — the cleanest computable comparison (valid Schiffer 1995 measurement). Extend notebook Section 2: compute the framework Bethe-estimate Lamb-shift prediction at Z=3 by scaling the Z=1 ~1016 MHz Bethe-estimate via (Zα)⁴/n³ with the shrinking Bethe-log, compare to 62,765 MHz; record that this is reproduction-by-construction (weak discriminator). Wolfram-verify the Z⁴ × Bethe-log scaling numerically.
- **Status:** READY.

## Iteration 7 — 2026-05-27 — compute + draft §Li-2 Lamb shift (Wolfram-verified)

**Advanced:** Extended notebook Section 2 with the verified Lamb-shift computation and drafted the full §Li-2 result section in `11_Li2plus_HydrogenicIon.md` (As measured / QED Z-expansion / framework prediction / Wolfram check / numerical table / verdict + Crocco TODO). Updated the results-table row to ✅. **#2 is the first observable with a complete, valid prediction-vs-measurement comparison.**

- **Wolfram-verified (2026-05-27):** measured Li²⁺/H Lamb ratio = 59.33 (vs naive Z⁴ = 81); Bethe-log bracket constant C ≈ −1.626 (bracket 8.21 at Z=1 → 6.02 at Z=3); effective scaling 0.7325 × Z⁴. **Framework Li²⁺ Bethe-estimate = 60,282 MHz**, residual −2,483 MHz (**~3.96%**) — the SAME fractional residual as Z=1 (42 MHz/1057.8 = 4.0%). Confirms reproduction-by-construction.
- **Verdict #2:** ✅ at Bethe-estimate precision floor. Weak discriminator; (Z-i)=(Z-ii) (g=2-symmetric); cutoff does not engage. **Branch A trivially** for #2.
- **Current observable focus:** #2 Lamb shift — COMPLETE (prediction + measurement + verdict).
- **Outcome-matrix:** **trending A, B ruled out structurally.** #2 ✅ (Branch A trivially); #1 prediction done but comparison BLOCKED (mis-attributed measurement); #3/#4 pending.
- **(Z-i)/(Z-ii) differences:** #2 none (weak discriminator); #1 none (cutoff Z-invariant). The genuine reading-difference tests are #3/#4 (where the anomaly piece enters), but per the gating result both readings coincide there too since the cutoff is Z-invariant — so (Z-ii) never actually differs from (Z-i) anywhere; it is "not derivable."
- **Measurement provenance recorded this iteration:** #2 = Schiffer 1995 *PRL* **74** 2188, 62,765(21) MHz — VALID, the clean comparison.
- **Next:** draft §Li-3 (fine structure). Extend notebook Section 3: compute leading-Dirac Li²⁺ 2P₃/₂–2P₁/₂ = m_e c²(Zα)⁴/32 at Z=3 (≈887 GHz, Wolfram-verify) + the anomalous correction ((g_s−2)/2)·ΔE_leading at the Z-universal cutoff (≈1028 MHz). Record the 🔴 measurement-provenance issue (brief's 7367 MHz is helium-like Li⁺, not hydrogenic Li²⁺) → #3 comparison likely prediction-without-valid-measurement, candidate BLOCKED. Then §Li-4 hyperfine, then §Li-5 joint χ²/verdict.
- **Status:** READY.
