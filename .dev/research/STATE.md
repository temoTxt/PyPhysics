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
