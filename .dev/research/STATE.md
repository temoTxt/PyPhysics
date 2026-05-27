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
