# Li²⁺ hyperfine — state log

## Iteration 0 — 2026-05-27 — initialized

- Branch `78-li2plus-hyperfine` checked out from main (post-#67 r_e thread closure).
- `.dev/research/brief.md` written (scope = observable #4 Li-7 1s hyperfine; $I=3/2$ nuclear-spin structure is the substantive AI extension of BS-§22's H $I=1/2$ apparatus).
- `.dev/research/loop_prompt.md` written.
- No prediction work yet.
- **Current observable focus:** #4 Li-7 1s hyperfine (the only observable owned by this branch).
- **Next:** read source-of-record §1 (`Bethe_Salpeter/06_Hyperfine.md` BS-§22 Fermi-contact apparatus) — record the Z=1 H 21-cm formula, the $(g_s/-2)^{n=1}$ scaling, and the $Z^3 g_s (m_e/m_p)$ leading dependence. Then look up the current best Li-7 1s hyperfine measurement (Beckmann 1974 / Schweikhard 1991 / modern refinement) and record DOI/year provenance.
- **Outcome-matrix:** not yet determinable.
- **Status:** READY.

## Iteration 1 — 2026-05-27 — read source-of-record §1 (BS-§22 Fermi-contact apparatus)

**What advanced:** Read `Bethe_Salpeter/06_Hyperfine.md` (BS-§22, the H 21-cm baseline). Recorded the Z=1 H Fermi-contact apparatus this branch must extend to Z=3 + I=3/2:

- **Fermi (1930) hyperfine Hamiltonian (hydrogenic 1S):**
  `H_HF = (8π/3) g_p μ_N g_s μ_B δ³(r) I·S` — leading term is *linear* in g_s (electron spin g-factor), confirming the n=1 r_e coupling.
- **1S splitting closed form:**
  `ΔE_HF(1S_{1/2}) = (4/3) g_p (m_e/M_p) α⁴ m_e c² (1 + QED corrections)`.
  Matrix element `⟨1S|δ³(r)|1S⟩ = 1/(π a_0³)`.
- **Z=1 H numerics (the baseline to scale):** leading Fermi at `g_s = -2` → **1418.4 MHz**; with textbook anomalous `(g_s/-2) = 1.00116` → 1420.04 MHz (triangulated r_e) vs measured **1420.405751768(2) MHz** (NIST 2020). Campaign verdict was ✅ at ~0.4 MHz precision floor, **explicitly "back-fit self-consistency, not independent corroboration"** (06_Hyperfine.md §22.1 back-fit caveat).
- **`(g_s/-2)^n` scaling:** n=1 for hyperfine (linear in g_s) — matches brief's claim; to be re-confirmed against `10_CrossComparison.md §2` next iter.
- **Triangulated cutoff:** `r_e/r_0 = 0.4994205099128317` (PR #62), giving `g_s = -2.00231930…`. This is the (Z-i) universal-cutoff value; (Z-ii) Z-scaled reading is still TBD.

**Angular-structure note (the substantive I=3/2 extension, queued):** H uses I=1/2 → F=0,1, splitting = A_H. Li-7 I=3/2 + S=1/2 → F=1,2; interval rule `ΔE(F→F-1) = A·F` gives top-interval splitting **2·A_{Li}** between F=2 and F=1. The `I·S` eigenvalues (F=2: +3/4, F=1: −5/4) are the substantive-AI move — to be derived/verified in the Wolfram notebook, not just asserted.

**Measurement provenance (still to confirm exact value next iter):** Beckmann 1974 *Z. Phys.* **270**, 173 (standard early Li²⁺ 1s HFS ref, ~12 732 MHz); refinements Schweikhard 1991 / Riis 1994. Li-7 nuclear: I=3/2, μ_I and g_I from CODATA/NUBASE — provenance to record when the value enters the prediction.

- **Current observable focus:** #4 Li-7 1s hyperfine.
- **Next:** read `10_CrossComparison.md §2` to confirm `n_HFS = 1` in the `(g_s/-2)^n` pattern and capture the six-observable back-fit framing (informs whether #4 lands in outcome-branch A vs C). Then look up + record the exact current-best Li²⁺ 1s HFS measurement (DOI/year) and Li-7 μ_I/g_I provenance.
- **Outcome-matrix:** not yet determinable (gathering apparatus; H baseline is back-fit self-consistency, which leans the universality question toward needing the Z=3 datapoint to discriminate A/B/C).
- **(Z-i)/(Z-ii):** universal-cutoff r_e/r_0 = 0.4994205099128317 recorded; Z-scaled (Z-ii) reading not yet constructed.
- **Status:** READY.

## Iteration 2 — 2026-05-27 — read source-of-record §2 (10_CrossComparison.md §2) — n=1 confirmed + out-of-sample insight

**What advanced:** Read `Bethe_Salpeter/10_CrossComparison.md` (campaign-closing cross-comparison). Two load-bearing facts captured:

1. **n=1 confirmed for hyperfine.** §2 states the scaling polynomial `f_i(g_s)` is "**linear for spin–orbit and Fermi contact contributions**; quadratic for two-fermion spin–spin contributions (positronium ortho-para)." Hyperfine (Fermi contact) → **n=1**, matching the brief and BS-§22.1. So Li²⁺ 1s HFS prediction = `(g_s/-2)^1 × (textbook Fermi-contact leading term at Z=3, I=3/2)`.

2. **The six-observable fit is back-fit self-consistency, NOT independent corroboration** (§1, §2 of CrossComparison, stated repeatedly). The triangulated `r_e/r_0 = 0.4994205099128317` is *defined* as the joint-best-fit value across six **Z=1-ish** g_s-dependent observables (H fine structure, H 21-cm, He ³P, He M1, positronium, muonium — all Z≤2 electron-bound-to-light-nucleus). `g_s = -2.00231930…` is substituted *by construction*.

**→ Substantive consequence for this branch (the discriminator):** Li²⁺ 1s HFS at **Z=3 is OUT-OF-SAMPLE** — it was *not* in the six-observable joint fit. Evaluating it at the *same universal* `r_e/r_0 = 0.4994205099128317` is therefore a genuine out-of-sample test, unlike the six in-sample observables which match by construction. This sharpens the outcome-matrix:
   - **(Z-i) universal cutoff:** g_s = -2.00231930 for all Z → Li²⁺ HFS = (g_s/-2)×textbook_Li. If this reproduces measurement at framework precision → **branch A**.
   - **(Z-ii) Z-scaled cutoff:** if the cutoff (hence g_s) runs with the Z=3 field strength, the multiplier differs → if a Z-scaled cutoff with framework-internal derivation is needed → **branch B**; if Li²⁺ needs its own per-Z back-fit → **branch C**.
   - Note: at n=1 the g_s correction is only +0.116%, so the dominant Z=3 physics is in the *textbook* Fermi-contact scaling (Z³, nuclear g_I, reduced mass). The r_e discriminating power lives in whether that +0.116% (and its possible Z-running) survives against the measurement precision.

- **Current observable focus:** #4 Li-7 1s hyperfine.
- **Next:** look up + record the exact current-best Li²⁺ 1s HFS measurement (value + DOI/year — Beckmann 1974 *Z. Phys.* 270, 173 baseline ~12 732 MHz; check for modern refinement) and the Li-7 nuclear moment μ_I / g_I with CODATA/NUBASE provenance. These are the inputs the Wolfram prediction cell needs. Then scaffold `r_e_Li2plus_hyperfine.wl`.
- **Outcome-matrix:** still not determinable (no numeric prediction yet), but now framed: Li²⁺ is the campaign's first **out-of-sample** g_s-dependent test → it can actually discriminate A vs B vs C, unlike the six in-sample observables.
- **(Z-i)/(Z-ii):** (Z-i) universal value fixed (g_s=-2.00231930). (Z-ii) Z-scaled reading: needs a model for how r_e runs with Z — flagged as a potential Tepper-input point (branch D risk) if the framework has no internal Z-running prescription.
- **Status:** READY.
