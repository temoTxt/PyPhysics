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

## Iteration 3 — 2026-05-27 — measurement + nuclear-moment provenance lookup (WebSearch)

**What advanced:** WebSearch for the Li²⁺ 1s HFS measurement + Li-7 nuclear moments. Two provenance findings, one of them materially affecting the verdict structure:

1. **⚠ No direct experimental ⁷Li²⁺ (hydrogenic, Z=3) 1s HFS measurement exists.** Per Pachucki, Patkóš, Yerokhin et al., *Hyperfine splitting in ⁶,⁷Li⁺* (arXiv:2309.00436, 2023): the Li²⁺ HFS "prediction is derived from the sum of theoretical [QED] contributions and the **experimental Li⁺** value; the uncertainty of the Li²⁺ prediction comes exclusively from the uncertainty of the experimental Li⁺ HFS." So the authoritative **comparator is a QED-theory value** (their Table VII), not a direct measurement. → This *changes the test*: the framework's verdict is "match to the QED-theory comparator" (weaker than a direct-measurement match, but still a real out-of-sample check at Z=3).
   - **Authoritative comparator source:** Pachucki et al., arXiv:2309.00436, Table VII (Li²⁺ 1s HFS). To extract the exact value next iter.
   - Historical/early ref: Beckmann, Boeklen, Elke 1974, *Z. Phys.* **270**, 173 (the brief's seed) — but this is the **Li / Li⁺ atomic-beam** HFS work, *not* the hydrogenic Li²⁺ 1s value.

2. **⚠ Brief's "~12 732 MHz" seed value is suspect.** A rough Fermi-contact scaling from H's 1420.4 MHz — `Δν_Li²⁺ ≈ Δν_H × Z³ × (g_I^{Li}/g_I^{p}) × (angular)`:
   - `Z³ = 27`; `g_I^{Li-7} = μ_I/I = 3.256427/1.5 = 2.1709 μ_N`; `g_I^{p} = 2.79285/0.5 = 5.5857 μ_N` → ratio 0.3887.
   - H splitting = A_H·(I+½) = A_H·1; Li-7 splitting (F=2↔1) = A_Li·(I+½) = A_Li·2 (interval rule).
   - → `A_Li/A_H = 0.3887 × 27 = 10.49`; `Δν_Li²⁺ ≈ 10.49 × 1420.4 × 2 ≈ 29 800 MHz ≈ 30 GHz` (pre reduced-mass/relativistic corrections).
   - This is ~2.3× the brief's 12 732 MHz. The 12 732 figure is most likely a *different* Li system (Li⁺ 1s2s³S₁ or neutral). **Resolve next iter** via Pachucki Table VII + the Wolfram computation.

3. **Li-7 nuclear data (provenance):** I=3/2; μ_I = +3.256427 μ_N (NUBASE/CODATA); g_I (=μ_I/I) = 2.1709. Confirm exact CODATA digits when value enters the notebook.

- **Current observable focus:** #4 Li-7 1s hyperfine.
- **Next:** scaffold `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_hyperfine.wl` with (a) the hydrogenic Fermi-contact closed form, (b) the H 1420.4 MHz reproduction as a sanity cell, (c) the Z=3 + I=3/2 scaling cell — which simultaneously *computes the textbook comparator* and resolves the 12 732-vs-30 000 MHz discrepancy. Header to carry the Crocco substantive-AI tag for the I=1/2→3/2 extension.
- **Outcome-matrix:** still not determinable, but the test is now reframed: **comparator is QED-theory (Pachucki 2023), not direct measurement** — record as a verdict caveat. Out-of-sample A/B/C discrimination still holds (Z=3 not in the six-observable fit).
- **(Z-i)/(Z-ii):** unchanged from iter 2.
- **Status:** READY.

## Iteration 4 — 2026-05-27 — scaffold `r_e_Li2plus_hyperfine.wl` (Fermi-contact apparatus, Z=3 + I=3/2)

**What advanced:** Created `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_hyperfine.wl` (new, distinct from other branches' notebooks). Four single-line MCP-ready cells (Wolfram gotchas obeyed — no bare `V`/`e`, descriptive symbols, `;`-joined):

- **Cell 1 — H 1s sanity:** `dnu_H = (4/3) gP (me/Mp) alpha^4 (me c²/h)`, gP=5.5856946; expect ~1418–1421 MHz (reproduces the 21-cm baseline).
- **Cell 2 — Li/H scaling ratio (the substantive core):** textbook splitting `dnu ~ mu_I(2I+1)/(2I) Z³`, so `dnu_Li/dnu_H = [muLi(2ILi+1)/(2ILi)]/[muP(2IP+1)/(2IP)] × Z³`. Hand-computed: nucFactorLi = 3.256427×(4/3) = 4.3419; nucFactorP = 2.79285×2 = 5.5857; ratio 0.7773 × 27 = **20.99**; × 1420.4058 MHz → **textbook Li²⁺ 1s HFS ≈ 29 810 MHz ≈ 29.8 GHz** (g_s=-2 baseline).
- **Cell 3 — framework (g_s/-2)¹ correction:** gsTri=-2.00231930436 at r_e/r_0=0.4994205099128317 → factor 1.001160 → **framework prediction ≈ 29 845 MHz**.
- **Cell 4 — I=3/2 angular structure (substantive-AI tag):** S=1/2⊗I=3/2 → F=1,2; I·S eigenvalues F=2:+3/4, F=1:−5/4; headline F=2↔1 interval = 2a = (I+½)a. This is the explicit I=1/2→3/2 extension of BS-§22's machinery.

**Discrepancy resolved (analytically, pending MCP confirm):** the brief's seed **12 732 MHz is NOT the hydrogenic Li²⁺ value** — the Fermi-contact scaling unambiguously gives ~29.8 GHz. 12 732 MHz is a different Li system (Li⁺ 1s2s³S₁ or neutral-Li manifold). Recorded.

- **Current observable focus:** #4 Li-7 1s hyperfine.
- **Next:** execute Cells 1–4 via Wolfram MCP to get exact numerics (confirm H sanity ~1420 MHz, textbook ~29.8 GHz, framework ~29.85 GHz, angular factor = 2a). Then fetch Pachucki 2023 Table VII for the exact QED-theory comparator to compute the residual.
- **Outcome-matrix:** still not finalized (need MCP exact value + Pachucki comparator), but provisional reading: since the (g_s/-2) correction is only +0.116% and the cutoff is held universal (Z-i), if the framework's 29.85 GHz matches Pachucki's QED comparator within the textbook-leading-g_s floor → **branch A** (universal cutoff survives the out-of-sample Z=3 test). Confirm next iters.
- **(Z-i)/(Z-ii):** (Z-i) prediction now numerically pinned (~29 845 MHz). (Z-ii) Z-scaled reading still needs a framework Z-running prescription — branch-D risk persists if absent.
- **Status:** READY.
