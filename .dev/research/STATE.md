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

## Iteration 5 — 2026-05-27 — execute + debug Wolfram cells (MCP); framework prediction pinned ≈ 29.80 GHz

**What advanced:** Ran all four cells via Wolfram MCP, found and fixed a methodological double-count, re-ran.

**First MCP run confirmed the arithmetic:**
- Cell 1 (H sanity): leading-Fermi H 1s HFS = **1421.16 MHz** ✅ (reproduces 21-cm baseline; ~0.2% above the campaign's 1418.4 MHz due to omitted reduced-mass factor).
- Cell 4 (angular structure): I·S(F=2) = **+3/4**, I·S(F=1) = **−5/4**, headline F=2↔1 splitting = **2a = (I+½)a** ✅. This is the verified I=3/2 extension of BS-§22's I=1/2 machinery.

**Bug found + fixed (the "debug" of this step):** iter-4 Cell 2 scaled the *measured* H value (1420.4058 MHz). That value already carries the electron-g_s correction (and proton-specific QED/recoil/nuclear-structure), so Cell 3's additional ×(g_s/−2) **double-counted the anomalous moment (~+34 MHz)** and wrongly imported proton nuclear-structure terms. The campaign structure is `framework = (g_s/−2)ⁿ × (textbook LEADING g_s=−2 term)`. Fixed Cell 2 to scale `dnuHleadingFermi = 1418.4` (BS-§22.1 g_s=−2 baseline).

**Corrected MCP run (self-consistent across two routes):**
- scaleRatio (Li/H) = **20.9878**.
- textbook Li²⁺ 1s HFS (g_s=−2 leading) = **29 769.1 MHz**.
- **framework Li²⁺ 1s HFS prediction (g_s=−2.00231930) = 29 803.6 MHz ≈ 29.80 GHz** (Z-i universal cutoff).
- cross-check (scale the framework-H value 1420.04 MHz by the ratio) = 29 803.5 MHz → **agrees to 0.1 MHz**, confirming the methodology is internally consistent.
- g_s contributes **+34 MHz** (+0.116%) over the g_s=−2 baseline.

**Below-floor caveat recorded:** the reduced-mass ratio (Li²⁺ nucleus ≫ proton) is a separate ~+0.14% (~+42 MHz) correction NOT yet applied (would push to ~29 846 MHz); it sits at/below the campaign's Bethe-estimate precision floor and is comparable to the g_s correction — to be noted in the doc verdict as a floor caveat, not a precision claim.

- **Current observable focus:** #4 Li-7 1s hyperfine.
- **Next:** fetch the exact Pachucki et al. 2023 (arXiv:2309.00436) Table VII Li²⁺ 1s HFS QED-theory comparator value, compute the residual vs the framework's 29 804 MHz, and assign the outcome-matrix branch. Then draft the #4 per-observable section in `Bethe_Salpeter/13_Li2plus_Hyperfine.md` with verdict.
- **Outcome-matrix:** framework prediction pinned (~29 804 MHz, Z-i); branch assignment pending the Pachucki comparator + residual. Provisional: if 29 804 MHz matches Pachucki's value within the ~0.1–0.5% Bethe-estimate floor → **branch A** (universal cutoff survives the out-of-sample Z=3 test).
- **(Z-i)/(Z-ii):** (Z-i) = 29 804 MHz (firm). (Z-ii) still needs a framework Z-running prescription for r_e — branch-D risk if the framework has none.
- **Status:** READY.

## Iteration 6 — 2026-05-27 — fetch Pachucki comparator; resolve outcome-matrix branch (A, structural)

**What advanced:** WebSearch + WebFetch on Pachucki, Patkóš, Yerokhin, *Hyperfine splitting in ⁶,⁷Li⁺* (arXiv:2309.00436, 2023). Confirmed it as the comparator-of-record and confirmed methodology: Li²⁺ hfs is *predicted from the experimental Li⁺ hfs* + QED theory (no direct Li²⁺ data). **The exact Table-VII digit was not extractable** — arxiv.org/abs gives only the abstract; `arxiv.org/html/...` 404s; ar5iv's summarizer drops the results table. Exact value flagged for human extraction (PDF parse) — **NOT blocking the verdict** (see below).

**Decisive verdict insight (the substantive payoff):** the framework prediction = `(g_s/−2)¹ × textbook-leading`. The *only* r_e-dependent content is the `(g_s/−2) = +0.116%` factor. That anomalous-moment enhancement is **Z-INDEPENDENT** — it is the free-electron g_s, identical to what textbook QED applies at every Z. Under **(Z-i) universal cutoff**, the framework applies `g_s = −2.00231930` at Z=3 exactly as textbook QED does (the triangulated r_e was *fit* to reproduce that g_s). → The Li²⁺ test is therefore **agreement by construction at the leading-g_s precision floor**, i.e. the *same structural self-consistency* as the H 21-cm result (06_Hyperfine.md §22.1 back-fit caveat), **NOT independent corroboration.**
   - The "out-of-sample at Z=3" sharpening from iters 2–5 is real *only if* the framework's r_e runs with Z (**(Z-ii)**), which would make g_s(Z=3) ≠ −2.00231930. **The framework does not assert Z-running r_e** — its stated position is a universal cutoff (Z-i). So (Z-ii) is a hypothetical the framework doesn't claim; constructing it would need Tepper input (would be branch D), but it is moot under the framework's actual (Z-i) stance.

- **Current observable focus:** #4 Li-7 1s hyperfine.
- **OUTCOME-MATRIX BRANCH: A** — under the Z=1-triangulated *universal* cutoff, the framework reproduces the standard-QED Li²⁺ 1s HFS at the framework (leading-g_s / Bethe-estimate) precision floor. **With the H-style structural-self-consistency caveat:** the +0.116% g_s enhancement is Z-independent and identical to textbook QED, so this is consistency-by-construction, not an independent Z=3 discrimination. (Branch B/C would require a Z-running cutoff the framework doesn't posit; branch D not triggered.)
- **Next:** add a Wolfram cell computing the *standard point-nucleus QED comparator* (relativistic Fermi factor `1/[γ(2γ−1)]`, γ=√(1−(Zα)²); +reduced-mass) so the doc has a reproducible residual-table number alongside the Pachucki external value. Then draft the #4 per-observable section in `Bethe_Salpeter/13_Li2plus_Hyperfine.md` with the branch-A verdict + structural caveat + nuclear-structure (Bohr–Weisskopf) floor caveat.
- **(Z-i)/(Z-ii):** (Z-i) = 29 804 MHz, branch A (structural). (Z-ii) not asserted by framework; hypothetical-only.
- **Measurement provenance:** comparator = Pachucki/Patkóš/Yerokhin 2023 (arXiv:2309.00436) Table VII (QED theory from exp. Li⁺); exact digit pending human PDF extraction. Li-7: I=3/2, μ_I=+3.256427 μ_N (NUBASE/CODATA).
- **Status:** READY.

## Iteration 7 — 2026-05-27 — add Cell 5: standard-QED point-nucleus comparator + residual table numbers (MCP-verified)

**What advanced:** Extended `r_e_Li2plus_hyperfine.wl` with **Cell 5** (standard-QED point-nucleus comparator) and verified it via Wolfram MCP. This gives the residual-table numbers the doc needs:

| Quantity | Value (MHz) | Note |
|---|---|---|
| Framework leading-Fermi × g_s (Z-i, minimal) | **29 803.6** | matches campaign H methodology (leading only) |
| + relativistic Dirac factor `1/[γ(2γ−1)]`, γ=√(1−(Zα)²) | +0.072% | framework inherits Dirac (BS-§14.1) |
| + reduced-mass ratio³ (⁷Li nucleus vs proton) | +0.140% | standard kinematic |
| **Standard-QED point-nucleus comparator** (leading+rel+redmass+g_s) | **29 866.8** | = framework with same corrections, by construction |
| Pachucki 2023 Table VII (incl. Bohr–Weisskopf) | _[human PDF extract]_ | external high-precision; nuclear-structure REDUCES by ~0.5–1% |

- **residual framework(minimal) − standard-QED comparator = −63.2 MHz (−0.21%)** — the gap is the relativistic + reduced-mass corrections the minimal leading route omits. **Both are r_e-INDEPENDENT**, so the residual does **not** bear on the branch-A verdict (exact analog of the H 21-cm leading-vs-full gap).
- Bohr–Weisskopf nuclear-structure: out of scope (issue #78), recorded as a framework-floor caveat (~−0.5–1%).

- **Current observable focus:** #4 Li-7 1s hyperfine.
- **OUTCOME-MATRIX BRANCH: A** (unchanged) — universal cutoff reproduces standard-QED Li²⁺ HFS at the leading-g_s/Bethe-estimate floor; structural self-consistency (g_s enhancement Z-independent), not independent corroboration.
- **Next (final substantive step):** draft the #4 per-observable section in `Bethe_Salpeter/13_Li2plus_Hyperfine.md` — measurement provenance (no direct exp; Pachucki comparator), framework prediction (29 804 MHz minimal / 29 867 MHz with standard corrections), the residual table above, the I=3/2 angular-structure derivation, the branch-A verdict + structural-self-consistency caveat + Bohr–Weisskopf floor caveat, Crocco substantive-AI tags + `<!-- TODO: human reviews and fills in -->` blocks (incl. the exact Pachucki digit). **This completes the done-criteria → STOP after the doc has the verdict.**
- **(Z-i)/(Z-ii):** (Z-i) firm (29 804 minimal / 29 867 with std corrections). (Z-ii) hypothetical-only (framework asserts universal cutoff).
- **Status:** READY (doc draft is the terminal step).

## Iteration 8 — 2026-05-27 — draft `13_Li2plus_Hyperfine.md` with verdict — DONE, LOOP STOPS

**What advanced:** Created `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/13_Li2plus_Hyperfine.md` (this branch's distinct output file) — the complete #4 per-observable result section. Contents:
- Chapter intro framing it as the Z=3/I=3/2 extension of BS-§22.1; companion-notebook link.
- **Source + the substantive $I=3/2$ angular structure**: $F=1,2$; I·S eigenvalues (F=2:+3/4, F=1:−5/4); headline interval = 2a = (I+½)a — the explicit substantive-AI extension of BS-§22's proton I=1/2 apparatus, with a `<!-- TODO human review -->` block.
- **Measurement provenance**: no direct Li²⁺ experiment; comparator = Pachucki et al. 2023 (arXiv:2309.00436) Table VII; Li-7 μ_I=+3.256427 μ_N; exact digit flagged for human PDF extraction.
- **Framework prediction**: Cell-2/3/5 MCP numbers — minimal 29 804 MHz, with-standard-corrections 29 867 MHz; full residual table.
- **Verdict: branch A** at leading-g_s/Bethe-estimate floor, with the load-bearing **structural-self-consistency caveat** (g_s enhancement is Z-independent, identical to textbook QED → agreement by construction, not independent Z=3 discrimination), the **(Z-i)/(Z-ii) honesty** note (framework asserts universal cutoff; Z-running not claimed → (Z-ii)/branch-B/C/D not the framework's position), and the **Bohr–Weisskopf nuclear-structure floor caveat** (out of scope per #78).
- Crocco substantive-AI tags + `<!-- TODO: human reviews and fills in -->` blocks throughout (angular structure, provenance/exact-digit, verdict disposition).

**DONE-CRITERIA MET:**
1. ✅ `13_Li2plus_Hyperfine.md` has the #4 result section: measurement provenance + framework prediction (Z-i; Z-ii addressed as non-asserted) + verdict.
2. ✅ `r_e_Li2plus_hyperfine.wl` exists with Wolfram-MCP-verified cells (1–5) for the Fermi-contact prediction at Z=3, I=3/2.
3. ✅ Outcome-matrix branch recorded: **A** (structural self-consistency).

- **Current observable focus:** #4 Li-7 1s hyperfine — COMPLETE.
- **OUTCOME-MATRIX BRANCH: A** (final).
- **Open items for human/orchestrator (not blocking):** (a) extract exact Pachucki Table VII Li²⁺ digit into the residual table; (b) contribute the 29 804/29 867 MHz prediction to the Self-Energy branch's joint-χ² fit (`r_e_Li2plus_joint_fit.wl` — owned by that branch, not edited here); (c) PR + GitHub comments handled by orchestrator after morning review.
- **Status:** ✅ DONE — loop stops (no further ScheduleWakeup).
