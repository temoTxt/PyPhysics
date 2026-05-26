# Follow-up note: empirical triangulation of $r_e/r_0$ — for Tepper Gill

**Date:** 2026-05-26.
**From:** Trey Morris (with Claude Opus 4.7).
**Re:** Author guidance recorded after the 2026-05-25 phone call following [PR #59](https://github.com/temoTxt/PyPhysics/pull/59) (the interim report); empirical triangulation closes [issue #61](https://github.com/temoTxt/PyPhysics/issues/61); first-principles rederivation continues in [issue #54](https://github.com/temoTxt/PyPhysics/issues/54).

<!-- TODO: human reviews and fills in — confirms the note's framing, the choice to send only the Pass B (= branch c) result as the headline, and the inclusion of the residual table without expanding the Pass A contrast in full. -->

---

Tepper — following your guidance that the branches (b) and (c) in [PR #59 §5](https://github.com/temoTxt/PyPhysics/pull/59) are bracketing guides from a uni-observable search rather than theoretical predictions, we have performed the joint fit across all six precision observables that you suggested. Brief summary of what we found; full notebook and per-observable diagnostics in [`Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl`](../Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl).

## Triangulated value

$$\boxed{\;r_e/r_0 \;=\; 0.499\,420\,509\,912\,831\,7 \;\pm\; 2.5\times 10^{-13}\;}$$

obtained from a precision-weighted $\chi^2$ joint fit. The uncertainty is the one-sigma value from the second derivative of $\chi^2$ at the optimum.

This value differs from the published $r_e/r_0 = 0.499\,857\,150\,068\,631$ in the fourth decimal place, and is **indistinguishable from the back-fit-to-$g_s$ value (branch (c)) to 16 significant figures**. The triangulation therefore makes quantitative what the Bethe–Salpeter campaign already documented structurally ([`Bethe_Salpeter/10_CrossComparison.md §2`](../Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md#2-the-r_e-back-fit-self-consistency-across-six-g_s-dependent-observables)): every $g_s$-dependent observable is $(g_s/-2)^n \times \text{textbook}$, so one back-fit applied to six observables yields one $r_e$ value.

## Per-observable residuals at the triangulated optimum

| Observable | Prediction | Measurement | Residual | Source of residual |
|---|---|---|---|---|
| Electron $g_s$ | $-2.00231930\ldots$ | $-2.00231930\ldots$ | $0$ | matches by construction |
| H $2P_{3/2}-2P_{1/2}$ | $10{,}962$ MHz | $10{,}969.13$ MHz | $-7.13$ MHz | Bethe-estimate floor (BS-§14.2) |
| H 1S hyperfine | $1{,}420.04$ MHz | $1{,}420.4058$ MHz | $-0.366$ MHz | sub-leading QED (BS-§22.1) |
| He ${}^3P_0-{}^3P_1$ | $29{,}616.95$ MHz | $29{,}616.952$ MHz | $-0.002$ MHz | kHz floor (BS-§72) |
| Positronium ortho-para | $203{,}389$ MHz | $203{,}389$ MHz | $0$ | matches by construction |
| Muonium hyperfine | $4{,}463.40$ MHz | $4{,}463.3028$ MHz | $+0.097$ MHz | sub-leading QED (BS-§80) |

Every non-zero residual is a documented framework-precision-floor effect (Bethe-estimate / sub-leading QED gap to measurement) and is *not* reducible by re-tuning $r_e/r_0$. There is no observable in tension with the triangulated value once the framework's known precision scope is respected.

## Honest scoping note on the objective function

One substantive-AI decision exposed by the joint fit. The literal "fit to all six measurements" reading — treating each measurement's stated uncertainty as the $\sigma_i$ in $\chi^2$ — gives a slightly different optimum ($r_e/r_0 = 0.499\,406\,1\ldots$, differing in the fifth decimal place). That optimum pulls $g_e$ off the measured value by $\sim 5.8\times10^{-5}$ trying to compensate the hyperfine and helium fine-structure residuals, and is an artifact of treating framework-precision-floor residuals as if they were measurement uncertainty.

The honest formulation adds a framework-floor noise term to each $\sigma_i$ (so $\sigma_i^2 = \sigma_{\text{meas}, i}^2 + \sigma_{\text{floor}, i}^2$), recognising that the leading-$(g_s/-2)^n$ model cannot fit Bethe-estimate-gap residuals by any choice of $r_e$ alone. Under that weighting the optimum is the value above ($0.499\,420\,509\,912\,831\,7$), matching branch (c) to 16 sig figs. Both passes are documented in the notebook; we are sending you the honest answer here, with the alternative recorded for transparency.

## Where this leaves the open question

The empirical triangulation has supplied a single concrete value for $r_e/r_0$. The remaining open question — tracked in [issue #54](https://github.com/temoTxt/PyPhysics/issues/54) — is the first-principles derivation: starting from the dual-Dirac renormalisation prescription (DRQM I §III.D), what derivation produces $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ at the precision of the cutoff identification? If the answer is "the working notebook computed exactly this value and the published digits are a transcription error," that closes Finding 2 cleanly. If the derivation gives a different value, the framework's $g_s$ prediction needs to be revisited. Either way, the triangulated value is now the empirical target the rederivation needs to reproduce.

## One specific question, if useful

We did not find any observable in tension with the triangulated value, so we don't have a stretched-fit signal to flag. The one question that would help us most is whether you can retrieve (or rederive) the original $r_e$ computation in DRQM I §III.D and confirm whether the intended digits are the triangulated $0.499\,420\,509\,912\,831\,7$ or something else. A one-line answer is sufficient.

Thanks again for the phone call and for the bracketing-guide guidance — it pointed exactly at the right next step. The campaigns continue if this disposition is useful.

— Trey
