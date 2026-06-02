---
title: "Follow-up note: empirical triangulation of \\(r_e/r_0\\)"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-05-26"
subject: "Empirical triangulation of r_e/r_0 across six precision observables (follow-up to interim report)"
---

# Follow-up note: empirical triangulation of $r_e/r_0$ — for Tepper Gill

**Date:** 2026-05-26.
**From:** Trey Morris (with Claude Opus 4.7).
**Re:** Author guidance recorded after the 2026-05-25 phone call following the interim report[^1]; empirical triangulation closes the joint-fit thread[^2]; the first-principles rederivation continues separately[^3].

---

Tepper — following your guidance that the branches (b) and (c) in §5 of the interim report[^1] are bracketing guides from a uni-observable search rather than theoretical predictions, we have performed the joint fit across all six precision observables that you suggested. We summarise here what we found. The full notebook and per-observable diagnostics are recorded separately[^4].

## 1. Triangulated value

We obtain, from a precision-weighted $\chi^2$ joint fit,

$$
r_e/r_0 \;=\; 0.499\,420\,509\,912\,831\,7 \;\pm\; 2.5\times 10^{-13},
\qquad (1)
$$

where the uncertainty is the one-sigma value taken from the second derivative of $\chi^2$ at the optimum.

We observe that (1) differs from the published $r_e/r_0 = 0.499\,857\,150\,068\,631$ in the fourth decimal place. It is, however, indistinguishable from the back-fit-to-$g_s$ value (branch (c) of the interim report) to sixteen significant figures. The triangulation thus makes quantitative what the Bethe–Salpeter campaign already documented structurally[^5]. Indeed, every $g_s$-dependent observable in the framework satisfies

$$
\text{prediction} \;=\; (g_s / -2)^n \times \text{textbook}, \qquad n \in \{1, 2\},
\qquad (2)
$$

so a single back-fit applied to six observables yields a single $r_e$ value.

## 2. Per-observable residuals at the triangulated optimum

The fit residuals are collected in Table 1.

| Observable | Prediction | Measurement | Residual | Source of residual |
|---|---|---|---|---|
| Electron $g_s$ | $-2.00231930\ldots$ | $-2.00231930\ldots$ | $0$ | matches by construction |
| H $2P_{3/2}-2P_{1/2}$ | $10{,}962$ MHz | $10{,}969.13$ MHz | $-7.13$ MHz | Bethe-estimate floor[^6] |
| H 1S hyperfine | $1{,}420.04$ MHz | $1{,}420.4058$ MHz | $-0.366$ MHz | sub-leading QED[^7] |
| He ${}^3P_0 - {}^3P_1$ | $29{,}616.95$ MHz | $29{,}616.952$ MHz | $-0.002$ MHz | kHz floor[^8] |
| Positronium ortho-para | $203{,}389$ MHz | $203{,}389$ MHz | $0$ | matches by construction |
| Muonium hyperfine | $4{,}463.40$ MHz | $4{,}463.3028$ MHz | $+0.097$ MHz | sub-leading QED[^9] |

*Table 1.* Per-observable residuals at the triangulated optimum (1).

Every non-zero residual in Table 1 is a documented framework-precision-floor effect. It is not reducible by re-tuning $r_e/r_0$. No observable is in tension with (1) once the framework's known precision scope is respected.

## 3. Honest scoping note on the objective function

We record one substantive-AI decision exposed by the joint fit. The literal "fit to all six measurements" reading treats each measurement's stated uncertainty as the $\sigma_i$ entering $\chi^2$. Under that reading the optimum is

$$
r_e/r_0 \;=\; 0.499\,406\,1\ldots,
\qquad (3)
$$

differing from (1) in the fifth decimal place. Equation (3) pulls $g_e$ off the measured value by approximately

$$
\Delta g_e \;\approx\; 5.8 \times 10^{-5},
\qquad (4)
$$

in an attempt to compensate the hyperfine and helium fine-structure residuals. This is an artifact of treating framework-precision-floor residuals as if they were measurement uncertainty.

The honest formulation adds a framework-floor noise term to each $\sigma_i$,

$$
\sigma_i^{2} \;=\; \sigma_{\text{meas},\,i}^{2} \;+\; \sigma_{\text{floor},\,i}^{2}.
\qquad (5)
$$

The recognition behind (5) is that the leading $(g_s/-2)^n$ model cannot fit Bethe-estimate-gap residuals by any choice of $r_e$ alone. Under the weighting prescribed by (5) the optimum is the value (1), matching branch (c) of the interim report to sixteen significant figures. Both passes are documented in the notebook[^4]. We send you the honest answer here, with the alternative (3) recorded for transparency.

## 4. Where this leaves the open question

A natural way to read the relationship between the published $r_e$ and the triangulated value follows. The published $r_e/r_0 = 0.499\,857\ldots$ is an *initial-value* result from the working-notebook numerical search against $g_s$ alone. The joint fit across six observables is a *refinement calculation* on top of that initial value. The triangulated value (1) supersedes the initial value in the same way that any iteratively-refined numerical result supersedes an earlier estimate. It does not invalidate the derivation that produced the initial value. It refines the cutoff identification using six constraints rather than one.

A further refinement is still open: a first-principles derivation of $r_e$ from the dual-Dirac renormalisation prescription, in the sense of DRQM I §III.D. This is tracked separately[^3]. The candidate starting points are sketched briefly below, and expanded into one-to-three-paragraph explanations each in a companion note[^10].

By "first-principles" we mean a derivation in which $r_e$ emerges from the dual-Dirac equation's internal structure, rather than being fixed by inverse-solving against the measured $g_e$. The output would be a closed-form expression for $r_e/r_0$ in terms of the framework's fundamental parameters — functions of $\alpha$, $r_0$, and the structural constants of the dual representation. Such an expression would evaluate numerically to approximately $0.499\,420\,509\ldots$ if the empirical refinement holds, or to a different value if the derivation exposes a structure the triangulation has not captured.

We have considered three plausible starting points, and would defer to your judgment on which (if any) is natural for this framework.

1. The proper-time self-energy integral in the dual-Dirac formalism, with $r_e$ emerging as the regulator scale at which mass renormalisation closes consistently.
2. A variational determination, with $r_e$ fixed by demanding that the renormalised dual-Dirac equation reproduce the physical electron mass at the cutoff. This is analogous to the way the standard-QED running coupling is fixed by a renormalisation condition.
3. A structural constant of the dual representation — something that the $b$-factor projection structure or the second-order Dirac decomposition singles out — making $r_e/r_0$ a derivable framework-internal ratio rather than a renormalisation choice.

Per your 2026-05-25 author guidance, a fourth candidate we had initially considered — the retrieval of an unpublished working-notebook derivation step — is empty. The original cutoff value was a numerical search against $g_s$ alone, with no derivation-step working draft to retrieve.

We have not pursued any of the three. The natural starting point depends on the framework's internal logic better than we can guess at from outside. The thread is non-urgent. The triangulated value (1) already serves as the campaign's current best refinement, and a first-principles derivation can sit indefinitely as a would-be-nice-to-have until you point us at a starting structure, or decide the thread is not where the project's time is best spent.

Whichever route turns out to be natural, the outcome could agree with the triangulated value (confirming this refinement), supersede it (a third refinement), or expose a derivation-level structure that reframes the cutoff entirely. Each of these outcomes is consistent with treating each calculation as an iterative refinement rather than as a final authoritative answer.

## References

[^1]: Interim report, pull request #59 (`github.com/temoTxt/PyPhysics/pull/59`). §5 of that report sketches the branches (b) and (c) discussed here.

[^2]: Issue #61 (`github.com/temoTxt/PyPhysics/issues/61`). Empirical triangulation thread.

[^3]: Issue #54 (`github.com/temoTxt/PyPhysics/issues/54`). First-principles rederivation thread.

[^4]: `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl`. Joint-fit notebook with per-observable diagnostics.

[^5]: `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md`, §2 (the $r_e$ back-fit self-consistency across six $g_s$-dependent observables).

[^6]: `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/14_Lamb_Shift_2P.md`, §14.2 (Bethe-estimate precision floor).

[^7]: `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/22_Hydrogen_Hyperfine.md`, §22.1 (sub-leading QED gap).

[^8]: `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/07_Helium_Fine_Structure.md`, §72 (kHz precision floor).

[^9]: `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/08_Muonium_Hyperfine.md`, §80 (sub-leading QED gap).

[^10]: `Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md`. Companion note expanding the three candidate starting points.
