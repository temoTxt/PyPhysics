---
title: "Three framework-specification questions to unblock first-principles \\(r_e\\) derivation"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-05-26"
subject: "Three framework-specification questions for hypothesis-(i) first-principles r_e derivation (follow-up to candidates note + overnight progress)"
---

# Three framework-specification questions for first-principles $r_e$ derivation — for Tepper Gill

**Date:** 2026-05-26 (PM, after the day's overnight progress runs).
**From:** Trey Morris (with Claude Opus 4.7).
**Re:** Follow-up to the [candidates note](2026-05_re_derivation_candidates_for_gill.md) and the [triangulation note](2026-05_re_triangulation_followup_for_gill.md). Tracks against [issue #75](https://github.com/temoTxt/PyPhysics/issues/75) under master [#67](https://github.com/temoTxt/PyPhysics/issues/67).

---

Tepper — following the candidates note we sent earlier today (the PM revision adds overnight progress from three parallel research branches), I want to consolidate the three framework-specification questions that surfaced as load-bearing across all three candidate branches. Each of these would unblock a first-principles derivation of $r_e/r_0$ that goes beyond the algebraic inheritance from QED's $a_e(\alpha)$ that the current closed-form result delivers.

## What the three branches converged on

- **Candidate 3 ([#66](https://github.com/temoTxt/PyPhysics/issues/66), [PR #70](https://github.com/temoTxt/PyPhysics/pull/70))** confirmed via a line-by-line read of §III.D Eqs. (III.18)–(III.23) that the as-published apparatus leaves $r_e/r_0$ as a free parameter — the algebraic-back-fit match to the closed-form $r_e/r_0 = (2 - \alpha/(2\pi))/(4 + \alpha/\pi)$ is necessary but not sufficient for intentional Schwinger encoding. Iter-5 of that branch added a cross-particle test against the FNAL Muon $g-2$ 2023 measurement: $r_\mu/r_0^\mu - r_e/r_0^e = 3.13\times 10^{-6}$ at $\sim 57{,}000\,\sigma_{a_\mu}$, which **empirically excludes any universal closed-form $r/r_0 = f(\alpha)$ in $\alpha$ alone**.

- **Candidate 2 ([#65](https://github.com/temoTxt/PyPhysics/issues/65), [PR #69](https://github.com/temoTxt/PyPhysics/pull/69))** found that the published expanded $K_D$ of (III.4) is structurally inadequate to pin $r_e$ via the variational principle: no scale exists where both (a) the NR expansion is valid AND (b) the cutoff couples meaningfully to the closure equation. At the electron-radius scale ($\hat{r}_e \sim 1$), kinetic energy dominates by four orders, so the NR expansion is invalid; at the Bohr scale where the expansion is valid, the cutoff $\hat{r}_e$ is suppressed by $\alpha^2$ and the closure matches plain hydrogen to 10 sig figs.

- **Candidate 1 ([#64](https://github.com/temoTxt/PyPhysics/issues/64), [PR #72](https://github.com/temoTxt/PyPhysics/pull/72))** delivered the closed form $r_e/r_0 = (2 - a_e)/(2(2 + a_e))$ matching the triangulated value to $3.45\times 10^{-13}$. The closed form is the algebraic inverse of the framework's existing g-formula (Eq. III.22) evaluated at QED's $a_e(\alpha)$ — i.e., *reproduction-by-inheritance* under what we have been calling hypothesis (ii) (photon propagator unchanged, all dual structure absorbed into the (II.3) Pauli source). The framework does not independently re-derive Schwinger's vertex result; it inherits it.

The convergent gap across all three is **a first-principles derivation of $a_e$ from framework-internal dynamics**, which requires three framework specifications the published apparatus does not currently supply. We are not asking you to do the derivation — we are asking whether the framework specifies, or has implicit choices for, these three pieces. Each is small enough to be answered in a few sentences.

## The three questions

### 1. The proper-time photon propagator $D_F^{(\tau)}(x-y)$

Two natural readings emerge from the framework's existing structure:

- **(i-a) Schwinger proper-time with $b$-dispersion**: $k^2 = (\omega/b)^2 - \mathbf{k}^2$, with $b = \sqrt{c^2 + \mathbf{u}^2}$. The photon dispersion inherits $b$-factors from the dual-current structure ($J^\mu_{\rm Gill} = (b\rho, \mathbf{J})$, Maxwell paper Eq. 12).
- **(i-b) Standard Feynman propagator with $b/c$ conversion absorbed into the source's coupling** rather than the propagator: $k^2 = (\omega/c)^2 - \mathbf{k}^2$, with the dual structure entering through modified vertex factors.

The two give different numerical $r_e/r_0$ at the same order in $\alpha$ — at the $\log(b/c)$ level, discriminable at the framework's precision. Which (if either) is the framework's intended choice?

### 2. The bound-state propagator $G_C(x,y;E)$ in the (II.3) "potential-in-the-mass" kernel

For the bound-electron one-loop self-energy, we need a propagator that handles the Coulomb potential non-perturbatively at the cutoff scale, where the NR expansion of (III.4) is invalid (per Candidate 2's iter-4 diagnostic above). Standard QED uses the Furry-picture Coulomb-bound propagator, summing the external Coulomb potential to all orders in $Z\alpha$ before perturbing in $\alpha$. Does the framework specify (or implicitly use) a dual analog of the Furry-picture construction, or some other treatment of the bound-state propagator under the (II.3) kernel?

### 3. The mass-renormalisation prescription at the cutoff

The natural framework-internal closure is

$$\langle K_D\rangle_{r_e} \;=\; m_e c^2 \;+\; \Delta E_{\rm bind}^{\rm framework} \;+\; \Delta E_{\rm SE}^{\rm framework}.$$

The textbook NR working hypothesis ($\Delta E_{\rm bind} = \langle V_0\rangle$, $\Delta E_{\rm SE} = 0$) gives a trivially-zero closure under the published expanded $K_D$ (per Candidate 2's iter-4 result above). What is the framework's intended form for $\Delta E_{\rm SE}^{\rm framework}(r_e)$ at the cutoff scale — the analog of QED's on-shell mass renormalisation condition $m_{\rm phys} = m + \Sigma(p)|_{p^2 = m^2c^2}$ in the dual framework's proper-time formulation?

## Posture

Any one of these would let us advance the corresponding candidate; all three would let us attempt the full hypothesis-(i) calculation. The thread is **non-urgent** in the sense that the triangulated $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ already serves as the campaign's current-best-refinement at the framework's precision floor and is consistent with the framework's apparatus at every observable in the joint fit. The substantive content of the questions is whether the framework's apparatus, as you intend it, *contains* the specifications that would close the gap from "empirical cutoff" to "framework-derived quantity" — or whether closing that gap is genuinely new theory work that the published apparatus does not yet support.

If the answer to any of (1), (2), (3) is "the framework does not yet specify this," that is itself a useful finding: it tells us the campaign's terminal verdict on $r_e/r_0$ is **Branch B with QED-inheritance, pending future framework development** rather than **Branch A from first principles**. Either outcome closes the question honestly; the questions here are about which honest framing is the right one.

We will record your responses on [issue #75](https://github.com/temoTxt/PyPhysics/issues/75) and pick up the derivation work in a follow-on branch once the specifications are in hand.

<!-- TODO: human reviews and fills in — Trey to verify the email's framing of (1)/(2)/(3) before sending; confirm the "non-urgent" posture is the intended tone; decide whether to include the specific cross-references (PR numbers, iter-N details) or to keep the email shorter and link out to the GitHub thread for those -->

— Trey
