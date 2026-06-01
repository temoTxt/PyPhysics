---
title: "Three framework-specification questions to unblock first-principles \\(r_e\\) derivation"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-05-26"
subject: "Three framework-specification questions for hypothesis-(i) first-principles r_e derivation (follow-up to candidates note + overnight progress)"
---

# Three framework-specification questions for first-principles $r_e$ derivation — for Tepper Gill

**Date:** 2026-05-26 (PM, after the day's overnight progress runs).
**From:** Trey Morris (with Claude Opus 4.7).
**Re:** Follow-up to the earlier candidates note[^1] and triangulation note[^2]; tracks against the framework-specifications issue[^3] under the master derivation thread[^4].

---

We consolidate here the three framework-specification questions that surfaced as load-bearing across the three parallel candidate branches for a first-principles derivation of $r_e/r_0$. Each question, if answered, would unblock a derivation that goes beyond the algebraic inheritance from QED's $a_e(\alpha)$ that the present closed-form result delivers.

## What the three branches converged on

Candidate 3 confirmed, by a line-by-line reading of §III.D Eqs. (III.18)–(III.23)[^5], that the as-published apparatus leaves $r_e/r_0$ as a free parameter. The algebraic-back-fit match to the closed form

$$
\frac{r_e}{r_0} \;=\; \frac{2 - \alpha/(2\pi)}{4 + \alpha/\pi}
\tag{1}
$$

is necessary but not sufficient for intentional Schwinger encoding. Iter-5 of the same branch added a cross-particle test against the FNAL Muon $g-2$ 2023 measurement, yielding

$$
\frac{r_\mu}{r_0^\mu} \;-\; \frac{r_e}{r_0^e} \;=\; 3.13 \times 10^{-6}
\tag{2}
$$

at roughly $57{,}000\,\sigma_{a_\mu}$. It follows that any universal closed form $r/r_0 = f(\alpha)$ in $\alpha$ alone is empirically excluded.

Candidate 2 found that the published expanded $K_D$ of Eq. (III.4) is structurally inadequate to pin $r_e$ via the variational principle[^6]. No scale exists where both (a) the NR expansion is valid and (b) the cutoff couples meaningfully to the closure equation. At the electron-radius scale ($\hat r_e \sim 1$) the kinetic-energy term dominates by four orders, so the NR expansion is invalid. At the Bohr scale the expansion is valid, but the cutoff $\hat r_e$ is suppressed by $\alpha^2$, and the closure matches plain hydrogen to ten significant figures.

Candidate 1 delivered the closed form

$$
\frac{r_e}{r_0} \;=\; \frac{2 - a_e}{2(2 + a_e)}
\tag{3}
$$

matching the triangulated value to $3.45 \times 10^{-13}$[^7]. Equation (3) is the algebraic inverse of the framework's existing $g$-formula (Eq. III.22) evaluated at QED's $a_e(\alpha)$ — that is, reproduction-by-inheritance under what we have been calling hypothesis (ii), in which the photon propagator is unchanged and all dual structure is absorbed into the (II.3) Pauli source. We observe that the framework does not independently re-derive Schwinger's vertex result; it inherits it.

The convergent gap across all three candidates is a first-principles derivation of $a_e$ from framework-internal dynamics. Such a derivation requires three framework specifications that the published apparatus does not currently supply. We are not asking you to perform the derivation; we are asking whether the framework, as you intend it, specifies (or implicitly fixes) the three pieces below. Each is small enough to answer in a few sentences.

## The three questions

### 1. The proper-time photon propagator $D_F^{(\tau)}(x-y)$

Two natural readings emerge from the framework's existing structure. The first is Schwinger proper-time with $b$-dispersion,

$$
k^2 \;=\; (\omega/b)^2 \;-\; \mathbf{k}^2,
\qquad
b \;=\; \sqrt{c^2 + \mathbf{u}^2},
\tag{4}
$$

in which the photon dispersion inherits $b$-factors from the dual-current structure $J^\mu_{\rm Gill} = (b\rho,\,\mathbf{J})$[^8]. The second is the standard Feynman propagator,

$$
k^2 \;=\; (\omega/c)^2 \;-\; \mathbf{k}^2,
\tag{5}
$$

with the $b/c$ conversion absorbed into the source's coupling rather than into the propagator, and the dual structure entering through modified vertex factors. The two readings give different numerical $r_e/r_0$ at the same order in $\alpha$. At the $\log(b/c)$ level they are discriminable at the framework's precision. We ask which (if either) is the framework's intended choice.

### 2. The bound-state propagator $G_C(x,y;E)$ in the (II.3) "potential-in-the-mass" kernel

For the bound-electron one-loop self-energy, we need a propagator that handles the Coulomb potential non-perturbatively at the cutoff scale, where the NR expansion of Eq. (III.4) is invalid (per Candidate 2's iter-4 diagnostic). Standard QED uses the Furry-picture Coulomb-bound propagator, summing the external Coulomb potential to all orders in $Z\alpha$ before perturbing in $\alpha$. We ask whether the framework specifies (or implicitly uses) a dual analog of the Furry-picture construction, or some other treatment of the bound-state propagator under the (II.3) kernel.

### 3. The mass-renormalisation prescription at the cutoff

The natural framework-internal closure is

$$
\langle K_D \rangle_{r_e}
\;=\; m_e c^2 \;+\; \Delta E_{\rm bind}^{\rm framework} \;+\; \Delta E_{\rm SE}^{\rm framework}.
\tag{6}
$$

The textbook NR working hypothesis ($\Delta E_{\rm bind} = \langle V_0\rangle$, $\Delta E_{\rm SE} = 0$) gives a trivially-zero closure under the published expanded $K_D$, per Candidate 2's iter-4 result. We ask what the framework's intended form is for $\Delta E_{\rm SE}^{\rm framework}(r_e)$ at the cutoff scale — the analog of QED's on-shell mass renormalisation condition $m_{\rm phys} = m + \Sigma(p)\big|_{p^2 = m^2 c^2}$ in the dual framework's proper-time formulation.

## Posture

Any one of these would advance the corresponding candidate; all three together would let us attempt the full hypothesis-(i) calculation. The thread is non-urgent, in the sense that the triangulated value

$$
r_e/r_0 \;=\; 0.499\,420\,509\,912\,831\,7
\tag{7}
$$

already serves as the campaign's current-best-refinement at the framework's precision floor[^2]. It is consistent with the framework's apparatus at every observable in the joint fit. The substantive content of the three questions is whether the framework, as you intend it, contains the specifications that would close the gap from empirical cutoff to framework-derived quantity, or whether closing that gap is genuinely new theory work the published apparatus does not yet support.

If the answer to any of (1), (2), (3) is that the framework does not yet specify the relevant piece, that is itself a useful finding. It tells us the campaign's terminal verdict on $r_e/r_0$ should be recorded as **Branch B, QED-inheritance, pending future framework development**, rather than **Branch A, from first principles**. Either outcome closes the question honestly; the questions here are about which honest framing is the right one.

We will record your responses on the framework-specifications issue[^3] and pick up the derivation work in a follow-on branch once the specifications are in hand.

— Trey

<!-- TODO: human reviews and fills in — Trey to verify the email's framing of (1)/(2)/(3) before sending; confirm the "non-urgent" posture is the intended tone; decide whether to include the specific cross-references (PR numbers, iter-N details) or to keep the email shorter and link out to the GitHub thread for those -->

## References

[^1]: Candidates note: `Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md`.
[^2]: Triangulation note: `Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md`.
[^3]: Issue #75 (framework-specifications), `github.com/temoTxt/PyPhysics/issues/75`.
[^4]: Master derivation thread, issue #67, `github.com/temoTxt/PyPhysics/issues/67`.
[^5]: Candidate 3: issue #66 and PR #70, `github.com/temoTxt/PyPhysics/issues/66`, `github.com/temoTxt/PyPhysics/pull/70`. Includes the line-by-line read of DRQM I §III.D Eqs. (III.18)–(III.23) and the FNAL Muon $g-2$ 2023 cross-particle test (iter-5).
[^6]: Candidate 2: issue #65 and PR #69, `github.com/temoTxt/PyPhysics/issues/65`, `github.com/temoTxt/PyPhysics/pull/69`. Includes the iter-4 diagnostic on the NR-expansion / cutoff-scale tension.
[^7]: Candidate 1: issue #64 and PR #72, `github.com/temoTxt/PyPhysics/issues/64`, `github.com/temoTxt/PyPhysics/pull/72`. Closed-form derivation of $r_e/r_0$ by algebraic inversion of Eq. (III.22) of DRQM I, evaluated at QED's $a_e(\alpha)$.
[^8]: Dual current $J^\mu_{\rm Gill} = (b\rho,\,\mathbf{J})$: Gill–Zachary, *Two Mathematically Equivalent Versions of Maxwell's Equations*, Eq. (12), in `Roadmapping/Converted_Markdown/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations/`.
