---
title: "Candidates for a first-principles derivation of \\(r_e/r_0\\)"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-05-26"
subject: "Four candidate starting points for deriving r_e from the dual-Dirac framework's internal structure (follow-up to the triangulation note)"
---

# Candidates for a first-principles derivation of $r_e/r_0$ — for Tepper Gill

**Date:** 2026-05-26 (revised 2026-05-26 with two waves of overnight progress on all three candidates).
**From:** Trey Morris (with Claude Opus 4.7).
**Re:** Follow-up to the 2026-05-26 triangulation note[^1]; Scope 1 of the master *r_e* tracker[^2]. The *Progress update* section reports findings from three parallel research branches[^3][^4][^5] under their common master[^6]. Two of the three branches have now reached definitive verdicts. **Candidate 2** is halted at Outcome D: the published §III.D non-relativistic expansion is structurally inadequate to pin $r_e$ at any scale where the expansion is valid. **Candidate 3** is halted at Outcome C-as-published: the closed-form Schwinger match is algebraically forced by the back-fit definition of the triangulated cutoff, not evidence of intentional encoding. Candidate 1 is live, scaffolded, and on a clear next-step plan.

---

## Context

The DRQM I §III.D value of $r_e$ was, originally, a cutoff parameter fixed by a numerical search against the measured electron $g_s$. The empirical refinement in a recent campaign pull request[^7] widened that search to a joint fit across six observables (electron $g_s$, hydrogen $2P_{3/2}-2P_{1/2}$, hydrogen $1S$ hyperfine, helium ${}^3P_0 - {}^3P_1$, positronium ortho-para, muonium hyperfine), and returned
$$
r_e/r_0 \;=\; 0.499\,420\,509\,912\,831\,7. \tag{1}
$$
The honest reading of (1), recorded in detail in the cross-comparison document[^8], is the following. The six observables are not six independent constraints on the cutoff. Each is computed as
$$
\mathcal{O}_n \;=\; \bigl( g_s/-2 \bigr)^n \,\times\, \mathcal{O}_n^{\text{textbook}}, \qquad n \in \{1,2\}, \tag{2}
$$
so all six share a single dependence on $g_s$, and therefore on $r_e$. The joint fit is one structural fact — the single-parameter cutoff prescription — applied to six manifestations, not six independent tests. What the multi-observable fit does confirm is that the six manifestations are mutually consistent at the framework's precision floor. The scaling (2) holds across all six. We read this as a self-consistency check on the framework's structural prescription, not as an independent corroboration of the cutoff's value.

There is, at present, no first-principles derivation of $r_e/r_0$ anywhere in the campaign, or in the published DRQM I material we have visibility into. The value originated as a cutoff search, and it remains a cutoff value. The triangulation confirms that the structure (2) is self-consistent across multiple observables under a single cutoff. It does not promote that cutoff to a derived quantity. A genuine first-principles derivation, in which $r_e$ emerges from the dual-Dirac equation's internal structure rather than being fit to data, would be a new piece of theory work. We note that Candidate 3 below offers a partial identification. The triangulated value already has a clean closed-form Schwinger reading, which recasts the derivation question into "why does the framework specify this cutoff prescription?" rather than "what numerical value should the cutoff have?"

This note collects the three candidate starting points we have considered for such a derivation. We give one to three paragraphs on each, describing what the route would look like, what its plausibility appears to be from outside the framework's internal logic, and what would be required to pursue it. A fourth candidate we had initially considered (retrieving a derivation step from the original DRQM I §III.D working notebook that did not make it into the published prose) has been ruled out by your 2026-05-25 author guidance. The original cutoff was a numerical search against $g_s$ alone, with no derivation-step working draft to retrieve.

We are not proposing any of these three. We defer to your judgment on which, if any, is natural for the framework. The thread is non-urgent. The triangulated value already serves as the campaign's current-best-refinement at the framework's precision floor. A first-principles derivation would refine the framework's standing, but it is not load-bearing for the campaign's current scope.

---

## Candidate 1 — Proper-time self-energy integral

The dual theory's distinguishing feature is the proper-time first-order Dirac structure. The standard Dirac equation
$$
(i\gamma^\mu \partial_\mu - m)\psi \;=\; 0 \tag{3}
$$
is replaced by a proper-time formulation in which the dynamical variable is $\partial_\tau$ rather than $\partial_t$. The two are mathematically equivalent but not physically equivalent; the loop-level content therefore need not coincide. The natural place to look for $r_e$ as a derived quantity is the framework's analog of the one-loop electron self-energy diagram. That is the diagram in which a single photon propagator and a single electron propagator close into a loop on the external electron line.

In standard QED, the one-loop self-energy $\Sigma(p)$ carries a UV divergence which is absorbed into mass renormalisation via a counterterm. The cutoff procedure (Pauli–Villars, dimensional regularisation, lattice cutoff) introduces a regulator scale, and in the renormalised expression that scale's dependence cancels modulo the running of $\alpha$. In the dual theory's proper-time formulation, the analog calculation would feature the proper-time propagator and a proper-time photon propagator. The UV behaviour might be qualitatively different from standard QED. It is not obvious that the integral is divergent at all; the proper-time integral's analytic properties could render it finite, or supply a different divergence structure.

Concretely, the route would set up a Schwinger-like proper-time integral representation of the dual-Dirac propagator, compute the one-loop self-energy with the proper-time photon propagator, and identify the natural scale at which the mass counterterm equals the physical-minus-bare mass difference. If the integral is finite, that scale is determined unambiguously. If divergent, it is set by the renormalisation condition. The output $r_e/r_0$ would then be a numerical function of $\alpha$ and the dual-Dirac structure constants ($b$, $c$, the projection-operator coefficients), evaluating to a definite numerical value comparable against the triangulated $0.499\,420\,509\ldots$ of (1).

---

## Candidate 2 — Variational determination via the renormalised dual-Dirac equation

In the dual-Dirac equation
$$
\hat{H}_{\text{dual}} \,\psi \;=\; E\,\psi, \tag{4}
$$
the cutoff $r_e$ appears as a parameter, either as a hard cutoff in the radial integration or as a regulator in the small-$r$ behaviour of the wavefunction. The eigenvalue $E$ depends on $r_e$, so we write $E = E(r_e)$. A variational determination would fix $r_e$ by demanding that the physical electron mass eigenvalue be recovered at the cutoff,
$$
E(r_e) \;=\; m_e c^2, \tag{5}
$$
including any binding-energy or self-energy contributions specified by the framework.

This is analogous to how the standard-QED running coupling $\alpha(\mu)$ is fixed by the renormalisation condition $\alpha(\mu_0) = \alpha_{\text{phys}}$ at a chosen reference scale. A framework parameter is fixed by a physical observable evaluated at the framework's natural scale. The advantage of this route over Candidate 1 is that it side-steps the UV-divergence question entirely. The variational principle operates within the framework's renormalised structure, not on the divergent bare-vertex calculation. The disadvantage is that the mass condition (5) alone may not uniquely determine $r_e$. Additional consistency conditions (gauge invariance of the photon propagator at the cutoff, current conservation at the radial boundary, perhaps the magnetic-moment relation closing the system) may be needed for closure.

If the right set of consistency conditions can be assembled, the variational route is the most direct path from the dual-Dirac equation's internal structure to a numerical $r_e/r_0$. The calculation would be a coupled eigenvalue and renormalisation-condition problem, solvable symbolically once the conditions are specified. The output is again comparable against the triangulated value of (1). The route's credibility, however, depends on whether the additional consistency conditions are framework-internally motivated or look ad hoc.

---

## Candidate 3 — The triangulated value already has a clean closed-form Schwinger reading

We observe that the triangulated value (1) is, to within $\sim 10^{-6}$, the closed-form value
$$
r_e/r_0 \;=\; \frac{2 - \alpha/(2\pi)}{4 + \alpha/\pi} \;=\; 0.499\,419\,632\,156\ldots, \tag{6}
$$
obtained by inverting
$$
g_r \;=\; 2\Bigl[\,1 \,-\, \frac{4\,r_0}{2r + r_0}\,\Bigr] \tag{7}
$$
against the Schwinger one-loop QED anomalous moment
$$
g_e^{(1\text{-loop})} \;=\; -2 \,-\, \alpha/\pi \;=\; -2.002\,322\,819\,465\,7\ldots \tag{8}
$$
The residual discrepancy between the closed-form $0.499\,419\,632$ of (6) and the triangulated $0.499\,420\,510$ of (1) tracks exactly the gap between (8) and the all-orders measured value
$$
g_e \;=\; -2.002\,319\,304\,362\,5\ldots \tag{9}
$$
The measured value (9) includes the Karplus–Kroll two-loop and higher corrections at $\sim 10^{-6}$ relative to (8). In other words, at one-loop QED precision, the cutoff prescription
$$
r_e \;=\; \frac{2 - \alpha/(2\pi)}{4 + \alpha/\pi}\,r_0 \tag{10}
$$
is what makes (7) reproduce the Schwinger result exactly.

This is, from outside the framework's internal logic, the cleanest and most natural reading of the cutoff identification. The DRQM I §III.D formula is, at the precision the framework's apparatus delivers, an algebraic re-encoding of the Schwinger one-loop anomalous moment through a particular cutoff substitution. The reason the triangulation across six observables returns essentially the same value as the uni-observable back-fit against $g_s$ alone is that every other observable's prediction follows the scaling (2). Substituting the measured $g_s$ is the same operation whether done once or six times. The measured $g_s$ is, to one-loop QED precision, exactly what the closed-form cutoff (10) is engineered to reproduce.

On this reading, the first-principles derivation question reduces to: why does the framework specify this particular cutoff prescription? Two sub-questions follow. First, is there a derivation in which the dual-Dirac equation's renormalisation structure produces (10) as the natural mass-renormalisation scale (Candidates 1 or 2 applied at one-loop precision)? Second, is the cutoff a structural constant of the dual representation itself, such that the combination $(2 - \alpha/(2\pi))/(4 + \alpha/\pi)$ comes out of the framework's algebra without reference to QED at all? The first is a one-loop calculation that should be tractable if the framework's renormalisation prescription is fully specified. The second requires identifying the framework-internal origin of an $\alpha$-dependent ratio, which is harder. If the first reading is the natural one, the campaign can attempt the calculation directly. If the second, we would need your guidance on which structural decomposition produces $\alpha$-dependent ratios in the framework.

We did not see this identification in the published DRQM I §III.D text. If it is a known property of the framework, please tell us; that would simplify the disposition considerably. If it is a coincidence at the campaign's precision (the triangulated value happens to be close to the Schwinger closed-form but is not engineered to be), the Karplus–Kroll-level residual is the test. With more precision-spectroscopy observables in the joint fit, the residual should not systematically point at the two-loop QED corrections unless the cutoff genuinely is encoding QED.

---

## Progress update (2026-05-26 PM) — three parallel research branches

Since the morning version of this note went out, three parallel research branches[^3][^4][^5] under their common master[^6] have begun attempting the three candidates. Per Crocco compliance, the work below is substantive AI, with per-paragraph human-review markers in the per-branch state logs[^9]. The headline findings are summarised here for your morning read.

### Candidate 1 (proper-time self-energy)[^3] — live; conceptual re-framing of $r_e$, baseline scaffolded

Four iterations: read DRQM I §II–III, *The Classical Electron Problem*, and the Bethe–Salpeter Lamb-shift document; scaffolded the dedicated notebook[^10] with a baseline cell. The baseline standard-QED on-shell shift,
$$
\frac{\delta m}{m} \;=\; \frac{3\alpha}{4\pi}\,\Bigl[\,\log\!\frac{\Lambda^2}{m^2} \,+\, \tfrac{1}{2}\,\Bigr], \tag{11}
$$
is the Bjorken–Drell Eq. (10.59) result, verified symbolically.

The headline finding is conceptual rather than computational. The existing Lamb-shift calculation uses the textbook non-relativistic Bethe (1947) UV cutoff $K \sim m c^2$, equivalently the Compton length $\lambda_C = \hbar/(mc)$. The triangulated $r_e$ sits well below that scale,
$$
r_e \;\sim\; r_0/2 \;=\; (\alpha/2)\,\lambda_C, \tag{12}
$$
so the parametric separation is
$$
r_e / \lambda_C \;\sim\; \alpha/2 \;\approx\; 3.65 \times 10^{-3}. \tag{13}
$$
We conclude that $r_e$ is not a UV-loop cutoff at all. It sits at the Coulomb-binding scale, not the Compton scale. Indeed, at $r = r_0/2$,
$$
\lvert V_0 \rvert \;=\; 2\,m c^2, \tag{14}
$$
which is exactly the pair-production threshold. The natural re-reading is that $r_e$ marks the radius inside which the bound-state wavefunction picks up virtual $e^+ e^-$ contributions, and inside which the §III.D small-component elimination
$$
\psi_2 \;=\; \frac{c\,(\boldsymbol{\sigma}\!\cdot\!\boldsymbol{\pi})\,\psi_1}{\lambda - V_0 + m c^2} \tag{15}
$$
stops being a valid approximation.

Technically, the (II.3) "potential-in-the-mass" form is the cleanest kernel. There are no $V \boldsymbol{\alpha}\!\cdot\!\boldsymbol{\pi}/(mc)$ or $\boldsymbol{\alpha}\!\cdot\!\boldsymbol{\nabla} V/(mc)$ pieces (those arise from operator non-commutativity in the (II.1) Dirac form), so the kernel reduces to the Pauli kinetic plus $V^2/(2 m c^2)$.

A heuristic sanity check at iteration 4 sharpens the discriminator. The naive identification $\Lambda = \hbar/(r_e c)$ at $r_e/r_0 = 0.5$ gives
$$
\log(\Lambda^2/m^2) \;=\; 11.23, \qquad \delta m/m \;=\; 2.04 \times 10^{-2}, \tag{16}
$$
versus the natural one-loop coupling $\alpha/(4\pi) = 5.81 \times 10^{-4}$ — a $35\times$ overestimate. Two distinguishable resolutions are queued for the next cells. First, the framework may supply $\Lambda \sim \hbar/(b\,r_e)$ with $b > c$ for bound states, suppressing the log. Second, the Bethe-log replacement $\log(K/\langle\Delta E\rangle)$ at fixed photon-loop measure (sum-over-states log) is parametrically smaller. The two are distinguishable. The first shifts $r_e/r_0$ via $\log(b/c)$ at fixed $\langle p^2\rangle$; the second shifts $r_e/r_0$ via the Bethe-log replacement at fixed photon-loop measure. **Status: live, on plan.** The outcome-matrix branch remains open. The next two cells are expected to discriminate the two resolutions and produce a numerical $r_e/r_0$.

### Candidate 2 (variational route)[^4] — halted at Outcome D; the published NR-expansion is structurally inadequate to pin $r_e$

Four iterations enumerated seven closure conditions. Only the mass-renormalisation condition (closure #7) was both framework-internal and potentially of sufficient precision. We pursued it under the textbook working hypothesis $\Delta E_{\rm bind}^{\rm framework} = \langle V_0\rangle$, $\Delta E_{\rm SE}^{\rm framework} = 0$, with trial $\psi_1 = N e^{-r/aa}$ on the cutoff-restricted domain $[r_e, \infty)$. The closure equation, in the dimensionless variables $\hat{a} = aa/r_0$ and $\hat{r}_e = r_e/r_0$, is
$$
E_{\rm dim}(\hat{a}, \hat{r}_e) \;=\; \frac{1}{2\alpha^2 \hat{a}^2} \;-\; \frac{\hat{a} + 2\hat{r}_e - 1}{\hat{a}^2 + 2\hat{a}\hat{r}_e + 2\hat{r}_e^2} \;=\; 0. \tag{17}
$$

The diagnostic table at $\alpha = 1/137.035999$ exposes the structural inadequacy.

| Regime | $\hat{a}$ | $\hat{r}_e$ | $E_{\rm dim}$ | Verdict |
|---|---|---|---|---|
| Electron-radius scale | $1$ | $0.5$ | $+9389.0$ | Kinetic dominates by four orders; **NR expansion invalid**. |
| Bohr scale | $1/\alpha^2 \approx 18\,778$ | $0.5$ | $-2.662\times 10^{-5}$ | Matches $-\alpha^2/2$; **cutoff invisible**. |
| Bohr scale, no cutoff | $1/\alpha^2$ | $0$ | $-2.662\times 10^{-5}$ | Identical to $\hat{r}_e = 0.5$ to ten significant figures. |

At the cutoff scale, both expansion parameters of (III.4), namely $V_0/(m c^2)$ and $\hbar/(m c\,r)$, are $O(1)$. The expansion is invalid there. At the Bohr scale, where the expansion is valid, the trial mass-density $r^2 e^{-2r/aa}$ peaks at $r \sim aa \sim 18\,000\,r_0$, suppressing the cutoff coupling by $\alpha^2$. It follows that no scale exists at which both (a) the NR expansion is valid and (b) the cutoff couples meaningfully to the closure. The published expanded $K_D$ cannot pin $r_e$. **Status: halted, Outcome D, blocked.** Two specific clarifications would unblock. First, a framework-internal $\Delta E_{\rm SE}^{\rm framework}(r_e)$ at the cutoff scale. Second, a redirect to the un-expanded full $H_D$ under a radial-cutoff regulator, which would be a five- to ten-iteration arc, distinct from any previously enumerated route, and not initiated without your endorsement. See Question 2 below.

### Candidate 3 (closed-form Schwinger)[^5] — halted at Outcome C-as-published; the KK + LR + KF residual match is algebraically forced

Four iterations executed the empirical-test path. At twenty-digit precision with $\alpha = 7.297\,352\,569\,3 \times 10^{-3}$, we record the comparison.

| Quantity | Value |
|---|---|
| $r_e/r_0$ closed-form $(2-\alpha/(2\pi))/(4+\alpha/\pi)$ | $0.499\,419\,632\,156\,99$ |
| $r_e/r_0$ triangulated (PR-refined) | $0.499\,420\,509\,912\,83$ |
| $\Delta g_e^{\rm obs} = g_e^{\rm meas} - g_e^{\rm Schwinger}$ | $+3.5151 \times 10^{-6}$ |
| Karplus–Kroll two-loop $+2 C_2 (\alpha/\pi)^2$, $C_2 = 0.328\,478\,965\,579$ | $+3.5446 \times 10^{-6}$ |
| Laporta–Remiddi three-loop $-2 C_3 (\alpha/\pi)^3$, $C_3 = 1.181\,241\,456\,587$ | $-2.961 \times 10^{-8}$ |
| Kinoshita–Fukuda four-loop $+2 C_4 (\alpha/\pi)^4$, $C_4 \approx 1.9106$ | $+1.1 \times 10^{-10}$ |
| Sum (KK + LR + KF) | $+3.5151 \times 10^{-6}$ |
| Residual: $\Delta g_e^{\rm obs} - \Delta g_e^{\rm pred,\,all}$ | $-9.7 \times 10^{-12}$ |

Iteration 3 read DRQM I §III.D Eqs. (III.18)–(III.23) line by line and confirmed that the as-published §III.D derivation does not produce $r_e/r_0$ as a closed form in $\alpha$. The numerical value at (III.22), $r_e = 0.499\,857\,150\,068\,631\,r_0$, is presented as an empirical fact (the uni-observable numerical search you confirmed on 2026-05-25), not derived. Consequently, the $10^{-11}$ agreement above is algebraically forced by the back-fit definition of the triangulated $r_e$. Any cutoff satisfying $g_r(r) = g_e^{\rm meas}$ will, when subtracted from the closed-form one-loop value, yield via $dg_r/dr$ propagation a $\Delta g_e$ identical to the all-orders-QED-beyond-one-loop content of measured $g_e$. That is, it must yield KK + LR + KF. The KK-consistency observation is necessary but not sufficient for intentional encoding. The as-published apparatus does not supply the sufficient piece.

**Status: halted, Outcome C-as-published. Finding 2 stays Marginal (characterised but not promoted).** The canonical record is now committed as a new subsection of the DRQM I verification document[^11] under §III.D ("Schwinger identification — empirical residual test"). The path to Outcome B (intentional Schwinger encoding implies Finding 2 promoted to **Pass**) now passes solely through a first-principles rederivation that produces the closed-form as a derived identity, that is, through Candidate 1 if it lands on the closed form, or through a sub-route of Candidate 2 if author input redirects to the un-expanded full $H_D$.

### Consolidated questions for you

After this overnight run, two questions remain load-bearing (the third was resolved by the iter-3 §III.D line-by-line read). Even one answer accelerates the remaining live branch.

1. **(Candidate 1, live branch.)** Does the framework specify a proper-time photon propagator form? Two natural candidates surface. First, a Schwinger proper-time form with $b$ replacing $c$ in the dispersion,
$$
k^2 \;=\; (\omega/b)^2 \,-\, \mathbf{k}^2. \tag{18}
$$
Second, the standard Feynman propagator with the $b/c$ conversion absorbed into the source coupling rather than the propagator. The two give different numerical $r_e/r_0$ at the same order in $\alpha$. The iter-4 sanity check (16) shows the difference is at the $\log(b/c)$ level, and hence discriminable.

2. **(Candidate 2, halted, unblocking.)** The published expanded $K_D$ of (III.4) is structurally inadequate to pin $r_e$, as the diagnostic table above shows. Two unblocking moves present themselves. **(2a)** Does the framework specify a $\Delta E_{\rm SE}^{\rm framework}(r_e)$ at the cutoff scale that we should have been carrying in $\langle K_D \rangle = m_e c^2 + \Delta E_{\rm bind} + \Delta E_{\rm SE}$? **(2b)** Is the variational determination intended to operate on the un-expanded full $H_D$ under a radial-cutoff regulator (radial-Dirac eigenvalue problem on $r \in [r_e, \infty)$) rather than the expanded $K_D$ of (III.4)? The un-expanded route is a five- to ten-iteration arc we have not begun without your endorsement.

The original Question 3 from earlier today (does §III.D derive $r_e/r_0$ as a closed form in $\alpha$?) was answered by the iter-3 line-by-line read. As published, §III.D does not derive a closed form; the numerical value is an empirical input from your uni-observable numerical search. The remaining open piece on the closed-form-encoding question is whether an in-progress or planned rederivation (per Candidate 1, or a redirected Candidate 2) would produce the closed form. That is now the substantive thread for the master tracker[^6].

---

## Closing — how this thread sits

We expect the most likely outcome, given your 2026-05-25 guidance and the campaign's posture, is that none of these three candidates rises above the triangulated cutoff in usefulness for the campaign's current scope. The triangulated value is empirically well-constrained at the framework's precision floor (the joint-fit consistency across six manifestations of the scaling (2)), agrees with measurement at that floor, and is sufficient for the campaign's downstream work without a first-principles derivation behind it. The Scope 1 thread is therefore non-urgent. It remains open as a refinement target for future framework development, not as a load-bearing item. Candidate 3's Schwinger closed-form identification, if intentional in the framework's construction, may already be the answer.

If one of the three candidates does look natural to you, we would be glad to pursue it via Mathematica symbolic calculation with your guidance on the framework's internal logic. If none of them looks natural, the triangulated value (1) stands as the campaign's $r_e$ disposition and Scope 1 can be closed without a derivation. Either way, our position is that the triangulation has confirmed the structure (2) is self-consistent under a single cutoff at the framework's precision floor; that is the campaign's current honest-scope position on the $r_e$ question. A first-principles derivation would be a refinement, not a prerequisite.

The Progress update above adds three substantive pieces of physical content beyond the morning version. First, the re-framing of $r_e$ from "UV cutoff" to "Coulomb-binding scale at the pair-production threshold" (Candidate 1) — a statement about what $r_e$ is in the framework's small-component-elimination apparatus, independent of which derivation route eventually settles its first-principles status, and worth your read even if Scope 1 closes without a derivation. Second, the structural inadequacy of the published expanded $K_D$ for variational $r_e$ determination (Candidate 2, Outcome D). At no scale do (a) the validity of the (III.4) NR expansion and (b) meaningful coupling of the cutoff to the closure equation coexist. Third, the algebraically-forced character of the Schwinger closed-form agreement (Candidate 3, Outcome C-as-published). The $10^{-11}$ KK + LR + KF residual match is necessary but not sufficient for intentional encoding, and the as-published §III.D does not provide the sufficient piece. The disposition of Scope 1 now turns on the live Candidate 1 branch (and on whether Candidate 2 redirects to the un-expanded full-$H_D$ arc with your endorsement).

— Trey

---

## References

[^1]: Triangulation note: `Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md`.
[^2]: Issue #54 (Scope 1 of the $r_e$ tracker), `github.com/temoTxt/PyPhysics/issues/54`.
[^3]: Issue #64 (Candidate 1, proper-time self-energy), `github.com/temoTxt/PyPhysics/issues/64`.
[^4]: Issue #65 (Candidate 2, variational route), `github.com/temoTxt/PyPhysics/issues/65`.
[^5]: Issue #66 (Candidate 3, closed-form Schwinger), `github.com/temoTxt/PyPhysics/issues/66`.
[^6]: Issue #67 (master tracker for #64–#66), `github.com/temoTxt/PyPhysics/issues/67`.
[^7]: Pull request #62 (empirical refinement of $r_e/r_0$ via joint fit), `github.com/temoTxt/PyPhysics/pull/62`.
[^8]: Cross-comparison document, §2 ("The $r_e$ self-consistency across six $g_s$-dependent observables at the triangulated value"): `Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md`.
[^9]: Per-branch state log on each Candidate's branch: `.dev/research/STATE.md`.
[^10]: Mathematica notebook for the Candidate 1 derivation: `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation_self_energy.wl`.
[^11]: DRQM I verification document, §III.D subsection ("Schwinger identification — empirical residual test"): `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md`.
