---
title: "Candidates for a first-principles derivation of \\(r_e/r_0\\)"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-05-26"
subject: "Four candidate starting points for deriving r_e from the dual-Dirac framework's internal structure (follow-up to the triangulation note)"
---

# Candidates for a first-principles derivation of $r_e/r_0$ — for Tepper Gill

**Date:** 2026-05-26.
**From:** Trey Morris (with Claude Opus 4.7).
**Re:** Follow-up to the 2026-05-26 triangulation note ([`Roadmapping/Author_Reports/2026-05_re_triangulation_followup_for_gill.md`](2026-05_re_triangulation_followup_for_gill.md)); Scope 1 of [issue #54](https://github.com/temoTxt/PyPhysics/issues/54).

---

## Context

The DRQM I §III.D $r_e$ value was originally a *cutoff parameter*, fixed by a numerical search against the measured electron $g_s$. The empirical refinement in [PR #62](https://github.com/temoTxt/PyPhysics/pull/62) widened that search from one observable to six (electron $g_s$, H $2P_{3/2}-2P_{1/2}$, H 1S hyperfine, He ${}^3P_0-{}^3P_1$, positronium ortho-para, muonium hyperfine), and returned $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$. Going from one observable to six does not change the *kind* of result — it is still a cutoff fit — but it does change the result's *credibility*: the cutoff is now constrained by six independent precision measurements rather than one, and the joint fit is consistent across all six at the framework's known precision floor.

There is, at present, **no first-principles derivation of $r_e/r_0$** anywhere in the campaign or in the published DRQM I material we have visibility into. The value originated as a cutoff search, and it remains a cutoff value. The triangulation makes that cutoff more credible by widening its constraint set; it does not promote it to a derived quantity. A genuine first-principles derivation — in which $r_e$ emerges from the dual-Dirac equation's internal structure rather than being fit to data — would be a genuinely new piece of theory work.

This note collects the four candidate starting points we have considered for such a derivation, with one to three paragraphs on each describing what the route would look like, what its plausibility appears to be from outside the framework's internal logic, and what would be required to pursue it. We are **not proposing any of these** — we defer to your judgment on which (if any) is natural for the framework. The thread is non-urgent: the triangulated value already serves as the campaign's current-best-refinement at six-observable credibility, and a first-principles derivation is a "would-be-nice-to-have" rather than a load-bearing item.

---

## Candidate 1 — Proper-time self-energy integral

The dual theory's distinguishing feature is the proper-time first-order Dirac structure: the standard Dirac equation $(i\gamma^\mu \partial_\mu - m)\psi = 0$ is replaced by a proper-time-formulated version in which the dynamical variable is $\partial_\tau$ rather than $\partial_t$. The natural place to look for $r_e$ as a derived quantity is the framework's analog of the one-loop electron self-energy diagram: the diagram in which a single photon propagator and a single electron propagator close into a loop on the external electron line.

In standard QED, the one-loop self-energy $\Sigma(p)$ has a UV divergence that is absorbed into mass renormalisation via a counterterm. The cutoff procedure (Pauli–Villars, dimensional regularisation, lattice cutoff) introduces a regulator scale, and in the renormalised expression that scale's dependence cancels modulo the running of $\alpha$. In the dual-theory's proper-time formulation, the analog calculation would feature the proper-time propagator and a proper-time photon propagator. The UV behavior might be qualitatively different from standard QED: possibly finite, possibly with a different divergence structure rooted in the proper-time integral's analytic properties.

Concretely, the route would set up a Schwinger-like proper-time integral representation of the dual-Dirac propagator, compute the one-loop self-energy with the proper-time photon propagator, and identify the natural scale at which the mass counterterm equals the physical-minus-bare mass difference. If the integral is finite, that scale is determined unambiguously; if divergent, it is set by the renormalisation condition. The output $r_e/r_0$ would be a numerical function of $\alpha$ and the dual-Dirac structure constants ($b$, $c$, the projection-operator coefficients), evaluating to a definite numerical value the campaign could compare against $0.499\,420\,509\ldots$

---

## Candidate 2 — Variational determination via the renormalised dual-Dirac equation

In the dual-Dirac equation $\hat{H}_{\text{dual}} \psi = E \psi$, the cutoff $r_e$ appears as a parameter — whether as a hard cutoff in the radial integration, or as a regulator in the small-$r$ behavior of the wavefunction. The eigenvalue $E$ depends on $r_e$: $E = E(r_e)$. A variational determination would fix $r_e$ by demanding that the physical electron mass eigenvalue be recovered at the cutoff, $E(r_e) = m_e c^2$ (including any binding-energy or self-energy contributions specified by the framework).

This is analogous to how the standard-QED running coupling $\alpha(\mu)$ is fixed by the renormalisation condition $\alpha(\mu_0) = \alpha_{\text{phys}}$ at a chosen reference scale: a framework parameter is fixed by a physical observable evaluated at the framework's natural scale. The advantage of this route over Candidate 1 is that it side-steps the UV-divergence question entirely — the variational principle operates within the framework's renormalised structure, not on the divergent bare-vertex calculation. The disadvantage is that mass-condition-alone might not uniquely determine $r_e$: additional consistency conditions (gauge invariance of the photon propagator at the cutoff, current conservation at the radial boundary, perhaps the magnetic-moment relation closing the system) may be needed for closure.

If the right set of consistency conditions can be assembled, the variational route is the most direct path from the dual-Dirac equation's internal structure to a numerical $r_e/r_0$. The calculation would be a coupled eigenvalue / renormalisation-condition problem solvable symbolically (Mathematica) once the conditions are specified. The output is again comparable against the triangulated $0.499\,420\,509\ldots$ — but the route's credibility depends on whether the additional consistency conditions are framework-internally motivated or look ad hoc.

---

## Candidate 3 — Structural constant of the dual representation

The dual theory's representation carries internal structural constants — the $b/c$ ratio between the proper-time and standard-time variables, the $b$-factor projection structure in the dual-Dirac decomposition, the second-order-Dirac decomposition into $(\sigma \cdot \pi)^2$ and Pauli/Darwin terms, the operator-ordering choices in the Foldy–Wouthuysen reduction. These are not free parameters: they are determined by the framework's algebraic construction once the dual representation is specified.

A *structural-constant* derivation of $r_e/r_0$ would identify some combination of these internal constants that evaluates numerically to $\approx 0.4994$. This is the cleanest possible outcome: $r_e/r_0$ would not be a renormalisation parameter at all, but a calculable framework constant — comparable to how $\alpha = e^2/(\hbar c) \approx 1/137.036$ is a fundamental dimensionless ratio of the universe rather than a fit parameter of a calculation. If the route works, the derivation might be short — once the relevant structural ratio is identified, the value follows directly from the algebra without solving any dynamical equation.

The plausibility of this route from outside the framework's internal logic is hard for us to assess. A value close to $1/2$ does suggest a possible origin in a simple algebraic decomposition: the projection-operator eigenvalues of the dual-Dirac decomposition are obvious candidates (operator $P$ with $P^2 = P$ has eigenvalues $\{0, 1\}$; combinations like $(1 - 4\delta)$ with $\delta$ small might generate the observed slight departure from $1/2$). But whether any such combination actually arises from the framework's construction at the relevant order — that we would need your input on. If this is the natural route, it would be by far the shortest derivation of the four; if it is not, the structural-constant decomposition's eigenvalues will not match.

---

## Candidate 4 — Working-notebook derivation step not in the published prose

The published DRQM I §III.D states the anomalous-$g$-factor formula $g_r = 2[1 - 4 r_0 / (2r + r_0)]$ and the cutoff value $r_e = 0.499\,857\,150\,068\,631\,r_0$, but does not show the derivation that produced the $r_e$ value. The natural reading of the published text — confirmed by your 2026-05-25 author guidance — is that $r_e$ was obtained by an inverse-solve numerical search against measured $g_e$. But the working notebook from which the DRQM I paper was prepared may contain calculation steps that did not make it into the published prose, possibly including a derivation of $r_e$ from a framework-internal route along the lines of Candidates 1, 2, or 3.

If such a step exists, retrieving and documenting it would close Scope 1 immediately, regardless of which of the three preceding candidates it actually instantiates. The cost of investigating this is minimal: it requires your recollection or access to your working notes from the DRQM I preparation period, not new theoretical work on our side. The retrieval itself would be Crocco-pragmatic (transcribing an existing calculation, not generating a new one); only the interpretation of the retrieved step in terms of the framework's current understanding would be substantive.

The honest downside: working notebooks are often partial, ambiguous, or hard to interpret out of the moment in which they were written. Even if a relevant derivation step is recovered, additional work may be needed to formalise it into a publishable derivation. But this remains the lowest-cost-to-investigate candidate, and it has the unique property of being unambiguously decidable by one party — you. If you remember the $r_e$ value being purely a numerical search and there was no derivation-step working draft, Candidate 4 is empty and the campaign continues with the triangulated value as the framework's empirical cutoff at six-observable credibility.

---

## Closing — how this thread sits

We expect the most likely outcome, given your 2026-05-25 guidance and the campaign's posture, is that none of these four candidates rises above the triangulated cutoff in usefulness for the campaign's current scope. The triangulated value is empirically well-constrained (six precision observables; $\sigma_r \sim 10^{-13}$), agrees with measurement at the framework's known precision floor, and is sufficient for the campaign's downstream work without a first-principles derivation behind it. The Scope 1 thread is therefore non-urgent — it remains open as a "would-be-nice-to-have" for future framework development, not as a load-bearing item.

If one of the four candidates does look natural to you, we would be glad to pursue it — either via Mathematica symbolic calculation (Candidates 1, 2, or 3 with your guidance on the framework's internal logic) or via documentation work (Candidate 4 if a working-notebook step exists). If none of them looks natural, the triangulated value stands as the campaign's $r_e$ disposition and Scope 1 can be closed without a derivation. Either way, our position is that the triangulation has made the cutoff more credible by widening its constraint set, and that is the campaign's current honest-scope position on the $r_e$ question.

— Trey
