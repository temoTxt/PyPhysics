---
title: "Research paths from the GR / Mercury-perihelion papers — a prioritized response"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-06-04"
subject: "Five candidate research paths opened by Zorn's 2026-06-04 papers, ordered by priority for the dual-theory framework"
---

# Research paths from the GR / Mercury-perihelion papers — for Zorn

**To:** Zorn (@Zorns-Lemmon).
**From:** Trey Morris (with Claude Opus 4.7).
**Date:** 2026-06-04.
**Re:** Your 2026-06-04 email[^1] on issue #139[^2]. Thank you for the weekend hunt. The papers you sent open five distinct threads. This note orders them by how directly they bear on the dual-theory framework, and gives a concrete next step for each.

---

## Cover note

This is a planning response, not a derivation. It ranks the threads your papers open and proposes a sequence. It does not yet contain results from those papers. The Bootello[^3] and Wex[^4] sources are not yet parsed into the repository, and the Edvardsson[^5] paper is access-blocked at the Lab. Where a paragraph interprets a source's relevance to the framework, that interpretation is flagged for your review.

<!-- TODO: human reviews and fills in -->

Per the lab's writing standard[^12] and AI-use compliance[^11], this report is substantive AI end-to-end. Every ranking and every framework-mapping claim is an interpretive move and is marked for your sign-off. No specific result, equation, or number is attributed to a paper we have not yet read. The standard textbook results quoted below as benchmarks are footnoted to the repository's existing perihelion documents, not to fabricated external citations.

---

## §1. The ordering, at a glance

You asked, in effect, which of these is worth pursuing. The order we propose is **1, 2, 4, 3, 5**, read as priority for the framework, not as the order you raised them.

| Rank | Path | Source | Why here | Status |
|---|---|---|---|---|
| 1 | Mercury perihelion via Bootello | Bootello[^3] | Plugs straight into the live perihelion thread[^6][^7] | **Lead** |
| 2 | Interaction-based gravity | Edvardsson[^5] | Deepest conceptual fit to the proper-time framing | **High-upside, access-blocked** |
| 3 | Pulsar-timing proper-time tests | Wex[^4] | A real experimental channel for a proper-time prediction | **Experimental** |
| 4 | Gravitoelectromagnetism | (your GEM reading) | Connects your weekend reading to our EM work | **Structural** |
| 5 | Light-bending consistency | shared | A cross-check gate, not a standalone program | **Validation** |

The reasoning for each rank follows.

<!-- TODO: human reviews and fills in -->

---

## §2. Path 1 — Mercury perihelion via Bootello (Lead)

This is the lead because it lands on top of work already in motion. The repository already carries a Mercury-perihelion derivation built on Corda's framework[^8], reviewed end to end[^6], and worked through Binet's orbit equation both the standard-GR way and the proper-time way[^7]. The benchmark is the general-relativistic anomalous advance per orbit,

$$
\Delta\phi_{\text{GR}} \;=\; \frac{6\pi G M}{c^{2}\,a\,(1-e^{2})}, \tag{1}
$$

which accounts for the observed Mercury residual of roughly forty-three arcseconds per century[^7].

The reason a second independent derivation matters is that our own thread is not yet closed. The framework's proper-time update sits in an unresolved adjudication. One route gives a Corda-style result near $+37.8$ arcseconds per century; the proper apsidal-angle route gives near $-7.2$ arcseconds per century[^6]. The structural identity behind that tension is

$$
\frac{\Delta\phi_{\text{framework}}}{\Delta\phi_{\text{GR}}} \;=\; -\frac{1}{6}, \tag{2}
$$

which falls out of the algebra as an identity, not an approximation[^7].

So Bootello is useful precisely as an outside check. If we parse it into the same Binet apparatus, we can ask whether his derivation corroborates one of our two framework numbers, reproduces the GR value (1), or stands apart from all three. A third independent perihelion calculation is the cheapest way to break the current tie.

<!-- TODO: human reviews and fills in -->

**Next step.** Download the published Bootello PDF[^3] into `External_Papers/`, parse it to markdown, and slot its force law into Binet's equation alongside the Corda and proper-time chains. One table, three derivations, same apparatus.

<!-- TODO: human reviews and fills in -->

---

## §3. Path 2 — Edvardsson's interaction-based gravity (High-upside, access-blocked)

This is ranked second because it is the deepest conceptual fit, gated only by access. You relayed that the abstract gives up the principle of equivalence and rebuilds gravity as a traditional interaction force inside special relativity, while claiming consistency with the Mercury precession and light bending[^5]. That is structurally close to the dual theory's stance. The dual theory's distinguishing feature is its proper-time, interaction-first formulation of dynamics. A gravity model that is interaction-first and equivalence-principle-free is the same shape of idea applied to gravity.

<!-- TODO: human reviews and fills in -->

The upside is high because, if the abstract holds, Edvardsson is a worked example of the framework's own philosophy in the one regime — gravitation — where the dual theory has so far said the least. The risk is that we cannot yet read it. It is in *Celestial Mechanics and Dynamical Astronomy*, the Lab has no subscription, and there is no arXiv copy[^5].

**Next step.** Try three acquisition routes in order: an author preprint or personal-page copy, an institutional or interlibrary-loan request, then a direct note to the author. Until we have the full text, we do not derive anything from it. We hold the conceptual parallel as a hypothesis, not a result.

<!-- TODO: human reviews and fills in -->

---

## §4. Path 4 — Pulsar-timing proper-time tests (Experimental)

You flagged this as fertile ground, and we agree it is the strongest *experimental* channel. Binary-pulsar timing measures the periastron advance directly, among other post-Keplerian effects, at high precision[^4]. The perihelion question of Path 1 and the precision-timing question here are the same physics observed in two systems: Mercury in the solar field, and a neutron star in a companion's field.

The framework's contribution would be a predicted fractional correction to the periastron advance,

$$
\dot{\omega}_{\text{obs}} \;=\; \dot{\omega}_{\text{GR}}\,\bigl(1 + \delta_{\text{framework}}\bigr), \tag{3}
$$

where $\delta_{\text{framework}}$ is the proper-time correction the dual theory predicts. The same $-1/6$-style structural factor that appears in (2) for Mercury should, if the framework is consistent, reappear in $\delta_{\text{framework}}$ here. A binary pulsar is where that prediction would actually be tested against data, since pulsar timing reaches a precision the solar-system perihelion does not.

<!-- TODO: human reviews and fills in -->

**Next step.** Pull the Wex review's LaTeX source[^4], identify the binary systems with the cleanest measured periastron advance, and fold a $\delta_{\text{framework}}$ forecast into the experimental-design plan[^9]. This is a longer thread than Paths 1 and 2, which is why it sits behind them.

<!-- TODO: human reviews and fills in -->

---

## §5. Path 3 — Gravitoelectromagnetism (Structural)

This connects your weekend reading to the repository's electromagnetism work[^10]. In the weak-field, slow-motion limit, linearized general relativity can be written in a Maxwell-analog form, with a gravito-electric field and a gravito-magnetic field obeying field equations of the schematic shape

$$
\nabla\cdot\mathbf{E}_g \;=\; -4\pi G\,\rho, \qquad
\nabla\times\mathbf{B}_g \;=\; -\frac{4\pi G}{c^{2}}\,\mathbf{j} + \frac{1}{c^{2}}\frac{\partial \mathbf{E}_g}{\partial t}. \tag{4}
$$

The convention-dependent numerical factors in (4) are left for the human pass to pin against a canonical reference, since they vary by a factor of up to four across the literature.

<!-- TODO: human reviews and fills in: attach the canonical GEM reference and fix the convention factors in (4) -->

The reason this is structural rather than a lead is that it is a reformulation route, not a new prediction by itself. The interesting question is whether the dual theory's existing EM structure maps onto the gravito-magnetic sector cleanly, and whether that map reproduces the perihelion advance of Path 1 by a different road. If it does, GEM becomes a unifying language for Paths 1 and 4. If it does not, it tells us where the framework's EM-gravity analogy breaks.

<!-- TODO: human reviews and fills in -->

**Next step.** Rewrite the GEM field equations (4) in the dual theory's own variables and check whether the gravito-magnetic term yields the same periastron advance the Binet route gives in Path 1.

<!-- TODO: human reviews and fills in -->

---

## §6. Path 5 — Light-bending consistency (Validation)

This is ranked last on purpose. It is a gate, not a program. Both Edvardsson[^5] and any GEM route[^10] claim to reproduce gravitational light deflection, so the standard benchmark

$$
\alpha \;=\; \frac{4 G M}{c^{2}\,b}, \tag{5}
$$

with $b$ the impact parameter, is the shared cross-check every candidate mechanism must pass. It is most useful applied to the outputs of Paths 2 and 3, rather than pursued on its own.

<!-- TODO: human reviews and fills in -->

**Next step.** Once any candidate mechanism produces a deflection prediction, compare it against (5). Treat agreement as a necessary condition, not as a finding in itself.

<!-- TODO: human reviews and fills in -->

---

## §7. Proposed sequence and what we need from you

The short version. Start Path 1 now, because it sharpens a tie we are already trying to break. Start the *acquisition* for Path 2 in parallel, because the bottleneck there is access, not effort. Treat Paths 4 and 3 as the next wave, and use Path 5 as the gate they must clear.

Two requests. First, the second arXiv link you mislaid[^1] — when you find it, drop it on issue #139 and we will fold it into the same parse. Second, if you have any back channel to the Edvardsson paper[^5], that unblocks the highest-upside thread.

<!-- TODO: human reviews and fills in -->

---

## References

[^1]: Email from Zorn (@Zorns-Lemmon), 2026-06-04, preserved verbatim in the first comment on issue #139.
[^2]: PyPhysics issue #139, "Parse Zorn's GR/Mercury papers and ship a research-paths response PDF to Zorn."
[^3]: A. Bootello, *Journal of Modern Physics* (SCIRP), paper id 145485. Published PDF and HTML at the journal site.
[^4]: N. Wex, *Testing Relativistic Gravity with Radio Pulsars*, arXiv:1402.5594.
[^5]: S. Edvardsson, *Relativistic gravitational force*, *Celestial Mechanics and Dynamical Astronomy* (no arXiv copy located; Lab subscription unavailable as of 2026-06-04).
[^6]: `Roadmapping/Author_Reports/2026-05_corda_perihelion_review_for_gill.md` — end-to-end review of the Corda framework; source of the $+37.8$ vs $-7.2$ arcsec/century adjudication.
[^7]: `Roadmapping/Author_Reports/2026-05_perihelion_classical_derivation_for_gill.md` — standard-GR and proper-time perihelion derivations via Binet's equation; source of the benchmark (1), the residual figure, and the identity (2).
[^8]: C. Corda, *The secret of planets' perihelion between Newton and Einstein*, *Physics of the Dark Universe* **32** (2021) 100834.
[^9]: `Roadmapping/Experimental_Design_Plan_Dual_Theory.md`.
[^10]: `Roadmapping/Electromagnetism/` — the framework's electromagnetism verification and derivation work.
[^11]: `Roadmapping/Tooling/CROCCO_COMPLIANCE.md`.
[^12]: `Roadmapping/Tooling/VOICE_MATCH_GILL.md`, §3.bis–§3.quinquies (equation, sentence, citation, and emoji discipline).
