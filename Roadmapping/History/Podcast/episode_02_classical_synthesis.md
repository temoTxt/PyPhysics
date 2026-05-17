---
episode: 02
title: "Classical Synthesis 1860-1900"
era: "1860-1900"
chapter: 02_classical_synthesis_1860_1900
speakers: [Historian, Physicist, Experimentalist]
target_runtime_min: 20
animations_cued:
  - hist_maxwell_synthesis
  - hist_michelson_morley
status: draft
---

# Episode 2 — Classical Synthesis (1860–1900)

> Companion dialogue script for [[02_classical_synthesis_1860_1900]]. First episode where the standard and dual framings *diverge* — and where the divergence has experimental teeth.

## Cold open

**Historian:** In November 1864, a 33-year-old Scottish physicist named James Clerk Maxwell stands before the Royal Society in London and reads a paper he calls *A Dynamical Theory of the Electromagnetic Field*. Twenty coupled differential equations. Five field quantities. A claim — almost a side remark — that the equations admit transverse wave solutions propagating at a speed that, within the precision of the constants then known, equals the speed of light.

**Experimentalist:** Maxwell does not stand at the Royal Society and shout "light is an electromagnetic wave!" That would be the cleaner story. He's careful. He writes — and I'm reading from the paper, paragraph 79 — that "we can scarcely avoid the inference that light consists in the transverse undulations of the same medium which is the cause of electric and magnetic phenomena." Avoid the inference. He's hedging, because in 1864 nobody has *seen* an electromagnetic wave other than light.

**Historian:** That confirmation takes another twenty-three years. Heinrich Hertz, in 1888, generates and detects them in his laboratory in Karlsruhe. We'll get to him.

**Physicist:** And the reminder I gave at the start of episode one — that this whole series is exploring a particular mathematical reframing of established physics, not proposing new physics — that reminder matters more this episode than last. Last episode was a null comparison. Everything held in standard and dual electromagnetism. *This* episode is where the dual framing first picks up real teeth.

**Experimentalist:** Specifically because of an experiment in Cleveland, Ohio, in the summer of 1887. Michelson and Morley try to detect the Earth's motion through the luminiferous aether and find — nothing. Forty times smaller than predicted. The standard resolution, after Einstein, is special relativity with length contraction. Tepper Gill's dual framework reaches the same null prediction by a *different* route. Same outcome. Different kinematics underneath.

**Physicist:** Both framings reproduce Michelson–Morley's null result. Where they diverge is in regimes where the source is itself moving at a non-negligible fraction of $c$. We'll spell that out at the end of the episode. For now: Maxwell, Hertz, Michelson, Lorentz, Thomson. Plus the first appearance of solid-state physics, with Edwin Hall. Forty-five years that close the classical era.

**Historian:** Let's start with Maxwell.

## Historical sweep

**Historian:** Maxwell reads [[faraday1832_experimental_researches_i|Faraday's *Experimental Researches*]] as a Cambridge undergraduate in the early 1850s. He's struck by the *lines of force* — Faraday's geometric, qualitative ontology. The Continental tradition — Ampère, Weber, Neumann — works in action-at-a-distance: charges exert forces on each other across empty space. Faraday's lines are different. They fill space continuously. They have *direction*. They have *tension*. Faraday believes they're real physical things, not just convenient diagrams.

**Physicist:** And the question Maxwell asks at 24 is: can we write that down in equations? His first paper, [[maxwell1855_lines_of_force]], builds an analogy with fluid mechanics. Imagine an incompressible fluid flowing along Faraday's lines. The velocity field of that fluid has certain divergence and curl properties — and those properties reproduce, exactly, the magnetostatic and electrostatic field equations as they were known in 1855. No new physics. But now there's a continuous field with a calculus around it.

**Historian:** Maxwell sits on the result for six years. Then in 1861 he publishes the strangest paper of his career: [[maxwell1861_physical_lines|*On Physical Lines of Force*]]. It proposes that space is filled with tiny mechanical *vortices* representing the magnetic field, separated by tiny mechanical *idle wheels* representing the electric current.

**Experimentalist:** Idle wheels. Like the things between gears that let two adjacent gears spin in the same direction.

**Historian:** Exactly. The whole construction is mechanical. Visualisable. And as Maxwell works through what it would take to make the model internally consistent — to satisfy conservation of charge, for instance — he is *forced* to add a term to Ampère's law that nobody had ever written down before: a term proportional to the time derivative of the electric displacement field $\partial_t \mathbf{D}$.

**Physicist:** The **displacement current**. This is the load-bearing addition. Ampère's law as it stood — curl of B equals current — is fine for steady currents, but it's incompatible with charge conservation when you have time-varying fields. Add Maxwell's displacement current and suddenly the four field equations close on themselves. They become *self-consistent*. And they admit transverse wave solutions.

**Experimentalist:** And the speed of the wave?

**Physicist:** Phase velocity comes out to $1/\sqrt{\mu_0 \varepsilon_0}$. Plug in the numbers. The product $\mu_0 \varepsilon_0$ has been measured — separately, electrostatically and magnetostatically, by [[weber_kohlrausch1856_elektrodynamische|Weber and Kohlrausch]] in 1856, using a current balance and a Leyden-jar capacitance comparison. The number that comes out: $3 \times 10^{10}$ centimeters per second. Within experimental error, that *is* the speed of light, which Fizeau and Foucault had measured separately on optical benches.

**Historian:** And in *On Physical Lines* Part III, 1862, Maxwell takes the hedge-and-step-back tone we'll see again three years later: "we can scarcely avoid the inference". Light is electromagnetic.

**Experimentalist:** And the [[maxwell1865_dynamical_theory|1865 paper]] drops the mechanical vortex model entirely.

**Historian:** Drops it cleanly. Twenty equations. Standalone. No idle wheels. Maxwell's later [[maxwell1873_treatise|*Treatise on Electricity and Magnetism*]] from 1873 systematises everything in two volumes, written in quaternion notation that Heaviside and Gibbs will replace with vector calculus in the 1880s. Most working physicists of the 1880s and 90s read the *Treatise*, not the 1865 paper.

`[cue: animation hist_maxwell_synthesis]`

**Physicist:** The animation we're cuing here walks through the assembly of the wave equation: write down Ampère's law with displacement current; take the curl; substitute Faraday's law; the result is $\Box \mathbf{E} = 0$ — a transverse wave equation in vacuum, propagating at exactly $c = 1/\sqrt{\mu_0 \varepsilon_0}$. Four steps. That's the whole unification.

**Historian:** And in 1888 — twenty-three years after Maxwell — Heinrich Hertz, at Karlsruhe, sets out to test the prediction. His apparatus: a transmitter that is essentially a spark gap on a dipole antenna, driven by a high-voltage induction coil. A receiver is a separate resonant loop with its own tiny adjustable spark gap. He places the receiver several meters from the transmitter. Each time the transmitter sparks, the receiver sparks too.

**Experimentalist:** Both sparks. Synchronous. And he can move the receiver around the room and the sparking depends on distance and orientation in exactly the way it should if the receiver were being driven by an electromagnetic wave from the transmitter.

**Historian:** [[hertz1888_funkenentladung]]. Then Hertz measures the wavelength by setting up a standing-wave pattern — reflects the radiation off a metal sheet and locates the nodes with the receiver. Wavelength times frequency gives propagation speed. Comes out, within experimental error, at $c$.

**Physicist:** And Hertz then shows the radiation reflects, refracts, polarises, absorbs — *behaves like light*. The classical synthesis is closed. Maxwell predicted it; Hertz saw it. The Continental and British traditions are reconciled in Maxwell's favour.

**Historian:** Now — and this is the dramatic part — one year *before* Hertz's experiments confirmed Maxwell, two experimentalists in Cleveland tried to push the framework further and broke it.

**Experimentalist:** Right. Albert Michelson and Edward Morley, at the Case School of Applied Science. They build an interferometer. Two perpendicular arms, about 11 meters each — they fold the path with mirrors to fit on a tabletop — mounted on a stone slab floating on mercury, so they can rotate the whole thing smoothly.

**Historian:** The reasoning: Maxwell's equations hold in some particular rest frame — the rest frame of the *aether*, the medium in which the electromagnetic waves propagate. Just like sound waves propagate in air at a definite speed relative to the air. The Earth orbits the Sun at about 30 kilometers per second. If the aether is at rest with respect to the Sun, then the Earth is moving through it at 30 km/s. And the speed of light, *measured from the moving Earth*, should depend on direction. Faster perpendicular to the motion, slower parallel.

**Physicist:** And the predicted fringe shift between two orthogonal arm orientations is roughly $L v^2 / \lambda c^2$. With $L = 11$ meters, $\lambda = 500$ nanometers, and $v = 30$ kilometers per second, that's about 0.4 fringes. Easily observable.

**Experimentalist:** And the observed upper bound: 0.01 fringes. Forty times smaller. They publish, [[michelson_morley1887_relative_motion]], and they're careful to say *they don't know what it means*. The result is incompatible with the Earth moving through a stationary aether at any plausible orbital velocity.

`[cue: animation hist_michelson_morley]`

**Historian:** Resolutions come quickly but uneasily. [[fitzgerald1889_aether|George FitzGerald, 1889]], proposes — in a single paragraph in *Science* — that matter must contract along the direction of motion by the factor $\sqrt{1 - v^2/c^2}$, which would cancel the predicted fringe shift exactly. He gives no mechanism. The proposal is *ad hoc*. But the algebra is right.

**Physicist:** [[lorentz1892_electron_theory|Hendrik Lorentz, 1892 and following]], arrives independently at the same factor via electron-theory motivations and works out the full *Lorentz transformations* — the coordinate change that leaves Maxwell's equations invariant between frames. Lorentz treats it as a mathematical device. Einstein in 1905 reinterprets it kinematically. That's next episode.

**Experimentalist:** And — preview — Gill's dual framework reaches the same null prediction without invoking length contraction at all. Let me hold on to that and come back to it in the proper-time interlude.

**Historian:** A few more developments to round out the era. In 1879, an American graduate student named [[hall1879_effect|Edwin Hall]], working under Henry Rowland at Johns Hopkins, discovers what we now call the Hall effect: a current through a thin gold leaf in a transverse magnetic field develops a transverse voltage perpendicular to both. It's the first major solid-state result with a clear electromagnetic mechanism — Lorentz force on the charge carriers, building up a static potential. The *sign* of the Hall voltage depends on whether the charge carriers are positive or negative — a fact that becomes diagnostic for band theory fifty years later. The solid-state thread sits mostly quiet in this chapter, but Hall's experiment opens the line that runs all the way to the transistor and the quantum Hall effect.

**Experimentalist:** And the last big event in the classical era — 1897 — is at the Cavendish in Cambridge. [[thomson1897_cathode_rays|J. J. Thomson]] resolves a fifteen-year German-versus-English argument about cathode rays — are they aether disturbances or particles? — by deflecting them in crossed electric and magnetic fields and measuring the ratio of charge to mass. The number comes out roughly 1,800 times larger than for hydrogen. Either the charge is enormous or the mass is tiny.

**Physicist:** Thomson concludes the mass is tiny. The rays are a stream of *corpuscles* — subatomic, common to all matter. He names them after no one, calling them simply "corpuscles", but the term *electron* — coined earlier by [[lorenz1867_identity_light_electricity|G. J. Stoney]] for the unit of charge — sticks within a few years. The electron exists. Lorentz's 1892 electron-theory programme has its substrate. Chapter 3, next episode, picks up from here.

## Physics deep dive

**Physicist:** Let me work through the two derivations the chapter's animations cover. First, the assembly of the wave equation. Start from Ampère's law with Maxwell's displacement-current correction:

$$\nabla \times \mathbf{B} = \mu_0 \mathbf{J} + \mu_0 \varepsilon_0 \, \partial_t \mathbf{E}.$$

In vacuum, $\mathbf{J} = 0$. Take the curl of both sides:

$$\nabla \times (\nabla \times \mathbf{B}) = \mu_0 \varepsilon_0 \, \partial_t (\nabla \times \mathbf{E}).$$

Use Faraday's law $\nabla \times \mathbf{E} = -\partial_t \mathbf{B}$ on the right side. And the vector identity on the left, with $\nabla \cdot \mathbf{B} = 0$:

$$-\nabla^2 \mathbf{B} = -\mu_0 \varepsilon_0 \, \partial_t^2 \mathbf{B}.$$

Rearranging:

$$\left( \nabla^2 - \mu_0 \varepsilon_0 \, \partial_t^2 \right) \mathbf{B} = 0.$$

That's a transverse wave equation with phase velocity $v = 1/\sqrt{\mu_0 \varepsilon_0}$. Plug in the numbers Weber and Kohlrausch had — that's $c$.

**Experimentalist:** And the displacement-current term is the load-bearing one. Without it, the right side of the original Ampère curl picks up a zero, the equation has no wave solution.

**Physicist:** Right. *Steady-state* magnetostatics doesn't need the displacement current. Time-varying electromagnetism *does*. It's also the term that makes electromagnetism consistent with charge conservation. Both are saying the same thing.

**Physicist:** Now Michelson–Morley. Standard derivation. Light travels at $c$ in the aether-rest frame. Lab frame moves through the aether at $v$. Along the arm parallel to the motion, light goes out at $c - v$ (relative to the lab mirror) and returns at $c + v$. Total round-trip time:

$$t_\parallel = \frac{L}{c - v} + \frac{L}{c + v} = \frac{2Lc}{c^2 - v^2} = \frac{2L}{c} \cdot \frac{1}{1 - v^2/c^2}.$$

Along the perpendicular arm, light travels at $c$ in the aether but the mirror has moved sideways while the light was in transit. The geometric path is longer by a factor of $1/\sqrt{1 - v^2/c^2}$. Round-trip time:

$$t_\perp = \frac{2L}{c} \cdot \frac{1}{\sqrt{1 - v^2/c^2}}.$$

The difference, $t_\parallel - t_\perp$, multiplied by the wave frequency, gives the predicted fringe shift. To leading order in $v^2/c^2$:

$$\Delta n \approx \frac{Lv^2}{\lambda c^2}.$$

That's 0.4 fringes for Michelson and Morley's parameters.

**Experimentalist:** And the observed shift is below 0.01 fringes. The aether is undetectable.

## Proper-time interlude

**Physicist:** Now the dual-framework reframing. In Tepper Gill's [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations|*Two Mathematically Equivalent Versions of Maxwell's Equations*]] — the load-bearing verification anchor for this whole chapter — the substitution $\partial_t = (b/c) \partial_\tau$ takes the standard Maxwell equations to a proper-time-equivalent form. The collaborative speed $b = \sqrt{c^2 + u^2}$ is a function of the source velocity $u$.

**Experimentalist:** Whose proper time?

**Physicist:** The source's. The source's clock. And here's the key for Michelson–Morley: the source — the lab apparatus, the mirrors, the photon paths — is *at rest in the lab frame*. So $u = 0$, $b = c$ exactly. The dual Maxwell equations collapse to the standard Maxwell equations identically. There is no asymmetry between the two interferometer arms. The predicted fringe shift is zero, derived without ever introducing length contraction.

**Experimentalist:** So the dual framework reaches the null result by a different route. What does that buy us?

**Physicist:** Two things. First, conceptually: the proper-time framework doesn't need length contraction as a separate kinematic postulate. It's not that lengths *don't* contract — it's that the velocity duality $\mathbf{w}/c = \mathbf{u}/b$, which is Eq. (1) of the Maxwell paper, already encodes whatever the relevant frame-dependent kinematics need to be. Length contraction in the standard Lorentz-transformation derivation is a *consequence* of the velocity-duality structure, not an additional postulate.

**Experimentalist:** And the second thing?

**Physicist:** Numerical agreement plus future testability. Both framings predict zero for Michelson–Morley. Both framings predict zero for the modern follow-up experiments — Brillet–Hall (1979), Müller-Peters–Chu (2003), and the like — at their stated precisions. Where the two diverge is in regimes where the source itself moves at a non-negligible fraction of $c$. For a GPS satellite at orbital velocity $u/c \sim 10^{-5}$, the collaborative speed is $b = c \sqrt{1 + 10^{-10}} \approx c (1 + 5 \times 10^{-11})$. The dual-framework correction is real but below current GPS clock precision. We work through that in [[07_PNT_GPS_SLR_QKD|Chapter 7]].

**Experimentalist:** So same null at experimental precision, different kinematics, with measurable divergence beyond a higher precision threshold than current GPS reaches.

**Physicist:** Exactly. The point of having the dual framework in our toolkit isn't to predict Michelson–Morley differently — Michelson–Morley is a *passed* test for both framings — but to set up the regime where the divergence is real. And, more importantly for this chapter, to provide an alternative *interpretation* of the null result that doesn't require the historical sequence of *ad hoc* length contraction (FitzGerald), *ad hoc* coordinate transformation (Lorentz), and *kinematic reinterpretation* (Einstein).

## Closing

**Historian:** So Chapter 2 in one episode. Maxwell from 1855 to 1873 — three foundational papers and a treatise. Hertz in 1888 confirming the waves exist. Michelson and Morley in 1887 — published a year earlier but interpreted later — breaking the aether picture. FitzGerald and Lorentz patching the breakage with length contraction. Hall opening the solid-state thread. Thomson identifying the electron and closing the classical era.

**Experimentalist:** And the load-bearing payoff: the verification campaign's anchor paper, Gill's *Two Mathematically Equivalent Versions of Maxwell's Equations*, lives in this chapter. Twenty-three of twenty-four equations verified by Wolfram in the existing campaign. The dual form of the equations and the dual derivation of Michelson–Morley both check out. There's one open finding — Eq. (24) has a missing $c$ factor and a missing potential term — flagged for author review in [[FINDINGS_for_author_review]] and addressed in [[05_QED_renormalization_solid_state_1948_1965|Chapter 5]] where the spin-1/2 fine-structure context makes the typo matter.

**Physicist:** And the framing-principle reminder: we're exploring how the same experimental record reads under a different mathematical convention. Michelson–Morley is a null result either way. Hertz sees electromagnetic waves either way. The dual framework reframes the *why* of the null — not the null itself.

**Historian:** Next episode: Planck, photoelectric, Bohr, de Broglie. The old quantum theory. Where the experimental record stops being something Maxwell's equations can carry alone. Thanks for listening.

`[cue: end card with bibliography wikilinks for show notes]`

## Show notes

Auto-generated from the chapter's bibliography by `lint_episode.py`; primary + retrospective sources cited above.
