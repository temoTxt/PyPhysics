---
chapter: 01
title: "Early Electromagnetism 1800-1860"
era: "1800-1860"
threads: [electromagnetism]
animations: [hist_faraday_induction]
verification_anchors: []
status: draft
---

# Chapter 1 — Early Electromagnetism (1800–1860)

> **Framing principle (load-bearing).** We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics. Every "Gill would predict X" claim below is qualified by what experimental precision would distinguish it from the standard prediction; the dual-theory framework reproduces all experimentally confirmed predictions of standard SR + QM + QED within current measurement precision.

## 1. Overview

The half-century from Volta's pile (1800) to Kirchhoff's telegrapher equation (1857) is the era in which electricity and magnetism, treated as separate phenomena since antiquity, are revealed to be two faces of a single subject. The arc begins with a brass-and-zinc stack of disks producing the first sustained electric current and ends with the empirical fact — half-noticed in 1857, fully assimilated by Maxwell a few years later — that the velocity of electromagnetic disturbance on a wire equals the speed of light. None of this era's physics is *relativistic*. Currents are quasi-static; sources are at rest in the lab frame; observed velocities `u ≪ c`. Under such conditions Gill's collaborative speed `b = √(c² + u²) ≈ c (1 + u²/2c²)` reduces to `c` to a precision orders of magnitude finer than any 19th-century experiment could detect. The chapter is therefore a **null comparison**: every load-bearing 1800–1860 result holds in standard and dual Maxwell with no distinguishable difference. The point of writing it is to establish the experimental record — the foundation on which both framings are erected — and to set up the genuine divergence in [[02_classical_synthesis_1860_1900]] when Maxwell's displacement current and the Michelson–Morley apparatus arrive together.

## 2. Historical narrative

### 2.1 1800 — Volta and the operational birth of current electricity

In a long letter to Sir Joseph Banks dated 20 March 1800, Alessandro Volta describes a column of alternating zinc and silver disks separated by brine-soaked cardboard. Touched between the top and bottom of the stack, the column delivers a continuous, mild electric shock — qualitatively different from the brief discharge of a Leyden jar. The "Voltaic pile" ([[volta1800_electricity_contact]]) is the first reliable source of sustained current. Within months, replicas are running across Europe.

Two distinct programmes follow from the pile. *Electrochemistry* — Nicholson and Carlisle decompose water by electrolysis in May 1800; Davy isolates sodium and potassium in 1807; Faraday formalises the laws of electrolysis in the 1830s — proceeds on its own track, mostly outside this chapter's scope. *Electromagnetism* — the subject of the chapter — begins as soon as anyone notices that a continuous current does something a static charge cannot.

### 2.2 1820 — Ørsted, Ampère, Biot–Savart: the link between electricity and magnetism

For two decades, nobody noticed. The pile was widely available but the connection to magnetism remained invisible until Hans Christian Ørsted's lecture demonstration in Copenhagen in spring 1820: a compass needle held parallel to a current-carrying wire swings *perpendicular* to the wire, not along it. Ørsted publishes the four-page Latin pamphlet *Experimenta circa effectum conflictus electrici in acum magneticam* ([[orsted1820_experimenta_circa]]) on 21 July 1820 and distributes it to learned societies across Europe.

The reaction is immediate. By the time the pamphlet reaches Paris in September, André-Marie Ampère is already at work. Within three weeks he has demonstrated that two parallel current-carrying wires attract when their currents flow in the same direction and repel when antiparallel. Reading on 18 and 25 September and 2 October 1820 to the Académie des Sciences, Ampère ([[ampere1820_action_mutuelle]]) gives the force law $dF \propto I_1 I_2 \, d\ell_1 \, d\ell_2 / r^2$ between current elements (the final form refined in his 1826 monograph). He coins the term *electrodynamics* and proposes that magnetism *is* the manifestation of microscopic, persistent molecular currents — the "Ampèrian current" model, the first attempt to reduce one phenomenon to the other.

Simultaneously, Jean-Baptiste Biot and Félix Savart ([[biot_savart1820_magnetisme_pile]]) measure the *field* of a long straight current-carrying wire and find that it scales as $1/r$. Refined with Laplace's analytical help in 1821, the result becomes the Biot–Savart law

$$d\mathbf{B} = \frac{\mu_0}{4\pi} \frac{I\,d\boldsymbol{\ell} \times \hat{\mathbf{r}}}{r^2}$$

— the field produced by a current element, paired with Ampère's force-on-a-current. The pair are the magnetostatic equivalents of Coulomb's law and the Lorentz force, written down forty years before they're unified by Maxwell.

### 2.3 1827 — Ohm's law and the legitimisation of quantitative electrical measurement

Georg Simon Ohm's monograph *Die galvanische Kette, mathematisch bearbeitet* ([[ohm1827_galvanische_kette]]) proposes — and verifies experimentally with thermocouple sources of known electromotive force — that for a metallic conductor

$$V = IR$$

with the resistance $R$ a property of the wire alone. The analogy is consciously Fourier-inspired: heat flux scales with the temperature gradient, current flux scales with the potential gradient. The early reception is hostile (Ohm is forced to resign his professorship for the "speculative" nature of the work), but by the 1840s the law has become indispensable to the burgeoning telegraph industry, and the Royal Society awards Ohm the Copley Medal in 1841.

### 2.4 1831–1834 — Faraday induction, Henry, and the Lenz sign

The decisive experimental discovery of the era. On 29 August 1831, Michael Faraday — working at the Royal Institution with two coils wound on opposite sides of a soft-iron ring — observes that closing the primary circuit produces a *momentary* deflection in a galvanometer attached to the secondary. Steady currents do nothing; only *change* in the primary current induces a current in the secondary. The phenomenon is the inverse of Ørsted's: *time-varying magnetic flux* generates an electric current. Published in November 1831 and printed in *Philosophical Transactions* the following March, the work ([[faraday1832_experimental_researches_i]]) is the first of Faraday's *Experimental Researches in Electricity*, a series that will run thirty volumes over the next thirty years.

Section 90 of the same paper describes the "Faraday disc": a copper disc rotating between the poles of a permanent magnet, with brushes at the axis and the rim. The disc generates a DC voltage proportional to its rotation rate — the first electromagnetic generator and a direct demonstration that *mechanical motion through a field* induces electromotive force.

The mathematical statement we now write as $\mathcal{E} = -d\Phi/dt$ is implicit in Faraday's prose but not formalised. Maxwell will derive it cleanly in *On Faraday's Lines of Force* (1855); we cover that step in [[02_classical_synthesis_1860_1900]].

Two near-contemporaneous results clean up the picture. **Joseph Henry** ([[henry1832_currents_and_sparks]]), working independently at the Albany Academy, discovers induction in summer 1831 — slightly before Faraday's published date — but doesn't publish until 1832 because of teaching obligations. Henry generously acknowledges Faraday's priority. His distinctive contribution is *self-inductance*: the spark and "kick" observed when breaking the circuit of a long coil. The SI unit of inductance is named for him.

**Heinrich Lenz** ([[lenz1834_richtung_strome]]), in St Petersburg, fixes the *direction* of the induced current: it is always such that its own magnetic field *opposes* the change that induced it. This is the negative sign in $\mathcal{E} = -d\Phi/dt$. Conceptually decisive: with the wrong sign, energy conservation fails. Lenz's law is effectively a statement of conservation of energy applied to electromagnetism, anticipating Helmholtz's general 1847 *Über die Erhaltung der Kraft* by thirteen years.

### 2.5 1832–1846 — Gauss, Weber, and absolute units

Side-by-side with the qualitative discoveries, a quantitative-metrological programme is building. Carl Friedrich Gauss's *Intensitas vis magneticae terrestris ad mensuram absolutam revocata* ([[gauss1832_intensitas]]) introduces *absolute* units for magnetic field strength: combining a swinging-magnetometer timing measurement with a coil-deflection experiment to disentangle the magnetic dipole moment from the field strength, Gauss expresses the Earth's magnetic field in cm, g, s alone. The Magnetic Union (Magnetischer Verein), founded jointly with Wilhelm Weber in 1834, coordinates simultaneous magnetic observations across European observatories using these absolute units. The methodological move — express electromagnetic quantities in pure mechanical units — is the conceptual ancestor of the modern SI dimensional system, and an essential step toward making electromagnetism *quantitative* rather than instrument-specific.

Weber's 1846 *Elektrodynamische Maassbestimmungen* extends this programme into a *theory* of electrodynamic forces between moving charges — a force law depending not only on distance but on relative velocity and acceleration. Weber's theory is action-at-a-distance, in the Continental tradition, and competes through the 1850s and 1860s with the British field theory of Faraday and Thomson; we return to that competition in [[02_classical_synthesis_1860_1900]].

### 2.6 1857 — Kirchhoff and the first appearance of `c`

Gustav Kirchhoff's analysis of telegraph-cable signal propagation ([[kirchhoff1857_bewegung_elektrizitat]]) uses the empirical ratio of electrostatic to electromagnetic units (measured by Weber and Kohlrausch in 1856) to deduce that an electric disturbance propagates along a wire at very nearly $3 \times 10^{10}$ cm/s — the speed of light in vacuum. Kirchhoff doesn't draw the optical conclusion; he is solving a telegraphy problem, and the coincidence is left as a striking but unexplained empirical fact.

Eight years later, Maxwell will pick up that fact and amplify it into the unification of electromagnetism with optics. That work belongs to [[02_classical_synthesis_1860_1900]]. But the empirical anchor is here, in 1857.

## 3. Proper-time commentary

### 3.1 What's directly verified

Nothing in this era's physics has been the subject of a Wolfram MCP verification under the Gill corpus. The verification campaign ([epic #3](https://github.com/temoTxt/PyPhysics/issues/3)) targets the proper-time substitutions in Maxwell's equations and downstream, none of which Gill applies to 1800–1860 sources. `#verified` claims are absent in this chapter by construction.

### 3.2 What's mechanically inferred

- **Faraday induction in the dual framework.** `#inferred` from [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (3′) — Proper-time-equivalent Maxwell's equations]]. The differential form $\mathcal{E} = -d\Phi/dt$ rewrites under the proper-time substitution $\partial_t = (b/c) \partial_\tau$ as $\mathcal{E} = -(c/b) d\Phi/d\tau$. **Experimental distinguishability:** at Faraday's apparatus velocities (rotor tip speed ≲ 10 m/s in the 1831 disc; primary-circuit make/break transient ≲ 10⁻³ s⁻¹), the correction factor $(c/b - 1) \sim u²/2c²$ is below 10⁻¹⁵, six orders of magnitude finer than mid-19th-century galvanometer precision. The two framings give numerically identical induced EMFs to any experiment Faraday or Henry could perform.

- **Telegrapher equation in the dual framework.** `#inferred` from [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (4) — Dual wave equations with dissipative term]]. Kirchhoff's 1857 derivation of $\partial_t^2 \psi - v_K^2 \partial_x^2 \psi = (\text{loss terms})$ on a transmission line carries forward under $t \to \tau$ into a structurally identical equation in proper time with $v_K \to v_K (c/b)$. **Experimental distinguishability:** telegraph cable signal-propagation velocities are ≲ 0.5 c (and that's an upper bound for late-19th-century low-loss lines; 1857 signals were much slower); the predicted dual-theory correction is below the noise floor of any Victorian galvanometer.

### 3.3 What Gill is silent on

- **Voltaic chemistry, Ohmic conduction in metals, magnetostatics of steady currents.** `#gill-silent`. The dual framework's distinctions sit downstream of Maxwell's equations in their dynamical (time-varying) form; for steady-state sources at rest in the lab, the two framings are notationally identical.

- **Ampèrian molecular currents as a model of magnetism.** `#gill-silent`. Ampère's model was overtaken by the electron theory (Lorentz, Zeeman) and quantum spin (Pauli, Dirac) long before the proper-time framework existed; Gill has no published commentary.

## 4. Key derivations worth animating

| Manim scene | Status | What it shows |
|---|---|---|
| [`hist_faraday_induction.py`](../Animations/manim_scenes/hist_faraday_induction.py) | rendered | Field lines through a moving conducting loop; induced EMF emerges from changing flux; Lenz-sign correctness visualised by counter-flux of the induced current. |

The Manim scene is pedagogical only — there is no proper-time content to animate in this era. The standard demonstration of $\mathcal{E} = -d\Phi/dt$ via a loop translating through a non-uniform field is what 1832 Faraday could *not* see analytically, and it sets up the formalism for the Maxwell-paper scenes in [[02_classical_synthesis_1860_1900]].

## 5. Primary sources cited

- [[volta1800_electricity_contact]] — the pile; first sustained DC source.
- [[orsted1820_experimenta_circa]] — current deflects compass needle.
- [[ampere1820_action_mutuelle]] — force law between parallel currents; the term *electrodynamics*.
- [[biot_savart1820_magnetisme_pile]] — field of a long straight current; the `1/r` law.
- [[ohm1827_galvanische_kette]] — `V = IR`.
- [[faraday1832_experimental_researches_i]] — electromagnetic induction; the rotating disc.
- [[henry1832_currents_and_sparks]] — independent discovery of induction; self-inductance.
- [[lenz1834_richtung_strome]] — sign convention; energy conservation.
- [[gauss1832_intensitas]] — absolute electromagnetic units; methodological precursor of SI.
- [[kirchhoff1857_bewegung_elektrizitat]] — first appearance of `c` in an EM context.

## 6. Retrospective reviews drawn on

- [[whittaker1910_aether_electricity]] — classic narrative; load-bearing for the 1820 Paris reception of Ørsted and the Faraday vs. Continental-tradition contrast.
- [[williams1965_faraday_biography]] — Faraday's lab notebooks and the "lines of force" ontology.
- [[darrigol2000_electrodynamics_ampere_einstein]] — modern academic synthesis of the two competing traditions.

## 7. Cross-references

- Forward: [[02_classical_synthesis_1860_1900]] — Maxwell, Hertz, Michelson–Morley.
- Verification anchors: [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] (forward reference; the dual-framework anchor for everything 1860 onward).
- Findings: none in this era (the corpus's three flagged findings — see [[FINDINGS_for_author_review]] — all sit in 1860s and later papers).
