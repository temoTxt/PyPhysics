---
chapter: 02
title: "Classical Synthesis 1860-1900"
era: "1860-1900"
threads: [electromagnetism, solid-state]
animations: [hist_maxwell_synthesis, hist_michelson_morley]
verification_anchors: ["Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations"]
status: draft
---

# Chapter 2 — Classical Synthesis (1860–1900)

> **Framing principle (load-bearing).** We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics. The dual framework reproduces every prediction of standard Maxwell + special relativity within current measurement precision. The contrast at the level of *interpretation* — most notably, the absorption of Michelson–Morley's null result into the proper-time formalism without invoking length contraction as a separate kinematic postulate — is the load-bearing payoff of this chapter.

## 1. Overview

The 1860s through 1890s are the era in which electromagnetism completes itself as a closed, mathematically self-consistent classical field theory and then runs into its own first inconsistency. The unification is Maxwell's: the [[maxwell1865_dynamical_theory]] gives twenty coupled differential equations from which all of [[01_early_electromagnetism_1800_1860|Chapter 1's]] discoveries fall out, plus the additional prediction that light *is* an electromagnetic wave propagating at $c = 1/\sqrt{\mu_0 \varepsilon_0}$. Twenty-three years later, [[hertz1888_funkenentladung|Hertz]] confirms experimentally that EM waves exist and behave like light. One year before Hertz's experiment, [[michelson_morley1887_relative_motion|Michelson and Morley]] try to measure the Earth's motion through the luminiferous aether — the postulated rest frame in which Maxwell's wave equation was supposed to hold — and find no motion at all. The aether's rest frame is undetectable. By 1897, [[thomson1897_cathode_rays|Thomson]] has identified the electron as the carrier of charge, providing the microscopic substrate Lorentz needs for the electron-theory programme. The chapter is anchored in the verification campaign's load-bearing Maxwell paper ([[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]) and is the first to exhibit a genuine — not just notational — difference between the standard and dual-theory framings of established physics. Next: [[03_old_quantum_theory_1900_1925]].

## 2. Historical narrative

### 2.1 1855–1865 — Maxwell formalises Faraday

Maxwell's path to the [[maxwell1865_dynamical_theory]] is staged across three papers spanning ten years.

**[[maxwell1855_lines_of_force]] — "On Faraday's Lines of Force".** A young Maxwell (24) reads Faraday's *Experimental Researches* and decides the field-line ontology can be made mathematical. He builds an analogy with incompressible-fluid streamlines: define a vector field whose flow lines coincide with Faraday's lines of force, and the divergence/curl operators reproduce Gauss's law and Ampère's law. No new physics yet — but the analytic language is now in place.

**[[maxwell1861_physical_lines]] — "On Physical Lines of Force".** Maxwell proposes a mechanical aether model: space is filled with vortices (representing magnetic field) separated by tiny "idle wheels" (representing electric current and displacement). Working through the internal consistency of this baroque construction, he discovers that Ampère's law $\nabla \times \mathbf{H} = \mathbf{J}$ must be modified — to be consistent with charge conservation $\partial_t \rho + \nabla \cdot \mathbf{J} = 0$ — by adding a term $\partial_t \mathbf{D}$. This is the **displacement current**. With the new term in place, the equations support transverse wave propagation at speed $1/\sqrt{\mu_0 \varepsilon_0}$, numerically very close to the speed of light measured by [[weber_kohlrausch1856_elektrodynamische|Weber and Kohlrausch]] in 1856. In Part III (1862), Maxwell concludes the wave *is* light.

**[[maxwell1865_dynamical_theory]] — "A Dynamical Theory of the Electromagnetic Field".** Maxwell drops the mechanical scaffolding and presents twenty coupled differential equations standalone. The equations govern E, B, current, charge, and the auxiliary potentials. Light is now a transverse EM wave; optics is a special case of electromagnetism. The [[maxwell1873_treatise|*Treatise*]] (1873) recasts these in (then-fashionable) quaternion notation; Heaviside and Gibbs in the 1880s give us the vector-calculus form physicists actually use.

Ludvig Lorenz in Copenhagen reaches the same conclusion — light is an EM disturbance — independently in [[lorenz1867_identity_light_electricity]], working in the framework of retarded action-at-a-distance and introducing what we now call the Lorenz gauge $\partial_\mu A^\mu = 0$.

### 2.2 1888 — Hertz: the waves exist

Heinrich Hertz, at the *Technische Hochschule Karlsruhe*, sets out to test Maxwell's prediction. His apparatus: a transmitter is a spark gap on a dipole antenna driven by a Ruhmkorff induction coil; the receiver is a separate resonant loop with its own micrometer-adjusted spark gap, placed several meters away. When the transmitter sparks, the receiver sparks too — sometimes only when correctly oriented, sometimes only at certain distances.

In [[hertz1888_funkenentladung]] and follow-on papers, Hertz measures the wavelength of his radiation by setting up standing waves (between the antenna and a reflecting metal sheet) and locating the nodes with the receiver. The wavelength, multiplied by the antenna's resonant frequency, gives a propagation speed within experimental error of $c$. Subsequent experiments show the waves can be reflected, refracted, polarised, and absorbed exactly like light. The British field-tradition and the Continental action-at-a-distance tradition reconcile in Hertz's favour of Maxwell.

### 2.3 1887 — Michelson and Morley's null result

One year *before* Hertz's spark detection — but with results that take a decade to interpret correctly — Albert Michelson and Edward Morley, at the Case School of Applied Science in Cleveland, attempt to measure the Earth's motion through the luminiferous aether.

The reasoning is straightforward in the late-19th-century picture. Maxwell's equations are presumed to hold in a particular rest frame — the rest frame of the aether — exactly the way the equations of acoustics hold in the rest frame of the air. If the Earth moves through that aether at orbital velocity $v \approx 30$ km/s, then light propagating parallel to that motion travels at $c - v$ relative to a lab-frame mirror, while light propagating perpendicular travels at $c$ in the aether but takes a slightly longer path in the lab. A two-arm interferometer with arms of length $L$, oriented one parallel and one perpendicular to the motion, should show a fringe shift of approximately

$$\Delta n \approx \frac{L v^2}{\lambda c^2}$$

between the two configurations. With $L \approx 11$ m (multiply folded by mirrors) and $\lambda \approx 500$ nm, the expected shift is roughly 0.4 fringes.

The measured upper bound: 0.01 fringes — *forty times smaller* than predicted. The result, [[michelson_morley1887_relative_motion]], is incompatible with the Earth moving through a stationary aether at any speed close to $v_\oplus$.

Two responses arrive within two years. [[fitzgerald1889_aether|FitzGerald (1889)]] proposes, in a one-paragraph letter to *Science*, that matter contracts along the direction of motion by a factor that just cancels the predicted shift — *ad hoc*, no mechanism, but algebraically correct. [[lorentz1892_electron_theory|Lorentz (1892, 1895, 1899, 1904)]] arrives independently at the same factor $\sqrt{1 - v^2/c^2}$ via a more detailed electron-theory motivation, and writes down the *Lorentz transformations* of space and time as the coordinate change that leaves Maxwell's equations invariant. Lorentz treats this as a mathematical device, not as a kinematic principle. Einstein (1905) will reinterpret it kinematically as Special Relativity — that's [[03_old_quantum_theory_1900_1925|Chapter 3]].

The **alternative resolution** that this campaign cares about: in Gill's dual-theory formulation, the null result is absorbed into the relation between observer time $t$ and proper time $\tau$ via the collaborative speed $b = \sqrt{c^2 + u^2}$ (where $u$ is the source velocity). Length contraction is not required as a separate kinematic postulate. We work through this in §3 below.

### 2.4 1879 — Hall, and the solid-state thread

Edwin Hall, a graduate student at Johns Hopkins, runs a current through a thin gold leaf in a transverse magnetic field and finds that a *transverse* potential difference develops perpendicular to both the current and the field — [[hall1879_effect]]. The **Hall effect** is conceptually simple (Lorentz force on the charge carriers, deflected sideways, building up a static potential to balance the force) but its sign in metals depends on whether the charge carriers are positive or negative — a discriminator that becomes diagnostic for band theory in the 1930s and beyond. The first major solid-state result tying directly to electromagnetism; the opening of the long line that runs through Bloch (1928), Wilson (1931), the transistor (1947), BCS (1957), and the quantum Hall effect (1980).

The solid-state thread sits mostly dormant in this chapter — quasi-static conduction in metals is uncomplicated — but is included because the *experimental capability* (precision galvanometry, thin-film deposition, controlled magnetic environments) developing in this era is what makes the [[05_QED_renormalization_solid_state_1948_1965|Chapter 5]] solid-state revolution possible.

### 2.5 1897 — Thomson identifies the electron

J. J. Thomson, at the Cavendish Laboratory in Cambridge, has been studying cathode rays — the glowing streams produced when a high voltage is applied across the electrodes of an evacuated tube. The German experimental tradition (Hertz himself, in 1883) had concluded the rays were aether disturbances; the English tradition (Crookes) had concluded they were charged particles. Thomson settles it by deflecting cathode rays in crossed electric and magnetic fields. Measuring the trajectory gives $e/m$. The value is roughly 1,800 times larger than that for the lightest known ion (hydrogen) — meaning either the carrier's charge is enormous (implausible) or its mass is tiny.

[[thomson1897_cathode_rays|Thomson's conclusion]] is that cathode rays are a stream of *corpuscles* — subatomic particles, common to all matter, with mass much smaller than the atom. The discovery confirms Lorentz's 1892 electron-theory hypothesis experimentally and provides the microscopic substrate every subsequent theory will rest on. The closing event of the classical era.

## 3. Proper-time commentary

This is the first chapter with `#verified` claims drawn from the verification campaign. The Maxwell equations of §2.1 are the anchor.

### 3.1 What's directly verified

**Maxwell's equations in dual form.** `#verified` from [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (3) — Standard (observer-time) Maxwell's equations]] and [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (3′) — Proper-time-equivalent Maxwell's equations]]. The standard Maxwell equations

$$\nabla \cdot \mathbf{E} = \rho/\varepsilon_0,\quad \nabla \cdot \mathbf{B} = 0,\quad \nabla \times \mathbf{E} = -\partial_t \mathbf{B},\quad \nabla \times \mathbf{B} = \mu_0 \mathbf{J} + \mu_0 \varepsilon_0 \partial_t \mathbf{E}$$

are *mathematically equivalent*, under the substitution $\partial_t = (b/c)\,\partial_\tau$ with $b = \sqrt{c^2 + u^2}$ and the [Eq. (1)] velocity duality $\mathbf{w}/c = \mathbf{u}/b$, to the proper-time form

$$\nabla \cdot \mathbf{E} = \rho/\varepsilon_0,\quad \nabla \cdot \mathbf{B} = 0,\quad \nabla \times \mathbf{E} = -(b/c)\partial_\tau \mathbf{B},\quad \nabla \times \mathbf{B} = \mu_0 \mathbf{J} + \mu_0\varepsilon_0 (b/c)\partial_\tau \mathbf{E}.$$

Both systems generate identical electromagnetic field configurations from identical source distributions. The duality is *not* an approximation; it is exact. The dual form is more convenient when the source is accelerating (the $\partial_\tau$ operator is the source's own clock), the standard form is more convenient when the field is being measured by a stationary lab apparatus. Wolfram MCP verification recorded in the linked anchor.

**Wave equation invariance and effective photon mass.** `#verified` from [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (4) — Dual wave equations with dissipative term]] and [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (6) — Effective photon mass]]. The vector-potential wave equation $\Box \mathbf{A} = \mu_0 \mathbf{J}$ in the dual framework picks up a dissipative term $-(\mathbf{u} \cdot \mathbf{a})/b^4 \cdot \partial_\tau \mathbf{E}$ proportional to source acceleration. Under the Liouville substitution $\psi = (b/c)^{1/2} \Psi_{\rm new}$ this term reorganises into an effective photon mass $\mu = $ (a function of $\mathbf{u}, \mathbf{a}$) — *non-zero only when the source accelerates*, vanishing in the inertial limit. For Hertz's spark-gap antenna (oscillating dipole moment), the source acceleration is non-zero — but at his radiation wavelengths and source velocities, the effective $\mu$ corresponds to a Compton wavelength orders of magnitude longer than any apparatus dimension and produces no measurable dispersion or attenuation.

### 3.2 What's mechanically inferred

**Michelson–Morley reframed via $b/c$.** `#inferred` from [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (10) — Proper-time boosts of position, velocity, acceleration]]. The standard derivation of the expected fringe shift assumes light propagation in the aether-rest frame at $c$, with the Earth moving through that aether at $v$, producing a frame-dependent "ether wind" $c \pm v$. The dual framework rejects this kinematic picture: there is no aether-rest frame, and the relevant velocity for a source at rest in the *lab* frame is $u = 0$, giving $b = c$ identically. The lab apparatus and the photons share a single proper-time parameterisation. The two interferometer arms experience the *same* propagation speed, and the predicted fringe shift is exactly zero. **Experimental distinguishability:** the standard SR (length-contraction-mediated) and dual-framework (no length contraction needed) derivations give the same null result for this experiment. They diverge only in regimes where the source velocity is itself relativistic — see [[07_PNT_GPS_SLR_QKD]] for the GPS satellite case at $u/c \sim 10^{-5}$ where the divergence is real but $10^{-10}$ in magnitude, well below current GPS clock precision.

**Lorentz transformations as a partial subset of $b/c$ rescaling.** `#inferred` from [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (9) — Mean-value form for $t$ in terms of $\tau$]]. The Lorentz time transformation $t' = \gamma(t - vx/c^2)$ corresponds, in the dual framework, to the $(c/b)$-weighted integral $t = (1/c) \int b(s) \, ds$ between the proper-time of an observer and the proper-time of a moving source. The Lorentz boost of *length* does not have a single dual-framework analog: in Gill's reformulation, lengths transform via the velocity Eq. (1) duality $\mathbf{w}/c = \mathbf{u}/b$ rather than via a separate kinematic postulate. The empirical predictions of standard SR and dual relativity coincide for any current experiment; the conceptual difference (whether length contraction is a fundamental kinematic effect or a derived consequence of velocity duality) sits in the interpretation, not the predictions.

### 3.3 What Gill is silent on

- **Hall effect, thermoelectricity, magnetoresistance** in metals. `#gill-silent` — quasi-static solid-state physics, dominated by carrier scattering and band structure, untouched by the proper-time formalism.

- **Thomson's $e/m$ measurement.** `#gill-silent` — lab-frame Lorentz force on slow electrons ($v \approx 0.1\,c$ in his apparatus, where $b \approx c \cdot 1.005$; the correction is well within his ~5% measurement uncertainty).

## 4. Key derivations worth animating

| Manim scene | Status | What it shows |
|---|---|---|
| [`hist_maxwell_synthesis.py`](../Animations/manim_scenes/hist_maxwell_synthesis.py) | rendered | Ampère's law + displacement current → wave equation $\Box \mathbf{E} = 0$ → propagation at $c = 1/\sqrt{\mu_0\varepsilon_0}$ → identification with light. The displacement-current step is the load-bearing one — without it, no wave; with it, optics is electromagnetism. |
| [`hist_michelson_morley.py`](../Animations/manim_scenes/hist_michelson_morley.py) | rendered | The interferometer geometry; the predicted fringe shift under the aether-wind model; the observed null result; and the dual-framework reframing (lab-rest source ⇒ $u=0$ ⇒ $b=c$ ⇒ no asymmetry). Closes on the fact that the two framings reproduce the same null but with different kinematics. |

Both scenes are non-trivial Manim work — Maxwell synthesis chains four equations through Transforms; Michelson–Morley overlays the apparatus schematic with the predicted-vs-observed fringe traces.

## 5. Primary sources cited

- [[maxwell1855_lines_of_force]] — first analytic translation of Faraday.
- [[maxwell1861_physical_lines]] — discovery of displacement current.
- [[maxwell1865_dynamical_theory]] — canonical equations + identification of light with EM. **Verification anchor for this chapter.**
- [[maxwell1873_treatise]] — systematic 2-volume treatise.
- [[weber_kohlrausch1856_elektrodynamische]] — ratio-of-units measurement that gives `c`.
- [[helmholtz1847_erhaltung_kraft]] — energy conservation as universal principle.
- [[lorenz1867_identity_light_electricity]] — independent derivation; Lorenz gauge.
- [[hertz1888_funkenentladung]] — experimental confirmation of EM waves.
- [[michelson_morley1887_relative_motion]] — null result; foundational for §3.
- [[fitzgerald1889_aether]] — ad hoc length-contraction proposal.
- [[lorentz1892_electron_theory]] — electron theory + Lorentz transformations.
- [[hall1879_effect]] — Hall effect; opens solid-state thread.
- [[thomson1897_cathode_rays]] — electron discovery; closes classical era.

## 6. Retrospective reviews drawn on

- [[whittaker1910_aether_electricity]] — classic narrative continuing from Ch 1; particularly load-bearing for the Lorenz/Maxwell priority point.
- [[buchwald1985_maxwell_microphysics]] — modern academic history of the post-Maxwell generation (Heaviside, Lodge, FitzGerald, Larmor).
- [[siegel1991_innovation_maxwell]] — close reading of the 1855–1865 Maxwell papers, defending the mechanical-aether model as a productive tool.
- [[darrigol2000_electrodynamics_ampere_einstein]] — modern synthesis spanning Ch 1 and Ch 2.
- [[jackson1998_classical_electrodynamics]] — graduate textbook; the standard-equation reference cited throughout the verification campaign.

## 7. Cross-references

- Backward: [[01_early_electromagnetism_1800_1860]].
- Forward: [[03_old_quantum_theory_1900_1925]] — Planck, photoelectric, Bohr.
- Verification anchor: [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]].
- Findings: [[FINDINGS_for_author_review]] — Finding 1 (Maxwell Eq. (24) typos) is in this paper.
