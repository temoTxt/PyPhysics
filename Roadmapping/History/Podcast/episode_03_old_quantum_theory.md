---
episode: 03
title: "Old Quantum Theory 1900-1925"
era: "1900-1925"
chapter: 03_old_quantum_theory_1900_1925
speakers: [Historian, Physicist, Experimentalist]
target_runtime_min: 22
animations_cued:
  - hist_bohr_proper_time
  - hist_compton_null
status: draft
---

# Episode 3 — Old Quantum Theory (1900–1925)

> Companion dialogue script for [[03_old_quantum_theory_1900_1925]].

## Cold open

**Historian:** December 14th, 1900. The Berlin Physical Society. Max Planck, 42 years old, gives a one-hour seminar. He has been working for six months on a single problem — the spectrum of light emitted by a *blackbody*, an idealised cavity that absorbs all radiation falling on it and re-emits it in thermal equilibrium. The experimental curves, measured at the Reichsanstalt by Lummer and Pringsheim, are extraordinarily precise. No classical theory fits them.

**Experimentalist:** The Rayleigh–Jeans formula works at low frequencies and diverges to infinity at high. The Wien formula works at high frequencies and fails at low. Neither matches across the whole spectrum.

**Historian:** Planck has a guess. He assumes — provisionally, with no physical motivation he is willing to defend in public — that the energy of each oscillator in the cavity walls is restricted to integer multiples of $h\nu$. With that single combinatorial assumption, the spectrum comes out right. Across the whole frequency range. To experimental precision.

**Physicist:** And Planck does *not* claim he's discovered something physical. He calls $h$ a "fortunate guess". He treats the quantisation as a mathematical artefact useful for the combinatorial argument. The physical interpretation — that energy is *actually* quantised, not just calculationally — has to wait for Einstein in 1905.

**Historian:** Five years later. Einstein's *annus mirabilis*. Four foundational papers in twelve months, two of which dominate the rest of this episode: the *light quantum* paper extending Planck's quantisation to light itself, and *special relativity* settling kinematics in moving frames once and for all. Plus the Mercury perihelion in 1915 and the GR foundational paper in 1916.

**Experimentalist:** And Bohr, in 1913, takes the quantum hypothesis and runs it through Rutherford's nuclear atom — quantised orbits, photon transitions between them, the Rydberg formula falling out exactly. Plus Sommerfeld in 1916 with the fine structure and the first appearance of the fine-structure constant. Plus Compton in 1923 confirming photons carry momentum, plus de Broglie in 1924 saying electrons carry wavelength.

**Physicist:** And the reminder I always lead with: this entire series is exploring how the same experimental record reads under a different mathematical convention — Tepper Gill's dual framework. For this episode specifically, the load-bearing connection is the Bohr-atom rederivation using the dual Hamiltonian $K = H^2/(2mc^2) + mc^2/2$ from Gill's *Foundations II* paper. Same Bohr energy levels at hydrogen velocities — agreement to about one part in $10^9$. No experimental contradiction. The conceptual win is that the *positive-definite* dual Hamiltonian gives the right answer without ever splitting into kinetic plus potential pieces — which dissolves several mathematical irritations that the old quantum theory inherited from classical mechanics and that proper quantum field theory inherits in turn.

**Historian:** Let's begin in Berlin in December 1900.

## Historical sweep

**Historian:** Planck's December 1900 paper, [[planck1900_blackbody]]. The setup: the energy density inside a cavity at temperature $T$, per unit frequency. The classical derivation says $u(\nu,T)$ is proportional to $\nu^2 k_B T$ — finite at any single frequency, but integrated over all frequencies the energy is infinite. That's the *ultraviolet catastrophe*.

**Experimentalist:** And it's clearly not what's measured. Lummer and Pringsheim's blackbody data, made with platinum cavities and bolometric detectors, gives a peaked spectrum. The peak frequency scales with temperature in the way Wien predicted. The low-frequency tail matches Rayleigh–Jeans. The high-frequency tail matches Wien. Nothing matches the whole curve.

**Historian:** Planck's quantisation. Energy per oscillator restricted to $E_n = nh\nu$. Plug into the Boltzmann combinatorial formula, and the average energy per oscillator becomes $\langle E \rangle = h\nu / (e^{h\nu/k_B T} - 1)$. Multiply by the density of cavity modes. Out comes

$$u(\nu, T) = \frac{8\pi h \nu^3}{c^3} \cdot \frac{1}{e^{h\nu/k_B T} - 1}.$$

Matches the experimental curve across the whole spectrum. Reduces to Rayleigh–Jeans at low frequency, to Wien at high.

**Physicist:** And Planck extracts $h \approx 6.55 \times 10^{-34}$ joule-seconds. Modern value is $6.626$. He gets it right to 1%, from blackbody data, in 1900.

**Historian:** Now five years jump. Einstein 1905. Four papers in the same year. The Brownian motion paper — [[einstein1905_brownian]] — predicts the random-walk dynamics of suspended particles, with mean-square displacement scaling linearly with time. Perrin verifies it in 1908. Decisive evidence for the atomic theory of matter against the lingering energeticist school.

**Experimentalist:** Then the photoelectric paper — [[einstein1905_photoelectric]]. Einstein argues from a statistical analogy with the Wien blackbody spectrum that light *itself* is composed of independent quanta $h\nu$, not just absorbed and emitted in quanta. The photoelectric effect — electrons ejected from a metal surface when illuminated — should have ejection energy $E_k = h\nu - W$, with $W$ the metal's work function. *Independent of intensity*. Linear in frequency. Threshold at $W$.

**Historian:** And the experimental verification takes ten years. [[millikan1916_photoelectric_verification|Millikan]], who is initially skeptical of the light-quantum hypothesis, spends a decade trying to disprove it and ends up confirming Einstein to 1% accuracy. Einstein's 1921 Nobel cites this paper, vindicated by Millikan.

**Physicist:** Then special relativity, [[einstein1905_specrel]]. Postulates: constancy of $c$, equivalence of inertial frames. Derives the Lorentz transformations as kinematic, not as a calculational device the way [[lorentz1892_electron_theory|Lorentz 1892]] had used them. Length contraction, time dilation, relativistic velocity addition, all follow.

**Experimentalist:** And [[einstein1905_emc2|the three-page sequel]] derives $E = mc^2$ from a thought experiment about an emitter radiating equal photons in opposite directions in two reference frames.

**Historian:** Three years later, 1908, Hermann Minkowski's lecture in Cologne — [[minkowski1908_raum_zeit]] — recasts everything in 4-dimensional geometric language. Spacetime as a pseudo-Riemannian manifold. Worldlines as 4-curves. The invariant interval $ds^2 = c^2 dt^2 - dx^2 - dy^2 - dz^2$. "Henceforth space by itself, and time by itself, are doomed to fade away into mere shadows."

**Physicist:** And here's a forward connection to the dual framework. Gill's [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL|*Foundations I: Mathematical*]] paper proves what they call the **Minkowski Incompatible Theorem** — Theorem 1.1, verified by Wolfram in the verification campaign. It shows that the standard Minkowski 4-velocity $u^\mu = dx^\mu/d\tau$ does not close on itself under Hamiltonian flow without an additional structural assumption. Gill's positive-definite dual Hamiltonian — the one we'll see used in the Bohr-atom rederivation — sidesteps the issue by squaring.

**Experimentalist:** Quick reality check. Minkowski geometry holds in standard SR. Holds in the dual framework. The two framings inherit the same Minkowski spacetime. What changes is how the *Hamiltonian dynamics* on that spacetime is written down.

**Physicist:** Right. Both framings have Minkowski. They differ in how they parameterise dynamics on it.

**Historian:** OK — back to 1911. Ernest Rutherford, in Manchester, interpreting the gold-foil scattering experiments of his colleagues Geiger and Marsden. A beam of alpha particles directed at thin gold foil. Most pass through almost undeflected. But a small fraction — about one in 8,000 — scatter back almost the way they came. As Rutherford put it: "almost as if you fired a 15-inch shell at a piece of tissue paper and it came back at you."

**Experimentalist:** Incompatible with Thomson's *plum pudding* model — a diffuse positive cloud. But perfectly consistent with a small, dense, positively charged nucleus surrounded by orbiting electrons. The nuclear atom. [[rutherford1911_alpha_scattering]].

**Historian:** But there's a problem. Classical electrodynamics demands that an orbiting — therefore accelerating — electron radiates. Energy loss should make every atom collapse into its nucleus on a microsecond timescale. Atoms shouldn't *exist* as stable structures.

**Physicist:** And Bohr in 1913 makes the conceptual leap. [[bohr1913_atom_constitution]]. Two postulates. First: electrons occupy discrete circular orbits in which — by fiat, by postulate — they do not radiate. The allowed orbits are those where angular momentum is an integer multiple of $\hbar$. Second: transitions between orbits emit or absorb photons with $h\nu = E_m - E_n$.

`[cue: animation hist_bohr_proper_time]`

**Physicist:** Apply Newton's second law plus Coulomb's law plus the angular-momentum quantisation $L = n\hbar$. You get $r_n = n^2 a_0$ where $a_0 = 4\pi\varepsilon_0\hbar^2 / (m_e e^2)$ is the Bohr radius. The total energy at orbit $n$ comes out as $E_n = -\frac{1}{2}\frac{e^2}{4\pi\varepsilon_0 a_0}\frac{1}{n^2}$. The Rydberg formula falls out. The numerical value of the Rydberg constant comes out in terms of $m_e$, $e$, $h$ to within experimental precision.

**Experimentalist:** And the animation we're cuing here walks through *two* derivations side by side. The standard Bohr derivation. And then the dual-framework redo using Gill's [[FoundationsII-Classical|Foundations II classical Hamiltonian]] $K = H^2/(2mc^2) + mc^2/2$. The dual derivation gives the same Bohr levels at hydrogen velocities — agreement to about one part in $10^9$.

**Physicist:** Which is to say: no experimental distinguishability for the Bohr spectrum itself. The conceptual win — and this is what makes the dual framework interesting at the *interpretive* level — is that the positive-definite Hamiltonian $K$ has no negative-energy modes. It's bounded below by $mc^2/2$ — half the rest energy. The full mathematical structure that makes [[Dual_Relativistic_Quantum_Mechanics_I|DRQM I]] work cleanly, with no negative-probability-density problem of the kind Klein–Gordon has, is in place from the classical level.

**Historian:** Sommerfeld in 1916 — [[sommerfeld1916_atombau]] — extends Bohr to elliptical orbits with a special-relativistic momentum correction. Out falls the fine structure. The first appearance of the fine-structure constant $\alpha = e^2/(\hbar c) \approx 1/137$ in spectroscopy.

**Experimentalist:** And the same constant recurs in DRQM I's dual-Dirac fine structure in Chapter 4.

**Historian:** Same year, 1916, Einstein completes general relativity. The [[einstein1915_perihelion|November 1915 Mercury perihelion paper]] derives Mercury's 43 arcseconds per century — Le Verrier's 1859 residual — exactly from the new field equations. The [[einstein1916_grundlage|1916 *Grundlage* paper]] is the systematic exposition. We'll re-derive Mercury at grad-student level in [[07_PNT_GPS_SLR_QKD|Chapter 7]] as the bridge to GPS clock corrections.

**Experimentalist:** Last big primary-source clusters before we get to the proper-time interlude. [[compton1923_xray_scattering|Compton's 1923 X-ray scattering paper]]: X-rays scattered off effectively free electrons shift wavelength by $\Delta\lambda = (h/m_e c)(1 - \cos\theta)$. Direct confirmation that photons carry momentum $p = h/\lambda$, vindicating Einstein's 1916 conjecture. And [[de_broglie1924_thesis|de Broglie's 1924 thesis]]: if photons have wavelength $\lambda = h/p$, electrons should too. Davisson–Germer confirm in 1927.

**Physicist:** Plus the solid-state thread quietly opens: [[onnes1908_helium_liquefaction|Onnes liquefies helium]] in 1908, then in 1911 discovers [[onnes1911_superconductivity|superconductivity in mercury at 4.19 K]]. Resistance drops abruptly to zero. *Nobody* knows why. Forty-six years until BCS theory explains it. We pick that up in Chapter 5.

## Physics deep dive

**Physicist:** Let me work through the Compton scattering geometry — that's the second animation. A photon of wavelength $\lambda_0$ scatters off an electron initially at rest. After the collision: photon at angle $\theta$ with wavelength $\lambda$; electron recoiling at angle $\phi$ with momentum $p_e$.

Energy conservation:

$$\frac{hc}{\lambda_0} + m_e c^2 \;=\; \frac{hc}{\lambda} + \sqrt{p_e^2 c^2 + m_e^2 c^4}.$$

Momentum conservation, components:

$$\frac{h}{\lambda_0} = \frac{h}{\lambda}\cos\theta + p_e \cos\phi, \qquad 0 = \frac{h}{\lambda}\sin\theta - p_e \sin\phi.$$

Square and add, eliminate $\phi$. After algebra:

$$\Delta\lambda \;\equiv\; \lambda - \lambda_0 \;=\; \frac{h}{m_e c}(1 - \cos\theta).$$

That's the Compton formula. The combination $h/(m_e c) \approx 2.43$ picometers is the Compton wavelength of the electron — the inherent length scale of the interaction.

`[cue: animation hist_compton_null]`

**Experimentalist:** And the experimental verification: vary the scattering angle, measure the wavelength shift, plot. Comes out linear in $1 - \cos\theta$ with the right slope. Bothe and Geiger in 1925 confirm momentum is conserved event-by-event, not statistically — killing the Bohr–Kramers–Slater attempt to keep classical EM and explain photoelectric/Compton via emission statistics.

## Proper-time interlude

**Physicist:** Now the dual-framework redo. Tepper Gill's [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations|Maxwell paper]] velocity duality is $\mathbf{w}/c = \mathbf{u}/b$ with $b = \sqrt{c^2 + u^2}$. For Compton scattering off an electron at rest, $u = 0$ identically, $b = c$, and the dual Compton formula reduces *exactly* to the standard one. No measurable deviation. Same shift, same angular dependence, same Compton wavelength.

**Experimentalist:** And for a moving target electron? Atomic orbital velocities are $u/c \sim Z\alpha \sim 10^{-2}$ for hydrogen, up to about $10^{-1}$ for inner-shell electrons in heavy atoms.

**Physicist:** So $b/c - 1 \approx u^2/(2c^2) \sim 10^{-4}$ at worst. Compton's apparatus precision in 1923 was about 5%. Modern high-precision atomic Compton-profile measurements are at the 0.1% level. Still well above the dual-framework correction. Compton scattering is *not* a regime where the standard and dual framings can be distinguished.

**Experimentalist:** And the Bohr atom? Same animation walks through that too.

**Physicist:** Right. The dual Hamiltonian $K = H^2/(2mc^2) + mc^2/2$ has the same stationary states as standard non-relativistic $H$ at non-relativistic velocities — the velocity dependence enters at $O(v^4/c^4)$, which is below 1913 Bohr-formula precision by something like five orders of magnitude. Same Rydberg formula. Same Lyman, Balmer, Paschen series.

**Experimentalist:** So what does the dual framework *buy* us in this era, if all the predictions coincide with standard QM?

**Physicist:** Three things, all interpretive. First, the positive-definite Hamiltonian — no negative-energy modes from the start. Klein–Gordon's famous negative-probability-density problem doesn't arise because the classical-level Hamiltonian is already positive-definite, and the same construction carries into the relativistic quantum case ([[Dual_Relativistic_Quantum_Mechanics_I|DRQM I]] in Chapter 4). Second, the *Minkowski Incompatible Theorem* — the standard 4-velocity doesn't close under Hamiltonian flow without an additional structural assumption; the dual framework's squaring sidesteps this. Third — and this is the most interesting one for this campaign — when we eventually get to QED in Chapter 5, the dual framework's reformulation gives a different account of the [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL|Dyson divergence conjecture]] that opens a different path to QED renormalisation.

**Experimentalist:** None of which changes the predicted spectrum lines of hydrogen.

**Physicist:** Correct. The predictions agree. The framings differ at the level of mathematical convention and conceptual interpretation. Same experimental record either way.

## Closing

**Historian:** Chapter 3 in one episode. Planck 1900 through de Broglie 1924. Old quantum theory: ad hoc quantisation rules added to classical mechanics, working well for some things and breaking down for others. Helium spectrum, anomalous Zeeman effect, molecular spectra all unaccounted for — those are the failures that motivate the 1925 revolution.

**Experimentalist:** Plus the solid-state thread quietly opening with Onnes's helium and the discovery of superconductivity — sitting unexplained for nearly half a century.

**Physicist:** And the dual-framework anchor for this chapter: Gill's [[FoundationsII-Classical|Foundations II classical Hamiltonian]], rederiving the Bohr atom with $K = H^2/(2mc^2) + mc^2/2$. Same hydrogen spectrum at agreement to $10^{-9}$. Conceptual win in the positive-definiteness of $K$, payoff deferred to Chapter 4's DRQM I and Chapter 5's QED.

**Historian:** Next episode: 1925 to 1948. Matrix mechanics, wave mechanics, the Schrödinger equation, the Dirac equation, the discovery of the positron, the Pauli exclusion principle, spin. The Bohr–Sommerfeld old quantum theory is replaced by genuine quantum mechanics. Thanks for listening.

`[cue: end card with bibliography wikilinks for show notes]`

## Show notes

Auto-generated from the chapter's bibliography by `lint_episode.py`.
