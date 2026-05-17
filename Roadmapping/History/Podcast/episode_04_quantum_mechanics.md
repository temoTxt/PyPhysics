---
episode: 04
title: "Quantum Mechanics 1925-1948"
era: "1925-1948"
chapter: 04_quantum_mechanics_1925_1948
speakers: [Historian, Physicist, Experimentalist]
target_runtime_min: 22
animations_cued:
  - hist_klein_gordon_vs_dual
  - hist_positron_isodual
status: draft
---

# Episode 4 — Quantum Mechanics (1925–1948)

> Companion dialogue script for [[04_quantum_mechanics_1925_1948]]. Second of the two load-bearing chapters: the dual framework's Dirac equation lives here, the Finding 2 g-factor flag is here, and the headline interpretive payoff — dissolution of Klein–Gordon's negative-probability problem — is here.

## Cold open

**Historian:** In June 1925, a 23-year-old physicist named Werner Heisenberg flees Göttingen for the rocky North Sea island of Helgoland. He has hay fever. He needs to escape the linden pollen of central Germany. He brings with him an unsolved problem he has been chewing on for two years: how to derive the spectral lines of hydrogen from first principles, without the *ad hoc* quantisation postulates that Bohr and Sommerfeld had imposed by fiat.

**Experimentalist:** The old quantum theory works for hydrogen. Doesn't work for helium. Doesn't work for the anomalous Zeeman effect. Doesn't work for molecular spectra. By 1924 everyone in the field knows the framework is incomplete. They're waiting for the right idea.

**Historian:** Heisenberg has the right idea on Helgoland in late June. Drop the orbital picture entirely. Replace position and momentum with arrays of *transition amplitudes* — one number per pair of atomic stationary states. Compose these arrays the way you'd multiply matrices. Discover they don't commute. Use the non-commutativity as the source of the spectral discreteness.

**Physicist:** [[heisenberg1925_quantum_mechanics|The Umdeutung paper]] is the genesis of matrix mechanics. Eight months later Schrödinger arrives at the same physics from the opposite direction — a continuous wave equation, [[schrodinger1926_quantisierung|*Quantisierung als Eigenwertproblem*]] — and the two formalisms turn out to be mathematically equivalent. Plus [[born1926_statistical_interpretation|Born's probability interpretation]] of $|\psi|^2$. Plus [[heisenberg1927_uncertainty|Heisenberg's uncertainty principle]]. Plus [[dirac1928_electron|Dirac's 1928 relativistic equation]] predicting spin, $g = 2$, antimatter. Plus the experimental confirmation of antimatter by [[anderson1932_positron|Anderson in 1932]]. Plus [[chadwick1932_neutron|Chadwick's neutron]] and [[fermi1934_beta_decay|Fermi's beta-decay theory]]. The solid-state thread accelerates with [[bloch1928_kristall|Bloch's band theory]] and the [[meissner_ochsenfeld1933_supraleiter|Meissner effect]]. Twenty-three years of revolution.

**Experimentalist:** And for this campaign — for the dual-theory framework we've been tracking through these episodes — this is *the second load-bearing chapter*. The first was Chapter 2 with Maxwell. This one is anchored on Tepper Gill's [[Dual_Relativistic_Quantum_Mechanics_I|*Dual Relativistic Quantum Mechanics I*]], which Trey co-authored, and which the verification campaign worked through equation by equation.

**Physicist:** The headline payoff is interpretive. The Klein–Gordon equation — the natural relativistic extension of Schrödinger — has a notorious negative-probability problem: the conserved current you'd want to interpret as probability density is *not positive-definite*. Schrödinger himself encountered the problem in early 1926, retreated to the non-relativistic limit, and Dirac eventually solved it by going to a first-order-in-time equation with a spinor structure. In Gill's dual framework, the negative-probability problem *doesn't arise at all*. The dual Hamiltonian is positive-definite by construction. We'll walk through that in the proper-time interlude.

**Historian:** Let's start with Heisenberg on Helgoland.

## Historical sweep

**Historian:** Heisenberg, on Helgoland, in late June 1925. The breakthrough. He spends nights walking the cliffs of the island, working through the algebra. Returns to Göttingen with a draft in early July. Shows it to Pauli, who is enthusiastic. Shows it to Born, who recognises within days that the arrays are *matrices* and that Heisenberg's product rule is matrix multiplication.

**Experimentalist:** Heisenberg himself doesn't know matrix theory. Born does. Born and Jordan write the formal paper. The famous *Dreimännerarbeit* — the "three-man paper" — comes out in early 1926 with Heisenberg, Born, and Jordan as co-authors. By autumn 1925, Pauli — without computing aid — has reproduced the entire hydrogen spectrum from matrix mechanics in a 30-page algebraic *tour de force*.

**Historian:** And while Heisenberg is doing matrix mechanics in Göttingen, Schrödinger in Zurich is doing something completely different. He's been thinking about [[de_broglie1924_thesis|de Broglie's matter-wave hypothesis]] from Chapter 3 — if electrons have wavelength, they should obey a wave equation. He writes one down. [[schrodinger1926_quantisierung|*Quantisierung als Eigenwertproblem*]] — quantisation as an eigenvalue problem. The time-independent Schrödinger equation, $H\psi = E\psi$. Plug in the Coulomb potential. Out comes the hydrogen spectrum.

**Physicist:** And here's an important subplot. Schrödinger initially derives a *relativistic* version of his equation. The one we now call Klein–Gordon — both [[klein_gordon1926|Klein and Gordon]] publish the same equation independently in early 1926, and Schrödinger's earlier draft never appears in print. The relativistic version gives *wrong* fine-structure predictions for hydrogen, so Schrödinger abandons it and publishes the non-relativistic Galilean version instead.

**Experimentalist:** So we have matrix mechanics and wave mechanics in 1926. They look completely different — discrete vs continuous, algebraic vs analytic, instantaneous-jumps vs smooth-evolution.

**Physicist:** By April 1926, Schrödinger himself proves the two formalisms mathematically equivalent. They're two representations of the same abstract Hilbert-space structure. [[von_neumann1932_grundlagen|Von Neumann's 1932 axiomatisation]] makes this rigorous. They are the same theory.

**Historian:** Then in mid-1926, [[born1926_statistical_interpretation|Born's *Quantum Mechanics of Collision Processes*]]. The probability interpretation. $|\psi(\mathbf{x})|^2$ is the probability density to find the particle at $\mathbf{x}$. Schrödinger resists this for years. He had hoped $\psi$ would describe a literal charge-density wave. Born's reading wins by 1927.

**Experimentalist:** And in [[heisenberg1927_uncertainty|early 1927, Heisenberg's uncertainty paper]]. Gamma-ray microscope thought experiment. Localising an electron to $\Delta x$ requires photons of wavelength $\lambda \leq \Delta x$, which carry momentum $h/\lambda$ and disturb the electron by $\Delta p \geq h/\Delta x$. So $\Delta x \cdot \Delta p \geq h$ — Robertson tightens this to $\hbar/2$ a couple of years later.

**Historian:** Bohr in September 1927 at Como crystallises *complementarity* — wave and particle descriptions are mutually exclusive but jointly necessary. The Copenhagen interpretation is born. Einstein objects, immediately and famously, at the 1927 Solvay Congress.

**Physicist:** And then in 1928, the most influential equation of the entire era. Paul Dirac, age 25, working at St John's College Cambridge. He wants a Lorentz-invariant wave equation that is *first-order* in time derivatives — to avoid the Klein–Gordon negative-probability problem. After some months he writes down

$$(i\gamma^\mu \partial_\mu - mc/\hbar)\psi = 0$$

with $\psi$ a four-component spinor and the $\gamma^\mu$ matrices satisfying $\{\gamma^\mu, \gamma^\nu\} = 2\eta^{\mu\nu}$. [[dirac1928_electron|The 1928 paper]] does four things at once. *First*: spin-1/2 emerges naturally from the 4-spinor structure. The Uhlenbeck–Goudsmit 1925 postulate is now a prediction. *Second*: the g-factor comes out as $g = 2$, matching the empirical anomaly. *Third*: the full fine structure of hydrogen including spin–orbit. *Fourth*: a positive-definite probability density $\psi^\dagger \psi$.

`[cue: animation hist_klein_gordon_vs_dual]`

**Physicist:** And the animation here puts the Klein–Gordon problem and the dual-framework resolution side by side. We'll come back to it in the proper-time interlude.

**Experimentalist:** But Dirac's equation also has *negative-energy* solutions. He initially treats them as a mathematical artefact. Then in 1930–31 he proposes the *hole theory*: the negative-energy continuum is filled by a Pauli-exclusion sea; a hole in the sea appears as a particle with the electron's mass but opposite charge.

**Historian:** And two years later, [[anderson1932_positron|Carl Anderson at Caltech]] is examining cosmic-ray cloud-chamber tracks. He sees a particle curving the wrong way in his magnetic field. Same curvature radius as an electron, opposite direction. He concludes it has the electron's mass but opposite charge.

**Physicist:** Antimatter. Predicted by Dirac in 1930–31. Observed by Anderson in 1932. Two years.

**Experimentalist:** Same year, in Cambridge, [[chadwick1932_neutron|James Chadwick]] discovers the neutron. Beryllium bombarded by polonium-emitted alphas produces an unidentified radiation. The radiation knocks recoil nuclei out of paraffin wax with energies that are inconsistent with a gamma-ray interpretation. Chadwick proposes — in a one-page *Nature* letter — that the radiation is a neutral particle with mass approximately one atomic unit. The neutron resolves the nitrogen-14 spin puzzle, completes the nucleus, and opens the experimental study of nuclear physics.

**Historian:** Two years after that, [[fermi1934_beta_decay|Fermi]] writes down the four-fermion theory of beta decay. $n \to p + e^- + \bar{\nu}_e$, with a contact interaction coupling $G_F$. The first quantitative theory of the *weak interaction*. Incorporates Pauli's 1930 neutrino hypothesis as a genuine particle. Will be replaced by the electroweak gauge theory in 1967–68 but remains the low-energy limit of weak processes.

**Experimentalist:** Solid-state runs alongside. [[bloch1928_kristall|Bloch's 1928 thesis]]: electrons in a perfectly periodic crystal don't scatter at all. Bloch waves. The mystery of why metals conduct — mean free paths of *centimeters*, not Ångströms — dissolves. Wilson in 1931 sorts metals, semiconductors, and insulators by band filling. Sommerfeld and Bethe's 1933 *Handbuch* article is the standard exposition.

**Historian:** And in 1933 — [[meissner_ochsenfeld1933_supraleiter|Meissner and Ochsenfeld]] discover that a superconductor cooled through its critical temperature in a magnetic field *expels* the field. Distinguishes superconductivity from a hypothetical "perfect conductor" (which would *trap* the field). Superconductivity is a thermodynamic phase, not just zero resistance.

**Physicist:** The microscopic explanation has to wait until 1957. That's Chapter 5.

## Physics deep dive

**Physicist:** The Klein–Gordon equation, with the negative-probability problem visible. Start from the relativistic energy-momentum relation $E^2 = p^2 c^2 + m^2 c^4$. Promote to operators: $E \to i\hbar \partial_t$, $\mathbf{p} \to -i\hbar \nabla$. Square both sides:

$$-\hbar^2 \partial_t^2 \psi \;=\; -\hbar^2 c^2 \nabla^2 \psi + m^2 c^4 \psi.$$

Rearrange:

$$\left(\Box + \frac{m^2 c^2}{\hbar^2}\right) \psi \;=\; 0,$$

with $\Box = (1/c^2)\partial_t^2 - \nabla^2$. The Klein–Gordon equation.

**Experimentalist:** And the conserved current?

**Physicist:** Construct it from the equation's Lagrangian. Comes out as

$$j^\mu \;=\; \frac{i\hbar}{2m}\left(\psi^* \partial^\mu \psi - \psi \, \partial^\mu \psi^*\right).$$

The zeroth component:

$$j^0 \;=\; \frac{i\hbar}{2m c^2}\left(\psi^* \partial_t \psi - \psi \, \partial_t \psi^*\right).$$

Not positive-definite. There exist initial wave functions for which $j^0 < 0$ in some spacetime region. Cannot be a probability density.

**Experimentalist:** And Dirac's first-order equation fixes this — gives $j^0 = \psi^\dagger \psi$, which is manifestly positive-definite.

**Physicist:** Right. The price Dirac pays is the four-spinor structure plus the negative-energy continuum that he then has to handle via hole theory.

## Proper-time interlude

**Physicist:** Now the dual-framework redo. Gill's [[Dual_Relativistic_Quantum_Mechanics_I|DRQM I]] paper writes down a positive-definite Hamiltonian by squaring:

$$K_D \;=\; \frac{(c\boldsymbol{\alpha}\cdot\mathbf{p} + \beta m c^2)^2}{2 m c^2} + \frac{m c^2}{2}.$$

This is the dual Dirac Hamiltonian. It's positive-definite by construction — the numerator is the square of a Hermitian operator, divided by a positive constant, plus a positive shift. Bounded below by $mc^2/2$, exactly half the rest energy.

**Experimentalist:** And the eigenvalues?

**Physicist:** Same eigenstates as the standard Dirac Hamiltonian. The spectrum is the squared spectrum of standard Dirac, scaled by $1/(2mc^2)$, shifted by $mc^2/2$. Numerically: for hydrogen, $K_D$ gives the same fine-structure splittings as standard Dirac. For the free particle: $K_D = (p^2 c^2 + m^2 c^4)/(2 m c^2) + mc^2/2$, which expands at low momentum to $mc^2/2 + p^2/(2m) + mc^2/2 = mc^2 + p^2/(2m)$ — the standard non-relativistic kinetic-plus-rest energy.

**Experimentalist:** And the probability density?

**Physicist:** Positive-definite. By construction. The dual square-root operator construction in [[Analytic_Representation_of_The_Square-Root_Operator]] gives a single-particle relativistic theory with positive-definite probability density. The Klein–Gordon problem doesn't arise.

**Experimentalist:** That's the headline payoff?

**Physicist:** That's the conceptual headline. No negative-energy modes. No hole theory required at the level of the equation. Antimatter is handled instead via Santilli's *isodual* mathematics — a separate companion equation for positrons, related to the dual Dirac by an isodual map that reverses metric sign, energy sign, and time direction simultaneously. Same experimental predictions as standard Dirac + hole theory. Different mathematical organisation.

`[cue: animation hist_positron_isodual]`

**Physicist:** The animation we're cuing here puts the two organisations side by side — standard Dirac with the negative-energy continuum, the Dirac sea, the hole reinterpretation; alongside dual Dirac + isodual companion equation = positron. Both produce the same observed positron cross-sections, decay widths, and g-factor.

**Experimentalist:** And the open finding from the verification campaign?

**Physicist:** Finding 2 in [[FINDINGS_for_author_review]]. The DRQM I paper gives the electron g-factor from the dual Dirac eigenvalue problem as $g_r = 2[1 - 4r_0/(2r + r_0)]$, where $r_0$ is the classical electron radius and $r$ is an internal parameter. The paper states $r = r_e \approx 0.499857150068631\,r_0$. Plug that in, you get $g \approx -2.0005714$. The experimental value is $g_e \approx -2.00231930436256$. The required $r$ to reproduce the experimental g-factor — what we computed via the verification campaign's Wolfram MCP work — is $r_e \approx 0.499420510\,r_0$.

**Experimentalist:** A numerical correction, not a structural change.

**Physicist:** Right. The dual Dirac equation's formula for $g$ is correct; the value of the parameter $r$ stated in the paper is off. Flagged for author review; the structural correctness of the equation isn't in question.

## Closing

**Historian:** Chapter 4 in one episode. 1925 through 1948. Matrix mechanics, wave mechanics, Born's probability rule, Heisenberg uncertainty, the Klein–Gordon and Dirac equations, antimatter discovered, the neutron, Fermi's weak interaction, Bloch's band theory, the Meissner effect. The chapter ends just before the QED renormalisation revolution.

**Experimentalist:** And the dual-framework payoff is real and substantive: the positive-definite dual Hamiltonian eliminates Klein–Gordon's negative-probability problem at the level of the *equation*, not just by retreating to a first-order spinor structure. Plus the verification anchor — DRQM I — has one open finding flagged for the authors: the g-factor parameter $r$ requires a numerical correction.

**Physicist:** Same predicted spectrum as standard QM at any current precision. Conceptual gain in the positive-definite Hilbert space and the isodual organisation of antimatter. One numerical finding for the authors of DRQM I to look at.

**Historian:** Next episode: 1948 to 1965. The Lamb shift breaks the Dirac fine structure. Schwinger, Tomonaga, Feynman, Dyson invent QED renormalisation independently. Wu measures parity violation in 1957. The transistor in 1947, BCS in 1957, Josephson in 1962. The chapter that closes the campaign's historical scope — and the chapter where the dual-framework payoff is biggest, with the *Dyson divergence conjecture* getting reframed through Gill's Foundations I and Feynman operator calculus. Thanks for listening.

`[cue: end card with bibliography wikilinks for show notes]`

## Show notes

Auto-generated from the chapter's bibliography by `lint_episode.py`.
