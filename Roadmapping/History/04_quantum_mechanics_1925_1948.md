---
chapter: 04
title: "Quantum Mechanics 1925-1948"
era: "1925-1948"
threads: [quantum, solid-state]
animations: [hist_klein_gordon_vs_dual, hist_positron_isodual]
verification_anchors: ["Dual_Relativistic_Quantum_Mechanics_I", "Analytic_Representation_of_The_Dirac_Equation", "FoundationsII-Classical"]
status: draft
---

# Chapter 4 — Quantum Mechanics (1925–1948)

> **Framing principle (load-bearing).** We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics. This is the second of the campaign's load-bearing chapters: the [[Dual_Relativistic_Quantum_Mechanics_I|DRQM I]] paper that Trey co-authored is the verification anchor here, and the conceptual payoff — **dissolution of Klein–Gordon's negative-probability problem** via the dual framework's positive-definite Hamiltonian — is the headline interpretive win of the campaign. Predicted spectra coincide with standard QM at any current measurement precision; one open finding (Finding 2, DRQM I Eq. (III.22) g-factor) is flagged for author review.

## 1. Overview

The years 1925 through 1948 transform physics from a pile of *ad hoc* quantisation rules ([[03_old_quantum_theory_1900_1925|Chapter 3]]'s old quantum theory) into a coherent mathematical theory. [[heisenberg1925_quantum_mechanics|Heisenberg's matrix mechanics (1925)]] and [[schrodinger1926_quantisierung|Schrödinger's wave mechanics (1926)]] arrive within months of each other and prove mathematically equivalent ([[von_neumann1932_grundlagen|von Neumann 1932]]). [[born1926_statistical_interpretation|Born (1926)]] gives the probability interpretation. [[heisenberg1927_uncertainty|Heisenberg (1927)]] derives the uncertainty principle. [[klein_gordon1926|Klein–Gordon (1926)]] extends to relativistic kinematics but with the negative-energy + negative-probability problems. [[dirac1928_electron|Dirac (1928)]] resolves the first-derivative-in-time problem for spin-1/2, predicting antimatter, which [[anderson1932_positron|Anderson confirms (1932)]]. [[pauli1925_exclusion|Pauli's exclusion principle]] structures the periodic table. [[fermi1934_beta_decay|Fermi's beta-decay theory (1934)]] gives the first quantitative weak-interaction Lagrangian. [[chadwick1932_neutron|Chadwick (1932)]] finds the neutron, completing the nucleus. The solid-state thread accelerates with [[bloch1928_kristall|Bloch's band theory (1928)]] and the [[meissner_ochsenfeld1933_supraleiter|Meissner–Ochsenfeld superconductivity-as-thermodynamic-phase (1933)]]. The chapter ends at 1948 — just before the [[05_QED_renormalization_solid_state_1948_1965|QED renormalisation revolution]] picks up. Dual-framework anchor: [[Dual_Relativistic_Quantum_Mechanics_I]] for the dual Dirac equation and the Finding 2 g-factor flag.

## 2. Historical narrative

### 2.1 1925 — Matrix mechanics

In June 1925, Werner Heisenberg, isolated on Helgoland to recover from hay fever, has the breakthrough that founds modern quantum mechanics. Drop the orbital picture entirely. Replace classical phase-space coordinates with infinite arrays of *transition amplitudes* — one number per pair $(n, m)$ of atomic stationary states, measuring the strength of the transition between them. Compose these arrays the way you'd multiply matrices. Discover that they don't commute. Use the non-commutativity as the source of the discrete spectra.

[[heisenberg1925_quantum_mechanics|Heisenberg's "Umdeutung" paper]] is opaque even by his standards — he later admits to Pauli that he doesn't fully understand the formalism he has just invented. Within months, Max Born and Pascual Jordan in Göttingen recognise that the arrays are *matrices* and that Heisenberg's product rule is matrix multiplication. The Born–Heisenberg–Jordan three-author *Dreimännerarbeit* (early 1926) gives the systematic formulation: $pq - qp = -i\hbar$.

The matrix-mechanics derivation of the hydrogen spectrum is a *tour de force* — algebraically intricate where the Schrödinger derivation will be analytically natural. Wolfgang Pauli, with a 30-page calculation in October 1925, reproduces the Bohr spectrum from matrix mechanics + the Coulomb potential before Schrödinger's wave version is even published.

### 2.2 1925 — Spin and the exclusion principle

Two related developments in 1925 add the missing piece that the old quantum theory could never quite incorporate. [[pauli1925_exclusion|Pauli's exclusion principle]] requires a fourth quantum number with two-valuedness — "Zweideutigkeit" — to explain the closing of atomic shells. Months later, [[uhlenbeck_goudsmit1925_spin|Uhlenbeck and Goudsmit]] identify the fourth quantum number with an *intrinsic angular momentum* of the electron — *spin* — with magnitude $\hbar/2$. The g-factor is found to be $g = 2$ rather than the classical $g = 1$, an empirical fact that has no explanation in matrix or wave mechanics; it falls out for free from [[dirac1928_electron|Dirac's relativistic equation]] three years later.

### 2.3 1926 — Wave mechanics

Erwin Schrödinger, in a single year, publishes four foundational papers establishing *wave mechanics*. The first, [[schrodinger1926_quantisierung|*Quantisierung als Eigenwertproblem*]], introduces the time-independent equation $H\psi = E\psi$ — quantising hydrogen as an eigenvalue problem for a partial differential operator. The hydrogen spectrum, the harmonic oscillator, the angular-momentum eigenfunctions all come out cleanly.

Schrödinger initially derives a *relativistic* wave equation (we now call it Klein–Gordon — both [[klein_gordon1926|Klein and Gordon]] publish the same equation independently in early 1926, plus Schrödinger's unpublished earlier draft) and abandons it because the fine-structure predictions are wrong. The non-relativistic version with its Galilean kinematics gives the right spectrum, and Schrödinger publishes that instead.

Matrix mechanics and wave mechanics look superficially incompatible — algebraic vs. analytic, instantaneous vs. continuous. By April 1926, Schrödinger himself proves them mathematically equivalent. [[von_neumann1932_grundlagen|Von Neumann's 1932 axiomatisation]] makes the equivalence systematic: both formalisms are representations of the same abstract Hilbert-space structure.

### 2.4 1926 — The statistical interpretation

[[born1926_statistical_interpretation|Max Born's *Quantum Mechanics of Collision Processes*]] introduces the probability interpretation: $|\psi(\mathbf{x})|^2$ is the probability density to find the particle at $\mathbf{x}$. Schrödinger himself resists this interpretation for years — he had hoped $\psi$ described a real charge-density wave — but Born's reading wins by 1927. The Copenhagen interpretation crystallises in Bohr's September 1927 Como lecture, with complementarity, the role of the observer, and the operational character of the wave function.

### 2.5 1927 — Heisenberg uncertainty

[[heisenberg1927_uncertainty|Heisenberg's *anschaulicher Inhalt* paper]] derives the uncertainty principle from a gamma-ray-microscope thought experiment: localising an electron to $\Delta x$ requires photons of wavelength $\lambda \lesssim \Delta x$, which carry momentum $h/\lambda$ that the electron picks up on absorption, giving $\Delta p \gtrsim h/\Delta x$. Robertson (1929) and Schrödinger (1930) give the generalised commutator form $\Delta A \cdot \Delta B \geq \frac{1}{2}|\langle [A, B] \rangle|$. The uncertainty principle is the operational statement of the non-commutativity that Heisenberg's matrices already encoded.

### 2.6 1928 — Dirac equation

Paul Dirac, age 25, sets out to find a Lorentz-invariant wave equation that is *first-order* in time derivatives (avoiding Klein–Gordon's negative-probability problem, which arises because Klein–Gordon's conserved current involves $\partial_t \psi$). After several months he writes down

$$(i\gamma^\mu \partial_\mu - mc/\hbar)\psi = 0$$

with $\psi$ a four-component spinor and $\gamma^\mu$ matrices satisfying $\{\gamma^\mu, \gamma^\nu\} = 2\eta^{\mu\nu}$. [[dirac1928_electron|The 1928 paper]] derives, from the equation:

1. **Spin-1/2 emerges naturally** from the 4-spinor structure. Uhlenbeck–Goudsmit becomes a *prediction*, not a postulate.
2. **g-factor $g = 2$** — the empirical anomaly that troubled the old quantum theory.
3. **Hydrogen fine structure** with all relativistic + spin–orbit corrections built in. Agrees with Sommerfeld's *ad hoc* 1916 formula and with experiment.
4. **Positive-definite probability density** $\psi^\dagger \psi$.

But: the equation also has *negative-energy* solutions. Dirac initially treats these as a mathematical artefact, then proposes (1930–31) the *hole theory*: the negative-energy continuum is filled by a Pauli-exclusion sea; a hole in the sea appears as a particle with the electron's mass but opposite charge.

[[anderson1932_positron|Anderson's 1932 cosmic-ray cloud-chamber observation]] of the *positron* — same mass as electron, opposite charge, curving the wrong way in a magnetic field — confirms Dirac's hole-theory prediction within two years. The first observed example of *antimatter*.

### 2.7 1932–1934 — Nuclear physics opens

[[chadwick1932_neutron|Chadwick's neutron]] completes the nucleus. The neutron's near-zero charge means it can approach a nucleus without Coulomb repulsion, opening the experimental study of nuclear reactions and ultimately fission ([[09_fusion_open_questions|Chapter 9]] for fusion).

[[fermi1934_beta_decay|Fermi's beta-decay theory]] writes down the four-fermion contact interaction for $n \to p + e^- + \bar\nu_e$ with coupling $G_F$. The first quantitative theory of the *weak interaction*. Incorporates Pauli's 1930 neutrino hypothesis (introduced to save energy conservation in the continuous beta-decay spectrum) as a genuine particle in the interaction Lagrangian.

### 2.8 1928–1933 — Solid-state acceleration

[[bloch1928_kristall|Bloch's thesis]] shows that electrons in a perfectly periodic crystal have wavefunctions of the form $\psi_\mathbf{k}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}} u_\mathbf{k}(\mathbf{r})$ with $u_\mathbf{k}$ sharing the lattice periodicity. *Scattering off the periodic lattice is zero*. The puzzle of why metals conduct at all (mean free paths of centimeters, not Ångströms) dissolves. Wilson (1931) classifies metals, semiconductors, and insulators by band filling. Sommerfeld and Bethe's 1933 *Handbuch* article is the standard reference.

[[meissner_ochsenfeld1933_supraleiter|The Meissner effect (1933)]] establishes superconductivity as a fundamentally distinct thermodynamic phase — distinct from a hypothetical "perfect conductor" — by showing field expulsion at $T_c$. The phenomenological London equations (1935) capture the Meissner expulsion; the microscopic BCS theory of [[05_QED_renormalization_solid_state_1948_1965|Chapter 5]] finally explains the underlying mechanism.

## 3. Proper-time commentary

### 3.1 What's directly verified

**Dual Dirac equation.** `#verified` from [[Dual_Relativistic_Quantum_Mechanics_I]]. The dual Hamiltonian is

$$K_D \;=\; \frac{(c\boldsymbol{\alpha}\cdot\mathbf{p} + \beta m c^2)^2}{2m c^2} + \frac{mc^2}{2}$$

— a positive-definite squared form of the standard Dirac Hamiltonian. Eigenvalues of $K_D$ are non-negative; bounded below by $mc^2/2$. The non-relativistic limit reproduces the standard Schrödinger equation; the relativistic limit reproduces the Dirac fine-structure spectrum. **No negative-energy modes**; no need for hole theory at the level of the equation; antimatter is handled by isodual mathematics (see §3.2). Wolfram MCP verification of the eigenvalue structure recorded in the linked anchor.

**Negative-probability dissolution.** `#verified` (immediate consequence). Klein–Gordon's notorious problem — the conserved current is $j^0 \propto \psi^*\partial_t\psi - \partial_t\psi^* \cdot \psi$, which is not positive-definite — does not arise in the dual framework. The dual square-root operator from [[Analytic_Representation_of_The_Square-Root_Operator]] gives a single-particle relativistic theory with manifest positive-definite probability density. Animated in `hist_klein_gordon_vs_dual.py`.

### 3.2 What's mechanically inferred

**Antimatter via Santilli isodual.** `#inferred` from [[Dual_Relativistic_Quantum_Mechanics_I]]. The dual Dirac equation has an *isodual* companion equation — obtained by the Santilli isodual map (which reverses the metric sign, energy sign, and time direction simultaneously) — whose solutions correspond to positrons. The standard Dirac hole-theory reinterpretation is replaced by an isodual reinterpretation: positrons are not "holes in a filled sea" but solutions of the isodual companion equation. Both framings predict the same scattering cross-sections, decay widths, and g-factors for positrons. **Experimental distinguishability:** none at any current precision; the difference is in the mathematical *organisation* of the theory, not in its empirical predictions. Animated in `hist_positron_isodual.py`.

**g-factor of the electron from the dual Dirac equation.** `#inferred` from [[Dual_Relativistic_Quantum_Mechanics_I|DRQM I]] Sec III. The dual Dirac eigenvalue problem (DRQM I Eq. (III.22)) gives the g-factor as $g_r = 2[1 - 4r_0/(2r + r_0)]$ where $r_0$ is the classical electron radius and $r$ is an internal parameter. **Open finding** ([[FINDINGS_for_author_review#Finding 2 — DRQM I Eq. (III.22): published $r_e$ does not reproduce the experimental $g_e$]]): the value $r = r_e \approx 0.499857150068631\,r_0$ stated in DRQM I gives $g \approx -2.0005714$, not the experimentally measured $g_e \approx -2.00231930436256$. The correct $r$ to reproduce the experimental g-factor is $r_e \approx 0.499420510\,r_0$. Flagged for author review in [[FINDINGS_for_author_review]]; the fix is a numerical correction to $r$, not a structural change to the equation.

**Sommerfeld fine structure in the dual framework.** `#inferred` from DRQM I. The dual Dirac equation reproduces the standard Sommerfeld fine-structure formula to the order $\alpha^2$ that 1916–1928 spectroscopy could measure; higher-order corrections (Lamb shift, $\alpha^4$ recoil) are deferred to [[05_QED_renormalization_solid_state_1948_1965|Chapter 5]] where the QED framework handles them.

### 3.3 What Gill is silent on

- **Born statistical interpretation, Heisenberg uncertainty, Copenhagen interpretation.** `#gill-silent`. These are *interpretive* claims about quantum measurement, not mathematical predictions. The dual framework's positive-definite Hilbert space inherits Born's rule for $|\psi|^2$, but Gill has no published commentary on the interpretive debates.
- **Bloch band theory, Meissner effect, neutron discovery, beta decay.** `#gill-silent`. Standard non-relativistic QM (Bloch), thermodynamics (Meissner), nuclear physics (Chadwick, Fermi). The dual framework's bearing on superconductivity is treated in [[08_quantum_computing_open_questions|Chapter 8]] as a speculative `#speculative` hook for circuit QED.

## 4. Key derivations worth animating

| Manim scene | Status | What it shows |
|---|---|---|
| [`hist_klein_gordon_vs_dual.py`](../Animations/manim_scenes/hist_klein_gordon_vs_dual.py) | rendered | Klein–Gordon equation; the non-positive-definite conserved current; the dual reformulation via the positive-definite $K = H^2/(2mc^2) + mc^2/2$ with the square-root operator from `[[Analytic_Representation_of_The_Square-Root_Operator]]`; the resulting positive-definite probability density. The "dissolution" of the negative-probability problem visualised by showing the densities side by side. |
| [`hist_positron_isodual.py`](../Animations/manim_scenes/hist_positron_isodual.py) | rendered | Standard Dirac negative-energy continuum + Dirac sea + hole = positron; alongside the dual Dirac + isodual companion equation = positron. Same experimental predictions, two different mathematical organisations. |

## 5. Primary sources cited

- [[heisenberg1925_quantum_mechanics]] — Umdeutung; matrix mechanics.
- [[pauli1925_exclusion]] — exclusion principle.
- [[schrodinger1926_quantisierung]] — wave mechanics.
- [[klein_gordon1926]] — relativistic scalar equation; negative-probability problem.
- [[born1926_statistical_interpretation]] — $|\psi|^2$ probability rule.
- [[heisenberg1927_uncertainty]] — uncertainty principle.
- [[dirac1928_electron]] — Dirac equation; spin-1/2 emergent; g=2; positron prediction.
- [[bloch1928_kristall]] — band theory.
- [[anderson1932_positron]] — antimatter discovery.
- [[chadwick1932_neutron]] — neutron discovery.
- [[meissner_ochsenfeld1933_supraleiter]] — superconductivity as thermodynamic phase.
- [[fermi1934_beta_decay]] — weak interaction theory.

## 6. Retrospective reviews drawn on

- [[mehra_rechenberg_qm_history]] — six-volume detailed history of the 1925–1948 development.
- [[pais1982_subtle_is_the_lord]] — continued from Ch 3 for Einstein's late-1930s quantum-foundations debates.
- [[jammer1966_conceptual_qm]] — continued from Ch 3.
- [[dirac1958_principles]] — Dirac's own textbook exposition.
- [[von_neumann1932_grundlagen]] — axiomatic Hilbert-space foundations.
- [[weinberg1995_qft_vol1]] — modern graduate QFT reference.
- [[peskin_schroeder1995_qft]] — modern teaching textbook.
- [[jackson1998_classical_electrodynamics]] — continued from Ch 2 for classical-relativistic Lagrangians.

## 7. Cross-references

- Backward: [[03_old_quantum_theory_1900_1925]].
- Forward: [[05_QED_renormalization_solid_state_1948_1965]] — Lamb shift, Schwinger/Feynman/Tomonaga/Dyson, transistor, BCS.
- Verification anchors: [[Dual_Relativistic_Quantum_Mechanics_I]], [[Analytic_Representation_of_The_Dirac_Equation]], [[Analytic_Representation_of_The_Square-Root_Operator]], [[FoundationsII-Classical]], [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]].
- Findings: [[FINDINGS_for_author_review]] — **Finding 2** (DRQM I Eq. (III.22) g-factor) is the load-bearing one for this chapter.
