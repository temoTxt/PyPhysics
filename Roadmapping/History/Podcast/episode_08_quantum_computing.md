---
episode: 08
title: "Quantum Computing: open questions and proper-time conjectures"
era: "forward"
chapter: 08_quantum_computing_open_questions
speakers: [Historian, Physicist]
target_runtime_min: 8
animations_cued: []
status: draft
---

# Episode 8 — Quantum Computing: open questions

> Companion dialogue script for [[08_quantum_computing_open_questions]]. First of the two *speculative* forward episodes — 2-voice (Historian + Physicist; Experimentalist drops out for the speculative chapters where there's less measured data to interrogate).

## Cold open

**Historian:** Quantum computing in 2026 sits in what John Preskill named the *NISQ era* — Noisy Intermediate-Scale Quantum. Hardware platforms with 100–1000 physical qubits, gate fidelities around 99.5–99.9%, coherence times of half a millisecond at best. No fault tolerance. A handful of contrived random-circuit-sampling demonstrations of quantum advantage, none of them on a useful problem.

**Physicist:** And the question for this chapter — a *speculative* forward chapter, not a derivation — is where Tepper Gill's dual-theory framework that we've been tracking through the campaign *might* bear on the open questions of quantum computing, when Gill hasn't published a thing on QC.

**Historian:** Five open questions in §4 of the chapter. Two of them tagged `#gill-silent` outright — topological qubits and quantum advantage at scale are entirely outside the dual framework's published reach. Three of them get speculative-but-explicit `#speculative` conjectures about how the dual framework's mathematical machinery could be applied.

**Physicist:** And the framing reminder, as always: this episode does not propose new physics. The conjectures are extrapolations from already-published Gill papers — Foundations I Mathematical for the KS-Hilbert space organisation, the Maxwell paper for the effective photon mass, DRQM I for the corrected g-factor — applied to questions Gill hasn't written about. We're explicit each time about the gap between "Gill has written X" and "extrapolating X to this QC question would suggest Y".

## Historical roots

**Historian:** The arc that lands us at quantum computing. Bohr's quantised atomic states in 1913. Heisenberg's matrix mechanics in 1925 — linear algebra of quantum mechanics, the language QC is written in. Dirac's spin-half equation in 1928 — spin-1/2 is the natural qubit. Bardeen, Cooper, Schrieffer in 1957 — superconductivity gives us the substrate for the dominant qubit modality. Josephson in 1962 — the nonlinear element that gives transmon qubits their anharmonicity. Lamb shift in 1947, cavity QED through the 1980s — the structures that become circuit QED in the 2000s.

**Physicist:** Then the explicitly-QC line. [[feynman1982_simulating_physics|Feynman 1982]] — *Simulating Physics with Computers* — the vision statement. Simulating quantum systems on classical computers requires exponential resources; a quantum computer should do it efficiently. [[deutsch1985_quantum_turing|Deutsch 1985]] — universal quantum computer formalised; quantum complexity theory begins. [[shor1994_factoring|Shor 1994]] — polynomial-time factoring of large integers. The killer app that motivated everything since. [[grover1996_search|Grover 1996]] — quadratic search speedup. [[kitaev2003_topological|Kitaev 2003]] — topological qubits via anyons, the speculative-but-actively-pursued alternative to surface-code fault tolerance.

## Five speculative bearings

**Physicist:** Five open questions in the chapter. Let me walk through how the dual framework *might* bear on three of them.

**Historian:** Pause first. The two it's silent on?

**Physicist:** Topological qubits — anyons in 2D fractional-quantum-Hall systems, Majorana zero modes in superconductor-semiconductor heterostructures. Anyon braiding statistics and 2D topological field theory are entirely outside what Gill has published. `#gill-silent`. And quantum advantage at scale — what's the first *useful* problem where a fault-tolerant quantum computer beats classical? That's an engineering and computer-science question; the dual-framework's bearing is at the level of QED corrections to the underlying physics, not on what algorithm to run.

**Historian:** OK — the three with `#speculative` bearing.

**Physicist:** First: decoherence channels in superconducting transmons. Standard QFT models transmon decoherence as a sum of Lindblad terms acting on the qubit's reduced density matrix, parameterised by phenomenological coupling constants to a bath. In Gill's [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL|Foundations I Mathematical]] Kuelbs–Steadman Hilbert space, the time-ordered operator algebra encodes bath coupling differently. *In principle* the FOC organisation could predict a different decoherence-channel decomposition. Whether it does — and whether the predicted decomposition matches experimental $T_1$ and $T_2$ measurements — is *not* something Gill has computed. The conjecture is "the KS-Hilbert reformulation might rearrange the decoherence-channel structure", not "Gill predicts X for transmon $T_1$".

**Historian:** And the second?

**Physicist:** Circuit-QED coherence and the effective photon mass. From [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] equation 6: the dual Maxwell equations have a dynamical effective photon mass $\mu$ that is non-zero whenever the source accelerates. In circuit-QED, the source is the transmon's own oscillating dipole; the effective acceleration is $\omega v / c \sim 10^{-5}$ at 5 GHz drive frequencies. The dual prediction for $\mu$ corresponds to a Compton wavelength orders of magnitude longer than the resonator size. No measurable dispersion at current circuit-QED precision. But the *form* of the dispersion relation is structurally different from standard EM. The speculative question is: at some future precision where this difference matters, does it predict a measurable departure from standard? Gill hasn't computed it for circuit-QED-relevant accelerations.

**Historian:** And the third?

**Physicist:** The g-factor as a QC test, via Finding 2 from [[FINDINGS_for_author_review]]. If a high-precision superconducting qubit array were designed to measure the electron g-factor in a novel regime — say, via a Penning-trap-style measurement embedded in a circuit-QED system — and reached precision beyond the current Harvard 13-digit measurement, the corrected $r_e$ value from Finding 2 could in principle be tested. By construction, the dual-framework $g_e$ with corrected $r_e$ matches the standard $g_e$ at the precision the formula was constructed to. Any differential test would require independent derivation of higher-order terms in the dual framework — not done.

## What this chapter is and isn't

**Historian:** A speculative roadmap. The conjectures are explicit about Gill having not published on the topics. Every claim tagged `#speculative` with an "extrapolated from `[[paper-anchor]]`" reference. Honest section §6 on experimental tests: no current experiment distinguishes standard from dual-framework predictions at QC-relevant precision.

**Physicist:** And the most concrete experimental hook — Penning-trap-style g-factor measurement in circuit-QED at 13+ digit precision — is on the *aspirational* side of current hardware capability. The bearing of the dual framework on quantum computing is conceptual rather than predictive. Alternative organisations of decoherence theory, alternative perturbative expansions that avoid the Dyson divergence asymptotic-series issue, alternative routes to circuit-QED coherence calculations.

## Closing

**Historian:** Chapter 8 — the first of the speculative forward chapters. Quantum computing is a consequence of the 1925–1965 historical arc. The dual framework's bearing on it sits at the level of mathematical organisation, with concrete hooks deferred to next-generation precision experiments.

**Physicist:** Next episode is Chapter 9 — fusion. Same format: speculative, roadmap rather than derivation, honest about how much is `#gill-silent`. The bearing of the dual framework on plasma physics is even thinner than on QC — fusion is mostly a `#gill-silent` chapter where the connection to the historical arc is context rather than prediction. Thanks for listening.

`[cue: end card with bibliography wikilinks for show notes]`

## Show notes

Auto-generated from the chapter's bibliography by `lint_episode.py`.
