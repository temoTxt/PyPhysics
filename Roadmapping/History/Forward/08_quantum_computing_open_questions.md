---
chapter: 08
title: "Quantum Computing: open questions and proper-time conjectures"
era: "forward"
variant: open-questions
threads: [quantum, solid-state]
animations: []
verification_anchors: ["FOUNDATIONS_FOR_QED_I_MATHEMATICAL", "Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations", "FINDINGS_for_author_review"]
status: draft
---

# Chapter 8 — Quantum Computing: open questions and proper-time conjectures

> **Framing principle (load-bearing).** We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics. This chapter is a *roadmap*, not a derivation: it identifies open questions in quantum computing where Gill has not published and where the dual-theory framework *might* bear, *without* completing the derivation. Every conjecture is tagged `#speculative` and is explicitly an "extrapolated from `[[paper-anchor]]`" claim, not a "Gill predicts" claim.

## 1. Overview

Quantum computing in 2026 is in the **NISQ era** — Noisy Intermediate-Scale Quantum, a term [[preskill2018_nisq|coined by Preskill]] for the ~50–1000 qubit regime without fault tolerance. The platforms that exist: superconducting transmons (IBM, Google, Rigetti, IQM), trapped ions (IonQ, Quantinuum), neutral atoms (Atom Computing, QuEra, Pasqal), and a handful of optical photonic systems. The platforms that are speculative but actively pursued: topological qubits (Microsoft's Majorana programme — see [[kitaev2003_topological|Kitaev 2003]]). The *consequences* of quantum computing — if scaled fault-tolerantly — include Shor's-algorithm-driven breaks on RSA cryptography, dramatic speedups in quantum simulation of chemistry and materials, and possibly novel optimization advantages via Grover-type and QAOA-type algorithms.

This chapter sits in the *speculative* forward bucket. Gill has not published on QC. The chapter identifies five places where the dual framework *might* bear, with explicit `#speculative` tags and explicit "Gill has not published on this — extrapolated from `[[paper-anchor]]`" disclaimers.

## 2. Historical roots

- **[[bohr1913_atom_constitution|Bohr atom]] (Ch 3):** quantisation of states. Qubit prototype.
- **[[heisenberg1925_quantum_mechanics|Heisenberg matrix mechanics]] (Ch 4):** linear algebra of QM. The lingua franca of QC.
- **[[dirac1928_electron|Dirac equation]] (Ch 4):** spin-1/2 systems. Natural qubit.
- **[[bcs1957_superconductivity|BCS Cooper pairs]] (Ch 5):** superconducting condensate. Substrate for superconducting qubits.
- **[[josephson1962_tunneling|Josephson junctions]] (Ch 5):** the nonlinear element that gives transmon qubits their anharmonicity.
- **[[lamb_retherford1947_microwave|Lamb shift]] (Ch 5):** cavity QED roots; eventually circuit QED.
- **[[feynman1982_simulating_physics|Feynman 1982]]:** the vision statement.
- **[[deutsch1985_quantum_turing|Deutsch 1985]]:** the universal quantum computer.
- **[[shor1994_factoring|Shor 1994]]:** the killer-app algorithm.

## 3. Current state (2026 perspective)

NISQ-class hardware operates with ~100–1000 physical qubits, gate fidelities of 99.5–99.9%, and coherence times of 50–500 μs depending on platform. A handful of demonstrations of quantum advantage on contrived random-circuit-sampling tasks ([[arute2019_quantum_supremacy|Google Sycamore 2019]], USTC photonic 2020, IBM 2023) — though the classical-vs-quantum advantage at any given depth is actively contested. Fault-tolerant quantum computation requires the surface code or an analog, with quoted overheads of $\sim 10^3$–$10^5$ physical qubits per logical qubit; not yet demonstrated end-to-end. Cryptographically relevant Shor's-algorithm-on-RSA-2048 estimates: ~20 million noisy physical qubits, or a few thousand logical qubits with $\sim 10^{-9}$ logical error rate per cycle. Optimistic timeline: 2035–2045.

## 4. Major open questions

Numbered list. Each: problem statement + one-sentence status.

1. **Decoherence channels in superconducting transmons.** What sets the $T_1$ time? Why hasn't transmon $T_1$ improved beyond ~500 μs in the past five years? *Current status:* a mix of two-level-system defects in the substrate/oxide, quasiparticle poisoning, and Purcell-type radiative loss; relative contributions unclear.
2. **Mechanism of error correlations.** Surface-code fault tolerance assumes errors are uncorrelated across qubits; recent experiments show otherwise (cosmic-ray events, correlated $T_1$ drops). *Current status:* mitigations being developed; theoretical framework for correlated noise immature.
3. **Threshold theorem in realistic hardware.** Surface-code threshold is ~1% gate error; current best is ~0.1%; but is the theorem's noise model adequate? *Current status:* threshold is robust to *some* noise correlations; not to others.
4. **Topological qubits.** Can Majorana zero modes (or anyons in 2D fractional-quantum-Hall systems) be braided non-Abelianly with low enough error rates to enable fault-tolerant QC without surface-code overhead? *Current status:* candidate signatures observed; no confirmed braiding operations as of 2025; Microsoft's 2023 Nature retraction set the field back.
5. **Quantum advantage at scale.** Beyond contrived random-circuit sampling, what's the first *useful* problem where a fault-tolerant quantum computer demonstrably outperforms classical? *Current status:* leading candidates are chemistry (FeMoco nitrogenase ground-state energy), materials (high-$T_c$ superconductor mechanism), and cryptanalysis (Shor on RSA).

## 5. Speculative proper-time implications

`#speculative` extrapolations. Each entry explicitly notes: Gill has not published on this; the conjecture is extrapolated from the cited verification anchor.

### 5.1 Question 1 (decoherence): does the dual framework's KS-Hilbert space give a different account?

- **Extrapolation source:** [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]] Sec 3–4 — the Kuelbs–Steadman Hilbert space + time-ordered Feynman Operator Calculus organisation of QFT.
- **Speculative bearing:** `#speculative`. Standard QFT models transmon decoherence as a sum of Lindblad terms acting on the qubit's reduced density matrix, parameterised by phenomenological coupling constants to a bath. In Gill's KS-Hilbert organisation, the time-ordered operator algebra encodes bath coupling differently — the time-ordering structure that the FOC papers organise on the KS space could in principle predict a different decoherence-channel decomposition. *Whether* it does, and whether the predicted decomposition matches experiment, is undertermined. The bath-coupling parameters would need to be re-fit; the structural prediction is the *shape* of the Lindblad operators.
- **Why this is not a prediction:** Gill has not done the explicit construction. The conjecture is "the KS-Hilbert reformulation might rearrange the decoherence-channel structure", not "Gill predicts X for transmon $T_1$".

### 5.2 Question 2 (effective photon mass + cavity QED): circuit-QED coherence

- **Extrapolation source:** [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] Eq. (6) — the dynamical effective photon mass $\mu$ that arises in the dual Maxwell equations when sources accelerate.
- **Speculative bearing:** `#speculative`. In circuit-QED, the qubit-resonator coupling depends on the photon field's dispersion relation. In standard EM, the photon is massless. In the dual EM framework, $\mu$ is non-zero whenever the source accelerates. For a transmon qubit driven by microwave fields at $\sim 5$ GHz, the source is the qubit's own oscillating dipole; the effective acceleration is $\sim \omega v / c \sim 10^{-5}$. The dual prediction for $\mu$ would correspond to a Compton wavelength orders of magnitude longer than the resonator size, producing no measurable dispersion at current circuit-QED precision. But the *form* of the dispersion is structurally different — the dual prediction includes a frequency-dependent term that standard EM doesn't.
- **Why this is not a prediction:** Gill has not computed $\mu$ for circuit-QED-relevant accelerations. The order-of-magnitude estimate ~$10^{-5}$ Compton-wavelength-correction is mechanical extrapolation, not a worked-out calculation.

### 5.3 Question 3 (Dyson divergence): does it affect high-precision QC test predictions?

- **Extrapolation source:** [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]] + [[Feynman_Operator_Calculus_Papers]].
- **Speculative bearing:** `#speculative`. QC verification at the next precision frontier (e.g., a 100-qubit transmon's full gate-error tomography at the $10^{-9}$ level) might require QED corrections at orders where the Dyson asymptotic-series truncation matters. The KS-Hilbert + FOC organisation gives the same numerical answers at current QED precision, but might give a *different* error-budget profile when extrapolated to next-generation precision tests. *Whether* this differs from standard QED's error budget is undetermined.
- **Why this is not a prediction:** Gill has not done the next-order calculation. The conjecture is "the dual organisation might give a different error budget", not "Gill predicts Y for QC verification".

### 5.4 Question 4 (g-factor as a QC test): could the Finding 2 correction be tested?

- **Extrapolation source:** [[FINDINGS_for_author_review|Finding 2 — DRQM I Eq. (III.22)]].
- **Speculative bearing:** `#speculative` — most concrete experimental hook. If a high-precision superconducting qubit array were designed to measure the electron g-factor in a novel regime — say, via a Penning-trap-style measurement embedded in a circuit-QED system — the precision could potentially distinguish standard QED's $g_e$ from a dual-framework prediction using the corrected $r_e$ from Finding 2. Standard $g_e \approx -2.00231930$. Dual $g_e^{\rm dual}$ with corrected $r_e$ matches by construction. The differential test would be at extreme precision, beyond what current QC hardware achieves.
- **Why this is not a prediction:** the dual-framework $g_e$ with corrected $r_e$ matches standard at the precision the formula was constructed to. Any test at a precision beyond that would require independent derivation of higher-order terms in the dual framework — not yet done.

### 5.5 Question 5 (topological QC): Gill-silent

- **Extrapolation source:** none.
- **Speculative bearing:** `#gill-silent`. Anyon braiding statistics and 2D topological field theory are entirely outside Gill's published scope. No commentary.

## 6. Experimental tests that could distinguish frameworks

Honest section. No experiment currently exists that would distinguish standard from dual-framework predictions at QC-relevant precision. The closest candidates:
- **Penning-trap-style g-factor measurement in circuit QED:** could potentially probe Finding 2's corrected $r_e$ if precision reaches 13+ digits (currently 13 at Harvard via Hanneke 2008).
- **Circuit-QED dispersion at $5$ GHz to $10^{-15}$ resonance fractional accuracy:** could in principle see the dual-framework effective photon mass; not currently achievable.

Most likely: *no* current experiment distinguishes. The framework's interesting bearing on QC is conceptual — alternative organisation of decoherence theory, alternative perturbative expansion that avoids Dyson divergence — rather than predictive.

## 7. Bibliography

### Modern reviews

- [[preskill2018_nisq]] — NISQ-era landscape.
- [[nielsen_chuang2000_quantum_computation]] — standard graduate textbook (from Ch 7).
- [[krantz2019_circuit_qed]] — superconducting-qubit engineering reference.
- [[pirandola2020_qkd_review]] — QKD landscape (from Ch 7), relevant for §5.2.

### Foundational sources

- [[feynman1982_simulating_physics]] — vision statement.
- [[deutsch1985_quantum_turing]] — universal quantum computer.
- [[shor1994_factoring]] — Shor's algorithm.
- [[grover1996_search]] — Grover's algorithm.
- [[kitaev2003_topological]] — topological QC.
- [[arute2019_quantum_supremacy]] — Google Sycamore demo.

### Historical roots (cross-references)

- [[bohr1913_atom_constitution]], [[heisenberg1925_quantum_mechanics]], [[dirac1928_electron]], [[bcs1957_superconductivity]], [[josephson1962_tunneling]], [[lamb_retherford1947_microwave]].

## 8. Cross-references

- Historical roots: Chapters [[03_old_quantum_theory_1900_1925]], [[04_quantum_mechanics_1925_1948]], [[05_QED_renormalization_solid_state_1948_1965]].
- Verification anchors: [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]], [[Feynman_Operator_Calculus_Papers]], [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], [[FINDINGS_for_author_review]].
- Companion forward chapters: [[07_PNT_GPS_SLR_QKD]], [[09_fusion_open_questions]].
