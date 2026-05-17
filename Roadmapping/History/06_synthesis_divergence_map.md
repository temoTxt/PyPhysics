---
chapter: 06
title: "Synthesis: where standard and proper-time agree, where they diverge"
era: "1800-1965"
threads: [electromagnetism, quantum, solid-state]
animations: [synthesis_tour]
verification_anchors: ["Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations", "Dual_Relativistic_Quantum_Mechanics_I", "FOUNDATIONS_FOR_QED_I_MATHEMATICAL", "Feynman_Operator_Calculus_Papers", "FINDINGS_for_author_review"]
status: draft
---

# Chapter 6 — Synthesis: where standard and proper-time agree, where they diverge

> **Framing principle (load-bearing).** We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics. This synthesis chapter draws the divergence map across all five historical eras: where the dual framework reproduces standard predictions identically (most cases), where it produces numerical corrections below current measurement precision (`#inferred` cases), where the dual framework simply hasn't published an alternative (`#gill-silent` cases), and where its reframing dissolves a known structural difficulty without changing the predictions (the three load-bearing interpretive wins).

## 1. Overview

The chapter is a **map**, not a discovery. Chapters 1–5 walked through 165 years of physics one era at a time, with each era carrying its own dual-framework annotations and confidence-tier tags. This chapter abstracts across all five, producing four classifications that future readers (and the original authors of the cited Gill papers) can use to navigate the campaign:

1. **Where standard and dual predictions coincide identically** — by construction. Most steady-state, low-velocity, lab-frame predictions fall in this category.
2. **Where standard and dual predictions diverge but the divergence is below current measurement precision** — `#inferred` claims with explicit precision thresholds.
3. **Where the dual framework has no published alternative** — `#gill-silent` topics, including all of GR and most of the strong + weak interactions.
4. **Where the dual framework reframes a known structural difficulty without changing predictions** — the three load-bearing interpretive wins (Maxwell's Michelson–Morley reframing; DRQM I's positive-definite Hamiltonian dissolving Klein–Gordon's negative probability; Foundations I + FOC's reframing of the Dyson divergence conjecture).

Forward chapters: [[07_PNT_GPS_SLR_QKD]] for the derivational PNT chapter; [[08_quantum_computing_open_questions]] and [[09_fusion_open_questions]] for the speculative forward chapters.

## 2. Where standard and dual predictions coincide identically

Most of the historical record. By construction: when the source is at rest in the lab frame, $u = 0$ and $b = \sqrt{c^2 + u^2} = c$ exactly, so the dual Maxwell, dual Dirac, and dual classical Hamiltonian *reduce to* their standard counterparts identically. The two framings give numerically indistinguishable predictions for:

| Chapter | Result | Why predictions coincide |
|---|---|---|
| 1 (1800–1860) | Volta, Ørsted, Ampère, Biot–Savart, Ohm, Faraday, Henry, Lenz, Gauss, Kirchhoff | All quasi-static; sources at rest; $u \ll c$; six orders of magnitude below 1850s-era detection thresholds. |
| 2 (1860–1900) | Maxwell vacuum wave propagation; Hertz wave detection; Thomson $e/m$ | Vacuum propagation: dual and standard wave equations are identical. Hertz: source acceleration negligible compared with $c$. Thomson: $v/c \lesssim 0.1$, dual correction within his ~5% precision. |
| 3 (1900–1925) | Bohr hydrogen spectrum; Compton (free electron) | Bohr orbits at $v/c \sim Z\alpha \approx 7\times 10^{-3}$; $b/c$ correction $\sim 10^{-9}$, well below 1913 spectroscopic precision. Free-electron Compton at rest: $u=0$, formulas identical. |
| 4 (1925–1948) | Schrödinger non-relativistic hydrogen; Dirac fine structure | Both equations have dual-framework analogs that reproduce identical spectra. |
| 5 (1948–1965) | Lamb shift; anomalous magnetic moment $\alpha/(2\pi)$ | Dual-framework one-loop QED reproduces the standard numerical results at the same order. |

## 3. Where standard and dual predictions diverge, below current precision

`#inferred` claims with explicit precision regimes. The dual framework predicts numerical departures from standard physics in these cases — but the departures are smaller than any current experiment can detect.

| Chapter | Claim | Standard prediction | Dual prediction | Divergence | Current precision | Distinguishable? |
|---|---|---|---|---|---|---|
| 2 | Michelson–Morley fringe shift | 0 (SR + length contraction) | 0 (velocity duality, no length contraction) | identical at any current precision | — | no |
| 3 | Bohr hydrogen spectrum | $E_n^{\rm std}$ | $E_n^{\rm dual} = E_n^{\rm std}[1 + O(v^4/c^4)]$ | ~$10^{-9}$ relative | $10^{-14}$ (modern hydrogen spectroscopy) | structurally yes; QED radiative corrections dominate at this precision in both framings |
| 3 | Compton scattering, bound electron at $u/c \sim 10^{-2}$ | $\Delta\lambda^{\rm std}$ | $\Delta\lambda^{\rm dual} = \Delta\lambda^{\rm std}\bigl[1 + O(u^2/c^2)\bigr]$ | ~$10^{-4}$ relative | $10^{-3}$ (Compton-profile measurements) | barely; below routine apparatus precision |
| 7 (forward) | GPS satellite clock correction | $-7\,\mu$s/day SR + $+45\,\mu$s/day GR | dual SR + GR (Gill-silent on GR) | $\sim 10^{-10}$ at orbital $u/c \sim 10^{-5}$ | $\sim 10^{-9}$ (current GPS clock precision) | not at routine GPS precision; possibly at future atomic-clock precision |
| 9 (forward) | Plasma cross-sections at $u/c \sim 10^{-2}$ | standard Coulomb | dual $b$-modified Coulomb | $\sim 10^{-4}$ | $\sim 10^{-2}$ (current ICF measurements) | not currently |

## 4. Where the dual framework has no published alternative — `#gill-silent`

Topics on which Gill has not published a proper-time / dual-framework treatment. The standard predictions are taken intact:

- **General relativity** (Mercury perihelion, gravitational time dilation, Schwarzschild, Shapiro delay, Hulse–Taylor binary, GR-based GPS clock correction, LIGO gravitational waves). `#gill-silent` for the entire GR programme. Most of [[07_PNT_GPS_SLR_QKD|Chapter 7's]] GPS treatment relies on standard GR results.

- **Quantum chromodynamics** (asymptotic freedom, confinement, hadron mass spectrum, lattice QCD). `#gill-silent`. The strong-interaction programme is entirely outside Gill's published scope.

- **Electroweak gauge theory** (W/Z bosons, Higgs mechanism, neutrino oscillations). `#gill-silent` for the gauge structure itself, though Fermi-type four-fermion contact interactions (the V–A weak Lagrangian) can be ported to the dual Dirac framework.

- **Cosmology** (Big Bang, inflation, CMB, dark matter, dark energy, structure formation). `#gill-silent`.

- **Solid-state mechanism** (Bloch band theory, BCS Cooper-pairing, Josephson junction phase coherence, transistor band-gap physics, magnetism). `#gill-silent`. [[08_quantum_computing_open_questions|Chapter 8]] proposes one *speculative* connection (effective photon mass and circuit-QED coherence) but tags it `#speculative` rather than `#inferred`.

- **Statistical mechanics** (Boltzmann, Gibbs, fluctuation-dissipation, phase transitions). `#gill-silent`, with one exception: [[Relativistic_Transformations_of_Thermodynamics]] does treat relativistic thermodynamics with dual-K specialisations.

## 5. The three load-bearing interpretive wins

The cases where the dual framework reframes a known structural difficulty *without changing the numerical predictions*. These are the campaign's central conceptual payoffs.

### 5.1 Michelson–Morley without length contraction

[[02_classical_synthesis_1860_1900|Chapter 2]] §3. The standard SR resolution of the 1887 null result requires length contraction as a kinematic postulate. The dual framework's velocity duality $\mathbf{w}/c = \mathbf{u}/b$ ([[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] Eq. (1)) reaches the same null prediction without a separate length-contraction postulate — the velocity-duality structure already encodes whatever frame-dependent kinematics the interferometer needs. **Predictions identical. Interpretation simpler.**

### 5.2 Klein–Gordon's negative-probability problem dissolved

[[04_quantum_mechanics_1925_1948|Chapter 4]] §3. The standard Klein–Gordon equation has a conserved current $j^0$ that is not positive-definite, making it unsuitable as a single-particle probability density. Standard QM resolves this by going to Dirac (a first-order-in-time equation with 4-spinor structure) and absorbing the negative-energy modes into hole theory + antimatter. The dual framework's positive-definite Hamiltonian $K = H^2/(2mc^2) + mc^2/2$ ([[FoundationsII-Classical]]) plus the square-root operator construction ([[Analytic_Representation_of_The_Square-Root_Operator]]) gives a single-particle relativistic theory with manifestly positive-definite probability density — without spinor doubling. Antimatter is then organised via Santilli isodual companion equations rather than as Pauli-filled negative-energy seas. **Predictions identical. Mathematical organisation cleaner.**

### 5.3 Dyson divergence conjecture reframed

[[05_QED_renormalization_solid_state_1948_1965|Chapter 5]] §3. [[dyson1952_divergence|Dyson 1952]] argues, on heuristic grounds, that the QED perturbation series in $\alpha$ is asymptotic rather than convergent. In practice this doesn't prevent 12-digit QED computations, but it does mean the series has no convergent resummation. Gill's [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL|Foundations I Mathematical]] organises the time-ordered operator algebra on a Kuelbs–Steadman Hilbert space; [[Feynman_Operator_Calculus_Papers|the Feynman Operator Calculus]] provides the perturbative scaffolding. In that representation the resummation problem is dissolved structurally; the Dyson divergence is a feature of one specific representation, not a barrier in principle. **Predictions identical (12+ digit agreement on $g-2$, Lamb shift, etc.). Structural soundness improved.**

## 6. The three flagged findings — [[FINDINGS_for_author_review]]

Three places where the verification campaign found *internal inconsistencies* — not predictions of new physics, but algebra that didn't reproduce — in cited Gill papers. The findings are reproducible from the recorded Wolfram MCP inputs, and they are the campaign's most concrete deliverable for the authors:

- **Finding 1.** [[02_classical_synthesis_1860_1900|Ch 2]] / [[maxwell1865_dynamical_theory|Maxwell paper anchor]]. *Eq. (24)* in Gill's *Two Mathematically Equivalent Versions of Maxwell's Equations* has a missing factor of $c$ in the term $e\hbar \boldsymbol{\Sigma}\cdot\mathbf{B}/(2m)$ (should be $/(2mc)$) and a missing $+V^2/(2mc^2)$ term. Spin-orbit and Pauli-relativistic-energy fix.

- **Finding 2** ⭐ **headline payoff.** [[04_quantum_mechanics_1925_1948|Ch 4]] / [[Dual_Relativistic_Quantum_Mechanics_I|DRQM I]]. *Eq. (III.22)*: the published $r_e \approx 0.499857150068631\,r_0$ gives $g \approx -2.0005714$, not the experimental $-2.00231930436256$. The correct $r_e$ is $\approx 0.499420510\,r_0$. Numerical correction; algebraic chain through Sec III is correct.

- **Finding 3.** [[05_QED_renormalization_solid_state_1948_1965|Ch 5]] / [[The_Classical_Electron_Problem|TCEP]]. *Eq. (4.16)*: the paper writes $v_g = v_g' - v$; algebra and the paper's own commentary give $v_g = v_g' + v$. Sign typo.

All three are flagged for author review (Trey is co-author of DRQM I — Finding 2 — and the verification work is framed accordingly).

## 7. Animations referenced

| Manim scene | Status | Role in the synthesis |
|---|---|---|
| [`synthesis_tour.py`](../Animations/manim_scenes/synthesis_tour.py) | rendered (PR #6) | The central-identity → three-pillars → three-findings tour that motivates the whole campaign. Extends naturally to span the full 165-year arc. |
| [`hist_faraday_induction.py`](../Animations/manim_scenes/hist_faraday_induction.py) through [`hist_lamb_shift_contrast.py`](../Animations/manim_scenes/hist_lamb_shift_contrast.py) | rendered (Chs 1–5) | The era-specific scenes; cross-referenced by chapter. |
| [`drqm_eq22_g_factor_finding.py`](../Animations/manim_scenes/drqm_eq22_g_factor_finding.py) | rendered (PR #6) | **Finding 2 visualised** — the campaign's headline payoff. |

## 8. Dataview indexes

This chapter ships with auto-generated index pages produced by `_tools/build_dataview_indexes.py`. The static markdown indexes provide Dataview-like queries without requiring the Obsidian Dataview plugin:

- [`_index_by_year.md`](_index_by_year.md) — every bibliography entry sorted chronologically.
- [`_index_by_tag.md`](_index_by_tag.md) — every entry grouped by tag (era, thread, framework-claim-type).
- [`_index_inferred_claims.md`](_index_inferred_claims.md) — every `#inferred` and `#speculative` claim across chapters, with anchors.

Regeneratable from the bibliography + chapter sources at any time via `uv run python Roadmapping/History/_tools/build_dataview_indexes.py`.

## 9. Cross-references

- Backward: [[01_early_electromagnetism_1800_1860]] through [[05_QED_renormalization_solid_state_1948_1965]].
- Forward: [[07_PNT_GPS_SLR_QKD]] (derivational), [[08_quantum_computing_open_questions]] + [[09_fusion_open_questions]] (speculative).
- Verification anchors: all four of the campaign's load-bearing physics papers.
- Findings: [[FINDINGS_for_author_review]].
