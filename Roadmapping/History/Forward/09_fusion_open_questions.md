---
chapter: 09
title: "Fusion: open questions and proper-time conjectures"
era: "forward"
variant: open-questions
threads: [quantum]
animations: []
verification_anchors: ["Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations", "Relativistic_Transformations_of_Thermodynamics"]
status: draft
---

# Chapter 9 — Fusion: open questions and proper-time conjectures

> **Framing principle (load-bearing).** We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics. This chapter is the most *thinly-connected* of the campaign: fusion is mostly classical plasma physics + nuclear cross-sections, both `#gill-silent`. The dual framework's bearing is marginal — quantitatively below current ICF and tokamak precision. The chapter's payoff is the *historical-arc contextualisation*: fusion as a direct descendant of [[einstein1905_emc2|Einstein 1905]] mass-energy equivalence + [[aston1922_mass_defect|Aston's binding-energy curve]] + [[bethe1939_ppchain|Bethe's stellar-fusion cycles]] + [[lawson1957_criterion|Lawson's confinement criterion]]. Tag-dominant `#gill-silent` with `#speculative` hooks where mechanically extrapolatable.

## 1. Overview

Fusion in 2026 is at an inflection point. [[nif2022_ignition|NIF demonstrated ignition in December 2022]] — fusion energy output exceeded laser-energy input to the target, the first laboratory demonstration of fusion gain > 1. The plug-in energy required to fire the lasers was 150× the fusion gain, so this is *not* practical power generation, but the underlying inertial-confinement physics is now demonstrated. ITER is under construction in France, targeting first plasma in 2034 and Q=10 D-T plasmas in 2039. Private fusion companies — Commonwealth Fusion Systems, TAE, Helion, Tokamak Energy, Zap Energy — have collectively raised >$6 billion and are targeting net-electrical-power demonstrations in the 2030s. The path from physics demonstrations to grid-scale electricity remains long, but the technical questions are increasingly well-defined.

This chapter is the *thinnest* of the campaign in dual-framework content. Most of fusion is `#gill-silent` — the cross-sections are nuclear-physics calculations (R-matrix, S-factor expansions), the plasma physics is MHD or kinetic theory, and the confinement engineering is materials science. The dual-framework's bearing is at the level of `#speculative` corrections to ion cross-sections in dense plasmas — see §5.

## 2. Historical roots

- **[[einstein1905_emc2|Einstein 1905]] $E = mc^2$ (Ch 3):** the fundamental energy-mass equivalence that *is* the source of fusion energy.
- **[[rutherford1911_alpha_scattering|Rutherford 1911]] (Ch 3):** atoms have a nucleus; the substrate for nuclear reactions.
- **[[aston1922_mass_defect|Aston 1920]] mass-spectrometry:** helium-4 is *less massive* than two protons + two neutrons. The mass defect $\Delta m \cdot c^2$ is the binding energy.
- **[[eddington1920_stellar_fusion|Eddington 1920]] conjecture:** stars are powered by hydrogen fusion via Aston's mass defect.
- **[[chadwick1932_neutron|Chadwick 1932]] neutron (Ch 4):** the neutral nuclear constituent; opens nuclear physics experimentally.
- **[[fermi1934_beta_decay|Fermi 1934]] weak interaction (Ch 4):** rate-limiting step of the pp-chain.
- **[[bethe1939_ppchain|Bethe 1939]] stellar-fusion cycles:** quantitative pp-chain and CNO-cycle nuclear-reaction theory; Nobel 1967.
- **Teller–Ulam thermonuclear (1951, classified):** military fusion, the H-bomb.
- **[[lawson1957_criterion|Lawson 1957]] criterion:** $n\tau \gtrsim 10^{20}$ m$^{-3}\cdot$s engineering target for fusion gain.

## 3. Current state (2026 perspective)

Three distinct confinement approaches:
- **Magnetic confinement (tokamaks, stellarators).** [[wesson2011_tokamaks|Wesson]] is the standard reference. JET (UK, decommissioned 2024) holds the magnetic-confinement fusion-energy record: 69 MJ in a 5-second D-T pulse at Q=0.33. ITER (France, first plasma ~2034) targets Q=10 in steady-state D-T operation. SPARC (Commonwealth Fusion Systems, MA) targets Q>2 by ~2027 with high-temperature-superconductor magnets.
- **Inertial confinement (laser-driven implosion).** [[atzeni_meyertervehn2004_inertial|Atzeni & Meyer-ter-Vehn]] is the standard reference. NIF (Livermore CA) achieved ignition in December 2022, target gain 1.5× now (2024); the engineering path to power-generation IFE is via repetitive driver technology not yet demonstrated.
- **Alternative concepts.** Z-pinch (Zap Energy), magnetised target (General Fusion, Helion), field-reversed configuration (TAE).

## 4. Major open questions

1. **The materials problem.** Reactor walls and breeding blankets in D-T fusion experience neutron fluxes of $\sim 10^{14}$ n/cm$^2$/s, fast (14 MeV) neutrons. No material has been tested in that environment at full reactor scale. *Current status:* the IFMIF-DONES neutron source under construction in Spain will provide reactor-relevant neutron exposure starting ~2030.
2. **Tritium breeding ratio.** Tritium is rare and short-lived ($t_{1/2} = 12.3$ years); a fusion reactor must breed its own tritium via Li-6(n, $\alpha$)T in the blanket. *Current status:* breeding ratios > 1 are predicted by simulation; not yet demonstrated experimentally.
3. **Disruption avoidance and mitigation in tokamaks.** A large-tokamak disruption can deposit megajoules of energy on millisecond timescales, damaging the wall. *Current status:* machine-learning disruption-prediction has reached useful accuracy; mitigation via massive gas injection is mature; predict-and-prevent at high duty cycle is open.
4. **Repetition-rate IFE.** NIF can fire one shot per day; an IFE power plant needs $\sim 10$ shots per second. Driver, target-injection, and chamber-clearing technology are all decades from product. *Current status:* the National Academies' 2024 IFE report identifies 8 critical R&D areas.
5. **Coulomb-barrier tunneling in dense ICF plasmas.** Standard nuclear physics treats the D-T fusion cross-section via R-matrix theory based on (mostly) low-density measurements. At ICF densities ($\sim 10^{32}$ ions/cm$^3$) at peak compression, screening effects modify the effective Coulomb barrier. *Current status:* plasma-screening corrections are calculated to first order in standard QED; second-order effects below current experimental sensitivity.

## 5. Speculative proper-time implications

`#speculative` extrapolations. Honest about how marginal the connection is.

### 5.1 Question 5 (Coulomb barrier in dense plasmas): does the dual framework's `b` modify the cross-section?

- **Extrapolation source:** [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] — the collaborative speed $b = \sqrt{c^2 + u^2}$ enters the EM-mediated Coulomb interaction.
- **Speculative bearing:** `#speculative`. For ICF plasmas at peak compression, ion velocities are $u \sim 10^6$ m/s, corresponding to $u/c \sim 3 \times 10^{-3}$. The dual-framework correction to the Coulomb barrier height enters as $b/c - 1 \sim u^2/(2c^2) \sim 5 \times 10^{-6}$. **Quantitatively below** current ICF capsule-diagnostics precision (~1% on areal density, gain). The conjecture is "the dual EM framework predicts a slightly modified Coulomb barrier at ICF densities"; the predicted modification is well below current measurement precision.
- **Why this is not a prediction:** Gill has not computed the plasma-screening correction in the dual framework. The order-of-magnitude estimate $\sim 5 \times 10^{-6}$ is mechanical extrapolation.

### 5.2 Modified Lorentz force in plasma

- **Extrapolation source:** [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] Eq. (18) — the dual Lorentz force has a $V/(mcb)$ correction term.
- **Speculative bearing:** `#speculative`. For typical fusion plasma ions with $V/(mc^2) \sim 10^{-6}$ (10 keV thermal energy vs MeV rest mass), the correction is below current MHD-modeling precision by many orders of magnitude. *Quantify exactly how much below:* $V/(mc^2) \sim 10^{-6}$, multiplied by typical $b/c$ correction, gives a relative force correction of order $10^{-8}$ — vs. MHD modeling errors of order $10^{-2}$ typically. Below detection.
- **Why this is not a prediction:** Gill has not done the dual-framework plasma MHD derivation.

### 5.3 Relativistic plasma thermodynamics

- **Extrapolation source:** [[Relativistic_Transformations_of_Thermodynamics]] — Nakamura RRT with dual-theory specialisations.
- **Speculative bearing:** `#inferred` (not `#speculative` — the relevant verification anchor exists). For ICF imploded cores ($T \sim 10$ keV, $kT/mc^2 \sim 10^{-2}$ for deuterium), relativistic thermodynamic corrections are at the percent level. The dual-K thermodynamics framework would give specific numerical predictions for equation-of-state and transport coefficients. *Not yet computed by the campaign* but the framework is in place; this is the most concrete dual-framework hook in the fusion chapter.

### 5.4 Coulomb-barrier tunneling and the strong interaction

- **Extrapolation source:** none.
- **Speculative bearing:** `#gill-silent`. The strong-interaction physics that determines D-T, D-D, T-T fusion cross-sections (R-matrix theory, optical models) is entirely outside Gill's published scope. The standard nuclear-physics treatment is taken intact.

## 6. Experimental tests that could distinguish frameworks

**None currently identified.** The dual-framework corrections to fusion cross-sections, Lorentz-force confinement dynamics, and plasma thermodynamics all sit *3–6 orders of magnitude below* current experimental precision. Future precision-cross-section measurements at e.g. IFMIF-DONES might reach the regime where second-order Coulomb-screening corrections are visible; the dual-framework prediction would enter alongside (and be dominated by) standard QED plasma-screening corrections at that precision.

The chapter is honest that fusion is *not* a productive testing ground for the dual framework. The connection is *historical-arc context* — fusion as a 165-year consequence of the kinematic and quantum physics tracked through Chapters 1–5 — not a prediction.

## 7. Bibliography

### Primary

- [[aston1922_mass_defect]] — mass defect.
- [[eddington1920_stellar_fusion]] — stars are powered by fusion.
- [[bethe1939_ppchain]] — pp-chain + CNO cycle.
- [[lawson1957_criterion]] — confinement engineering criterion.
- [[nif2022_ignition]] — NIF ignition demonstration.

### Retrospective

- [[wesson2011_tokamaks]] — standard tokamak reference.
- [[atzeni_meyertervehn2004_inertial]] — standard ICF reference.
- [[jackson1998_classical_electrodynamics]] — Jackson Ch 8–12 for classical-plasma EM.

### Historical roots (cross-references)

- [[einstein1905_emc2]], [[rutherford1911_alpha_scattering]], [[chadwick1932_neutron]], [[fermi1934_beta_decay]].

## 8. Cross-references

- Historical roots: Chapters [[03_old_quantum_theory_1900_1925]], [[04_quantum_mechanics_1925_1948]].
- Verification anchors: [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] (Eq. 18 modified Lorentz force; collaborative speed $b$ in plasma), [[Relativistic_Transformations_of_Thermodynamics]] (dual-K thermodynamics — relevant for ICF imploded cores).
- Companion forward chapters: [[07_PNT_GPS_SLR_QKD]], [[08_quantum_computing_open_questions]].
