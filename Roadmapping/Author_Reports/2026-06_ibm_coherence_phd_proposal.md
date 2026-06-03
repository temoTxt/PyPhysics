---
title: "A proper-time dissipative term in driven superconducting qubits: a dissertation research proposal"
author: "Trey Morris"
date: "2026-06-03"
subject: "PhD dissertation research proposal expanding the IBM-quantum-coherence white paper (issue #103) into a falsifiable, multi-year program testing whether Gill's dual / proper-time electrodynamics predicts a measurable angle- and shape-dependent dephasing signature in driven transmon qubits."
---

# A proper-time dissipative term in driven superconducting qubits

## Abstract

Decoherence sets the present ceiling on superconducting quantum computation. The dominant theoretical account of transmon dephasing rests on the standard Maxwell field coupled to environmental noise channels. We propose to test an alternative. Gill's dual, proper-time reformulation of electrodynamics is mathematically equivalent to Maxwell's theory under a change of variables, yet it carries one structural term that the textbook theory does not[^1]. The term is a first-derivative dissipative contribution whose coefficient is the proper-time scalar $(\mathbf{u}\!\cdot\!\mathbf{a})/b^4$. It vanishes for an unaccelerated source, and it inherits an angular dependence through $\mathbf{u}\!\cdot\!\mathbf{a}=|\mathbf{u}||\mathbf{a}|\cos\theta$. A transmon qubit driven by a shaped microwave pulse is an accelerated source by construction. We therefore ask whether the dual term contributes a dephasing channel with a drive-angle and pulse-shape signature distinct from standard radiation reaction. The proposal develops the test in three stages: a derivation of the term's contribution to the qubit master equation, a concrete prediction for the angle and shape dependence of $T_2$, and an experiment on IBM Quantum hardware. We state at the outset that the most likely outcome is a null result, and we treat the honest specification of a falsifiable test as the program's primary contribution.

<!-- TODO: human reviews and fills in — Trey to confirm the abstract's framing, in particular that "null result most likely" is the intended public posture for a dissertation proposal (substantive AI framing choice; the white paper, issue #103, adopted the same posture for internal planning). -->

## 1. Background and motivation

The transmon is the workhorse of contemporary superconducting quantum processors[^2]. Its coherence time $T_2$ limits the depth of any circuit executed before error correction is required. Decoherence in these devices is well characterised experimentally and is attributed to a small set of channels: flux noise, charge noise, dielectric loss, and quasiparticle tunnelling[^3]. Each channel is modelled within standard electrodynamics coupled to an environmental spectral density. The standard account is successful, and we do not dispute it.

We observe, however, that the standard account assumes the textbook Maxwell field. Gill's proper-time formulation is an alternative to that field theory[^1]. The two formulations are mathematically equivalent under the substitutions $c\to b$, $t\to\tau$, and $\mathbf{w}\to\mathbf{u}$, but they are not, in the sense Gill makes precise, physically equivalent. The proper-time wave equation carries one term with no textbook counterpart. To the extent that this term couples to the driven condensate of a transmon, it predicts a dephasing contribution that the standard account omits. The motivation for this proposal is that the framework has never been tested against a condensed-matter device, and a driven qubit is the cleanest accelerated electromagnetic source available in the laboratory.

It is worth being precise about where a deviation could hide. The standard channels are characterised in the idle qubit, where the source is unaccelerated[^3]. Gate operations are different. During a gate the drive accelerates the condensate on a timescale of tens of nanoseconds, and it is exactly the accelerated regime in which the term of Eq. (2.2) is nonzero. The standard noise budget is therefore measured in the regime where the dual term vanishes by construction, and the gate regime where it does not vanish is the regime least constrained by existing characterisation. This is the gap the proposal targets.

## 2. Problem statement and research questions

The dual wave equation, written for the electric field in the proper-time formulation, takes the form

$$
\nabla^2 \mathbf{E} \;-\; \frac{1}{b^2}\,\partial_\tau^2 \mathbf{E} \;-\; \frac{\mathbf{u}\!\cdot\!\mathbf{a}}{b^4}\,\partial_\tau \mathbf{E} \;=\; \mathcal{S}, \qquad (2.1)
$$

where $\mathcal{S}$ is the source term and $b$ is the proper-time speed parameter[^1]. The first two terms reproduce the textbook wave equation under the substitutions above. The third term is the structural deviation. Its coefficient is the dissipative scalar

$$
\Gamma_{\!d} \;\equiv\; \frac{\mathbf{u}\!\cdot\!\mathbf{a}}{b^4} \;=\; \frac{|\mathbf{u}|\,|\mathbf{a}|\cos\theta}{b^4}, \qquad (2.2)
$$

in which $\theta$ is the angle between the source velocity and its acceleration. Equation (2.2) is the object of this proposal. It vanishes when $\mathbf{a}=\mathbf{0}$, and it carries an explicit angular dependence through $\cos\theta$.

The research questions follow directly. First, what kinematic quantities of the driven Cooper-pair condensate play the role of $\mathbf{u}$ and $\mathbf{a}$ in Eq. (2.2)? Second, does the term contribute a dephasing rate large enough to measure at present $T_2$ values? Third, can the contribution's angular and shape dependence be separated experimentally from the standard radiation-reaction term, which depends on $\dot{\mathbf{a}}$ rather than on $\theta$? These three questions define the three phases of the program.

## 3. Hypotheses and predictions

We advance one principal hypothesis and one falsifier. The principal hypothesis is that the dissipative term of Eq. (2.2) contributes an additional pure-dephasing channel to the transmon master equation during gate operations. Written as a Lindblad dephasing rate, the proposed contribution is

$$
\Gamma_\varphi^{(d)}(\theta, A) \;=\; \kappa \int_0^{t_g} \frac{\mathbf{u}(\tau)\!\cdot\!\mathbf{a}(\tau)}{b^4}\;\big|A(\tau)\big|^2\,d\tau, \qquad (3.1)
$$

where $A(\tau)$ is the drive envelope, $t_g$ is the gate time, and $\kappa$ is a coupling constant to be fixed in Phase A. Equation (3.1) is a proposed form, not a derived result. The derivation of $\kappa$ and of the correct integrand is the work of Phase A.

The falsifier is sharper. Standard Abraham–Lorentz–Dirac radiation reaction depends on $\dot{\mathbf{a}}$ and carries no explicit $\cos\theta$ dependence at fixed $\dot{\mathbf{a}}$[^4]. The dual term depends on $\mathbf{u}\!\cdot\!\mathbf{a}$. We therefore predict that varying the drive angle $\theta$ at fixed $|\dot{\mathbf{a}}|$, fixed amplitude, and fixed gate time produces a $T_2$ variation of the form

$$
\frac{1}{T_2(\theta)} \;=\; \frac{1}{T_2^{(0)}} \;+\; C\,\cos\theta, \qquad (3.2)
$$

with $C$ a constant whose sign and magnitude Phase B will predict. A measured angular dependence of the form (3.2), not attributable to the standard channels at fixed $\dot{\mathbf{a}}$, supports the framework. Its absence at the predicted sensitivity refutes the framework for this device class.

<!-- TODO: human reviews and fills in — the Lindblad form (3.1) and the linear-in-cos(theta) prediction (3.2) are substantive AI proposals, not derivations. Trey or a superconducting-qubit theorist to confirm (a) that a pure-dephasing channel is the right master-equation slot for the term, and (b) that the cos(theta) separation from ALD at fixed adot is physically realisable in the IQ-modulation plane. -->

## 4. Identification of the condensate kinematics

The principal unknown is the identification of $\mathbf{u}$ and $\mathbf{a}$ for a Cooper-pair condensate. We propose three candidate identifications and will adjudicate among them in Phase A. The first is the supercurrent velocity set by the phase gradient,

$$
\mathbf{u}_s \;=\; \frac{\hbar}{2 m_e}\,\nabla\varphi, \qquad (4.1)
$$

where $\varphi$ is the macroscopic condensate phase[^5]. The second is the AC-Josephson velocity at the drive frequency, fixed by the Josephson relation $\partial_t\varphi = 2eV/\hbar$ and therefore tunable through the drive amplitude. The third is the time derivative of the Cooper-pair dipole moment in the transmon charge basis.

The three identifications give different effect sizes. We summarise the order-of-magnitude estimates, all taken at $b\approx c$ so that $\beta=|\mathbf{u}|/b$:

| Identification | $\beta^4$ scale | Status against present $T_2$ floor |
|---|---|---|
| Cooper-pair drift velocity | $\lesssim 10^{-30}$ | Unmeasurable. |
| Fermi velocity | $\sim 10^{-12}$ | Below floor by roughly six orders. |
| AC-Josephson velocity | drive-amplitude dependent | The only candidate that may reach the floor. |

We observe that only the AC-Josephson identification offers a route to a measurable effect, and then only at large drive amplitude. This observation already constrains the program. The honest reading is that two of the three candidates are excluded by their own effect-size estimates before any derivation begins.

<!-- TODO: human reviews and fills in — the three candidate identifications and the order-of-magnitude table are carried over from the white paper (issue #103) and are substantive AI moves. The table assumes b ~ c; if the framework prescribes a different b in a superconducting medium, the estimates shift. A condensed-matter theorist should vet which condensate velocity is the correct kinematic variable in the transmon charge basis. -->

## 5. Prior work

The transmon and its decoherence are documented in a mature literature. Koch and co-workers introduced the transmon as a charge-noise-insensitive design derived from the Cooper-pair box[^2]. Ithier and co-workers characterised decoherence in a superconducting qubit circuit and established the channel decomposition still in use[^3]. Krantz and co-workers later consolidated the engineering account in a review that we take as the standard reference for the device physics[^6]. The Josephson effect that underlies the AC-Josephson identification of Section 4 is due to Josephson[^5].

On the framework side, the dual Maxwell wave equation and its dissipative third term are recorded in the verification of Gill's electrodynamics paper[^1]. The proper-time radiation-reaction structure, which the falsifier of Section 3 must be separated from, is treated in the verification of Gill's classical-electron work[^7]. The dual relativistic quantum mechanics that supplies the proper-time Hamiltonian is recorded in its own verification document[^8]. We note that none of these framework documents addresses a condensed-matter device; the conceptual extension to a transmon is the novel content of this proposal.

The methodological precedent for the test is the strong-field radiation-reaction program. There, the experimental difficulty is not detecting energy loss but separating a velocity-and-acceleration-dependent reaction force from competing effects of similar magnitude. Recent laser–electron-beam experiments attacked exactly this separation problem[^10]. Our falsifier inherits the same logic in a quantum-device setting: the dual term and standard radiation reaction are of comparable structure, and the experimental content lies in separating them through their differing dependence on the kinematic angle. We regard the radiation-reaction experiments as the closest existing analogue to the measurement proposed here, and we will draw on their analysis methods for baseline subtraction.

The general question behind the proposal is how two field theories that are mathematically equivalent, yet physically distinct, are to be told apart in the laboratory. Gill's phrase for the relation is that the two are mathematically equivalent but not physically equivalent. An equivalence under a change of variables does not guarantee an equivalence of predictions once the variables are tied to laboratory clocks and laboratory sources. The transmon supplies a setting in which the tie is concrete and the residual term is, in principle, switched on by the drive.

## 6. Proposed methodology

The program proceeds in three phases, each with an explicit acceptance criterion and an explicit closure condition.

**Phase A — apparatus mapping.** We will derive the contribution of the dissipative term to the transmon effective master equation, starting from Eq. (2.1) and committing to one of the three identifications of Section 4. The derivation proceeds in two steps. We first reduce the dual field at the qubit through the proper-time Hamiltonian of the dual relativistic quantum mechanics[^8], so that the field term of Eq. (2.2) enters as a perturbation of the drive Hamiltonian. We then trace out the field to obtain the induced dissipator, identifying the coefficient $\kappa$ of Eq. (3.1). Each algebraic step is checked symbolically in the Wolfram language, following the campaign's verification discipline; the checks are written as single-line evaluations so that the transport does not silently break a multi-line expression. The numerical baseline will be computed with `qiskit-dynamics`, whose `Solver` integrates the Lindblad equation under a time-dependent drive envelope entered as a `Signal`. Acceptance requires a closed-form or tractable expression for $\Gamma_\varphi^{(d)}(\theta, A)$ with explicit $b\to c$ and $\beta\to 0$ limits, together with a numerical $\Delta T_2/T_2$ at IBM Heron nominal drive parameters. The closure condition is honest and explicit. If $\Delta T_2/T_2 < 10^{-3}$ at the highest computable order, or if the contribution factorises into the standard radiation-reaction account with no residual angle dependence separable from $\dot{\mathbf{a}}$, the program closes at Phase A.

**Phase B — concrete predictions.** Conditional on a nonzero Phase-A residual, we will specify drive angles and pulse-shape envelopes that maximise and minimise the dissipative contribution at fixed drive amplitude. Acceptance requires at least one prediction in the form of Eq. (3.2), with a definite sign and magnitude for $C$, compared against state-of-the-art transmon $T_2$ of order $100$ to $300~\mu\text{s}$ and against the angular resolution of present IBM pulse control. The predictions will be generated and stress-tested in simulation before any hardware time is requested.

**Phase C — IBM Quantum test.** We will pre-flight the experiment with `qiskit-aer`, augmenting a backend noise model with a custom Lindblad channel carrying the Phase-A operator. On-hardware execution will use a single high-coherence transmon on Eagle- or Heron-class hardware through `qiskit-ibm-runtime`, under the `T2Hahn` and `T2Ramsey` protocols of `qiskit-experiments`. The discriminating observable is the angular dependence at fixed dynamical-decoupling sequence, measured as a residual after subtracting matched-baseline decoupling runs. We adopt this residual construction because CPMG and Uhrig sequences suppress low-frequency dephasing and would otherwise mask the framework's contribution. The verdict scheme is threefold. **Pass** — the framework predicts the measured angle and shape residual within the campaign's precision. **Marginal** — the framework matches at the precision floor only, consistent but not discriminating. **Refuted** — the predicted signature is absent at the predicted sensitivity.

The angular sweep fixes the statistical design. We will sample the drive angle at a set of values spanning $0$ to $\pi$, so that a $\cos\theta$ residual of the form (3.2) is fit against the constant and against any even-in-$\theta$ systematic. The detectable coefficient $C$ scales with the inverse square root of the total shot count at each angle, and with the inverse of the per-point $T_2$ uncertainty. The Phase-B prediction for $C$ therefore sets the required shot budget before any hardware time is requested, and a $C$ below the budgeted sensitivity returns a Marginal verdict rather than a Refutation. We will pre-register the angle set, the shot budget, and the residual-subtraction procedure, so that the verdict is fixed by the prediction and not by the data.

<!-- TODO: human reviews and fills in — the three-phase granularity, the qiskit toolchain choices, and the residual-after-matched-decoupling observable are substantive AI methodology choices inherited from the white paper. A superconducting-qubit experimentalist should confirm the T2Hahn/T2Ramsey + matched-DD residual is the right protocol, and whether pulse-level access (qiskit.pulse) is required or gate-level DRAG sweeps suffice. -->

## 7. Preliminary results

The white paper that seeds this proposal establishes three results on which the program builds[^9]. First, the dissipative term of Eq. (2.2) is a verified structural feature of the dual formulation, not a conjecture; it is recorded as Eq. (4) of the Maxwell verification and recurs as the third term of the modified Liénard–Wiechert fields[^1]. Second, the order-of-magnitude estimates of Section 4 are in hand, and they already exclude two of the three candidate identifications. Third, the simulation toolchain is identified end-to-end, so that every phase has a laptop-scale analogue before any IBM credit is consumed.

We are careful about what these preliminary results do and do not establish. They establish that the term exists and that a test is specifiable. They do not establish that the term contributes measurably to transmon dephasing; that is precisely the open question of Phase A.

The precision-atomic campaign supplies the reason to look at a condensed-matter device at all. In the bound-electron g-factor and the hydrogenic spectra, the framework reproduces the textbook quantum-electrodynamic result by construction; the cutoff carries no nuclear-charge dependence, and the framework inherits the bound-state structure of standard theory rather than predicting a deviation from it[^11]. A test that the framework passes by construction cannot discriminate. The driven condensate is attractive precisely because the dissipative term is switched on by acceleration, a regime the atomic campaign never probes. The proposal is therefore the campaign's first attempt at a setting where the framework could, in principle, depart from standard theory rather than reproduce it.

## 8. Timeline and milestones

We propose a four-year program. Year one is devoted to Phase A: the symbolic derivation, the Wolfram verification, and the `qiskit-dynamics` baseline, ending in the Phase-A acceptance or closure decision. Year two, conditional on a nonzero residual, develops the Phase-B predictions and the full `qiskit-aer` simulation sweep over drive angle and pulse shape. Year three executes Phase C on IBM Quantum hardware and performs the residual analysis against matched-decoupling baselines. Year four is reserved for analysis, replication, and dissertation writing. We note that the Phase-A closure condition may terminate the program at the end of year one; in that case the dissertation documents the closure as the framework's terminal verdict for this device class, and the remaining years redirect to a fallback question identified in Section 10.

The milestones and their deliverables are as follows.

| Year | Phase | Milestone deliverable | Gate |
|---|---|---|---|
| 1 | A | Wolfram-verified expression for $\Gamma_\varphi^{(d)}(\theta,A)$; `qiskit-dynamics` baseline $\Delta T_2/T_2$. | Accept if residual $\geq 10^{-3}$; else close. |
| 2 | B | Pre-registered $\cos\theta$ prediction (3.2) with signed $C$; full `qiskit-aer` angle/shape sweep. | Accept if $C$ exceeds budgeted sensitivity. |
| 3 | C | Hardware angle sweep on Eagle/Heron; residual after matched decoupling. | Pass / Marginal / Refuted per Section 6. |
| 4 | — | Replication, analysis, dissertation. | Defence. |

The gate column makes each year's continuation conditional on the prior year's result, so that a null finding closes the program cleanly rather than leaving it open-ended.

## 9. Expected contributions and broader impact

The program's primary contribution is the specification of a falsifiable test that does not presently exist in the framework's repertoire. At present no condensed-matter prediction of the dual framework is on the books. We regard the test as worth specifying whether or not the framework passes it. A secondary contribution, independent of the framework's fate, is a reusable `qiskit-dynamics` and `qiskit-experiments` pipeline for inserting a custom Lindblad channel and measuring its angular signature as a decoupling residual; this pipeline is of use to any program that posits a non-standard dephasing channel. The broader impact is methodological. A clean experimental adjudication of a framework that is mathematically equivalent to Maxwell's theory, yet physically distinct, would sharpen the general question of how such equivalences are to be tested.

<!-- TODO: human reviews and fills in — the contribution claims are substantive. Trey to confirm the "reusable pipeline" secondary contribution is real and not overstated, and that the broader-impact paragraph does not drift into hyperbole forbidden by the voice standard. -->

## 10. Risks, feasibility, and alternatives

The principal risk is that the effect is unmeasurably small. The order-of-magnitude table of Section 4 makes this risk explicit, and the Phase-A closure condition converts it into a clean decision rather than an open-ended search. A second risk is that the dual contribution is not separable from standard radiation reaction; the falsifier of Eq. (3.2) is constructed precisely to test separability, and a negative separability result is itself a publishable closure. A third risk is access: pulse-level control on IBM hardware may not be provisioned, in which case gate-level DRAG-coefficient sweeps provide a coarser but adequate angular handle.

The fallback question, should Phase A close the primary program, is whether the same dissipative term contributes to the coherence of trapped-ion or neutral-atom qubits, where the accelerated source is a single charged particle rather than a condensate and the identification of $\mathbf{u}$ and $\mathbf{a}$ is unambiguous. This fallback inherits the entire methodology of Sections 6 through 8 and requires only a change of platform. We regard the fallback as a strength of the proposal, because the central derivation is platform-independent.

<!-- TODO: human reviews and fills in — the trapped-ion / neutral-atom fallback is a substantive AI proposal added for this dissertation-scale document (it is not in the white paper). Trey to confirm the fallback is well-posed and that the dissertation scope should include it. -->

## 11. Resources required

The program is modest in physical resources and concentrated in analytical effort. Years one and two require only commodity computing for the Wolfram verification and the `qiskit-dynamics` and `qiskit-aer` simulations, both of which run at laptop scale. Year three requires IBM Quantum Network access to a single Eagle- or Heron-class device, with pulse-level control through `qiskit.pulse` if it can be provisioned and gate-level DRAG-coefficient sweeps as the fallback. The shot budget is fixed by the Phase-B prediction for the coefficient $C$, and credit is requested only after that prediction is in hand. The advisory requirement is the one genuine dependency: the apparatus-mapping of Section 4 needs a superconducting-qubit theorist to confirm the condensate kinematics, and the framework-internal derivation of Phase A benefits from contact with the framework's authors. We treat both as collaborations to be secured in year one, not as open risks to be discovered later.

<!-- TODO: human reviews and fills in — Trey to confirm the resource scope, in particular whether IBM Quantum Network access is already available through an existing affiliation, and which superconducting-qubit theorist is the intended Phase-A collaborator. -->

## References

[^1]: `Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`. Verification of Gill's dual Maxwell formulation; the dissipative third term with coefficient $(\mathbf{u}\!\cdot\!\mathbf{a})/b^4$ is recorded as Eq. (4) and recurs in the modified Liénard–Wiechert fields.
[^2]: `koch2007_transmon`. Koch et al., "Charge-insensitive qubit design derived from the Cooper pair box," Phys. Rev. A 76, 042319 (2007).
[^3]: `ithier2005_decoherence`. Ithier et al., "Decoherence in a superconducting quantum bit circuit," Phys. Rev. B 72, 134519 (2005).
[^4]: `Roadmapping/Equation_Verification/The_Classical_Electron_Problem.md`. Verification of Gill, Zachary & Lindesay; proper-time treatment of radiation reaction, against which the dual term must be separated.
[^5]: `josephson1962_tunneling`. Josephson, "Possible new effects in superconductive tunnelling," Phys. Lett. 1, 251 (1962).
[^6]: `krantz2019_circuit_qed`. Krantz et al., "A quantum engineer's guide to superconducting qubits," Appl. Phys. Rev. 6, 021318 (2019).
[^7]: `Roadmapping/Equation_Verification/The_Classical_Electron_Problem.md`. As [^4]; the proper-time radiation-reaction structure.
[^8]: `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md`. Verification of DRQM I; the proper-time Hamiltonian and dual Dirac apparatus.
[^9]: `Roadmapping/Author_Reports/2026-05_quantum_coherence_whitepaper.md`. The seed white paper; issue #103, github.com/temoTxt/PyPhysics/issues/103.
[^10]: Issues #43 and #48, the proper-time radiation-reaction and bremsstrahlung prediction threads (github.com/temoTxt/PyPhysics/issues/43, github.com/temoTxt/PyPhysics/issues/48), which collect the strong-field radiation-reaction experiments and the framework's treatment of them.
[^11]: Issues #78 and #82, the Li$^{2+}$ Z-extension and hydrogenic g-factor Z-scan (github.com/temoTxt/PyPhysics/issues/78, github.com/temoTxt/PyPhysics/issues/82); the Z-scan returned Outcome C, in which the cutoff carries no Z-dependence and the framework inherits the bound-state quantum-electrodynamic structure.
