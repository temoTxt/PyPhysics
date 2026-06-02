---
title: "GPS satellite relativistic corrections — standard derivation, proper-time companion, and open questions for a curved-spacetime extension"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-05-27"
---

# §1 Cover note

This is a short follow-up to the interim report of 2026-05-25[^1]. We narrow it to one applied-relativity question that opened in the public repository the same week, namely the relativistic clock-correction problem for GPS satellites. The question was prompted by a popular-science article on the $\pm 38\,\mu$s/day offset[^2]. Two campaigns followed. The first[^3] is a standard-method derivation reproducing Ashby's *Living Reviews* treatment[^4] end-to-end. The second[^5] is a proper-time companion that applies the Gill–Zachary substitution rules[^6] to each of the nine effects. The companion campaign surfaced four open questions about how the framework should extend to curved spacetime. The framework as-published does not yet answer them. This report consolidates the two campaigns and the four questions in a form short enough to read alongside the broader interim report.

<!-- TODO: human reviews and fills in — confirms the cover-note framing as a focused follow-up (not solicited, not presuming engagement beyond the prior interim report). -->

## What this report is and is not

We intend this as a *consolidation* of two campaigns whose underlying work is committed to the public repository in full. It is not a re-derivation. It is not a journal-style summary. It is not a request that you read the per-effect documents to interpret the questions in §5. Each question is self-contained. The load-bearing section is §5. Everything else is context to make §5 readable in one sitting.

# §2 The GPS question

GPS satellites carry atomic clocks. Once in orbit, those clocks tick faster than ground clocks by

$$+38.57\,\mu\text{s/day}, \qquad (1)$$

a combined gain decomposing into a gravitational-redshift contribution of $+45.65\,\mu$s/day and a special-relativistic velocity-dilation contribution of $-7.11\,\mu$s/day. The engineering response is to slow the satellite clock by exactly this fractional rate before launch. We observe that the fundamental clock frequency of $10.230\,000\,000\,00$ MHz is built at $10.229\,999\,995\,43$ MHz, a fractional offset

$$\Delta f / f = -4.4647 \times 10^{-10}. \qquad (2)$$

Once in orbit, the relativistic gain exactly cancels the pre-launch slowdown. The satellite clock then matches a ground clock at the geoid potential. This is the most stringently-tested applied-relativity result available. Indeed, the satellite constellation has operated since 1977 with the relativistic correction validated daily against the $\sim 10^{-10}$ precision of cesium and rubidium clocks. Ashby (2003) is the canonical reference[^4]; the IAU 2000 geoid potential supplies the ground-side reference, and the operational pre-launch offset is recorded in the GPS Interface Control Document[^7].

The proper-time framework's perspective is operationally aligned with how GPS actually works. The satellite's atomic clock *is* its proper time $\tau$. The question "how does $\tau_{\rm sat}$ relate to $\tau_{\rm grnd}$?" is exactly the kinematic question the framework's apparatus

$$b = \sqrt{c^2 + \mathbf{u}^2}, \qquad \frac{dt}{d\tau} = \frac{b}{c} \qquad (3)$$

is designed to answer. The campaign set out to determine, per effect, whether the framework reproduces standard SR + GR at GPS-relevant precision. The expected outcome, given the framework's claim of mathematical equivalence[^6], is that it does. The campaign also sought to surface any structural difference that would be operationally distinguishable at the next-generation optical-clock precision floor (of order $10^{-17}$), or in regimes the framework as-published does not cover.

<!-- TODO: human reviews and fills in — confirms the framing of GPS as operationally aligned with the proper-time apparatus (satellite clock IS τ); confirms the expected-outcome framing of the campaign goal so the §4 "no observable deviation" finding reads as the positive consistency check it is rather than a null result. -->

# §3 Standard derivation — what the first campaign reproduces

The standard campaign[^3] covered nine separately-named effects. Each is documented in its own per-effect file with a Wolfram-MCP numerical check and a comparison against Ashby (2003)[^4]. The combined decomposition reproduces Ashby end-to-end at the precision the published parameters can sustain.

| Effect | Standard-method magnitude | Ashby cross-reference |
|---|---|---|
| 02 Gravitational time dilation | $+45.65\,\mu$s/day | §3.1 Eq. (19) |
| 03 Velocity time dilation | $-7.11\,\mu$s/day | §3.2 Eq. (20) |
| 04 Combined offset including $J_2$ | $+38.57\,\mu$s/day | §3.3 Eq. (39) |
| Pre-launch frequency offset | $\Delta f/f = -4.4647 \times 10^{-10}$ | §3.3 Eq. (39) |
| 05 Eccentricity periodic ($e = 0.02$) | $\pm 46$ ns peak | §4 Eq. (43) |
| 06 Sagnac maximum | $\sim 137$ ns per signal | §3.4 Eq. (35) |
| 07 Shapiro maximum | $\sim 62$ ps per signal | §6 Eq. (53) |
| 08 $J_2$ contribution to geoid | $+0.033\,\mu$s/day | §3.5 Eq. (24) |

The campaign's value to the repository is *consolidation*. We collect Ashby's per-effect derivations into one place with Wolfram-MCP verification of every quoted number against named parameter values[^8]. This is a reproduction-of-standard exercise. No novel claim is made. We flag the campaign here because §4's proper-time companion is conditional on the §3 baseline being right.

<!-- TODO: human reviews and fills in — confirms the "consolidation, not novel" framing for §3 is the right honest scope; the value to the repository is the per-effect Wolfram-MCP verification, not the combined offset itself. -->

# §4 Proper-time companion — what the second campaign found

The companion campaign[^5] re-derived the same nine effects under the Gill–Zachary substitution rules[^6],

$$\frac{\mathbf{w}}{c} = \frac{\mathbf{u}}{b}, \qquad \frac{1}{c}\frac{\partial}{\partial t} = \frac{1}{b}\frac{\partial}{\partial \tau}, \qquad b^2 = c^2 + \mathbf{u}^2. \qquad (4)$$

We recorded, per effect, whether the framework applies as-published or requires a speculative extension. The nine per-effect documents split into three categories.

**Pass — direct framework application (pt_03, pt_06).** Two effects admit the verified substitution rules without modification. The most direct demonstration is pt_03 velocity time dilation. We observe that the framework's identity

$$\frac{d\tau}{dt} = \frac{c}{b} = \sqrt{1 - \mathbf{w}^2/c^2} = \frac{1}{\gamma} \qquad (5)$$

follows by algebraic identity from the velocity-duality substitution of (4). The numerical value $-7.11\,\mu$s/day is identical to the standard prediction at all decimal places, not at any precision floor. For pt_06 Sagnac, the rotating-frame photon-null condition runs through with $c \to b$ substitution. The result is

$$\Delta t = -\frac{2\boldsymbol{\omega} \cdot \mathbf{A}}{b^2}, \qquad (6)$$

which reduces to the standard $-2\boldsymbol{\omega} \cdot \mathbf{A} / c^2$ at $b \approx c$ for an Earth-rotating receiver. The residual sits at $\sim 10^{-10}$ ns, invisible to current GPS clocks and marginal at next-generation optical-clock precision.

**Marginal — kinematic Pass plus speculative GR extension (pt_01, pt_04, pt_05).** Three effects admit the kinematic piece as direct framework application but require a curved-spacetime extension that the framework as-published does not provide. The campaign adopts a *minimal-extension hypothesis*: replace $c$ by the local $b$ in the Schwarzschild line element. Under that hypothesis, the predictions reduce to standard at GPS-relevant velocities. The combined offset of (1) is reproduced. The operational pre-launch offset of (2) is unchanged. The pt_05 eccentricity case additionally surfaces a framework-specific third-term contribution,

$$\frac{\mathbf{u} \cdot \mathbf{a}}{b^4} \, \tau_{\rm signal}, \qquad (7)$$

arising from the wave-equation dissipative coefficient in Maxwell-paper Eq. (4)[^6]. For GPS circular orbits this vanishes ($\mathbf{u} \perp \mathbf{a}$). For elliptical orbits at $e = 0.02$ it is of order $10^{-34}$ s per signal, operationally invisible. We note that the same coefficient appears in the radiation-reaction sub-investigation[^9], where it contributes at order unity. The framework is thus internally consistent on which structures appear where.

**Out of scope — framework as-published does not extend (pt_02, pt_07, pt_08).** Three effects are purely gravitational: the Schwarzschild redshift, the Shapiro propagation delay, and the $J_2$ multipole. The framework's substitution rules govern SR Maxwell and flat-space dynamics; they do not extend to curved spacetime in any verified paper. Under the same minimal $c \to b$ Schwarzschild extension, each effect reduces to standard at GPS precision. However, the speculative extension is the load-bearing assumption, not a derived result. We therefore record these three as out-of-scope for the as-published framework, not as refuted by it. Every document in this category carries a substantive-AI TODO block on its §0 "Framework applicability" paragraph.

**Operational summary.** At GPS precision ($\mathbf{u}^2/c^2 \approx 1.7 \times 10^{-10}$; gravitational $GM/(rc^2) \approx 10^{-10}$; clock stability $\sim 10^{-14}$), the framework reproduces the standard prediction for every effect. The pre-launch frequency offset of (2) is unchanged; the operational engineering response is unchanged; no observable deviation was found. This is the expected outcome of mathematical equivalence in the appropriate limit; it is consistent with, but does not independently validate, either formulation. We observe that any disagreement with measurement at GPS precision would have been a serious problem for the framework; none was found.

<!-- TODO: human reviews and fills in — confirms the three-category classification of the nine effects is the right honest framing, particularly the "speculative extension" labelling of pt_02 / pt_07 / pt_08. The minimal-extension hypothesis (c → b in the Schwarzschild metric) is a substantive-AI proposal carried in this PR, not a derived consequence of any verified Gill paper; the steel-man pass on PR #73 made this explicit at every per-effect doc. -->

# §5 Four open questions for a curved-spacetime extension

This is the load-bearing section of the report. The four questions all surface from the proper-time companion's Marginal and Out-of-scope per-effect documents. Each is operationally invisible at current-generation GPS precision. Each is structurally load-bearing for the framework's extension to regimes where the differences would matter. Each question is one short paragraph by design. The per-effect documents under the GPS-relativity folder[^10] carry the supporting context.

**Q1. What is the correct curved-spacetime extension of the framework?** The minimal-extension hypothesis used in pt_01, pt_02, pt_04, pt_07, and pt_08 is $c \to b$ everywhere in the Schwarzschild line element, with $b \to c$ at low velocity. This is one possible extension. We considered three others without resolution:

- (a) the substitution $c \to b$ applied only in the time-time metric component $g_{00}$;
- (b) treating null and timelike geodesics differently, since the photon's local frame has $\mathbf{u} = 0$ and so $b = c$, but the worldline along which the metric is evaluated does not;
- (c) replacing the metric structure entirely with a $b$-dependent connection that reduces to the Schwarzschild connection in the $b = c$ limit.

For GPS at the $10^{-10}$ precision floor, all four hypotheses reduce to the same operational answer. The question matters for next-generation optical-clock GPS, where $10^{-17}$ stability is at the edge of distinguishing them, and for strong-field regimes where $GM/(rc^2)$ is no longer of order $10^{-10}$. The campaign cannot resolve this from the verified corpus. We flag it as the principal Tepper-input ask of this report.

**Q2. How does $\mathbf{u} \cdot \mathbf{a}$ transform under the proper-time boost?** The boost in question is the one from Maxwell-paper Eq. (11)[^6]. This is the same open question recorded in the proper-time cheat-sheet's open-questions list[^11], flagged for the radiation-reaction sub-investigation[^9]. It surfaces again in pt_05 of the GPS campaign. We observe that the third-term contribution of (7) depends on $\mathbf{u} \cdot \mathbf{a}$ having a frame-independent operational meaning. The covariance of this 3-vector dot product under the proper-time boost has not been derived in any verified paper. For GPS the question is moot ($\mathbf{u} \cdot \mathbf{a} = 0$ for circular orbits; of order $10^{-34}$ s for elliptical). However, it bears on whether the radiation-reaction predictions of the sub-investigation[^9] are frame-independent at the level the experiments need. Even a one-line confirmation that $\mathbf{u} \cdot \mathbf{a}$ transforms as a Lorentz scalar (or, conversely, that it does not) would resolve this for both campaigns.

**Q3. Does the framework predict a satellite-side Sagnac contribution?** The standard Sagnac correction $\Delta t = -2\boldsymbol{\omega} \cdot \mathbf{A} / c^2$ carries $c$ because it is derived in the rotating frame of the *receiver*. Under the proper-time framework, the receiver's effective light-speed is

$$b_{\rm recv} = \sqrt{c^2 + u_{\rm recv}^2}, \qquad (8)$$

where $u_{\rm recv} \leq 465$ m/s is the ground rotation, giving $b_{\rm recv}^2/c^2 - 1 \approx 1.2 \times 10^{-12}$. The receiver-side derivation in pt_06 uses $b_{\rm recv}$ and recovers the standard prediction at this precision. A symmetric reading of the framework would have the satellite's

$$b_{\rm sat} = \sqrt{c^2 + u_{\rm sat}^2} \qquad (9)$$

enter the Sagnac formula on the source side as well, with $u_{\rm sat} \approx 3874\,\text{m/s}$ giving $b_{\rm sat}^2/c^2 - 1 \approx 1.7 \times 10^{-10}$. That reading would predict an additional satellite-side correction at the $\sim 10^{-10}$ level of the standard Sagnac. This sits below current GPS precision, but at the edge of next-generation optical-clock receivers. We ask whether the framework intends the source-side and receiver-side $b$ to enter symmetrically, or whether the receiver-frame derivation is the complete prediction.

**Q4. Are there framework-specific $J_2$-like corrections beyond the speculative extension?** Earth's quadrupole moment $J_2 = 1.083 \times 10^{-3}$ modifies the gravitational potential through a $P_2(\cos\theta)$ multipole term. Under the speculative $c \to b$ extension, the $J_2$ correction enters via the same $1/b^2$ factor as the leading point-mass term. It reproduces the standard $+0.033\,\mu$s/day contribution to the combined offset. A first-principles derivation of multipole-induced clock corrections from the framework's underlying dynamics would clarify whether the standard contribution is the leading-order term, or whether additional $b$-dependent contributions arise at higher multipole order. The Lense–Thirring (frame-dragging) contribution from Earth's rotation sits in the same category. It is at the $\sim 10^{-16}$ level for current GPS, but within reach of ACES on the ISS. The framework's prediction here would need a curved-spacetime gravitomagnetic extension that we have not derived. For GPS the Q4 effect is sub-picosecond per day. For next-generation Earth-orbiting clock missions it is the leading post-Schwarzschild correction.

<!-- TODO: human reviews and fills in — confirms (a) Q1 captures the right load-bearing question (curved-spacetime extension); (b) Q2 is framed consistently with the way it appears in _proper_time_cheatsheet.md §4a and in issue #43; (c) Q3's "satellite-side Sagnac" reading is a sensible interpretation rather than a misreading of the framework; (d) Q4's J_2 + Lense-Thirring framing is the right scope (GPS-relevant only as future-mission context, not as a current-precision ask). -->

# §6 Suggested next steps

We order these by what would unblock the most downstream work.

1. **Resolve Q1.** Of the four questions, Q1 is the one that determines whether the others have well-defined answers; Q3 and Q4 are partially dependent on it, and Q2 is independent and propagates to the radiation-reaction sub-investigation[^9] regardless. Even a directional answer would suffice. Such an answer might take the form "the right extension is $c \to b$ everywhere", or "it is the $c \to b$ only in $g_{00}$ variant", or "the framework's intended GR extension is not yet written down, but the right answer is the cleanest minimal coupling consistent with the Maxwell paper's substitution rules". Any of these would let the GPS campaign convert the affected Marginal and Out-of-scope per-effect documents from speculative to derived.

2. **Address Q2 across two campaigns.** A $\mathbf{u} \cdot \mathbf{a}$ transformation rule would close out the third-term contribution of (7) in both pt_05 (GPS) and the radiation-reaction sub-investigation[^9]. The radiation-reaction sub-investigation is the regime where $\mathbf{u} \cdot \mathbf{a}$ is operationally observable; the GPS context is incidental. The cross-campaign consistency is the load-bearing reason to write the rule down.

3. **Defer Q3 and Q4 to next-generation precision.** Neither question affects current GPS operations. We flag them here for completeness and for future-mission scope, namely optical-clock GPS, ACES, and mission concepts at the $10^{-17}$–$10^{-18}$ clock-stability floor. A directional answer on either would be sufficient to scope follow-on threads.

A natural next campaign, should the disposition of Q1 warrant it, is a derived curved-spacetime extension of the framework — pursued either through the dual-Dirac equation or the proper-time photon propagator. We would scope such a campaign as a follow-on issue, keyed off the direction given. If, alternatively, the as-published framework's silence on GR is intentional — the framework is SR-and-flat-space by design, leaving GR to the standard GR apparatus — confirming that would let us close the four open questions in this report with that disposition recorded.

The prior interim report's §6 $r_e$ question[^1] remains the highest-load ask across all campaigns. This report's Q1 is the highest-load ask within the GPS scope; it is, however, a much smaller decision, since GPS is operationally non-falsifying at current precision.

<!-- TODO: human reviews and fills in — confirms the §6 ordering is the right disposition (Q1 first, Q2 cross-campaign, Q3/Q4 deferred); confirms the closing note about the prior interim report's §6 remaining the highest-load ask is the right framing relative to this report's narrower scope. -->

## References

[^1]: Interim report of 2026-05-25. `Roadmapping/Author_Reports/2026-05_interim_for_gill.md`.
[^2]: Issue #57 — GPS relativistic clock-correction question. `github.com/temoTxt/PyPhysics/issues/57`.
[^3]: Standard-method derivation campaign. PR #71, `github.com/temoTxt/PyPhysics/pull/71`.
[^4]: Ashby, N. (2003). *Relativity in the Global Positioning System.* Living Reviews in Relativity.
[^5]: Proper-time companion campaign. PR #73, `github.com/temoTxt/PyPhysics/pull/73`.
[^6]: Gill, T. L. and Zachary, W. W. *Two Mathematically Equivalent Versions of Maxwell's Equations.* See `Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md` for the per-equation verification of the substitution rules used here.
[^7]: GPS Interface Control Document ICD-200. Operational pre-launch frequency offset; IERS / WGS-84 / EGM2008 supply geoid and Earth-model parameters.
[^8]: Parameter sources: IERS conventions, WGS-84, EGM2008, GPS ICD-200.
[^9]: Radiation-reaction sub-investigation. Issue #43, `github.com/temoTxt/PyPhysics/issues/43`.
[^10]: Per-effect documents. `Roadmapping/GPS_Relativity/`.
[^11]: Proper-time cheat-sheet, open-questions list §4a. `Roadmapping/Electromagnetism/_proper_time_cheatsheet.md`.
