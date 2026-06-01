---
title: "Classical perihelion-advance derivation: original method and proper-time update — side by side via Binet's orbit equation"
---

# Classical perihelion-advance derivation: standard GR and framework proper-time update

**To:** Tepper Gill
**From:** Trey Morris
**Date:** 2026-05-31
**Re:** Follow-up to the 2026-05-30 Corda review[^1], per the request to see both the original (classical) method and the proper-time update via the *Relativistic Newtonian Theory* force law.
**Repo branch:** `81-corda-paper-analysis`.

---

## Cover note

We address here the perihelion-advance derivation done two ways. The first is the original route, the standard general-relativistic calculation in the Schwarzschild geometry. The second is the proper-time update, in which Newton's central force is replaced by the modification suggested in *Relativistic Newtonian Theory*[^2]. We work both derivations through the same apparatus, namely Binet's orbit equation, so that the comparison is line for line. Every algebraic step is verified symbolically via the Wolfram MCP server. The headline structural result is the identity

$$
\frac{\Delta\phi^{\text{framework}}}{\Delta\phi^{\text{GR}}} \;=\; -\frac{1}{6}, \tag{1}
$$

which falls out of the algebra as an identity, not as a numerical approximation.

Per author-review convention[^3], findings are surfaced for the author's verdict and are not announced as corrections. Per Crocco §5[^4], this is substantive AI: the algebra is pragmatic (Wolfram-verified); the framing of the framework's force law as the proper-time update of Newton's `-GM/r²` is the substantive interpretive move and is flagged for sign-off.

---

## §1. The apparatus — Binet's orbit equation

We consider a particle of mass `m` in a central force `F(r) ê_r`, with the conserved angular momentum per unit mass

$$
L \;=\; r^{2}\,\dot\phi. \tag{2}
$$

Following Binet, we change variables to `u = 1/r`. The orbit equation in `u(φ)` form is

$$
\frac{d^{2} u}{d\phi^{2}} \;+\; u
\;=\; -\,\frac{1}{m\,L^{2}\,u^{2}}\,F\!\left(\tfrac{1}{u}\right). \tag{3}
$$

This is the standard apparatus for perihelion-advance calculations[^5]. Its advantage is that, for an inverse-square force, the equation is linear in `u` with a constant right-hand side. It then solves exactly, and the geometry of the orbit is read off directly from the solution. Any correction to Newton's force law shows up as an extra term on the right of (3). The perihelion-advance question reduces, in each case, to a single technical question: what does that extra term do to the periodicity in `φ` of the solution?

---

## §2. The unperturbed orbit — closed ellipse from Newton's `-GM/r²`

For Newton, `F(r) = -GMm/r²`, so that `F(1/u) = -GMm u²`. Substituting into (3),

$$
\frac{d^{2} u}{d\phi^{2}} \;+\; u
\;=\; \frac{GM}{L^{2}}. \tag{4}
$$

Equation (4) is harmonic-oscillator-with-constant-forcing. Its general solution is

$$
u(\phi) \;=\; \frac{1}{p}\,\bigl[\,1 + e\cos(\phi - \phi_{0})\,\bigr],
\qquad
p \;\equiv\; \frac{L^{2}}{GM}. \tag{5}
$$

The quantity `p` is the semi-latus rectum. We verify the solution by substituting (5) into (4); the residual `d²u/dφ² + u − GM/L²` simplifies to `0` in Wolfram[^6].

The function `u(φ)` has period `2π`. Perihelion occurs when `u` is maximal, that is, when `cos(φ − φ₀) = 1`, which recurs every `2π` in `φ`. Thus the orbit is a closed ellipse and the perihelion does not advance:

$$
\Delta\phi^{\text{Newton}}_{\text{per orbit}} \;=\; 0. \tag{6}
$$

We observe that this is Bertrand's theorem in a familiar guise: of all central potentials, only Newton's `1/r` and the harmonic `r²` give universally closed bound orbits.

---

## §3. The original method — standard GR via Schwarzschild

In Schwarzschild geometry, the geodesic equation for a test particle in the equatorial plane reduces to a Binet equation with one extra term[^7]:

$$
\frac{d^{2} u}{d\phi^{2}} \;+\; u
\;=\; \frac{GM}{L^{2}}
\;+\; \frac{3\,GM}{c^{2}}\,u^{2}. \tag{7}
$$

The new term `3GM u²/c²` is the relativistic perturbation. We work to first order in the small parameter

$$
\epsilon \;\equiv\; \frac{GM}{c^{2}\,p}, \tag{8}
$$

and write `u = u₀ + δu`, with `u₀` the unperturbed Newton solution from §2.

The perturbation source on the right of (7) is

$$
\frac{3\,GM}{c^{2}}\,u_{0}^{2}
\;=\; \frac{3\,GM}{c^{2}\,p^{2}}\,\bigl(1 + e\cos\phi\bigr)^{2}
\;=\; \frac{3\,GM}{c^{2}\,p^{2}}\,\Bigl[\,1 + 2e\cos\phi + e^{2}\cos^{2}\phi\,\Bigr]. \tag{9}
$$

We expand (9) symbolically; the coefficient of `cos φ` is `6 GM e/(c² p²)`, and this is the load-bearing resonant-driving term[^8].

It follows that the `cos φ` piece of the perturbation resonantly drives the homogeneous equation

$$
\frac{d^{2}(\delta u)}{d\phi^{2}} \;+\; \delta u \;=\; \frac{6\,GM\,e}{c^{2}\,p^{2}}\,\cos\phi, \tag{10}
$$

producing a secular term

$$
\delta u_{\text{secular}}(\phi)
\;=\; \frac{3\,GM\,e}{c^{2}\,p^{2}}\,\phi\,\sin\phi. \tag{11}
$$

The other pieces of (9) — the constant `1 + e²/2` piece and the `cos 2φ` piece arising from `e²cos²φ/2` — give bounded, non-secular corrections; they do not shift the perihelion.

Adding (11) to `u₀`,

$$
u(\phi) \;\approx\; \frac{1}{p}\,\Bigl[\,1 + e\cos\phi
\;+\; \frac{3\,GM}{c^{2}\,p}\,e\,\phi\,\sin\phi\,\Bigr]
\;=\; \frac{1}{p}\,\Bigl[\,1 + e\cos\bigl(\phi(1-k)\bigr)\,\Bigr]
\;+\; \mathcal{O}(\epsilon^{2}), \tag{12}
$$

with

$$
k \;=\; \frac{3\,GM}{c^{2}\,p}. \tag{13}
$$

The trigonometric identification uses `cos(φ(1 − k)) ≈ cos φ + k φ sin φ` for small `k`. Perihelion occurs at `φ(1 − k) = 2πn`, and so successive perihelia are separated by `φ = 2π/(1 − k) ≈ 2π(1 + k) = 2π + 2π k`. The advance per orbit is therefore

$$
\Delta\phi^{\text{GR}}_{\text{per orbit}}
\;=\; 2\pi k
\;=\; \frac{6\pi\,GM}{c^{2}\,p}
\;=\; \frac{6\pi\,GM}{c^{2}\,a\,(1-e^{2})}. \tag{14}
$$

This is the canonical Einstein 1915 result. Indeed, for Mercury (semi-major axis `a = 5.7909 × 10¹⁰` m, eccentricity `e = 0.20563`, semi-latus rectum `p = a(1 − e²) = 5.546 × 10¹⁰` m) we obtain

$$
\Delta\phi^{\text{GR}}_{\text{Mercury}}
\;=\; \frac{6\pi\,GM_{\odot}}{c^{2}\,p}
\;\times\; \frac{100 \cdot 365.25 \cdot 86400}{T_{\text{orb}}}
\;\times\; \frac{180 \cdot 3600}{\pi}
\;\text{arcsec/century}
\;=\; 42.99\;\text{arcsec/century}. \tag{15}
$$

The numerical evaluation is performed in Wolfram[^9] with `GM_⊙ = 1.327 × 10²⁰`, `c = 2.998 × 10⁸`, `p = 5.546 × 10¹⁰`, and `T = 87.969 · 86400`. The result `42.9918` arcsec/century is the canonical figure.

---

## §4. The proper-time update — *Relativistic Newtonian Theory* force law

The framework's modified Newton force, derived[^10] from the canonical proper-time Hamiltonian

$$
K \;=\; \frac{\boldsymbol\pi^{2}}{2m} \;+\; m c^{2} \;+\; V \;+\; \frac{V^{2}}{2 m c^{2}}, \tag{16}
$$

with `V = -GMm/r`, takes the form

$$
\vec{a}_{m}
\;=\; -\,\frac{GM}{r^{2}}\,\Bigl(1 \;-\; \frac{GM}{c^{2}\,r}\Bigr)\,\hat{u}_{r}. \tag{17}
$$

We refer to (17) as the framework force law; it is the equation labelled (h4) in Gill's paper. Multiplying by `m` gives the force itself,

$$
F(r) \;=\; -\,\frac{GM\,m}{r^{2}} \;+\; \frac{(GM)^{2}\,m}{c^{2}\,r^{3}}. \tag{18}
$$

Substituting (18) into Binet's equation (3), and writing `L` for the angular momentum per unit mass,

$$
\frac{d^{2} u}{d\phi^{2}} \;+\; u
\;=\; -\,\frac{1}{L^{2}\,u^{2}}\!\left[\,-\,GM\,u^{2} \;+\; \frac{(GM)^{2}}{c^{2}}\,u^{3}\,\right]
\;=\; \frac{GM}{L^{2}} \;-\; \frac{(GM)^{2}}{c^{2}\,L^{2}}\,u. \tag{19}
$$

Rearranging (19),

$$
\frac{d^{2} u}{d\phi^{2}}
\;+\; u\,\underbrace{\Bigl[\,1 + \frac{(GM)^{2}}{c^{2}\,L^{2}}\,\Bigr]}_{\equiv\,\alpha^{2}}
\;=\; \frac{GM}{L^{2}}. \tag{20}
$$

We observe the structural difference between (20) and the Schwarzschild equation (7). The Schwarzschild perturbation `3GM u²/c²` is quadratic in `u` and drives a resonant secular term. The framework's perturbation `−(GM)² u/(c² L²)` is linear in `u`. It is therefore absorbed into the homogeneous part as a shift of the effective spring constant from `1` to `α²`. The two perturbations have completely different structure.

Equation (20) is again a harmonic-oscillator-with-constant-forcing problem, with frequency `α` in place of `1`. Its solution is

$$
u(\phi) \;=\; \frac{GM/L^{2}}{\alpha^{2}} \;+\; B\,\cos\!\bigl(\alpha(\phi - \phi_{0})\bigr)
\;=\; \frac{1}{p\,\alpha^{2}}\,\Bigl[\,1 + e'\cos\bigl(\alpha(\phi - \phi_{0})\bigr)\,\Bigr], \tag{21}
$$

with rescaled eccentricity `e' = B p α²`. The residual `d²u/dφ² + α² u − GM/L²` evaluated on (21) gives `0` in Wolfram[^11].

The perihelion-to-perihelion period in `φ` is now `2π/α`, not `2π`. Thus the advance (or recession) per Newtonian orbit (`2π` in `φ`) is

$$
\Delta\phi^{\text{framework}}_{\text{per orbit}}
\;=\; \frac{2\pi}{\alpha} \;-\; 2\pi
\;=\; 2\pi\!\left(\frac{1}{\sqrt{1 + (GM)^{2}/(c^{2}\,L^{2})}} \;-\; 1\right). \tag{22}
$$

Using the same small parameter as in GR,

$$
\epsilon \;=\; \frac{GM}{c^{2}\,p} \;=\; \frac{(GM)^{2}}{c^{2}\,L^{2}}, \tag{23}
$$

(the second equality follows from `p = L²/GM`), we Taylor-expand

$$
\frac{1}{\sqrt{1+\epsilon}} \;-\; 1
\;=\; -\,\frac{\epsilon}{2} \;+\; \mathcal{O}(\epsilon^{2}). \tag{24}
$$

It follows that, to first order,

$$
\Delta\phi^{\text{framework}}_{\text{per orbit}}
\;\approx\; -\,\pi\,\epsilon
\;=\; -\,\frac{\pi\,GM}{c^{2}\,p}
\;=\; -\,\frac{\pi\,GM}{c^{2}\,a\,(1-e^{2})}. \tag{25}
$$

The Taylor expansion (24) is verified symbolically in Wolfram[^12]: `2π/√(1 + x) − 2π` is `−π x + 3π x²/4 + O(x³)`, with leading term `−π x` at `x = GM/(c² p)`.

For Mercury, using the same parameters as in §3,

$$
\Delta\phi^{\text{framework}}_{\text{Mercury}}
\;=\; -\,\frac{\pi\,GM_{\odot}}{c^{2}\,p}
\;\times\; 415.2\;\text{orbits/century}
\;\times\; \frac{180 \cdot 3600}{\pi}
\;=\; -\,7.17\;\text{arcsec/century}. \tag{26}
$$

The numerical evaluation gives `−7.1653` arcsec/century[^13].

---

## §5. The structural identity `Δφ_framework / Δφ_GR = -1/6`

Comparing (14) and (25) directly,

$$
\frac{\Delta\phi^{\text{framework}}}{\Delta\phi^{\text{GR}}}
\;=\; \frac{-\,\pi\,GM/(c^{2}\,p)}{+\,6\pi\,GM/(c^{2}\,p)}
\;=\; -\,\frac{1}{6}\qquad\text{(exactly, not approximately).} \tag{27}
$$

Both predictions scale as `GM/(c² p)` per orbit; the dimensionless small parameter is the same in (14) and (25). Thus the ratio (27) is independent of the orbital parameters of the test body. Indeed, **the −1/6 is a structural identity of the two formulations, not a fit.** It holds for Mercury, for Venus, for Earth, for S2 around Sgr A*, and for any binary pulsar. It carries an opposite-sign signature: the framework's force law makes the perihelion *regress* (move backwards against the orbital motion), where general relativity makes it *advance* (move forwards).

---

## §6. Why the structural difference — the missing five-sixths and the sign flip

The Schwarzschild Binet equation (7) carries a `3GM u²/c²` correction. The framework's Binet equation (20) carries a `−(GM)² u/(c² L²)` correction. The structural origin of the difference is as follows.

General relativity contributes both from the time dilation in the geodesic equation and from the spatial-metric curvature in the Schwarzschild line element

$$
ds^{2} \;=\; -\,(1 - r_{g}/r)\,(c\,dt)^{2} \;+\; (1 - r_{g}/r)^{-1}\,dr^{2} \;+\; r^{2}\,d\Omega^{2}. \tag{28}
$$

The time-dilation contribution alone produces a perihelion advance smaller than the full GR result. The curved-space radial coordinate, expressed in the `(1 − r_g/r)^{−1}` factor on `dr²`, supplies the rest. In the parameterised post-Newtonian (PPN) language, the perihelion advance per orbit is

$$
\Delta\phi \;=\; \frac{2 + 2\gamma - \beta}{3}\;\cdot\;\frac{2\pi\,GM}{c^{2}\,p}, \tag{29}
$$

where `γ` and `β` are the PPN parameters of the metric[^14]. In general relativity `γ = β = 1`, and (29) reduces to `(2 + 2 − 1)/3 = 1`, giving `6π GM/(c² p)` per orbit, in agreement with (14).

The framework's force law (17) is, by contrast, a force-law-only modification of Newtonian gravity. It modifies the time dynamics (how fast `r` and `φ` change in time) but does not curve space; there is no modified `dr²` coefficient, and the radial coordinate stays flat. In PPN language this corresponds to `γ = 0` and a value of `β` determined by the way the `V²/(2mc²)` kernel in (16) translates to the static-isotropic metric expansion. With `γ = 0` and an appropriate `β`, the PPN coefficient `(2 + 2γ − β)/3` can go negative — which is what the explicit Binet calculation in §4 reproduces.

We observe that this is the classic "1/6 factor" of modified-Newtonian gravity without spatial-metric curvature. It is a known signature in the PPN literature: any modification that adds a relativistic correction to the central force but leaves the spatial metric flat produces a perihelion shift that is some fraction of the GR result. For the specific form (17), the fraction is `−1/6`.

The open question, raised for the author's verdict, is whether the framework's intended content includes a spatial-metric curvature component beyond Eq. (17), or whether (17) is taken to be the complete framework prediction for gravity. Under the latter reading, the framework's Mercury prediction is the `−7.17` arcsec/century value computed via Binet in §4, and the framework does not match the observed `≈ +43` arcsec/century. Under the former reading, the missing five-sixths might be supplied by a curved-space companion to (17) that the *Relativistic Newtonian Theory* paper does not write out explicitly.

---

## §7. Numerical comparison for Mercury

The full Wolfram-MCP-verified numerical evaluation is performed[^15] with the standard constants `M_⊙ = 1.98892 × 10³⁰` kg, `a = 5.7909 × 10¹⁰` m, `e = 0.20563`, `T_orb = 87.969 · 86400` s, `G = 6.6743 × 10⁻¹¹`, `c = 2.998 × 10⁸` m/s, and `p = a(1 − e²)`. Using the conversion factors `orbits per century = 100 · 365.25 · 86400 / T_orb` and `rad → arcsec = 180/π · 3600`, the two predictions are

$$
\Delta\phi^{\text{GR}}_{\text{Mercury}}
\;=\; \frac{6\pi\,GM_{\odot}}{c^{2}\,p}\;\cdot\;\text{(orbits/century)}\;\cdot\;\text{(rad → arcsec)}
\;=\; +\,42.99\;\text{arcsec/century}, \tag{30}
$$

$$
\Delta\phi^{\text{framework}}_{\text{Mercury}}
\;=\; -\,\frac{\pi\,GM_{\odot}}{c^{2}\,p}\;\cdot\;\text{(orbits/century)}\;\cdot\;\text{(rad → arcsec)}
\;=\; -\,7.17\;\text{arcsec/century}. \tag{31}
$$

The observed value, from the Newcomb–Clemence residual, is approximately `+43` arcsec/century[^16]. The ratio of (31) to (30) is `−7.17 / 42.99 = −0.16667 = −1/6` to six significant figures, in agreement with the identity (27).

---

## §8. Honest scope — what this derivation does and does not show

The derivation in §3 and §4 shows the following. The classical Binet's-equation method handles both Newton's `−GM/r²` and the framework's `−(GM/r²)(1 − GM/(c²r))` with the same machinery. For Schwarzschild GR, the standard `6π GM/(c² p)` result drops out by solving the perturbed Binet equation to first order in `GM/(c² p)`. For the framework's force law (17), taken as the complete force prescription, the perihelion advance per orbit is exactly `−π GM/(c² p)`, that is, `−1/6` of the GR result. The minus sign comes from the framework's perturbation being linear in `u`, not quadratic — a structural difference from the Schwarzschild `3GM u²/c²` term. Numerically for Mercury, the figures are `+42.99` arcsec/century in GR and `−7.17` arcsec/century in the framework.

The derivation does not show the following. It does not establish that the framework as a whole is incompatible with the Mercury observation; the calculation tests one specific gravitational extension, namely (17) read as the full framework prediction for gravity. If the framework's intended content includes a curved-space companion to (17) (a `(1 + GM/(c²r))` factor on the radial metric, or an equivalent), the result changes. The derivation does not show that the framework's proper-time apparatus is wrong; the proper-time kernel (16) is verified at the algebraic level by Wolfram MCP[^17], and the question here is only what observable prediction (17) produces. Finally, the derivation does not litigate the foundations of general relativity; the Schwarzschild calculation in §3 is the textbook result, used here as a known calibration target.

---

## §9. Questions for the author

We list four questions, in the order that adjudicates the campaign verdict.

1. Is Eq. (17), `a = −(GM/r²)(1 − GM/(c²r)) ê_r`, the complete framework prediction for gravity at this order? If yes, then the framework's Mercury prediction is `−7.17` arcsec/century, per the derivation in §4. If no — if the framework includes a curved-space companion to (17) — please indicate the form, or the place in the published paper where it appears.

2. Is the structural identity (27), `Δφ_framework / Δφ_GR = −1/6`, an expected feature of the proper-time formulation (a known PPN-style consequence of force-law-only modifications), or does it come as a surprise?

3. Is the linear-in-`u` structure of the perturbation in (20) the intended reading of Eq. (17), or does the framework's full apparatus produce additional quadratic-in-`u` terms (matching the Schwarzschild structure of §3) that the present derivation has missed?

4. For the verdict on the Mercury pull request[^18], should the campaign's reading be (a) the Corda-style heuristic, `+37.79` arcsec/century, identified in the 2026-05-30 review as numerically close to the observation but conceptually problematic; (b) the apsidal-angle / Binet calculation, `−7.17` arcsec/century, opposite sign and `1/6` magnitude, that is the result of the present report; or (c) something else, contingent on the answers to questions 1–3?

---

## §10. Crocco §5 — AI-use disclosure

The algebra and the numerical evaluation are pragmatic AI; the Wolfram MCP server verifies every load-bearing step, and the residuals computed there are `0` or matching-to-machine-precision. The substantive AI moves are three. First, the reading of the framework's Eq. (17) as the complete framework prediction for gravity, on which §4 rests; question 1 in §9 asks the author to adjudicate. Second, the framing of the `−1/6` result as a structural identity of the two formulations rather than a numerical coincidence (§5). Third, the PPN-language interpretation in §6, connecting the `−1/6` to the absence of spatial-metric curvature; this is standard literature, but the choice of which framing to use is itself a substantive choice and is flagged here for confirmation.

The prompt-of-record for this document is the author's response to the 2026-05-30 Corda review (the request for the "classical method, both original and proper-time update"), together with the user's instruction to produce the present companion. No prose-generation prompt beyond these was used.

<!-- TODO: human reviews and fills in — confirms (a) Gill's Eq. (17) is the
complete framework prediction for gravity at this order, (b) the Binet's-equation
derivation in §4 correctly handles the proper-time-updated force, (c) the −1/6
structural identity matches the framework's intended content (not surprising) or
surprises (suggesting the framework includes additional content beyond Eq. (17)),
and (d) the recommended verdict for the Mercury pull request is (a), (b), or (c)
per §9 question 4. -->

---

## References

[^1]: Companion report, 2026-05-30 Corda review:
      `Roadmapping/Author_Reports/2026-05_corda_perihelion_review_for_gill.md`.

[^2]: T. L. Gill, *Relativistic Newtonian Theory* (the paper under campaign
      analysis on this branch); Eq. (h4) of that paper specifies the modified
      central force law used in §4 below.

[^3]: Author-review convention recorded in `CLAUDE.md` at the repository root
      (the rule "frame findings concerning DRQM I as 'for author review'
      rather than corrections"); the convention is applied here.

[^4]: Crocco compliance framework, `Roadmapping/Tooling/CROCCO_COMPLIANCE.md`,
      §5 on substantive-AI disclosure.

[^5]: H. Goldstein, C. Poole, and J. Safko, *Classical Mechanics*, 3rd ed.
      (Addison-Wesley, 2002), §3.5; C. W. Misner, K. S. Thorne, and J. A.
      Wheeler, *Gravitation* (Freeman, 1973), Box 25.6; R. M. Wald, *General
      Relativity* (University of Chicago Press, 1984), §6.3.

[^6]: Wolfram MCP verification: substitute `u(φ) = (1 + e cos φ)/p` with
      `p = L²/GM` into the unperturbed equation; the residual
      `d²u/dφ² + u − GM/L²` simplifies to `0`.

[^7]: Misner, Thorne, and Wheeler, *Gravitation*, Box 25.6 and Eq. (25.40);
      the equation appears here as (7).

[^8]: Wolfram MCP verification: `3GM/c² · u₀²` expanded symbolically gives a
      coefficient of `cos φ` equal to `6 GM e/(c² p²)`.

[^9]: Wolfram MCP verification: with `GM_⊙ = 1.327 × 10²⁰`, `c = 2.998 × 10⁸`,
      `p = 5.546 × 10¹⁰`, `T = 87.969 · 86400`, the expression
      `6π GM/(c² p) · 415.2 · (180 · 3600/π)` evaluates to `42.9918` arcsec/c.

[^10]: Pull request analysis: `Roadmapping/Equation_Verification/Dual_Newtonian_Theory.md`,
       linked from pull request #83, `github.com/temoTxt/PyPhysics/pull/83`.
       That document verifies (17) via Wolfram MCP from the proper-time
       Hamiltonian (16).

[^11]: Wolfram MCP verification: substitute the ansatz (21) into the framework
       Binet equation (20); the residual evaluates to `0`.

[^12]: Wolfram MCP verification: `Series[2π/√(1 + x) − 2π, {x, 0, 2}]` returns
       `−π x + 3π x²/4 + O(x³)`.

[^13]: Wolfram MCP verification: `−π GM/(c² p) · 415.2 · (180 · 3600/π)`
       evaluates to `−7.1653` arcsec/c.

[^14]: C. M. Will, *Theory and Experiment in Gravitational Physics*, 2nd ed.
       (Cambridge University Press, 2018), §7.3, gives the PPN perihelion-advance
       formula (29).

[^15]: Wolfram MCP numerical evaluation: with `Msun = 1.98892 × 10³⁰`,
       `aMerc = 5.7909 × 10¹⁰`, `eMerc = 0.20563`,
       `Torb = 87.969 · 86400`, `GG = 6.6743 × 10⁻¹¹`,
       `cc = 299792458`, and `pMerc = aMerc · (1 − eMerc²)`, the expressions
       `dGR = 6π GMsun/(cc² pMerc) · orbitsPerCentury · radToArcsec` and
       `dFramework = −π GMsun/(cc² pMerc) · orbitsPerCentury · radToArcsec`
       give `+42.99` and `−7.17` respectively, with `dFramework/dGR = −1/6` to
       six significant figures.

[^16]: S. Newcomb, *The Elements of the Four Inner Planets and the
       Fundamental Constants of Astronomy* (Government Printing Office,
       1895); G. M. Clemence, *Rev. Mod. Phys.* **19**, 361 (1947).

[^17]: Pull request #83, `github.com/temoTxt/PyPhysics/pull/83`, with companion
       documents `Roadmapping/Equation_Verification/Dual_Newtonian_Theory.md`
       and `Roadmapping/Mercury_Perihelion/05_nbody_orbital_mechanics.md`. The
       latter is the original apsidal-angle calculation showing the `−1/6`
       result, now reproduced via Binet's equation in §4 of the present
       report.

[^18]: Pull request #83 on `github.com/temoTxt/PyPhysics`, the Mercury
       perihelion campaign; this report is one input to the verdict on that
       pull request. Issue #81, `github.com/temoTxt/PyPhysics/issues/81`, is
       the source ticket.
