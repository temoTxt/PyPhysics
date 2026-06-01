---
title: "Corda (2021) on the perihelion of planets — a math and physical-intuition check for the *Relativistic Newtonian Theory* foundations"
---

# Author-review report: Corda (2021), *The secret of planets' perihelion between Newton and Einstein* — a math and physical-intuition check

**To:** Tepper Gill (DRQM I and *Relativistic Newtonian Theory* senior author).
**From:** Trey Morris.
**Date:** 2026-05-30.
**Re:** Issue #81[^1], the basis for the Mercury perihelion derivation in *Relativistic Newtonian Theory*.

**Source paper under review:** Christian Corda, *Phys. Dark Univ.* **32** (2021) 100834[^2]. The pdf in the repository is at `Roadmapping/External_Papers/Secret_Corda.pdf`[^3]; the marker-pdf conversion is at `Roadmapping/Converted_Markdown/Secret_Corda/Secret_Corda.md`[^4].

---

## Cover note (for author review, not for journal submission)

The author indicated that the perihelion derivation in *Relativistic Newtonian Theory* takes Corda's framework as its starting point. The present report checks Corda's algebra and physical intuition end-to-end. Its purpose is to let the author adjudicate whether the framework's intended Mercury prediction is the Corda-style `πm/M + relativistic correction` template (the chain analysed in pull request #83[^5], giving `+37.79` arcsec/century) or the proper apsidal-angle calculation of the modified force law (document 05 in that pull request, giving `−7.17` arcsec/century).

Per author-review convention[^6], the findings below are surfaced for the author's verdict; they are not announced as corrections. Every algebraic claim is verified symbolically via the Wolfram MCP server. Every interpretive claim is flagged as substantive AI.

Per Crocco §5[^7] compliance, this is substantive AI. The prompt-of-record is the user-instruction sequence on issue #81 and the user's 2026-05-30 instruction to "write the complete analysis in a markdown file for further review by Tepper." No prose-generation prompt beyond the user's instructions was used.

---

## §1. The structure of Corda's paper (eight pages, seventy-seven numbered equations)

We list the principal load-bearing equations of Corda's paper and the result each section produces.

| § | Topic | Headline equation | Result |
|---|---|---|---|
| §1 | Introduction, GR perihelion | Eq. (1): `Δφ_GR = 24π³ a²/[T₀² c² (1−e²)] = 3π r_g/[a(1−e²)]` | `42.99` arcsec/century for Mercury |
| §2 | "Circular orbit" Newtonian back-reaction | Eq. (20): `Δφ = π m/M` | `44.66` arcsec/century (Approach 1) |
| §3 | Harmonic-oscillator / apsidal-angle | Eq. (31): `φ₀ = 2π [3 + F'(r)r/F(r)]^(−1/2)` | the same `44.66` arcsec/century (Approach 2) |
| §4 | Reduced-mass Kepler | Eq. (57): `T = T₀ (1 + m/M)^(−1/2)` | the same `44.66` arcsec/century (Approach 3) |
| §5 | Breakdown of `πm/M` for Venus and Earth | (numerical only) | `30×` and `51×` too large |
| §6 | Gravitational time dilation | Eq. (67): `Δφ_F = (5π/2) r_g/r₀` | `34.31` arcsec/century (TD only) |
| §7 | Rotational (Langevin) time dilation | Eq. (72): `τ ≈ (1 − ½ r²ω²/c²) t` | (intermediate) |
| §8 | Combined gravitational and rotational TD | Eq. (77): `Δφ_F = 3π r_g/r₀` | `41.17` arcsec/century (matches `e = 0` GR) |
| §9 | Conclusion: "Newtonian plus corrections matches GR" | — | — |

The paper's two thesis claims appear verbatim in §9. The first is that it is not correct that Newtonian theory cannot predict the anomalous rate of precession of the perihelion of planets' orbits; the real problem is instead that a pure Newtonian prediction is too large. The second is that the perihelion's precession can be achieved with the same precision as general relativity by extending Newtonian gravity through the inclusion of gravitational and rotational time-dilation effects.

The math and physical-intuition check below addresses each thesis separately.

---

## §2. Every load-bearing equation surfaced and verified

We state Corda's equations verbatim (translated from the marker-pdf conversion), record the Wolfram-MCP-verified algebra, and flag the physical-intuition status. To preserve readability across the seventeen equations of Corda's §§2–8, the numbering below is Corda's original; our own derived expressions carry our document's numbering and are placed in display form.

### §2.1 Eqs. (3)–(5) — circular-orbit Newton baseline (no back-reaction)

$$
\frac{G M m}{r_{0}^{2}} \;=\; \frac{m\,v_{0}^{2}}{r_{0}}, \tag{3}
$$

$$
v_{0} \;=\; \sqrt{G M / r_{0}}, \tag{4}
$$

$$
T_{0} \;=\; \frac{2\pi\,r_{0}}{v_{0}} \;=\; \frac{2\pi\,r_{0}^{3/2}}{\sqrt{G M}}. \tag{5}
$$

These are standard textbook expressions (Newton with centripetal balance). The algebra is correct. The physical content is that Mercury orbits a stationary Sun. **Pass.**

### §2.2 Eqs. (6)–(7) — Newtonian angular velocity and unperturbed angle

$$
\omega_{0} \;=\; 2\pi / T_{0}, \tag{6}
$$

$$
\varphi_{0} \;=\; \omega_{0}\,T_{0} \;=\; 2\pi. \tag{7}
$$

These are definitionally correct. **Pass.**

### §2.3 Eqs. (8)–(13) — the substitution that drives Corda's §2 result

Corda's core substitution is to replace `M` in (3) with `M + m`, on the ground that the Sun-centred observer sees Mercury's reaction force on the Sun:

$$
\frac{G(M+m)\,m}{r_{0}^{2}} \;=\; \frac{m\,v^{2}}{r_{0}}. \tag{8}
$$

He justifies the substitution by writing the Newton's-third-law pair:

$$
\vec F_{G} \;=\; \frac{G M m}{r_{0}^{2}}\,\hat u_{r}, \tag{9}
$$

$$
M\,a_{s}\,\hat u_{r} \;=\; \frac{G M m}{r_{0}^{2}}\,\hat u_{r}
\;\;\Rightarrow\;\;
a_{s} \;=\; \frac{G m}{r_{0}^{2}}, \tag{10}
$$

$$
m\,a_{m}\,\hat u_{r} \;=\; -\,\frac{G M m}{r_{0}^{2}}\,\hat u_{r}
\;\;\Rightarrow\;\;
a_{m} \;=\; -\,\frac{G M}{r_{0}^{2}}, \tag{11}
$$

$$
a\,\hat u_{r} \;\equiv\; a_{m} - a_{s} \;=\; -\,\frac{G(M+m)}{r_{0}^{2}}\,\hat u_{r}, \tag{12}
$$

$$
F\,\hat u_{r} \;=\; -\,\frac{G(M+m)\,m}{r_{0}^{2}}\,\hat u_{r}. \tag{13}
$$

The algebra of (8)–(13) is correct. It is the standard two-body relative-acceleration result. The same `G(M + m)` appears in textbook reduced-mass treatments (Goldstein §3.1, MTW §25). What `(M + m)` captures is the kinematic fact that the *relative* acceleration of Mercury with respect to the Sun is bigger than `GM/r²` by the factor `(M + m)/M = 1 + m/M`. **Pass on algebra.**

**Physical-intuition status, for the author's verdict.** What is being computed is the *relative* acceleration in the centre-of-mass frame; it is not a Sun-fixed force on Mercury. The standard reduced-mass two-body decomposition replaces the two-body problem with a one-body problem for a reduced-mass particle `μ = Mm/(M + m)` orbiting a fixed centre at `M + m`. Corda's Eq. (12) is correct. The interpretation in §2, in which the substitution is read as a "perturbation" to Mercury's force, muddles the standard story.

### §2.4 Eqs. (14)–(17) — perturbed velocity and period

$$
v \;=\; \sqrt{G(M+m)/r_{0}}, \tag{14}
$$

$$
T \;=\; \frac{2\pi\,r_{0}}{v} \;=\; \frac{2\pi\,r_{0}^{3/2}}{\sqrt{G(M+m)}}, \tag{15}
$$

$$
(M+m)^{-1/2} \;=\; M^{-1/2}\,(1 + m/M)^{-1/2}, \tag{16}
$$

$$
T \;=\; T_{0}\,(1 + m/M)^{-1/2}. \tag{17}
$$

The algebra of (17) is verified by Wolfram MCP[^8]: the residual `T − T₀/√(1 + m/M)` simplifies to `0`. The period decreases slightly because the effective gravitational constant in the one-body problem is `G(M + m)`, not `GM`. This is a real kinematic effect on the sidereal period. **Pass on algebra.**

### §2.5 Eqs. (18)–(20) — angular velocity and "perihelion advance"

$$
\omega \;=\; \omega_{0}\,(1 + m/M)^{1/2}, \tag{18}
$$

$$
\varphi \;=\; \omega\,T_{0} \;=\; 2\pi\,(1 + m/M)^{1/2}
\;\simeq\; 2\pi\,(1 + m/(2M)), \tag{19}
$$

$$
\Delta\varphi \;=\; \varphi - \varphi_{0} \;\simeq\; \pi\,m/M. \tag{20}
$$

The Taylor expansion is verified by Wolfram MCP[^9]: `Series[2π √(1 + x) − 2π, {x, 0, 2}]` returns `π x − π x²/4 + O(x³)`, with leading term `π m/M`. Numerically for Mercury (`m/M = 1.66 × 10⁻⁷`), one obtains

$$
\Delta\varphi \;=\; 5.21 \times 10^{-7}\,\text{rad/orbit}
\;\times\; 415.2\,\text{orbits/century}
\;\times\; \tfrac{180 \cdot 3600}{\pi}
\;=\; 44.66\,\text{arcsec/century}. \tag{R1}
$$

The arithmetic is Wolfram-verified. **Pass on algebra.**

**Physical-intuition status, for the author's verdict — the load-bearing critique.** What Eq. (20) computes is the angle Mercury sweeps in the *unperturbed* period `T_0` at the *perturbed* angular velocity `ω`. This quantity is not a perihelion advance. The perihelion advance is the angular rotation of the orbit's major axis between successive perihelia. It is a property of the *shape* of the orbit, not of the period or angular velocity.

For an inverse-square force, the orbit is a closed ellipse (Bertrand's theorem). The angle between successive perihelia is therefore exactly `2π`, irrespective of `T` and `ω`. Changing `G` from `GM` to `G(M + m)` rescales the orbital speed and the period, but it leaves the orbit a closed ellipse. The `πm/M` quantity computed in (20) is a real kinematic effect; it is not the perihelion-advance quantity. The §3 derivation makes this explicit; we return to it in §2.7.

### §2.6 Eqs. (21)–(31) — Corda's own apsidal-angle formalism (the paper's §3)

Following Price and Rush (1979), Corda writes the radial equation, imposes angular-momentum conservation, Taylor-expands around `r_0`, and arrives at the classical Bertrand-style apsidal angle:

$$
F_{c0}(r) \;=\; m\,\bigl(\ddot r - \omega_{0}^{2}\,r\bigr), \tag{21}
$$

$$
J_{0} \;=\; m\,r^{2}\,\omega_{0}, \tag{22}
$$

$$
F_{c0}(r) \;=\; m\,\Bigl(\ddot r \;-\; \frac{J_{0}^{2}}{m^{2}\,r^{3}}\Bigr), \tag{23}
$$

$$
F_{c0}(r_{0}) \;=\; -\,\frac{J_{0}^{2}}{m\,r_{0}^{3}}. \tag{24}
$$

After perturbing `x = r − r_0` and series-expanding,

$$
\ddot x \;+\; m^{-1}\!\left[\,-\,\frac{3\,F_{c0}(r_{0})}{r_{0}} \;-\; F'_{c0}(r_{0})\,\right] x \;=\; 0, \tag{27}
$$

$$
T_{0} \;=\; 2\pi\!\left(\frac{m}{-\,3\,F_{c0}(r_{0})/r_{0} \;-\; F'_{c0}(r_{0})}\right)^{1/2}, \tag{28}
$$

$$
\frac{\varphi_{0}}{2}
\;=\; \pi\!\left(\frac{m}{-\,3\,F_{c0}/r_{0} \;-\; F'_{c0}}\right)^{1/2}\,\frac{J_{0}}{m\,r_{0}^{2}}. \tag{29}
$$

Using (24) in (29):

$$
\varphi_{0}
\;=\; 2\pi\!\left[\,3 \;+\; \frac{F'_{c0}(r_{0})\,r_{0}}{F_{c0}(r_{0})}\,\right]^{-1/2}. \tag{31}
$$

We verify (31) for `F = −k/r²` symbolically in Wolfram[^10]: the ratio `F'(r) r / F(r)` evaluates to `−2`, and the apsidal angle becomes `φ_0 = 2π/√(3 + (−2)) = 2π`. Indeed, for an inverse-square force the apsidal angle is `2π`, the orbit is closed, and there is no precession. This is precisely Bertrand's theorem at work.

The algebra of (21)–(31) is correct, and this conclusion contradicts the §2 reading. Corda derives the apsidal-angle formula from scratch in his Eq. (31), and specialises immediately afterwards to the inverse-square case `F_{c0} = F_G`. He records there that the apsidal angle reduces to `φ_0 = 2π`, the same value as his own Eq. (7), and notes explicitly that this is the closed-orbit result. **Pass on algebra; the conflict is with the §2 interpretation.**

### §2.7 Eqs. (32)–(41) — Corda's "promotion" to the `(1 + m/M)` force

To extend the apsidal-angle calculation to the reduced-mass case, Corda substitutes the central force law:

$$
F_{c0}(r) \;\to\; F_{c}(r) \;=\; \Bigl(1 + \frac{m}{M}\Bigr)\,F_{c0}(r), \tag{32}
$$

$$
\omega_{0} \;\to\; \omega, \tag{33}
$$

$$
J_{0} \;\to\; J \;=\; m\,r^{2}\,\omega. \tag{34}
$$

This is meant to mirror the §2 substitution `M → M + m` inside the central-force law. Then

$$
F'_{c}(r_{0}) \;=\; \Bigl(1 + \frac{m}{M}\Bigr)\,F'_{c0}(r_{0}), \tag{36}
$$

$$
\frac{F'_{c}(r_{0})}{F_{c}(r_{0})}
\;=\; \frac{F'_{c0}(r_{0})}{F_{c0}(r_{0})}. \tag{37}
$$

Equation (37) is verified by Wolfram[^11]: the residual `F'_c/F_c − F'_{c0}/F_{c0}` simplifies to `0`, since the `(1 + m/M)` factor cancels in `F'/F`. **Pass on algebra.**

Here Corda's own apparatus contradicts his own §2 conclusion. The apsidal-angle formula (31) depends *only* on the ratio `F'/F`. Equation (37) says that constant rescaling of `F_{c0}` leaves `F'/F` unchanged. Therefore the apsidal angle for the rescaled force `F_c = (1 + m/M)\,F_{c0}` is the *same* as for the unscaled force `F_{c0}`, which is `2π` (from Eq. 31 with `F_{c0} = F_G = −k/r²`). The "perturbed" force still gives a closed orbit. There is no precession.

Corda's subsequent Eqs. (38)–(41) derive `T = T_0 (1 + m/M)^(−1/2)` again from the period formula (28), then `ω = ω_0 (1 + m/M)^(1/2)`, then `φ = ωT_0 ≈ 2π(1 + m/(2M))` — the same `πm/M` as §2. This is exactly the §2 mistake repeated: a period change rephrased as a perihelion advance. Corda's own Eq. (31), applied to his own Eq. (32) force, contradicts the §3 conclusion. **Algebra pass; interpretation refuted by the apparatus.**

### §2.8 Eqs. (42)–(57) — reduced-mass Kepler derivation

Corda's §4 ("Third Kepler's law") works through the standard centre-of-mass plus reduced-mass two-body decomposition:

$$
\vec R \;=\; \frac{m\,\vec r_{m} + M\,\vec r_{M}}{M + m},
\qquad
\vec r \;=\; \vec r_{m} - \vec r_{M}, \tag{43}
$$

$$
J \;=\; \mu\,r^{2}\,\omega \;=\; 2\mu\,\dot A, \tag{45}
$$

$$
T \;=\; 2\pi\,\sqrt{\frac{a^{3}}{G\,M_{T}}},
\qquad
M_{T} \;=\; M + m, \tag{53}
$$

$$
\frac{a^{3}}{T^{2}} \;=\; \frac{G(M+m)}{4\pi^{2}}, \tag{54}
$$

$$
T \;=\; T_{0}\,(1 + m/M)^{-1/2}. \tag{57}
$$

The algebra is correct — standard textbook treatment (Goldstein §3.7, MTW §25). It lands at the same `T = T_0 (1 + m/M)^(−1/2)` as in (17) and as in §2.7's (39). **Pass.**

**Physical-intuition status, for the author's verdict.** The §4 derivation is impeccable as a derivation of the period correction due to the two-body kinematics. That is what the reduced-mass treatment gives. Corda then writes (59), `φ = ω T_0 ≈ 2π(1 + m/(2M))`, and calls the excess a "perihelion advance" — the same conflation as in §2.7.

### §2.9 Numerical summary: the three approaches in Corda §§2–4 all give the same `Δφ`

The numerical evaluation is performed[^12] with `G = 6.6743 × 10⁻¹¹`, `M_⊙ = 1.98892 × 10³⁰` kg, `m_⊙^{Mercury} = 3.3011 × 10²³` kg, `a_Mercury = 5.7909 × 10¹⁰` m, and `T_orb = 87.969 · 86400` s. With `orbits/century = 100 · 365.25 · 86400 / T_orb` and `rad → arcsec = 180/π · 3600`, the quantity `Δφ = π m/M` per orbit becomes

$$
\Delta\varphi \;=\; \pi\,m/M
\;=\; 5.21424 \times 10^{-7}\,\text{rad}, \tag{R2}
$$

and per century becomes

$$
\Delta\varphi
\;=\; \pi\,(m/M) \cdot \text{(orbits/century)} \cdot \text{(rad → arcsec)}
\;=\; 44.66\,\text{arcsec/century}. \tag{R3}
$$

The arithmetic in (R3) matches Corda's reported `44.39` arcsec/century to within rounding[^12]. The numerical match to the observed `43` arcsec/century is correct arithmetic on a quantity that, per §2.6 and §2.7, is not a perihelion advance.

### §2.10 Eq. §5 — breakdown for Venus and Earth, the paper's self-refutation

Corda himself shows that the `πm/M` formula fails for non-Mercury planets. The numerical evaluation[^13] gives, with `m_Venus = 4.87 × 10²⁴`, `T_Venus = 224.701 · 86400`, `m_Earth = 5.97 × 10²⁴`, and `T_Earth = 365.256 · 86400`:

$$
\Delta\varphi^{\text{Venus}}_{\text{Corda}}
\;=\; \pi\,(m_{\text{Venus}}/M_{\odot})
\cdot \text{(orbits/century)}_{\text{Venus}}
\cdot \text{(rad → arcsec)}
\;=\; 257.91\,\text{arcsec/century}, \tag{R4}
$$

$$
\Delta\varphi^{\text{Earth}}_{\text{Corda}}
\;=\; \pi\,(m_{\text{Earth}}/M_{\odot})
\cdot \text{(orbits/century)}_{\text{Earth}}
\cdot \text{(rad → arcsec)}
\;=\; 194.50\,\text{arcsec/century}. \tag{R5}
$$

In summary:

| Planet | `πm/M` (Corda's formula) | Observed residual | Discrepancy |
|---|---|---|---|
| Mercury | `44.66` arcsec/c | `≈ 43` arcsec/c | `1.04×` (the "agreement") |
| Venus | `257.91` arcsec/c | `≈ 8.62` arcsec/c | `30×` too big |
| Earth | `194.50` arcsec/c | `≈ 3.83` arcsec/c | `51×` too big |

A formula that is right by `4%` for one planet and wrong by factors of `30` and `50` for the others is not a physics result. It is a numerical coincidence at one particular value of `m/M`. Corda acknowledges the breakdown (§5) but rationalises it with n-body effects ("Venus has Mercury as an interior planet"). The rationalisation does not explain why a strictly two-body formula would be valid for one two-body system and invalid for others. The honest reading is that `πm/M` is not the perihelion advance for any planet, including Mercury.

### §2.11 Eqs. (60)–(67) — gravitational time dilation

Corda's §6 applies Schwarzschild weak-field corrections. The standard relations are

$$
t_{g} \;=\; \sqrt{g_{00}(r_{0})}\,t_{l}, \tag{60}
$$

with the isotropic-coordinate weak-field Schwarzschild line element

$$
ds^{2} \;=\; \Bigl(1 - \frac{2GM}{r\,c^{2}}\Bigr)(c\,dt)^{2}
\;-\; \Bigl(1 + \frac{2GM}{r\,c^{2}}\Bigr)(dr^{2} + r^{2}\,d\theta^{2} + r^{2}\,\sin^{2}\theta\,d\phi^{2}), \tag{61}
$$

$$
t_{g}
\;=\; \sqrt{1 - r_{g}/r_{0}}\,t_{l}
\;\approx\; \Bigl(1 - \tfrac{1}{2}\,\tfrac{r_{g}}{r_{0}}\Bigr)\,t_{l}, \tag{62}
$$

where `r_g = 2GM/c²` is the Sun's gravitational radius (`≈ 2.95` km). Corda also asserts a corresponding radial-distance contraction `r_0 → r_0 (1 − ½ r_g/r_0)` (Eq. 63), arguing from `r = c t`.

**Physical-intuition status, for the author's verdict — marginal.** The time-dilation factor (62) is standard GR weak-field. The radial-distance contraction in (63) is unusual: in isotropic coordinates, the radial coordinate `r` of a point at Schwarzschild radius `R` is `r ≈ R(1 − GM/(c² R))` (a small isotropic-to-areal coordinate shift), so there is a `½ r_g/r_0` correction when comparing isotropic and areal radial coordinates; treating this as "the comoving Newtonian observer sees a different distance because `r = c t`" is non-standard. Worth flagging. The final-result coefficient in (67) is, however, unaffected to first order.

Substituting (63) into (17),

$$
T_{F}
\;=\; \frac{2\pi\,r_{0}^{3/2}\,\bigl(1 - \tfrac{1}{2}\,\tfrac{r_{g}}{r_{0}}\bigr)^{3/2}\,(1 + m/M)^{-1/2}}{\sqrt{G M}}. \tag{64}
$$

Computing `ω_F = 2π/T_F`,

$$
\omega_{F}
\;\simeq\; \omega\,\bigl(1 - \tfrac{1}{2}\,\tfrac{r_{g}}{r_{0}}\bigr)^{-3/2}\,
\bigl(1 + \tfrac{1}{2}\,\tfrac{r_{g}}{r_{0}}\bigr), \tag{65}
$$

$$
\varphi_{F}
\;\simeq\; 2\pi\,\bigl(1 + \tfrac{5}{4}\,\tfrac{r_{g}}{r_{0}}\bigr), \tag{66}
$$

$$
\Delta\varphi_{F}
\;\simeq\; \tfrac{5\pi}{2}\,\tfrac{r_{g}}{r_{0}}. \tag{67}
$$

The leading-order expansion in (66) is verified by Wolfram[^14]: `Series[2π (1 − x/2)^(−3/2) (1 + x/2), {x, 0, 1}]` returns `2π + (5π/2) x + O(x²)`, with leading coefficient `(5/2) π r_g/r_0`. The algebra is correct, conditional on the substitution (63). The result (67) is roughly two-thirds of the GR answer.

### §2.12 Eqs. (68)–(72) — rotational time dilation (Langevin)

The standard SR derivation of the Langevin line element for a rotating frame:

$$
ds^{2} \;=\; c^{2}\,dt^{2} \;-\; dr^{2} \;-\; r^{2}\,d\phi^{2} \;-\; dz^{2}, \tag{68}
$$

$$
\phi \;=\; \phi' + \omega\,t', \tag{69}
$$

$$
ds^{2} \;=\; \Bigl(1 - \tfrac{r^{2}\,\omega^{2}}{c^{2}}\Bigr)\,c^{2}\,dt^{2}
\;-\; 2\omega\,r^{2}\,d\phi'\,dt'
\;-\; dr^{2} \;-\; r^{2}\,d\phi'^{2}
\;-\; dz'^{2}, \tag{70}
$$

$$
d\tau
\;=\; \sqrt{1 - r^{2}\,\omega^{2}/c^{2}}\,dt
\;\simeq\; \Bigl(1 - \tfrac{1}{2}\,\tfrac{r^{2}\,\omega^{2}}{c^{2}}\Bigr)\,dt. \tag{71}
$$

At `r = r_0`, `ω = ω_0`,

$$
\tau
\;\simeq\; \Bigl(1 - \tfrac{1}{2}\,\tfrac{r_{0}^{2}\,\omega_{0}^{2}}{c^{2}}\Bigr)\,t. \tag{72}
$$

The algebra is textbook special relativity. **Pass.**

**Physical-intuition status, for the author's verdict.** Corda invokes Einstein's Equivalence Principle to argue that the centrifugal field in the rotating frame can be treated "as if" it were a gravitational field. This is the conceptual move that permits the additive combination of gravitational and rotational time dilation. The application of the EEP here is standard. The relevant question is whether the Langevin correction, applied to a Keplerian orbital period, is the right way to recover the GR perihelion advance. The result in §8 suggests that, in the circular limit, it is — and this matches the standard PPN result for the gravitational-redshift plus Doppler combination.

### §2.13 Eqs. (73)–(77) — the combined correction

Using (5)–(6) to express the rotational time-dilation in terms of `r_g`:

$$
\tau
\;=\; \sqrt{1 - \tfrac{1}{2}\,\tfrac{r_{g}}{r_{0}}}\,t
\;\simeq\; \Bigl(1 - \tfrac{1}{4}\,\tfrac{r_{g}}{r_{0}}\Bigr)\,t. \tag{73}
$$

The comoving Newtonian observer replaces `T_F → T_T = T_F (1 − ½ r_g/r_0)^{1/2}`:

$$
T_{T} \;=\; T_{0}\,(1 + m/M)^{-1/2}\,\bigl(1 - \tfrac{1}{2}\,r_{g}/r_{0}\bigr)^{2}, \tag{74}
$$

$$
\omega_{F}
\;=\; 2\pi/T_{T}
\;\simeq\; \omega\,\bigl(1 - \tfrac{1}{2}\,r_{g}/r_{0}\bigr)^{-2}
\;\simeq\; \omega\,\bigl(1 + \tfrac{3}{2}\,r_{g}/r_{0}\bigr), \tag{75}
$$

$$
\varphi_{F}
\;=\; \omega_{F}\,T_{0}\,(1 - m/(2M))
\;\simeq\; 2\pi\,\Bigl(1 + \tfrac{3}{2}\,\tfrac{r_{g}}{r_{0}}\Bigr), \tag{76}
$$

$$
\Delta\varphi_{F}
\;\simeq\; 3\pi\,\tfrac{r_{g}}{r_{0}}. \tag{77}
$$

We flag an OCR / marker-pdf ambiguity. The conversion of Eq. (75) shows the last rotational factor as `(1 − r_g/(2 r_0))^{+1/2}`, which would give a leading coefficient of `1`, not `3/2`. Reading instead as `(1 − r_g/(2 r_0))^{−1/2}`, consistent with `ω = 2π/T` and the `T_T = T_F (1 − r_g/(2 r_0))^{1/2}` substitution, recovers the `3/2`. Either (75) in the paper has a typographical error or the converted markdown does. The final equation (77) is correct and matches the standard PPN gravitational-redshift plus Doppler combination.

The leading-order expansion in (76) is verified by Wolfram[^15] under the `−1/2` reading: `Series[(1 − x/2)^(−3/2) (1 + x/2) (1 − x/2)^(−1/2), {x, 0, 1}]` returns `1 + 3x/2 + O(x²)`, with leading coefficient `(3/2) r_g/r_0`. **Pass on algebra, modulo the suspected typographical error in (75).**

### §2.14 Comparison of Corda's Eq. (77) to standard GR

The standard GR result, from the Schwarzschild geodesic[^16], is

$$
\Delta\varphi_{\text{GR}}
\;=\; \frac{6\pi\,GM}{c^{2}\,a\,(1-e^{2})}
\;=\; \frac{3\pi\,r_{g}}{a\,(1-e^{2})}. \tag{R6}
$$

For Mercury, with `e = 0.20563` and `1 − e² ≈ 0.9577`, the standard result is `Δφ_GR = 42.99` arcsec/century. For the circular limit (`e = 0`, `a = r_0`) the standard result reduces to `3π r_g/r_0`. We observe that Corda's Eq. (77) reproduces the circular-orbit limit of GR. The numerical evaluation[^17] gives

$$
\Delta\varphi_{\text{Corda §8, circular}}
\;=\; 3\pi\,r_{g}/a_{\text{Mercury}}\;\cdot\;\text{(orbits/century)}\;\cdot\;\text{(rad → arcsec)}
\;=\; 41.17\,\text{arcsec/century}, \tag{R7}
$$

$$
\Delta\varphi_{\text{GR, Mercury}}
\;=\; 3\pi\,r_{g}/[\,a_{\text{Mercury}}\,(1 - e_{\text{Mercury}}^{2})\,]\;\cdot\;\text{(orbits/century)}\;\cdot\;\text{(rad → arcsec)}
\;=\; 42.99\,\text{arcsec/century}. \tag{R8}
$$

Corda's §§6–8 derivation is missing the `1/(1 − e²)` eccentricity factor. For Mercury this is a `4.4%` shortfall (`−1.82` arcsec/century).

**For the author's verdict.** This discrepancy is significant for Mercury specifically because of its high eccentricity. For Venus and Earth (more circular), the discrepancy shrinks. A complete derivation of the GR result needs the eccentricity from the genuine orbit, for example via the apsidal-angle formula or via the relativistic correction to the orbit equation.

---

## §3. The load-bearing critique — three findings

We surface three findings, in order of weight.

### Finding 1 (algebra-driven). Bertrand's theorem rules out Corda §§2–4

Corda's apsidal-angle formula (Eq. 31), applied to *any* inverse-square force — including the `(1 + m/M)`-rescaled version of Eq. (32) — gives `φ_0 = 2π` exactly. Indeed, for `F = −k/r²` one has `F'(r) r / F(r) = −2`, and so `φ_0 = 2π/√(3 − 2) = 2π`. For `F = (1 + m/M)(−k/r²)` one has the same `F'/F = −2/r`, and so the same `φ_0 = 2π`.

This is Bertrand's theorem in action. Of all central potentials, only `1/r` (Kepler, Newton) and `r²` (harmonic) give universally closed bound orbits. Inverse-square forces (Newton's gravity, Eq. 9) produce closed ellipses with no precession, at any energy or angular momentum. Constant rescaling of an inverse-square force does not change this; it just changes the orbital period.

It follows that Corda's Eq. (20), `Δφ = πm/M`, is not an apsidal precession. It is the angle by which Mercury's position differs from `2π` when measured at time `T_0` (the unperturbed period) under the perturbed angular velocity `ω`. After the true perturbed period `T = T_0(1 + m/M)^(−1/2)` has elapsed, Mercury has swept exactly `2π` again and has returned to perihelion. **No precession. Refuted.**

### Finding 2 (empirics-driven). The §5 self-refutation

Applying Corda's own formula `Δφ = π m/M × (orbits/century)` to Venus and Earth gives the numbers tabulated in §2.10 above; from (R4) and (R5), the formula predicts `257.91` arcsec/century for Venus against the observed `8.62`, and `194.50` arcsec/century for Earth against the observed `3.83`. The Mercury "agreement" of `44.66` against `43` is the lone case in which the formula and the data line up.

Corda's response, that n-body effects are responsible and that Mercury has no interior planets, does not work as a defence. The formula `Δφ = πm/M` makes no reference to n-body effects. The derivation in §§2–4 is strictly two-body throughout, with no interior planets entering. If the formula were the correct two-body result, it would apply to any two-body sub-problem, and the n-body corrections would be smaller perturbations rather than discrepancies by a factor of `30` to `50`.

**For the author's verdict.** The strongest available evidence that `πm/M` is not the correct two-body perihelion advance is that the formula is empirically wrong for every planet other than Mercury. The Mercury "agreement" is coincidence: the value `m/M ≈ 1.66 × 10⁻⁷`, multiplied by Mercury's `415.2` orbits/century, happens to give a number near `43` arcsec/century. Venus and Earth have different `(m/M, orbits/century)` products and give wildly wrong predictions. **Refuted.**

### Finding 3 (formal-derivation-driven). §§6–8 reproduce only the `e = 0` GR limit

Corda's Eq. (77), `Δφ_F = 3π r_g/r_0`, is correct algebra (modulo the §2.13 typographical error) but recovers only the circular-orbit limit of the true GR result `3π r_g/[a (1 − e²)]`. For Mercury (`e = 0.21`) this is a `4.4%` shortfall (`−1.82` arcsec/century).

**For the author's verdict.** §§6–8 is a real first-order GR derivation route. What it gets right is the "PPN sum of gravitational and Doppler time dilation" result, well known in GR pedagogy. What it does not capture is the eccentricity dependence, which makes it numerically equivalent to GR only for nearly circular orbits. For Mercury the eccentricity correction matters; for Venus and Earth (smaller `e`), the expression `3π r_g/a` is much closer to the full GR answer.

The two §§6–8 results work in their domain. The trouble is purely the combination with §§2–4. Corda writes in §9 that "Newtonian plus corrections matches GR." The honest reading is that §§6–8 *alone* matches GR (in the circular limit); §§2–4 contributes a fictitious `πm/M` that does not correspond to a real precession; and the two together do not add up to GR — the `πm/M` is simply spurious. **Marginal on §§6–8 by themselves; refuted on the §9 combination claim.**

---

## §4. Implication for *Relativistic Newtonian Theory*

Gill's *Relativistic Newtonian Theory* extends Corda's framework by adding the dual-theory `(1 + V/(m c²))` factor to the force law:

$$
\vec a_{m}
\;=\; -\,\frac{GM}{r^{2}}\,\Bigl(1 - \frac{GM}{c^{2}\,r}\Bigr)\,\hat u_{r}
\qquad\text{(framework, Eq. (h4)).}\tag{R9}
$$

The expression (R9) is transcribed from the pull request #83 verification document[^18]. All algebra in that document is Wolfram-MCP verified.

### §4.1 What pull request #83 computed (Corda-style template)

Document 03 of the pull request[^19] followed Corda's heuristic and obtained

$$
\Delta\varphi_{d}
\;=\; 2\pi\!\left[\Bigl(1 + \frac{m}{M}\Bigr)
\Bigl(1 - \frac{G(M^{2} + m^{2})}{(M + m)\,c^{2}\,r_{0}}\Bigr)\right]^{1/2} - 2\pi, \tag{R10}
$$

with first-order Taylor expansion

$$
\Delta\varphi_{d_{1}}
\;\simeq\; 2\pi\!\left[\frac{m}{2 M} - \frac{G M\,(1 + m/M)}{2\,c^{2}\,r_{0}}\right]. \tag{R11}
$$

The decomposition is

| Contribution | Value | Source |
|---|---|---|
| Reduced-mass `π m/M` | `+44.66` arcsec/century | Corda §§2–4 (not a real precession) |
| Framework relativistic correction `−π GM/(c² r_0)` | `−6.87` arcsec/century | wrong sign vs. the GR `+42.99` arcsec/century |
| Net `Δφ_d` | `+37.79` arcsec/century | sum, vs. observed `≈ 43` arcsec/century |

The `12%` shortfall (`+37.79` vs. `+43`) was interpreted in pull request #83 as a near miss. Under the present analysis, this is a coincidence on top of a coincidence: `πm/M` is fictitious, and the relativistic correction has the wrong sign.

### §4.2 What pull request #83 document 05 computed (apsidal angle, the genuine answer)

Document 05 of the pull request[^20] applies Corda's *own* apsidal-angle formula (Eq. 31) to the framework's modified force law,

$$
F(r) \;=\; -\,\frac{G M m}{r^{2}}\,\Bigl(1 - \frac{GM}{c^{2}\,r}\Bigr), \tag{R12}
$$

and obtains

$$
\Delta\varphi^{\text{genuine}}_{\text{framework}}
\;=\; -\,\frac{1}{6}\,\Delta\varphi_{\text{GR}}
\;=\; -\,7.17\,\text{arcsec/century}. \tag{R13}
$$

The opposite sign and the one-sixth magnitude of GR are exact, not approximate. This is the well-known `1/6` factor of a force-law-only modification of Newtonian gravity that lacks GR's spatial-metric curvature[^21].

The n-body refinement cannot close this gap. Planetary perturbations contribute framework relativistic corrections at the level of `10⁻⁴` to `10⁻⁵` of the Sun–Mercury term, since `GM/(c² r)` scales with the attracting mass and the other planets are light. The `−7.17` arcsec/century is the framework's structural prediction for the eccentric two-body case, computed via the same apparatus (apsidal angle) Corda himself uses.

### §4.3 The recommendation, for the author

The Corda-style `πm/M + relativistic correction` template that *Relativistic Newtonian Theory* inherits should be replaced by the apsidal-angle calculation as the framework's genuine perihelion prediction. The relevant numbers are

| Calculation | Result | Status |
|---|---|---|
| Corda-style heuristic (PR #83 doc 03) | `+37.79` arcsec/century | numerical near miss; structurally fictitious |
| Apsidal angle (PR #83 doc 05) | `−7.17` arcsec/century | genuine precession of the modified force law |
| Standard GR | `+42.99` arcsec/century | observed `≈ 43` arcsec/century |

If the framework's intended Mercury prediction is the apsidal-angle result, the verdict reads: wrong sign and one-sixth of the GR magnitude. This is the structural signature of any modified-Newtonian gravity that adds a `(1 + V/(m c²))`-style factor to the force without modifying the spatial metric. The `12%` near miss recorded in pull request #83 document 03 is an artefact of taking Corda's `πm/M` heuristic as the baseline.

---

## §5. Specific questions for the author

We list six questions, drawn from the analysis above, in the order they fall.

1. **Bertrand's theorem applicability.** Do you agree that Eq. (31), applied to an inverse-square force (including the `(1 + m/M)`-rescaled version of Eq. 32), gives `φ_0 = 2π`, that is, no apsidal precession? This is the load-bearing critique of §§2–4. If you read Bertrand differently, or if you identify a step in the derivation of Eq. (31) that does not apply, please flag it.

2. **Period versus precession.** Do you agree that the `πm/M` quantity in Corda's Eqs. (19)–(20) is the angle Mercury sweeps in the *unperturbed* period `T_0` at the *perturbed* angular velocity `ω`, and not the rotation of the orbit's major axis between successive perihelia? This is the conceptual distinction that determines whether *Relativistic Newtonian Theory* inherits a fictitious baseline.

3. **The §5 self-refutation.** The same `πm/M` formula gives Venus `258` arcsec/century and Earth `195` arcsec/century, against the observed `8.6` and `3.8` arcsec/century. Is Corda's "n-body effects" rationalisation sufficient, or does the empirical failure mean that the formula is not a perihelion advance? Pull request #83's reading is the latter. We would like your reading recorded.

4. **§§6–8 and eccentricity.** Corda's `3π r_g/r_0` reproduces the `e = 0` GR limit and misses `1/(1 − e²)`. For Mercury this is a `4.4%` shortfall. Do you read §§6–8 as a valid pedagogical derivation of the circular-orbit GR result (the PPN sum of gravitational and Doppler time dilation), or as a complete GR derivation? The distinction matters for whether the eccentricity correction needs to be added to the framework's prediction.

5. **The apsidal-angle calculation for Gill's force law.** Document 05 of pull request #83 returns `−7.17` arcsec/century, that is, `−1/6` of the GR value, for the framework's force law `a = −(GM/r²)(1 − GM/(c² r)) ê_r`. This is the proper answer to "what does Gill's modified Newton predict for Mercury?". Do you read this as the framework's intended prediction, or is the Corda-style heuristic (`+37.79` arcsec/century) intended?

6. **The `1/6` structural signature.** Pull request #83 document 05 identifies `−1/6` of GR as the classic signature of a force-law-only modification of Newton without spatial-metric curvature. Does this match your understanding of what *Relativistic Newtonian Theory* delivers, or does the framework include a spatial-metric component that pull request #83 has missed?

The answers to these six adjudicate whether *Relativistic Newtonian Theory*'s Mercury prediction is `+37.79` arcsec/century (Corda template), `−7.17` arcsec/century (apsidal angle), or something else.

---

## §6. Honest framing — what this report claims and does not claim

### What this report claims

Corda's algebra is mostly correct; in particular, equations (17), (20), (31), (37), (67), and (77) are all Wolfram-MCP verified. Corda's §§2–4 result `Δφ = πm/M` is not a perihelion advance; it is a period-and-angular-velocity rescaling that does not change the apsidal angle of a closed inverse-square orbit. This is verified via Corda's own Eq. (31). The §5 numerical breakdown for Venus and Earth is empirical evidence that the formula is not a real perihelion advance. The §§6–8 derivation reproduces the GR result in the circular-orbit limit only. *Relativistic Newtonian Theory*'s `+37.79` arcsec/century prediction inherits the Corda framing's fictitious `πm/M` baseline together with a wrong-sign relativistic correction. The framework's genuine perihelion prediction, computed via the apsidal-angle formula on the modified force law, is `−7.17` arcsec/century, that is, `−1/6` of the GR value.

### What this report does not claim

The report does not claim that Corda's paper is not a valid published work; it is peer-reviewed in *Physics of the Dark Universe* (a respectable Elsevier journal), and the findings here concern the *interpretation* of correctly derived algebra. The report does not claim that *Relativistic Newtonian Theory* as a whole is wrong; the Mercury perihelion is one test, and other tests (the four candidate curved-spacetime extensions tracked in the GPS author-report Q1[^22]) remain open. The report does not claim that the framework cannot deliver a Mercury-matching prediction; the framework may have alternative extensions to gravity beyond the one analysed in *Relativistic Newtonian Theory*, and this report only addresses the specific extension that paper proposes. The report does not claim that the apsidal-angle critique is novel; it is straightforward textbook orbital mechanics. The novelty here is applying it explicitly to Corda's framework and to *Relativistic Newtonian Theory* as an inheritance check.

### What the report asks of the author

The report asks for six adjudications (§5 above). The author's verdict on these determines what the campaign's Mercury verdict should read; the campaign's current pull request #83 documents both the Corda-style heuristic and the apsidal-angle critique, and is held pending the author's reading.

---

## §7. Crocco §5 — human-acceptance section

<!-- TODO: human reviews and fills in — confirms that (a) Bertrand's theorem
correctly rules out perihelion advance from a magnitude rescaling of an
inverse-square force; (b) the §3 reading of Corda's Eq. (31) as "Corda's own
apparatus contradicts his §2 conclusion" is the right framing rather than a
pedantic reading; (c) the "period change does not equal perihelion advance"
distinction in Finding 1 is the load-bearing critique for the author to
adjudicate; (d) the §5 self-refutation reading (Venus 258 arcsec/century,
Earth 195 arcsec/century) is correct empirical evidence rather than a
mis-application of the formula; (e) the §§6–8 status as a circular-orbit GR
derivation route is correctly characterised; (f) the §4 implication for
*Relativistic Newtonian Theory* — that the campaign should adopt the
apsidal-angle calculation as the framework's genuine Mercury prediction in
place of the Corda-style heuristic — is the right direction for the campaign;
and (g) the six questions in §5 are the right ones to ask in order to close
the Mercury perihelion verdict. -->

---

## References

[^1]: Issue #81, `github.com/temoTxt/PyPhysics/issues/81` — the source ticket for the Mercury perihelion campaign and for the present report.

[^2]: C. Corda, *The secret of planets' perihelion between Newton and Einstein*, *Phys. Dark Univ.* **32** (2021) 100834. `doi:10.1016/j.dark.2021.100834`.

[^3]: Repository pdf at `Roadmapping/External_Papers/Secret_Corda.pdf`.

[^4]: Marker-pdf conversion at `Roadmapping/Converted_Markdown/Secret_Corda/Secret_Corda.md`.

[^5]: Pull request #83 on `github.com/temoTxt/PyPhysics`, the *Relativistic Newtonian Theory* / Mercury perihelion analysis branch. The two relevant per-document references are footnotes [^19] and [^20] below.

[^6]: Author-review convention recorded in `CLAUDE.md` at the repository root (the rule "frame findings concerning DRQM I as 'for author review' rather than corrections"); the convention is applied here.

[^7]: Crocco compliance framework, `Roadmapping/Tooling/CROCCO_COMPLIANCE.md`, §5 on substantive-AI disclosure.

[^8]: Wolfram MCP verification: the residual `T − T₀/√(1 + m/M)` for `T = 2π r₀/v`, `v = √(G(M + m)/r₀)`, `T₀ = 2π r₀^(3/2)/√(GM)` simplifies to `0`.

[^9]: Wolfram MCP verification: `Series[2π √(1 + x) − 2π, {x, 0, 2}]` returns `π x − π x²/4 + O(x³)`; the leading term is `π m/M`.

[^10]: Wolfram MCP verification: with `F[r] := −kk/r²`, the ratio `F'(r) r / F(r)` at `r = r₀` is `−2`, and `φ_0 = 2π/√(3 + (−2)) = 2π`. The inverse-square force gives no precession.

[^11]: Wolfram MCP verification: with `Fc[r] := (1 + mMerc/Msun) Fc0[r]`, the expression `Fc'(r₀)/Fc(r₀) − Fc0'(r₀)/Fc0(r₀)` simplifies to `0`.

[^12]: Wolfram MCP numerical evaluation: with `GG = 6.6743 × 10⁻¹¹`, `Msun = 1.98892 × 10³⁰`, `mMerc = 3.3011 × 10²³`, `aMerc = 5.7909 × 10¹⁰`, `Torb = 87.969 · 86400`, `orbitsPerCentury = 100 · 365.25 · 86400 / Torb`, and `radToArcsec = 180/π · 3600`, the per-orbit value `π mMerc/Msun = 5.21424 × 10⁻⁷` rad, and the per-century value `π (mMerc/Msun) · orbitsPerCentury · radToArcsec = 44.6557` arcsec/century. Corda's paper reports `44.39` arcsec/century, within rounding.

[^13]: Wolfram MCP numerical evaluation: with `mVenus = 4.87 × 10²⁴`, `TVenus = 224.701 · 86400`, `mEarth = 5.97 × 10²⁴`, `TEarth = 365.256 · 86400`, the Corda-formula values are `Venus 257.91` arcsec/century and `Earth 194.50` arcsec/century, against the observed `Venus 8.62` arcsec/century and `Earth 3.83` arcsec/century.

[^14]: Wolfram MCP verification: `Series[2π (1 − x/2)^(−3/2) (1 + x/2), {x, 0, 1}]` returns `2π + (5π/2) x + O(x²)`; the leading coefficient is `(5/2) π r_g/r_0`.

[^15]: Wolfram MCP verification under the `−1/2` reading: `Series[(1 − x/2)^(−3/2) (1 + x/2) (1 − x/2)^(−1/2), {x, 0, 1}]` returns `1 + 3 x/2 + O(x²)`; the leading coefficient is `(3/2) r_g/r_0`.

[^16]: C. W. Misner, K. S. Thorne, and J. A. Wheeler, *Gravitation* (Freeman, 1973), §40.5 (the standard Schwarzschild-geodesic derivation of the perihelion advance).

[^17]: Wolfram MCP numerical evaluation: with `aMerc = 5.7909 × 10¹⁰`, `eMerc = 0.20563`, `cc = 299792458`, `rg = 2 GG Msun/cc²`, `orbitsPerCentury = 100 · 365.25 · 86400 / (87.969 · 86400)`, and `radToArcsec = 180/π · 3600`, the circular-limit value `3π rg/aMerc · orbitsPerCentury · radToArcsec = 41.17` arcsec/century, and the full GR value `3π rg/[aMerc (1 − eMerc²)] · orbitsPerCentury · radToArcsec = 42.99` arcsec/century.

[^18]: Pull request #83 verification document `Roadmapping/Equation_Verification/Dual_Newtonian_Theory.md`, on `github.com/temoTxt/PyPhysics/pull/83`. All algebra in that document is Wolfram-MCP verified.

[^19]: Pull request #83 document 03, `Roadmapping/Mercury_Perihelion/03_numerical_predictions.md`, on `github.com/temoTxt/PyPhysics/pull/83`.

[^20]: Pull request #83 document 05, `Roadmapping/Mercury_Perihelion/05_nbody_orbital_mechanics.md`, on `github.com/temoTxt/PyPhysics/pull/83`. The apsidal-angle calculation in that document yields the `−1/6` of GR result.

[^21]: C. M. Will, *Theory and Experiment in Gravitational Physics*, 2nd ed. (Cambridge University Press, 2018), §7.3, on the PPN parameter analysis; see also C. W. Misner, K. S. Thorne, and J. A. Wheeler, *Gravitation* (Freeman, 1973), §40.6.

[^22]: GPS author-report on the curved-spacetime-extension question: `Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.md`, Q1.
