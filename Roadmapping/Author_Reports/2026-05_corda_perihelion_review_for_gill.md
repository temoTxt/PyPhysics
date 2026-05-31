---
title: "Corda (2021) 'The secret of planets' perihelion between Newton and Einstein' — math + physical-intuition check for the *Relativistic Newtonian Theory* foundations"
---

# Author-review report: Corda (2021) "The secret of planets' perihelion between Newton and Einstein" — math + physical-intuition check

**To:** Tepper Gill (DRQM I + *Relativistic Newtonian Theory* senior author)
**From:** Trey Morris
**Date:** 2026-05-30
**Re:** Issue [#81](https://github.com/temoTxt/PyPhysics/issues/81), basis for Tepper's Mercury perihelion derivation
**Source paper under review:**
Christian Corda, *Phys. Dark Univ.* **32** (2021) 100834,
[doi:10.1016/j.dark.2021.100834](https://doi.org/10.1016/j.dark.2021.100834).
PDF in repo at [`Roadmapping/External_Papers/Secret_Corda.pdf`](../External_Papers/Secret_Corda.pdf);
marker-pdf conversion at
[`Roadmapping/Converted_Markdown/Secret_Corda/Secret_Corda.md`](../Converted_Markdown/Secret_Corda/Secret_Corda.md).

---

## Cover note (for author review, not for journal submission)

You indicated that the perihelion derivation in your *Relativistic Newtonian
Theory* paper takes Corda's framework as its starting point. This report
checks Corda's math and physical intuition end-to-end so you can adjudicate
whether the framework's intended Mercury prediction is the Corda-style
`πm/M + relativistic correction` template (the chain analysed in PR
[#83](https://github.com/temoTxt/PyPhysics/pull/83) giving `+37.79″/century`) or the proper apsidal-angle
calculation of the modified force law (PR #83 doc 05 giving `−7.17″/century`).

Per CLAUDE.md author-review convention, the findings below are surfaced for
your verdict — not announced as corrections. Every algebraic claim is
verified symbolically via Wolfram MCP; every interpretive claim is flagged
as substantive AI and tagged with author-review markers `[for Tepper:]`.

Per Crocco §5 compliance, this is substantive AI. The prompt-of-record is
the user instruction sequence on issue #81 and the user's 2026-05-30 ask
"write the complete analysis in a markdown file for further review by
Tepper." No prose generation prompt beyond the user's instructions was used.

---

## §1. Corda's paper — structure (8 pages, 77 numbered equations)

| § | Topic | Headline equation | Result |
|---|---|---|---|
| §1 | Introduction, GR perihelion | Eq. (1): `Δφ_GR = 24π³ a² / [T₀² c² (1−e²)] = 3π r_g / [a(1−e²)]` | `42.99″/c` for Mercury |
| §2 | "Circular orbit" Newtonian back-reaction | Eq. (20): `Δφ = πm/M` | `44.66″/c` (Approach 1) |
| §3 | Harmonic-oscillator / apsidal-angle | Eq. (31): `φ₀ = 2π[3 + F'(r)r/F(r)]^(−1/2)` | same `44.66″/c` (Approach 2) |
| §4 | Reduced-mass Kepler | Eq. (57): `T = T₀ (1+m/M)^(−1/2)` | same `44.66″/c` (Approach 3) |
| §5 | Breakdown of `πm/M` for Venus / Earth | (numerical only) | `30×` and `51×` too large |
| §6 | Gravitational time dilation | Eq. (67): `Δφ_F = (5π/2) r_g/r₀` | `34.31″/c` (TD-only) |
| §7 | Rotational (Langevin) time dilation | Eq. (72): `τ ≈ (1 − ½ r²ω²/c²) t` | (intermediate) |
| §8 | Combined gravitational + rotational TD | Eq. (77): `Δφ_F = 3π r_g/r₀` | `41.17″/c` (matches e=0 GR) |
| §9 | Conclusion: "Newtonian + corrections matches GR" | — | — |

**The paper's two thesis claims (§9 verbatim):**

> (i) It is not correct that Newtonian theory cannot predict the anomalous
> rate of precession of the perihelion of planets' orbit. The real problem
> is instead that a pure Newtonian prediction is too large.
>
> (ii) Perihelion's precession can be achieved with the same precision of
> general relativity by extending Newtonian gravity through the inclusion of
> gravitational and rotational time dilation effects.

The math + physical-intuition check below addresses each thesis separately.

---

## §2. Every load-bearing equation surfaced and verified

Each subsection below states Corda's equation verbatim (translated from the
marker-pdf conversion), notes the Wolfram-MCP-verified algebra, and flags
the physical-intuition status.

### §2.1 Eqs. (3)–(5) — circular-orbit Newton baseline (no back-reaction)

$$
\frac{GMm}{r_0^2} = \frac{m v_0^2}{r_0}     \qquad (3)
$$
$$
v_0 = \sqrt{GM/r_0}                          \qquad (4)
$$
$$
T_0 = \frac{2\pi r_0}{v_0} = \frac{2\pi r_0^{3/2}}{\sqrt{GM}}  \qquad (5)
$$
Standard textbook (Newton + centripetal). ✅ Algebraically correct.
Physical content: Mercury orbits a stationary Sun. Standard.

### §2.2 Eq. (6)–(7) — Newtonian angular velocity and unperturbed angle

$$
\omega_0 = 2\pi / T_0                         \qquad (6)
$$
$$
\varphi_0 = \omega_0 T_0 = 2\pi               \qquad (7)
$$
Definitionally correct. ✅

### §2.3 Eqs. (8)–(13) — the **substitution** that drives Corda's §2 result

Corda's *core substitution*: replace `M` in Eq. (3) with `M + m` because the
Sun-centred observer sees Mercury's reaction force on the Sun too.

$$
\frac{G(M+m)m}{r_0^2} = \frac{m v^2}{r_0}     \qquad (8)
$$
He justifies this by writing the Newton's-third-law pair (Eqs. 9–13):

$$
\vec F_G = \frac{GMm}{r_0^2}\hat u_r           \qquad (9)
$$
$$
M a_s\hat u_r = \frac{GMm}{r_0^2}\hat u_r \;\Rightarrow\; a_s = \frac{Gm}{r_0^2}  \qquad (10)
$$
$$
m a_m\hat u_r = -\frac{GMm}{r_0^2}\hat u_r \;\Rightarrow\; a_m = -\frac{GM}{r_0^2}  \qquad (11)
$$
$$
a\hat u_r \equiv a_m - a_s = -\frac{G(M+m)}{r_0^2}\hat u_r  \qquad (12)
$$
$$
F\hat u_r = -\frac{G(M+m)m}{r_0^2}\hat u_r     \qquad (13)
$$
**Algebra ✅** — this is the standard two-body relative-acceleration result.
The same `G(M+m)` shows up in textbook reduced-mass treatments (Goldstein
§3.1, MTW §25). What `(M+m)` captures is the kinematic fact that the
*relative* acceleration of Mercury w.r.t. the Sun is bigger than `GM/r²` by
the factor `(M+m)/M = 1 + m/M`.

**Physical intuition status [for Tepper:]** what is being computed here is
the **relative** acceleration in the centre-of-mass frame, not a "Sun-fixed"
force on Mercury. The standard reduced-mass two-body decomposition replaces
the two-body problem with a one-body problem for a reduced-mass particle
`μ = Mm/(M+m)` orbiting a fixed centre at `M + m`. Corda's Eq. (12) is
correct, but the interpretation in §2 (that this constitutes a "perturbation"
to Mercury's force) muddles the standard story.

### §2.4 Eqs. (14)–(17) — perturbed velocity and period

$$
v = \sqrt{G(M+m)/r_0}                          \qquad (14)
$$
$$
T = \frac{2\pi r_0}{v} = \frac{2\pi r_0^{3/2}}{\sqrt{G(M+m)}}  \qquad (15)
$$
$$
(M+m)^{-1/2} = M^{-1/2}\,(1 + m/M)^{-1/2}      \qquad (16)
$$
$$
T = T_0\,(1 + m/M)^{-1/2}                      \qquad (17)
$$

**Wolfram MCP check:**
```mathematica
ClearAll[GG, Msun, mMerc, r0, T, T0];
v   = Sqrt[GG (Msun + mMerc)/r0];
T   = 2 Pi r0 / v;
T0  = 2 Pi r0^(3/2) / Sqrt[GG Msun];
FullSimplify[T - T0 / Sqrt[1 + mMerc/Msun]]
(* Result: 0  ✅ *)
```

**Algebra ✅.** Period decreases slightly because the effective gravitational
constant in the one-body problem is `G(M+m)` not `GM`. This is a real
kinematic effect on the **sidereal period**.

### §2.5 Eqs. (18)–(20) — angular velocity and "perihelion advance"

$$
\omega = \omega_0\,(1 + m/M)^{1/2}             \qquad (18)
$$
$$
\varphi = \omega T_0 = 2\pi\,(1 + m/M)^{1/2} \simeq 2\pi\,(1 + m/(2M))  \qquad (19)
$$
$$
\boxed{\;\Delta\varphi = \varphi - \varphi_0 \simeq \pi m / M\;}  \qquad (20)
$$

**Wolfram MCP Taylor check:**
```mathematica
Series[2 Pi Sqrt[1 + x] - 2 Pi, {x, 0, 2}]
(* Result: Pi x - Pi x^2/4 + O[x]^3  ✅ leading term is π m/M *)
```

**Algebra ✅.** Numerically for Mercury (`m/M = 1.66 × 10⁻⁷`):
`Δφ = 5.21 × 10⁻⁷ rad/orbit × 415.2 orbits/century × (180·3600/π)`
`= 44.66″/century`. Wolfram-verified.

**🔴 Physical-intuition status [for Tepper:]** here is the load-bearing
critique. **What Eq. (20) computes is the angle Mercury sweeps in the
unperturbed period `T_0` at the perturbed angular velocity `ω`.** That is
not a perihelion advance. The perihelion advance is the angular rotation of
the orbit's major axis between successive perihelia — a property of the
shape of the orbit, not of the period or angular velocity.

For an inverse-square force the orbit is a closed ellipse (Bertrand's
theorem), so the angle between successive perihelia is exactly `2π`,
**regardless of `T` and `ω`**. Changing `G` from `GM` to `G(M+m)` rescales
the orbital speed and period but leaves the orbit a closed ellipse. The
`πm/M` quantity Corda computes is a real kinematic effect, but it is not
the perihelion-advance quantity. The §3 derivation makes this completely
explicit — see §2.7 below.

### §2.6 Eqs. (21)–(31) — Corda's own apsidal-angle formalism (§3 of paper)

Following Price & Rush (1979), Corda writes the radial equation,
imposes angular-momentum conservation, Taylor-expands around `r_0`, and
arrives at the **classical Bertrand-style apsidal angle**:

$$
F_{c0}(r) = m\left(\ddot r - \omega_0^2 r\right)                  \qquad (21)
$$
$$
J_0 = m r^2 \omega_0                                              \qquad (22)
$$
$$
F_{c0}(r) = m\left(\ddot r - \frac{J_0^2}{m^2 r^3}\right)         \qquad (23)
$$
$$
F_{c0}(r_0) = -\frac{J_0^2}{m r_0^3}                              \qquad (24)
$$
After perturbing `x = r − r_0` and series-expanding:
$$
\ddot x + m^{-1}\!\left[-\frac{3 F_{c0}(r_0)}{r_0} - F'_{c0}(r_0)\right] x = 0  \qquad (27)
$$
$$
T_0 = 2\pi\!\left(\frac{m}{-\,3 F_{c0}(r_0)/r_0 - F'_{c0}(r_0)}\right)^{1/2}  \qquad (28)
$$
$$
\frac{\varphi_0}{2} = \pi\!\left(\frac{m}{-3F_{c0}/r_0 - F'_{c0}}\right)^{1/2}\!\frac{J_0}{m r_0^2}  \qquad (29)
$$
Using Eqs. (24) and (29):
$$
\boxed{\;\varphi_0 = 2\pi\left[\,3 + \frac{F'_{c0}(r_0)\,r_0}{F_{c0}(r_0)}\,\right]^{-1/2}\;}  \qquad (31)
$$

**Wolfram MCP for `F = −k/r²`:**
```mathematica
ClearAll[Fc, r, r0, kk];
Fc[r_] := -kk/r^2;
ratio = (Fc'[r] r/Fc[r]) /. r -> r0;
phi0  = 2 Pi / Sqrt[3 + ratio];
{Simplify[ratio], Simplify[phi0]}
(* Result: {-2, 2 Pi}  ✅ For inverse-square: F'r/F = −2, φ_0 = 2π, NO PRECESSION *)
```

**Algebra ✅ — and this contradicts §2.** Corda derives the apsidal-angle
formula from scratch (Eq. 31), specialises to inverse-square in his own text
(line right after Eq. 31: *"by setting F_{c0} = F_G in Eq. (31)… one finds
φ_0 = 2π, which is exactly Eq. (7)"*), and explicitly notes that this is the
closed-orbit result. **For an inverse-square force, the apsidal angle is 2π,
the orbit is closed, and there is no precession.** Verified by Wolfram MCP
([Bertrand's theorem](https://en.wikipedia.org/wiki/Bertrand%27s_theorem) at
work).

### §2.7 Eqs. (32)–(41) — Corda's "promotion" to the (1+m/M) force

To extend Eq. (31)'s apsidal-angle calculation to the reduced-mass case,
Corda makes the substitution:
$$
F_{c0}(r) \;\to\; F_c(r) = \left(1 + \frac{m}{M}\right) F_{c0}(r)  \qquad (32)
$$
$$
\omega_0 \to \omega                                                \qquad (33)
$$
$$
J_0 \to J = m r^2 \omega                                           \qquad (34)
$$
This is meant to mirror the §2 substitution `M \to M+m` inside the central
force law. He then writes:
$$
F'_c(r_0) = \left(1 + \frac{m}{M}\right) F'_{c0}(r_0)              \qquad (36)
$$
$$
\boxed{\;\frac{F'_c(r_0)}{F_c(r_0)} \;=\; \frac{F'_{c0}(r_0)}{F_{c0}(r_0)}\;}  \qquad (37)
$$

**Wolfram MCP check:**
```mathematica
ClearAll[Fc, Fc0, r, r0, mMerc, Msun];
Fc[r_] := (1 + mMerc/Msun) Fc0[r];
FullSimplify[Fc'[r0]/Fc[r0] - Fc0'[r0]/Fc0[r0]]
(* Result: 0  ✅ — the (1+m/M) factor cancels in F'/F *)
```

**Algebra ✅, and 🔴 here is where Corda's apparatus contradicts his own
§2 conclusion.** The apsidal-angle formula Eq. (31) depends *only* on the
ratio `F'/F`. The Eq. (37) Wolfram-verified identity says that constant
rescaling of `F_{c0}` leaves `F'/F` unchanged. Therefore the apsidal angle
for the rescaled force `F_c = (1+m/M) F_{c0}` is *the same* as for the
unscaled force `F_{c0}` — which is `2π` (Eq. 31 with `F_{c0} = F_G = −k/r²`).
**The "perturbed" force still gives a closed orbit. No precession.**

Corda's subsequent Eqs. (38)–(41) derive `T = T_0(1+m/M)^(−1/2)` again from
the period formula (28), then `ω = ω_0(1+m/M)^(1/2)`, then
`φ = ω T_0 ≈ 2π(1 + m/(2M))` — the same `πm/M` as §2. **But this is exactly
the §2 mistake repeated: a period change rephrased as a "perihelion
advance."** Corda's own Eq. (31), applied to his own Eq. (32) force,
*contradicts* the §3 conclusion.

### §2.8 Eqs. (42)–(57) — reduced-mass Kepler derivation

Corda's §4 ("Third Kepler's law") works out the standard centre-of-mass +
reduced-mass two-body decomposition:
$$
\vec R = \frac{m\vec r_m + M\vec r_M}{M + m},\qquad
\vec r = \vec r_m - \vec r_M                                       \qquad (43)
$$
$$
J = \mu r^2 \omega = 2\mu\,\dot A                                  \qquad (45)
$$
$$
T = 2\pi\!\sqrt{\frac{a^3}{G M_T}},\quad M_T = M + m              \qquad (53)
$$
$$
\frac{a^3}{T^2} = \frac{G(M+m)}{4\pi^2}                            \qquad (54)
$$
$$
T = T_0 (1 + m/M)^{-1/2}                                           \qquad (57)
$$

**Algebra ✅** — standard textbook treatment (Goldstein §3.7, MTW §25).
Lands at the same `T = T_0(1+m/M)^(−1/2)` as Eqs. (17) and (39).

**Physical-intuition status [for Tepper:]** §4's derivation is impeccable as
a derivation of the *period correction* due to the two-body kinematics.
That's what the reduced-mass treatment gives. But Corda then writes Eq. (59)
`φ = ωT_0 ≈ 2π(1 + m/(2M))` and calls the excess a "perihelion advance"
— same conflation as §2.7.

### §2.9 Numerical summary of Corda's three approaches (§2/§3/§4 all give the same Δφ)

```mathematica
ClearAll[GG, Msun, mMerc, aMerc, Torb, orbitsPerCentury, radToArcsec, dphiNewton];
GG = 6.6743*^-11; Msun = 1.98892*^30; mMerc = 3.3011*^23;
aMerc = 5.7909*^10; Torb = 87.969*86400.;
orbitsPerCentury = 100 * 365.25 * 86400. / Torb;
radToArcsec = 180/Pi * 3600.;
dphiNewton = Pi mMerc / Msun;
Print["Δφ = πm/M per orbit  = ", N[dphiNewton], " rad"];
Print["Δφ = πm/M per century = ", N[dphiNewton * orbitsPerCentury * radToArcsec], " arcsec"];
(* Result:                                                                       *)
(*   Δφ = πm/M per orbit  = 5.21424×10⁻⁷ rad                                     *)
(*   Δφ = πm/M per century = 44.6557 arcsec  ← Corda's reported 44.39″/c        *)
```

**The numerical match to observed `43″/century` is correct arithmetic on a
quantity that — per §2.6 and §2.7 above — is **not** a perihelion advance.**

### §2.10 §5 — Breakdown for Venus and Earth (the paper's self-refutation)

Corda himself shows the `πm/M` formula fails for non-Mercury planets:

```mathematica
ClearAll[mVenus, mEarth, TVenus, TEarth, orbitsVenus, orbitsEarth];
mVenus = 4.87*^24; TVenus = 224.701*86400.;
mEarth = 5.97*^24; TEarth = 365.256*86400.;
orbitsVenus = 100*365.25*86400./TVenus;
orbitsEarth = 100*365.25*86400./TEarth;
{"Venus",  N[Pi mVenus/Msun * orbitsVenus * radToArcsec], "vs observed 8.62″/c"}
{"Earth",  N[Pi mEarth/Msun * orbitsEarth * radToArcsec], "vs observed 3.83″/c"}
(* Results:                                                                      *)
(*   Venus  formula prediction: 257.91 ″/c   vs observed  8.62″/c  →  30× too big *)
(*   Earth  formula prediction: 194.50 ″/c   vs observed  3.83″/c  →  51× too big *)
```

| Planet | `πm/M` (Corda's formula) | Observed residual | Discrepancy |
|---|---|---|---|
| Mercury | `44.66″/c` | `≈ 43″/c` | ✅ 1.04× ← "agreement" |
| Venus | `257.91″/c` | `≈ 8.62″/c` | 🔴 **30× too big** |
| Earth | `194.50″/c` | `≈ 3.83″/c` | 🔴 **51× too big** |

**A formula that's right by 4% for one planet and wrong by factors of 30–50
for others is not a physics result — it is a numerical coincidence at one
particular `m/M` value.** Corda himself acknowledges the breakdown (§5)
but rationalises it with n-body effects ("Venus has Mercury as an interior
planet"), which does not explain why a strictly two-body formula would be
valid for one two-body system and invalid for others. **The honest reading:
`πm/M` is not the perihelion advance for any planet, including Mercury.**

### §2.11 Eqs. (60)–(67) — gravitational time dilation

Corda's §6 applies Schwarzschild weak-field corrections.
Standard relations:
$$
t_g = \sqrt{g_{00}(r_0)}\, t_l                                     \qquad (60)
$$
Isotropic-coordinate weak-field Schwarzschild line element:
$$
ds^2 = \left(1 - \frac{2GM}{r c^2}\right)(c\,dt)^2
       - \left(1 + \frac{2GM}{r c^2}\right)(dr^2 + r^2 d\theta^2 + r^2 \sin^2\theta\,d\phi^2)
                                                                   \qquad (61)
$$
$$
t_g = \sqrt{1 - r_g/r_0}\, t_l \;\approx\; \left(1 - \tfrac{1}{2}\tfrac{r_g}{r_0}\right) t_l   \qquad (62)
$$
where `r_g = 2GM/c²` is the Sun's gravitational radius (`≈ 2.95 km`).
Corda also asserts a corresponding radial-distance contraction:
`r₀ → r₀(1 − ½ r_g/r₀)` (Eq. 63), arguing from `r = ct`.

**🟡 Physical-intuition status [for Tepper:]** the time-dilation factor
(Eq. 62) is standard GR weak-field. The radial-distance contraction in
Eq. (63) is unusual — in isotropic coordinates the radial coordinate `r` of
a point at Schwarzschild radius `R` is `r ≈ R(1 − GM/(c²R))` (a small
isotropic-to-areal coordinate shift), so there *is* a `½ r_g/r₀` correction
when comparing isotropic to areal radial coordinates, but treating this as
"the CNO sees a different distance because `r = ct`" is non-standard.
Worth flagging but the final-result coefficient in Eq. (67) is unaffected
to first order.

Substituting Eq. (63) into Eq. (17):
$$
T_F = \frac{2\pi r_0^{3/2}\left(1 - \tfrac{1}{2}\tfrac{r_g}{r_0}\right)^{3/2}(1 + m/M)^{-1/2}}{\sqrt{GM}}  \qquad (64)
$$
Computing `ω_F = 2π/T_F`:
$$
\omega_F \simeq \omega\left(1 - \tfrac{1}{2}\tfrac{r_g}{r_0}\right)^{-3/2}\!\left(1 + \tfrac{1}{2}\tfrac{r_g}{r_0}\right)   \qquad (65)
$$
$$
\varphi_F \simeq 2\pi\left(1 + \tfrac{5}{4}\tfrac{r_g}{r_0}\right)   \qquad (66)
$$
$$
\boxed{\;\Delta\varphi_F \;\simeq\; \tfrac{5\pi}{2}\,\tfrac{r_g}{r_0}\;}  \qquad (67)
$$

**Wolfram MCP Taylor check:**
```mathematica
Series[2 Pi (1 - x/2)^(-3/2) (1 + x/2), {x, 0, 1}]
(* Result: 2 Pi + (5 Pi/2) x + O[x]^2   ✅ leading is (5/2)π r_g/r_0 *)
```

**Algebra ✅** given the substitution in Eq. (63). This is roughly
two-thirds of the GR answer.

### §2.12 Eqs. (68)–(72) — rotational time dilation (Langevin)

Standard SR derivation of the Langevin line element for a rotating frame:
$$
ds^2 = c^2 dt^2 - dr^2 - r^2 d\phi^2 - dz^2                       \qquad (68)
$$
$$
\phi = \phi' + \omega t'                                           \qquad (69)
$$
$$
ds^2 = \left(1 - \tfrac{r^2\omega^2}{c^2}\right)c^2 dt^2 - 2\omega r^2 d\phi' dt' - dr^2 - r^2 d\phi'^2 - dz'^2   \qquad (70)
$$
$$
d\tau = \sqrt{1 - r^2\omega^2/c^2}\,dt \;\simeq\; \left(1 - \tfrac{1}{2}\tfrac{r^2\omega^2}{c^2}\right) dt   \qquad (71)
$$
At `r = r_0`, `ω = ω_0`:
$$
\tau \simeq \left(1 - \tfrac{1}{2}\tfrac{r_0^2\omega_0^2}{c^2}\right) t   \qquad (72)
$$

**Algebra ✅** — textbook special relativity.
**Physical-intuition status:** Corda invokes Einstein's Equivalence Principle
to argue that the centrifugal field in the rotating frame can be treated
"as if" it were a gravitational field. This is the conceptual move that
permits combining gravitational + rotational TD additively. **[for Tepper:]**
the EEP application here is standard; the question is whether the Langevin
correction, applied to a Keplerian orbital period, is the right way to
recover the GR perihelion advance. The §8 result suggests yes (in the
circular limit), and this matches the standard PPN result for the
gravitational-redshift + Doppler combination.

### §2.13 Eqs. (73)–(77) — combined correction

Using Eqs. (5)–(6) to express the rotational TD in terms of `r_g`:
$$
\tau = \sqrt{1 - \tfrac{1}{2}\tfrac{r_g}{r_0}}\, t
     \simeq \left(1 - \tfrac{1}{4}\tfrac{r_g}{r_0}\right) t        \qquad (73)
$$
The CNO replaces `T_F \to T_T = T_F (1 - \tfrac{1}{2} r_g/r_0)^{1/2}` and:
$$
T_T = T_0 (1 + m/M)^{-1/2}\,(1 - \tfrac{1}{2} r_g/r_0)^2           \qquad (74)
$$
$$
\omega_F = 2\pi/T_T \simeq \omega\,(1 - \tfrac{1}{2} r_g/r_0)^{-2} \simeq \omega\,(1 + \tfrac{3}{2}\, r_g/r_0)   \qquad (75)
$$
$$
\varphi_F = \omega_F T_0 (1 - m/(2M)) \simeq 2\pi\left(1 + \tfrac{3}{2}\,\tfrac{r_g}{r_0}\right)  \qquad (76)
$$
$$
\boxed{\;\Delta\varphi_F \;\simeq\; 3\pi\,\tfrac{r_g}{r_0}\;}  \qquad (77)
$$

> **OCR/marker-pdf note.** The conversion of Eq. (75) shows the last
> rotational factor as `(1 - r_g/(2r_0))^{+1/2}`, which would give a leading
> coefficient of `1` (not `3/2`). Reading instead as `(1 - r_g/(2r_0))^{-1/2}`
> (consistent with `ω = 2π/T` and the `T_T = T_F (1 - r_g/(2r_0))^{1/2}`
> substitution) recovers the `3/2`. Either Eq. (75) in the paper has a typo
> or the converted markdown does. The final Eq. (77) is correct and matches
> the standard PPN gravitational-redshift + Doppler combination.

**Wolfram MCP Taylor check (assuming `−1/2`):**
```mathematica
Series[(1 - x/2)^(-3/2) (1 + x/2) (1 - x/2)^(-1/2), {x, 0, 1}]
(* Result: 1 + 3x/2 + O[x]^2   ✅ leading is (3/2) r_g/r_0 *)
```

**Algebra ✅** (modulo the suspected typo in Eq. 75).

### §2.14 Comparison of Corda's Eq. (77) to standard GR

Standard GR (Schwarzschild geodesic, MTW §40.5):
$$
\Delta\varphi_{\rm GR} = \frac{6\pi GM}{c^2 a (1-e^2)} = \frac{3\pi r_g}{a(1-e^2)}
$$
For Mercury (`e = 0.20563`, `1 - e^2 ≈ 0.9577`): `Δφ_GR = 42.99″/century`.

For the **circular** limit (`e = 0`, `a = r_0`): `Δφ_GR = 3π r_g/r_0`.
**Corda's Eq. (77) reproduces the circular-orbit limit of GR.**

```mathematica
ClearAll[GG, Msun, aMerc, eMerc, cc, rg, orbitsPerCentury, radToArcsec];
GG = 6.6743*^-11; Msun = 1.98892*^30; aMerc = 5.7909*^10;
eMerc = 0.20563; cc = 299792458.; rg = 2 GG Msun / cc^2;
orbitsPerCentury = 100 * 365.25 * 86400. / (87.969 * 86400.);
radToArcsec = 180/Pi * 3600.;
{"Corda §8 (3π r_g/a, circular)",  3 Pi rg/aMerc * orbitsPerCentury * radToArcsec,
 "GR Mercury  (3π r_g/[a(1-e²)])", 3 Pi rg/(aMerc (1 - eMerc^2)) * orbitsPerCentury * radToArcsec}
(* Results:                                                                  *)
(*   Corda §8 circular:  41.17″/c                                            *)
(*   GR Mercury (e=0.21): 42.99″/c                                           *)
```

Corda's §§6–8 derivation is missing the `1/(1-e²)` eccentricity factor. For
Mercury this is a `4.4%` shortfall (`−1.82″/century`). **[for Tepper:]**
this discrepancy is significant for Mercury specifically because of its
high eccentricity; for Venus and Earth (more circular) the discrepancy
shrinks. A true derivation of the GR result needs the eccentricity from the
genuine orbit (e.g.~via the apsidal-angle formula or the relativistic
correction to the orbit equation).

---

## §3. The load-bearing critique — three findings

### Finding 1 (algebra-driven). Bertrand's theorem rules out §2–§4

Corda's apsidal-angle formula (Eq. 31), applied to *any* inverse-square
force — including the `(1+m/M)`-rescaled version of Eq. (32) — gives
`φ_0 = 2π` exactly:
- `F = −k/r²` ⇒ `F'(r) r / F(r) = −2` ⇒ `φ_0 = 2π/√(3 − 2) = 2π`. ✅
- `F = (1+m/M)(−k/r²)` ⇒ same `F'/F = −2/r` ⇒ same `φ_0 = 2π`. ✅

This is Bertrand's theorem: of all central potentials,
*only* `1/r` (Kepler / Newton) and `r²` (harmonic) give closed bound orbits
universally. Inverse-square forces (Newton's gravity, Eq. 9) produce closed
ellipses with no precession at any energy or angular momentum. Constant
rescaling of an inverse-square force does not change this — it just changes
the orbital period.

**Therefore Eq. (20)'s `Δφ = πm/M` is not an apsidal precession.** It is
the angle by which Mercury's *position* differs from `2π` when measured at
time `T_0` (the unperturbed period) under the perturbed angular velocity
`ω`. After the true perturbed period `T = T_0(1+m/M)^(−1/2)` has elapsed,
Mercury has swept exactly `2π` again and returned to perihelion. **No
precession.**

### Finding 2 (empirics-driven). The §5 self-refutation

Applying Corda's own formula `Δφ = πm/M × (orbits/century)` to Venus and
Earth:

| Planet | `m × orbits/c × radToArcsec / Msun × π` | Observed | Discrepancy |
|---|---|---|---|
| Mercury | `44.66″/c` | `43″/c` | ✅ |
| Venus | `257.91″/c` | `8.62″/c` | 🔴 **`30×`** |
| Earth | `194.50″/c` | `3.83″/c` | 🔴 **`51×`** |

Corda's response ("n-body effects, Mercury has no interior planets") doesn't
work as a defence:
- The formula `Δφ = πm/M` makes no reference to n-body effects.
- The derivation in §§2/3/4 is strictly two-body throughout (no interior
  planets enter).
- If the formula were the correct two-body result, it would apply to any
  two-body sub-problem and the n-body corrections would be smaller
  perturbations — not corrections of `factor 30–50`.

**[for Tepper:]** the strongest available evidence that `πm/M` is not the
correct two-body perihelion advance is that the formula is empirically
wrong for every planet other than Mercury. The Mercury "agreement" is
coincidence at `m/M ≈ 1.66 × 10⁻⁷` happening to give a number near `43″/c`
once multiplied by Mercury's `415.2 orbits/century`. Venus and Earth have
different `(m/M, orbits/c)` products and give wildly wrong predictions.

### Finding 3 (formal-derivation-driven). §§6–8 reproduce only the e=0 GR limit

Corda's Eq. (77) `Δφ_F = 3π r_g/r_0` is correct algebra (modulo the §2.13
typo) but recovers only the circular-orbit limit of the true GR result
`3π r_g/[a(1−e²)]`. For Mercury (e = 0.21) this is a 4.4% shortfall
(`−1.82″/century`).

**[for Tepper:]** §§6–8 is a real first-order GR derivation route; what it
gets right is the "PPN sum of gravitational + Doppler time dilation"
result, which is well-known in GR pedagogy. What it doesn't capture is the
eccentricity dependence, which makes it numerically equivalent to GR only
for nearly-circular orbits. For Mercury, the eccentricity correction
matters; for Venus and Earth (smaller `e`), `3π r_g/a` is much closer to
the full GR answer.

**The two §§6–8 results work in their domain.** The trouble is purely the
combination with §§2–4. Corda writes (§9) that "Newtonian + corrections
matches GR." The honest reading: §§6–8 *alone* matches GR (in the circular
limit); §§2–4 contribute a fictitious `πm/M` that doesn't correspond to a
real precession; the two together don't add up to GR — the `πm/M` is just
spurious.

---

## §4. Implication for Tepper's *Relativistic Newtonian Theory*

Gill's *Relativistic Newtonian Theory* extends Corda's framework by adding
the dual-theory `(1 + V/(mc²))` factor to the force law:

$$
a_m = -\frac{GM}{r^2}\left(1 - \frac{GM}{c^2 r}\right)\hat u_r   \qquad (\text{Eq. h4})
$$

(transcribed from PR [#83](https://github.com/temoTxt/PyPhysics/pull/83)'s
[`Equation_Verification/Dual_Newtonian_Theory.md`](../Equation_Verification/Dual_Newtonian_Theory.md);
all algebra Wolfram-MCP verified).

### §4.1 What PR #83 computed (Corda-style template)

PR #83 [`03_numerical_predictions.md`](../Mercury_Perihelion/03_numerical_predictions.md)
followed Corda's heuristic:

$$
\Delta\varphi_d = 2\pi\!\left[\left(1 + \tfrac{m}{M}\right)
                  \left(1 - \tfrac{G(M^2 + m^2)}{(M+m) c^2 r_0}\right)\right]^{1/2} - 2\pi
$$

with first-order Taylor expansion:

$$
\Delta\varphi_{d_1} \;\simeq\; 2\pi\!\left[\frac{m}{2M} - \frac{GM(1+m/M)}{2 c^2 r_0}\right]
$$

This decomposes as:

| Contribution | Value | Source |
|---|---|---|
| Reduced-mass `πm/M` | `+44.66″/c` | Corda §§2–4 (**not a real precession**) |
| Framework relativistic correction `−πGM/(c² r_0)` | **`−6.87″/c`** | wrong sign vs. GR's `+42.99″/c` |
| Net `Δφ_d` | **`+37.79″/c`** | sum, vs. observed `≈ 43″/c` |

The 12% shortfall (`+37.79″` vs `+43″`) was interpreted in PR #83 as a
near-miss. **Under the present analysis, this is a coincidence on top of a
coincidence**: `πm/M` is fictitious, and the relativistic correction has
the wrong sign.

### §4.2 What PR #83 doc 05 computed (apsidal-angle, the genuine answer)

[`Mercury_Perihelion/05_nbody_orbital_mechanics.md`](../Mercury_Perihelion/05_nbody_orbital_mechanics.md)
applied Corda's *own* apsidal-angle formula (Eq. 31) to Gill's modified
force law:

$$
F(r) = -\frac{GMm}{r^2}\left(1 - \frac{GM}{c^2 r}\right)
$$

The apsidal-angle formula gives a precession per orbit that is:

$$
\boxed{\;\Delta\varphi_{\rm framework}^{\rm genuine}
       \;=\; -\frac{1}{6}\,\Delta\varphi_{\rm GR}
       \;=\; -7.17\,\text{″/century}\;}
$$

**Opposite sign, one-sixth magnitude of GR exactly.** This is the famous
"1/6 factor" of a force-law-only modification of Newtonian gravity that
lacks GR's spatial-metric curvature. (See e.g. Misner-Thorne-Wheeler
*Gravitation* §40.6 for the analogous PPN parameter analysis.)

The N-body refinement cannot close this gap: planetary perturbations
contribute framework relativistic corrections at `~10^{-4}–10^{-5}` of the
Sun–Mercury term (since `GM/(c²r)` scales with the attracting mass and the
other planets are light). The `−7.17″/century` is the framework's structural
prediction for the eccentric two-body case, computed via the same
apparatus (apsidal angle) Corda himself uses.

### §4.3 The recommendation [for Tepper]

The Corda-style `πm/M + relativistic correction` template that *Relativistic
Newtonian Theory* inherits should be replaced by the apsidal-angle
calculation as the framework's genuine perihelion prediction. The relevant
numbers:

| Calculation | Result | Status |
|---|---|---|
| Corda-style heuristic (PR #83 doc 03) | `+37.79″/century` | numerical near-miss; structurally fictitious |
| Apsidal-angle (PR #83 doc 05) | `−7.17″/century` | genuine precession of modified force law |
| Standard GR | `+42.99″/century` | observed `~43″/c` |

If the framework's intended Mercury prediction is the apsidal-angle result,
the verdict is: **wrong sign, `1/6` magnitude**. This is the structural
signature of any modified-Newtonian gravity that adds a `(1 + V/(mc²))`-style
factor to the force without modifying the spatial metric. The 12% near-miss
in PR #83 doc 03 is an artefact of taking Corda's `πm/M` heuristic as the
baseline.

---

## §5. Specific questions for Tepper

Six questions for your verdict, drawn from the analysis above:

1. **Bertrand's theorem applicability.** Do you agree that Eq. (31) applied
   to an inverse-square force (including the `(1+m/M)`-rescaled version of
   Eq. 32) gives `φ_0 = 2π`, i.e.~no apsidal precession? This is the
   load-bearing critique of §§2–4. If you read Bertrand differently or
   identify a step in Eq. (31)'s derivation that doesn't apply, please flag.

2. **Period vs. precession.** Do you agree that the `πm/M` in Corda's
   Eqs. (19)–(20) is computing the angle Mercury sweeps in the *unperturbed*
   period `T_0` at the *perturbed* angular velocity `ω` — not the rotation
   of the orbit's major axis between successive perihelia? This is the
   conceptual distinction that determines whether *Relativistic Newtonian
   Theory* inherits a fictitious baseline.

3. **§5 self-refutation.** The same `πm/M` formula gives `Venus 258″/c` and
   `Earth 195″/c` vs. observed `8.6″/c` and `3.8″/c`. Is Corda's "n-body
   effects" rationalisation sufficient, or does the empirical failure mean
   the formula is not a perihelion advance? (PR #83's reading is the
   latter; we'd like your reading recorded.)

4. **§§6–8 vs. eccentricity.** Corda's `3π r_g/r_0` reproduces the e=0 GR
   limit and misses `1/(1−e²)`. For Mercury this is a 4.4% shortfall. Do you
   read §§6–8 as a valid pedagogical derivation of the circular-orbit GR
   result (PPN sum of gravitational + Doppler TD), or as a complete GR
   derivation? The distinction matters for whether the eccentricity
   correction needs to be added to the framework's prediction.

5. **The apsidal-angle calculation for Gill's force law.** PR #83 doc 05's
   apsidal-angle result `−7.17″/century` (`−1/6` of GR exactly) for the
   framework's force law `a = −(GM/r²)(1 − GM/(c²r))·ê_r` is the proper
   answer to "what does Gill's modified Newton predict for Mercury?". Do
   you read this as the framework's intended prediction, or is the
   Corda-style heuristic (`+37.79″/century`) intended?

6. **The `1/6` structural signature.** PR #83 doc 05 identifies `−1/6` of
   GR as the classic signature of a force-law-only modification of Newton
   without spatial-metric curvature. Does this match your understanding of
   what *Relativistic Newtonian Theory* delivers, or does the framework
   include a spatial-metric component that PR #83 missed?

The answers to these six adjudicate whether *Relativistic Newtonian
Theory*'s Mercury prediction is `+37.79″/century` (Corda template),
`−7.17″/century` (apsidal-angle), or something else.

---

## §6. Honest framing — what this report does and does not claim

### What this report claims

- Corda's algebra is mostly correct (specifically: Eqs. 17, 20, 31, 37, 67,
  77 are all Wolfram-MCP verified).
- Corda's §§2–4 result `Δφ = πm/M` is *not* a perihelion advance — it is a
  period-and-angular-velocity rescaling that does not change the apsidal
  angle of a closed inverse-square orbit. This is verified via Corda's own
  Eq. (31).
- The §5 numerical breakdown for Venus and Earth is empirical evidence that
  the formula is not a real perihelion advance.
- §§6–8 reproduce the GR result in the circular-orbit limit only.
- *Relativistic Newtonian Theory*'s `+37.79″/century` prediction inherits
  the Corda framing's fictitious `πm/M` baseline plus a wrong-sign
  relativistic correction.
- The framework's *genuine* perihelion prediction, computed via the
  apsidal-angle formula on the modified force law, is `−7.17″/century`
  (`−1/6` of GR).

### What this report does not claim

- That Corda's paper is not a valid published work. It is a peer-reviewed
  paper in *Physics of the Dark Universe* (a respectable Elsevier journal).
  The findings here are about the *interpretation* of correctly-derived
  algebra.
- That *Relativistic Newtonian Theory* as a whole is wrong. The Mercury
  perihelion is one test; other tests (the four candidate curved-spacetime
  extensions tracked in the GPS author report's Q1) remain open.
- That the framework cannot deliver a Mercury-matching prediction. The
  framework may have alternative extensions to gravity beyond the one
  analysed in *Relativistic Newtonian Theory*; this report only addresses
  the specific extension that paper proposes.
- That the apsidal-angle critique is novel. It is straightforward
  textbook orbital mechanics. The novelty here is applying it explicitly
  to Corda's framework and to *Relativistic Newtonian Theory* as an
  inheritance check.

### What the report is asking from you

Six adjudications (§5 above). Your verdict on these determines what the
campaign's Mercury verdict should read; the campaign's current PR #83
documents both the Corda-style heuristic and the apsidal-angle critique
and is held pending your reading.

---

## §7. Cross-references

- **This document**: [`Roadmapping/Author_Reports/2026-05_corda_perihelion_review_for_gill.md`](2026-05_corda_perihelion_review_for_gill.md)
- **Corda PDF**: [`Roadmapping/External_Papers/Secret_Corda.pdf`](../External_Papers/Secret_Corda.pdf)
- **Corda markdown** (marker-pdf): [`Roadmapping/Converted_Markdown/Secret_Corda/Secret_Corda.md`](../Converted_Markdown/Secret_Corda/Secret_Corda.md)
- **Concise equation-verification doc**: [`Roadmapping/Equation_Verification/Secret_Corda_Perihelion.md`](../Equation_Verification/Secret_Corda_Perihelion.md)
- **Gill's *Relativistic Newtonian Theory* analysis (PR #83)**:
  - PR: https://github.com/temoTxt/PyPhysics/pull/83
  - [`Roadmapping/Equation_Verification/Dual_Newtonian_Theory.md`](../Equation_Verification/Dual_Newtonian_Theory.md)
  - [`Roadmapping/Mercury_Perihelion/README.md`](../Mercury_Perihelion/README.md)
  - [`Roadmapping/Mercury_Perihelion/05_nbody_orbital_mechanics.md`](../Mercury_Perihelion/05_nbody_orbital_mechanics.md) (apsidal-angle calculation; the `−1/6` of GR result)
- **GPS author report Q1 (curved-spacetime extension question)**: [`Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.md`](2026-05_gps_relativity_summary_for_gill.md)
- **Issue #81**: https://github.com/temoTxt/PyPhysics/issues/81

## §8. Human-acceptance section (Crocco §5 substantive AI)

<!-- TODO: human reviews and fills in — confirms that (a) Bertrand's theorem
correctly rules out perihelion advance from a magnitude rescaling of an
inverse-square force; (b) the §3 reading of Corda's Eq. (31) as "Corda's own
apparatus contradicts his §2 conclusion" is the right framing rather than a
pedantic reading; (c) the "period change ≠ perihelion advance" distinction
in Finding 1 is the load-bearing critique for Tepper to adjudicate; (d) the
§5 self-refutation reading (Venus 258″/c, Earth 195″/c) is correct empirical
evidence rather than a mis-application of the formula; (e) the §§6–8 status
as a circular-orbit GR derivation route is correctly characterised; (f) the
§4 implication for *Relativistic Newtonian Theory* — that the campaign should
adopt the apsidal-angle calculation as the framework's genuine Mercury
prediction in place of the Corda-style heuristic — is the right direction
for the campaign; and (g) the six questions in §5 are the right ones to ask
to close the Mercury perihelion verdict. -->
