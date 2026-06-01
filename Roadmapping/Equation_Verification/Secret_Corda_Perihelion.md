# Equation verification — Corda (2021), "The secret of planets' perihelion between Newton and Einstein"

**Source:** Christian Corda, *Phys. Dark Univ.* **32** (2021) 100834,
[doi:10.1016/j.dark.2021.100834](https://doi.org/10.1016/j.dark.2021.100834).
Attached to [issue #81](https://github.com/temoTxt/PyPhysics/issues/81); converted via
[`parse_papers.py`](../parse_papers.py) (marker-pdf) into
[`Roadmapping/Converted_Markdown/Secret_Corda/Secret_Corda.md`](../Converted_Markdown/Secret_Corda/Secret_Corda.md).

**Why this matters for the campaign.** Per Tepper Gill's note on issue #81, his
*Relativistic Newtonian Theory* (Gill's *Dual Newtonian Theory*, the paper
analysed in PR [#83](https://github.com/temoTxt/PyPhysics/pull/83)) takes Corda's framework as
its starting point. The Mercury perihelion derivation in
[`Roadmapping/Mercury_Perihelion/`](../Mercury_Perihelion/) inherits Corda's
two key moves: (a) identifying `πm/M` as a "Newtonian precession" from
back-reaction (Corda §§2–4), and (b) adding gravitational + rotational time
dilation to recover the GR result (Corda §§6–8). If Corda's framing has
errors, Tepper's derivation inherits them.

This document checks every load-bearing equation symbolically and physically.

> **Honest framing.** Corda is a published GR theorist (the paper appeared in a
> peer-reviewed Elsevier journal). The criticisms below are not about
> computational errors in the algebra — most of Corda's algebraic
> manipulations are correct as stated. The critique is that **the algebra
> proves a different thing than Corda claims**. The numerical coincidence at
> Mercury (`πm/M ≈ 43″/century`) is exactly that — a coincidence — and is
> explicitly self-refuted by §5 of the paper itself.

---

## §1. Paper structure (8 pages, 77 equations)

| Section | Content | Headline equation |
|---|---|---|
| §1 Introduction | GR perihelion formula, Mercury background | Eq. (1): `Δφ_GR = 24π³a²/(T₀²c²(1−e²)) = 3π r_g/[a(1−e²)]` |
| §2 Circular orbit | "Newtonian back-reaction" via `F = G(M+m)m/r²` | Eq. (20): `Δφ ≈ πm/M` |
| §3 Harmonic oscillator | Apsidal-angle formalism, same result | Eq. (31): `φ₀ = 2π[3 + F'r/F]^(−1/2)` |
| §4 Kepler's third law | Reduced-mass treatment, same result | Eq. (57): `T = T₀(1+m/M)^(−1/2)` |
| §5 Breakdown | Formula gives `~30×` too large for Venus, `~50×` for Earth | (numerical refutation) |
| §6 Gravitational TD | Schwarzschild weak-field correction | Eq. (67): `Δφ_F = (5π/2) r_g/r₀` |
| §7 Rotational TD | Langevin metric correction | Eq. (72): `τ ≈ (1 − ½r²ω²/c²) t` |
| §8 Combined | Sum of gravitational + rotational TD | Eq. (77): `Δφ_F = 3π r_g/r₀` |
| §9 Conclusion | "Newtonian + corrections matches GR" | — |

## §2. Verification of Corda's algebra (symbolic, Wolfram-MCP)

**The algebra is, for the most part, correct.** This section confirms each
manipulation Corda makes is mathematically valid — the trouble is in §3 and
§4 (interpretation, next section).

### 2.1 Eq. (17) — Period under reduced-mass force

Corda's Eq. (8): `F = G(M+m)m/r²` (force as seen by a Sun-centred observer
treating Sun as fixed but using `M+m` in the gravitational constant). Equating
to centripetal gives `v = √(G(M+m)/r₀)` (Eq. 14) and
`T = T₀(1+m/M)^(−1/2)` (Eq. 17).

```mathematica
ClearAll[Msun, mMerc, r0, v, T, T0, GG];
v = Sqrt[GG (Msun + mMerc)/r0];
T = 2 Pi r0 / v;
T0 = 2 Pi r0^(3/2) / Sqrt[GG Msun];
FullSimplify[T - T0 / Sqrt[1 + mMerc/Msun]]
(* Result: 0  ✅ *)
```

✅ Algebraically correct.

### 2.2 Eq. (20) — Taylor expansion to `πm/M`

`Δφ = 2π[(1+m/M)^(1/2) − 1] ≈ π m/M − (π/4)(m/M)² + …`

```mathematica
Series[2 Pi Sqrt[1 + x] - 2 Pi, {x, 0, 2}]
(* Result: Pi x - Pi x^2/4 + O[x]^3  ✅ leading term is π m/M *)
```

✅ Algebraically correct.

### 2.3 Eqs. (21)–(31) — Apsidal-angle formalism

Corda derives the apsidal angle from the radial equation
`F_c0(r) = m(r̈ − ω₀²r)` (Eq. 21), uses angular-momentum conservation
`J₀ = mr²ω₀` (Eq. 22), Taylor-expands around `r = r₀` (Eqs. 25–27), and lands
at the classical **Bertrand apsidal-angle formula** (Eq. 31):
```math
φ₀ = 2π [3 + F'(r₀) r₀ / F(r₀)]^(−1/2)
```
**For an inverse-square force `F = −k/r²`:** `F'(r)·r/F(r) = −2`, so
`φ₀ = 2π/√(3−2) = 2π` — a closed orbit, **no precession**. This is
[Bertrand's theorem](https://en.wikipedia.org/wiki/Bertrand%27s_theorem) at
work.

```mathematica
ClearAll[Fc, r, r0, kk];
Fc[r_] := -kk/r^2;
ratio = (Fc'[r] r/Fc[r]) /. r -> r0;
phi0 = 2 Pi / Sqrt[3 + ratio];
{Simplify[ratio], Simplify[phi0]}
(* Result: {-2, 2 Pi}  ✅ apsidal angle = 2π, NO precession *)
```

✅ Algebraically correct, and it says **the orbit does not precess**.

### 2.4 Eq. (37) — `(1+m/M)` cancels in F'/F

Corda's "more precise" §3 treatment promotes `F_c0 → F_c = (1+m/M) F_c0`
(Eq. 32). The ratio `F'_c/F_c = F'_c0/F_c0` because the prefactor cancels.

```mathematica
ClearAll[Fc, Fc0, r, r0, mMerc, Msun];
Fc[r_] := (1 + mMerc/Msun) Fc0[r];
FullSimplify[Fc'[r0]/Fc[r0] - Fc0'[r0]/Fc0[r0]]
(* Result: 0  ✅ — the (1+m/M) factor scales magnitude only, not the apsidal
   structure *)
```

✅ Algebraically correct — **and this is exactly why it predicts no precession**.

### 2.5 Eq. (52)–(57) — Reduced-mass Kepler

Standard centre-of-mass + reduced-mass two-body decomposition. The relations
`J = √[μak(1−e²)]`, `T = 2π√(a³μ/k) = 2π√(a³/(GM_T))`,
`a³/T² = G(M+m)/(4π²)` are all textbook ([MTW](https://archive.org/details/gravitation00misn) §25;
Goldstein §3.7).

✅ Algebraically correct.

### 2.6 Eqs. (60)–(67) — Gravitational time dilation

Corda's `t_g ≈ (1 − ½ r_g/r₀) t_l` (Eq. 62) is the standard Schwarzschild
weak-field proper-time relation. He then substitutes the dilated `r₀` and `T`
into Eq. (17) to get `T_F` (Eq. 64), Taylor-expands `ω_F = 2π/T_F`, and lands
at `Δφ_F ≈ (5π/2) r_g/r₀` (Eq. 67).

```mathematica
Series[2 Pi (1 - x/2)^(-3/2) (1 + x/2), {x, 0, 1}]
(* Result: 2 Pi + (5 Pi/2) x + O[x]^2  ✅ leading is (5/2)π r_g/r_0 *)
```

✅ Algebraically correct given the substitution in Eq. (63).

### 2.7 Eq. (75)–(77) — Combined gravitational + rotational TD

Corda's combined factor `(1−x/2)^(−3/2)(1+x/2)(1−x/2)^(−1/2)` should give
`1 + (3/2)x + O(x²)` (and Corda's final Eq. (77) reads `Δφ_F = 3π r_g/r₀`).

> **Possible OCR/typo artifact in Eq. (75).** The marker-pdf rendering shows
> the last factor as `(1−r_g/(2r₀))^(+1/2)`, which would give a leading
> coefficient of `1` (not `3/2`). Reading instead as `(1−r_g/(2r₀))^(−1/2)`
> (consistent with `ω = 2π/T` and `T_T = T_F · (1 − r_g/(2r₀))^(1/2)`) recovers
> the `3/2`. Either the paper has a typo or the converted markdown does.
> Eq. (77)'s stated final `3π r_g/r₀` is the version that's internally
> consistent with the rotational TD derivation in §7.

```mathematica
Series[(1 - x/2)^(-3/2) (1 + x/2) (1 - x/2)^(-1/2), {x, 0, 1}]
(* Result: 1 + 3x/2 + O[x]^2  ✅ leading is (3/2) r_g/r_0 *)
```

Assuming the `−1/2` reading, ✅ Eq. (77) is algebraically correct.

---

## §3. Where the physics goes wrong — the load-bearing critique

The algebra is mostly correct; the **interpretation** is not. Specifically,
three independent problems with the §§2–4 result `Δφ = πm/M`:

### 3.1 Bertrand's theorem says no precession

The apsidal-angle formula Corda derives (Eq. 31) applied to **any**
inverse-square force gives `φ₀ = 2π` exactly. This is Bertrand's theorem: the
only central forces with universally closed bound orbits are `1/r²` and `r`
(harmonic). Promoting `F_c0 → (1+m/M) F_c0` rescales the magnitude of an
inverse-square force but leaves it inverse-square; the apsidal angle stays
`2π`, the orbit stays closed, **no precession**.

Corda's own equation says so (verified in §2.3 above). His subsequent claim
that `Δφ = πm/M` is a precession comes from a **different path** — not from the
apsidal angle in (31) but from the period mismatch in (39) and the angular
velocity in (40). Those are the same `T = T₀(1+m/M)^(−1/2)` and
`ω = ω₀(1+m/M)^(1/2)` from §2. **Period and angular velocity are not the same
thing as apsidal precession.**

### 3.2 Conflating "period change" with "perihelion advance"

Corda's argument structure is:
1. Without back-reaction: `T_0`, `ω_0`, Mercury sweeps `2π` per period.
2. With back-reaction: `T < T_0`, `ω > ω_0`.
3. "Therefore" Mercury sweeps `ω · T_0 > 2π` per (the original) period.
4. The excess `ω T_0 − 2π ≈ πm/M` is the "perihelion advance."

Step 4 is the mistake. A perihelion **advance** is the angular rotation of
the orbit's major axis between successive perihelia — i.e.~the difference
between the angle swept between two perihelia and `2π`. The period over
which this is measured is the **anomalistic period** (perihelion to perihelion),
*not* the unperturbed period `T_0` of some other calculation.

When the orbit is a closed ellipse (inverse-square law, Bertrand), the
anomalistic period equals the sidereal period, and the angle swept between
consecutive perihelia is exactly `2π` — **no matter what `T` and `ω` are**.
Changing the orbital period because the force magnitude changes (here, by
`(1+m/M)`) just means the orbit goes faster or slower; it does *not* make
the closed ellipse start precessing.

The relevant question is: **does the apsidal angle change?** For an
inverse-square force, the answer is **no** (Bertrand). Corda's Eq. (31)
confirms this.

### 3.3 The §5 self-refutation

The strongest evidence that the `πm/M` formula is not a real precession is
that **the paper itself shows it fails for every planet other than Mercury**:

```mathematica
(* Apply Corda's πm/M formula to all three inner planets, with NASA fact-sheet
   masses, against the modern observed perihelion-advance residuals.        *)
ClearAll[Msun, mMerc, mVenus, mEarth, TMerc, TVenus, TEarth,
         orbitsPerCentury, radToArcsec];
Msun = 1.98892*^30; mMerc = 3.3011*^23; mVenus = 4.87*^24; mEarth = 5.97*^24;
TMerc = 87.969*86400.; TVenus = 224.701*86400.; TEarth = 365.256*86400.;
radToArcsec = 180/Pi*3600.;
Table[{planet, Pi mass/Msun * (100*365.25*86400./T) * radToArcsec, observed},
      {{planet, mass, T, observed},
       {{"Mercury", mMerc, TMerc, "≈ 43″/c"},
        {"Venus",   mVenus, TVenus, "≈ 8.6″/c"},
        {"Earth",   mEarth, TEarth, "≈ 3.8″/c"}}}]
```

| Planet | `πm/M` (Corda formula) | Observed residual | Ratio |
|---|---|---|---|
| Mercury | **44.66″/century** | ≈ 43″/c | **1.04×** (claimed success) |
| Venus | **257.91″/century** | ≈ 8.6″/c | **30× too big** |
| Earth | **194.50″/century** | ≈ 3.8″/c | **51× too big** |

The formula's mass-dependence is `Δφ ∝ m · (orbits/c)`. For Mercury this
gives a number that happens to match GR; for Venus and Earth it gives wildly
wrong numbers. **A formula that's right by a factor of 1.04 for one planet
and wrong by factors of 30–50 for the others is not physics; it is a
coincidence.**

Corda's response in §5 is to invoke n-body effects ("Venus has Mercury as an
interior planet"). But the `πm/M` formula does not contain n-body content;
it is a strictly two-body result. The post-hoc rationalisation does not
explain why a strictly two-body formula should be valid for one two-body
system (Sun–Mercury) and invalid for others (Sun–Venus, Sun–Earth). The
honest reading: the formula is not a perihelion advance.

### 3.4 What the `πm/M` result actually is

It is the **fractional change in orbital angular frequency** when one moves
from the "Sun fixed at origin" approximation to the "centre-of-mass at
origin" two-body treatment:

```
ω/ω₀ = √(1 + m/M) ≈ 1 + m/(2M).
```

This is a real effect on the **sidereal period** (Mercury takes very slightly
less time to go around because the gravitational constant in the effective
one-body problem is `G(M+m)`, not `GM`). It does **not** correspond to a
perihelion advance because the orbit is still a closed ellipse —
Bertrand's theorem.

---

## §4. §§6–8 — Gravitational + rotational time dilation

This is a separate calculation from §§2–4 and has different status.

### 4.1 What §§6–8 actually compute

Corda applies first-order GR corrections (gravitational time dilation,
isotropic-coordinate radial-coordinate shift, Langevin rotational time
dilation) to the orbital period `T₀` and angular velocity `ω₀`. The result
`Δφ_F = 3π r_g/r₀` (Eq. 77) **reproduces the e=0 (circular-orbit) limit of
the GR perihelion advance**.

```mathematica
ClearAll[GG, Msun, aMerc, eMerc, cc, rg, orbitsPerCentury, radToArcsec];
GG = 6.6743*^-11; Msun = 1.98892*^30; aMerc = 5.7909*^10;
eMerc = 0.20563; cc = 299792458.;
rg = 2 GG Msun / cc^2;
orbitsPerCentury = 100 * 365.25 * 86400. / (87.969 * 86400.);
radToArcsec = 180/Pi * 3600.;
{"Corda §8 (3π r_g/a, circular)",   3 Pi rg/aMerc * orbitsPerCentury * radToArcsec,
 "GR (3π r_g/[a(1−e²)], Mercury)",  3 Pi rg/(aMerc (1 - eMerc^2)) * orbitsPerCentury * radToArcsec}
(* Result:                                                                      *)
(*   Corda §8 (3π r_g/a, circular):           41.17″/century                    *)
(*   GR (3π r_g/[a(1−e²)], Mercury e=0.21):   42.99″/century                    *)
```

For Mercury (`e = 0.21`, `1−e² ≈ 0.958`), GR predicts
`3π r_g/[a(1−e²)] = 42.99″/century`, **the famous Einstein 1915 result**.
Corda's `3π r_g/a` gives `41.17″/century` — short by the eccentricity factor
`1/(1−e²) ≈ 1.044` (≈ 4.4 %, or `−1.82″/century`).

### 4.2 Is this an alternative derivation of GR's result?

**Mostly yes, in the circular limit.** The combination of gravitational
proper-time dilation + kinematic (rotational) time dilation, computed to
first order in `r_g/r`, is a well-established way to recover the GR
perihelion advance in the weak-field, slow-motion limit. Corda's calculation
is a particular route through this — substituting the dilated `r` and `T`
into the Newtonian Kepler relation and reading off the angular-frequency
shift. The same answer falls out of more standard treatments (e.g.~PPN
analysis, Will *Theory and Experiment in Gravitational Physics* §6).

What Corda's derivation is **missing** is the eccentricity dependence. The
full GR result `3π r_g/[a(1−e²)]` carries the `1/(1−e²)` factor because the
perihelion advance per orbit depends on how close the orbit comes to the
massive body — at perihelion, `r = a(1−e)`, and the relativistic effect is
strongest there. Corda's circular-orbit treatment evaluates the correction at
`r = r₀ = a`, missing this.

For Mercury this is a 4.4 % discrepancy. For more circular orbits (Venus,
Earth) the discrepancy shrinks. For Mercury, the discrepancy *happens to
roughly cancel against* the `πm/M` "reduced-mass advance" Corda adds from
§§2–4 (`44.66 + 0 = 44.66` vs `41.17` from §8 alone — and somehow these are
both presented as agreeing with the observed `43″`). The agreement is not
robust under variation of parameters.

### 4.3 The structural worry

Corda's framework is internally inconsistent: §§2–4 give `πm/M` (which §3.1
above shows is **not** a real perihelion advance) and §§6–8 give `3π r_g/r₀`
(which **is** a real first-order GR correction, but only in the e=0 limit).
The paper presents them as additive corrections that together give the GR
result. The honest read: only the §§6–8 correction is a precession; the
§§2–4 `πm/M` is a period change that is not a precession; and even the §§6–8
result misses the `1/(1−e²)` eccentricity factor.

---

## §5. Implication for Tepper's *Relativistic Newtonian Theory*

Gill's *Relativistic Newtonian Theory* paper (analysed in PR #83) extends
Corda's framework by adding the dual-theory `(1 + V/(mc²))` factor to the
force law:
```
a_m = −(GM/r²)(1 − GM/(c²r)) · ê_r       (Eq. h4 in Gill's paper)
```
and computes the perihelion advance using Corda's `πm/M + relativistic correction`
template. The headline result in PR #83's
[`03_numerical_predictions.md`](../Mercury_Perihelion/03_numerical_predictions.md):

| Prediction | `Δφ` per century |
|---|---|
| Corda `πm/M` (reduced-mass; **not a real precession**) | `+44.66″` |
| Gill framework relativistic correction `−πGM/(c²r₀)` (**opposite sign from GR**) | `−6.87″` |
| Gill net `Δφ_d` (circular-orbit heuristic) | `+37.79″` |
| Standard GR `6πGM/[c²a(1−e²)]` | `+42.99″` |
| Observed (Newcomb–Clemence residual) | `≈ 43″` |

The chain of inheritance:
1. Gill's "perihelion advance from back-reaction" = Corda's `πm/M` →
   **not a perihelion advance** (§3.1, §3.2 above).
2. Gill's relativistic correction has the **opposite sign** of GR's (and
   PR #83 doc 05's proper apsidal-angle calc shows it gives **`−1/6` of GR
   exactly** — `−7.17″/century`).
3. Gill's framework prediction `+37.79″` is therefore the sum of a fictitious
   precession (`+44.66″`) and a wrong-sign relativistic correction
   (`−6.87″`). The numerical near-miss to observed `43″` is the same kind of
   coincidence as Corda's near-miss for Mercury.

**The load-bearing recommendation for Tepper:** the `πm/M` baseline is not a
perihelion advance. PR #83's
[`05_nbody_orbital_mechanics.md`](../Mercury_Perihelion/05_nbody_orbital_mechanics.md)
computes the framework's genuine perihelion advance from the apsidal angle
(the formula Corda himself derives in Eq. 31, applied to Gill's modified
force law): **`−7.17″/century` — opposite sign and one-sixth the magnitude
of GR**. This is the structural mismatch the framework needs to address;
the `+37.79″` "agreement" from Corda's framing is a coincidence on top of a
sign error.

---

## §6. Summary verdict

| Claim | Algebra | Physics | Verdict |
|---|---|---|---|
| Eq. (17) `T = T₀(1+m/M)^(−1/2)` | ✅ correct | ✅ real period change | ✅ |
| Eq. (20) `Δφ = πm/M` (as a **period change**) | ✅ correct | ✅ real period change | ✅ |
| Eq. (20) `Δφ = πm/M` (as a **perihelion advance**) | ✅ correct | 🔴 conflates period with precession | 🔴 |
| Eq. (31) `φ₀ = 2π[3 + F'r/F]^(−1/2)` | ✅ correct | ✅ for inverse-square, gives 2π (no precession) | ✅ (and contradicts §2 reading) |
| Eq. (32) `F_c = (1+m/M) F_c0` (as a "perturbation") | ✅ algebraically | 🔴 inverse-square stays inverse-square; no precession | 🔴 |
| Eqs. (60)–(67) gravitational TD `Δφ_F = (5π/2) r_g/r₀` | ✅ correct | ✅ first-order Schwarzschild | ✅ |
| Eq. (77) combined TD `Δφ_F = 3π r_g/r₀` | ✅ correct (modulo §2.7 typo) | ⚠ matches **e=0** GR; misses `1/(1−e²)` for Mercury | ⚠ |
| §9 "Newtonian matches GR" | — | 🔴 §5 self-refutation; conflation of §§2–4 with §§6–8 | 🔴 |

**Bottom line for Tepper's derivation.** The `πm/M` baseline that
*Relativistic Newtonian Theory* inherits from Corda is not a perihelion
advance. The genuine perihelion advance of Gill's modified force law,
computed via the apsidal-angle formula Corda himself derives (Eq. 31),
is `−1/6` of GR — opposite sign, one-sixth magnitude — independent of the
`πm/M` Corda-style framing. The numerical near-miss to `43″/century` from
Corda's framing (and from Gill's `37.79″`) is the same kind of coincidence
as the Mercury-only success of `πm/M`; it does not survive scrutiny via the
genuine perihelion-advance formula or via the same-formula test on Venus
and Earth.

If the framework's intended Mercury prediction is the genuine perihelion
advance of the modified force law `a = −(GM/r²)(1 − GM/(c²r)) · ê_r`, that
prediction is **`−7.17″/century`**, not `+44.39″` or `+37.79″`. The sign and
magnitude mismatch with observed `+43″/c` is a structural feature of any
force-law-only modification that lacks GR's spatial-metric curvature — the
classic "1/6 factor" of modified Newtonian gravity without a curved space
metric.

---

## §7. Cross-references

- **Source PDF**:
  [`Roadmapping/External_Papers/Secret_Corda.pdf`](../External_Papers/Secret_Corda.pdf)
  (Corda, *Phys. Dark Univ.* 32 (2021) 100834).
- **Converted markdown**:
  [`Roadmapping/Converted_Markdown/Secret_Corda/Secret_Corda.md`](../Converted_Markdown/Secret_Corda/Secret_Corda.md).
- **Companion Gill paper analysis (PR #83)**:
  [`Equation_Verification/Dual_Newtonian_Theory.md`](Dual_Newtonian_Theory.md) and
  [`Mercury_Perihelion/`](../Mercury_Perihelion/).
- **Genuine perihelion-advance calculation for Gill's force law**:
  [`Mercury_Perihelion/05_nbody_orbital_mechanics.md`](../Mercury_Perihelion/05_nbody_orbital_mechanics.md)
  (apsidal-angle treatment giving `−1/6` of GR exactly).
- **GPS author report cross-reference**:
  [`Author_Reports/2026-05_gps_relativity_summary_for_gill.md`](../Author_Reports/2026-05_gps_relativity_summary_for_gill.md)
  §5 Q1 (curved-spacetime extension of the framework).
- **Issue #81**: original ticket.

## §8. Crocco compliance (Substantive AI)

This analysis is substantive AI. The algebra checks are pragmatic (Wolfram MCP
runs); the physical-intuition critiques (Bertrand's theorem applicability,
"period vs precession" framing, identification of `πm/M` as a coincidence
rather than a derivation) are substantive interpretive moves. The
`<!-- TODO: human reviews and fills in -->` block below is left for the human
reviewer.

<!-- TODO: human reviews and fills in — confirms that (a) Bertrand's theorem
correctly rules out perihelion advance from a magnitude rescaling of an
inverse-square force, (b) the §3.1 reading of Corda's Eq. (31) as
"Corda's own apparatus contradicts his §2 conclusion" is the right framing
rather than a pedantic reading, (c) the §3.2 distinction between "period
change" and "perihelion advance" is the load-bearing critique for Tepper to
adjudicate, (d) the §5 implication that Gill's *Relativistic Newtonian
Theory* inherits this issue is the right conclusion to communicate, and
(e) the recommended next step — adopting the apsidal-angle calculation
(PR #83 doc 05) as the framework's genuine Mercury prediction in place of
the Corda-style `πm/M + relativistic correction` heuristic — is the right
direction. -->
