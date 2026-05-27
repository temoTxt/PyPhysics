# Effect 03 — Velocity (kinematic) time dilation

**One-line summary:** The GPS satellite moves at `v ≈ 3 874 m/s` in the non-rotating geocentric inertial (ECI) frame, faster than a ground station's rotational speed of `≤ 465 m/s` at the equator. The resulting time dilation makes the satellite clock tick *slow* by `−7.11 μs/day` — the smaller, opposite-sign partner to effect 02.

**Status:** ✅ drafted.
**Ashby (2003) cross-reference:** §3.2, Eq. (20).

## 1. Effect statement

Special-relativistic time dilation acts on any clock moving relative to the inertial frame in which the comparison is made. For GPS, the natural inertial frame is the **Earth-Centered Inertial (ECI) frame** — origin at the Earth's center, axes non-rotating with respect to distant stars. In ECI, the GPS satellite circles at speed `v_sat = √(GM⊕/a_GPS) ≈ 3 874 m/s` (the circular-orbit relation); the ground clock orbits at speed `v_grnd = ω⊕ R⊕ cos λ ≤ 465 m/s`, where the upper bound applies at the equator.

The fractional rate offset is `−(v² − v_eq²)/(2c²)` — the velocity term from the master equation of effect 01, applied differentially between the two clocks. Both have the same algebraic sign (both make their clock tick slow); the satellite's larger speed wins, and the *difference* makes the satellite tick slow relative to the ground.

## 2. Setup

| Symbol | Value | Meaning |
|---|---|---|
| `v_sat` | `3 873.96 m/s` | circular-orbit speed in ECI, `√(GM/a)` |
| `v_grnd` (equator) | `465.10 m/s` | rotational speed at equator, `ω⊕ R⊕` |
| `v_grnd` (latitude λ) | `465.10 cos λ m/s` | rotational speed at latitude λ |

The satellite's orbital speed is fixed by the orbit's semi-major axis through the circular-orbit relation `v² = GM/a` (Kepler's third law in differential form). The ground station's speed depends on latitude. For the headline calculation we use the equator (the largest ground speed) — this is the conservative (most generous to the satellite) choice. The latitude dependence is sub-percent of the total velocity correction and is folded into the operational geoid value used in effect 02's alternative decomposition.

## 3. Derivation

From effect 01's master equation, the velocity term gives:

$$\left.\frac{d\tau}{dt}\right|_{\rm vel} \;\approx\; 1 \;-\; \frac{v^2}{2c^2}.$$

For two clocks, satellite (`v_sat`) and ground (`v_grnd`), the differential rate from velocity alone is:

$$\left.\frac{d\tau_{\rm sat}}{dt}\right|_{\rm vel} - \left.\frac{d\tau_{\rm grnd}}{dt}\right|_{\rm vel} \;=\; -\,\frac{v_{\rm sat}^2 - v_{\rm grnd}^2}{2c^2}.$$

Substituting `v_sat² = GM/a_GPS` (circular orbit) and `v_grnd² = ω⊕² R⊕² cos²λ` (Earth's rotation), evaluated at the equator where `cos λ = 1`:

$$\boxed{\;\left.\frac{\Delta f}{f}\right|_{\rm vel} \;=\; -\,\frac{1}{2c^2}\left(\frac{GM_\oplus}{a_{\rm GPS}} \;-\; \omega_\oplus^2 R_\oplus^2\right).\;}$$

The sign is negative: the satellite ticks slower than the ground clock by this fractional rate.

## 4. Numerical evaluation

Substituting:

$$\frac{v_{\rm sat}^2}{2c^2} = \frac{GM_\oplus}{2\,a_{\rm GPS}\, c^2} = 8.349\,071\,2 \times 10^{-11},$$

$$\frac{v_{\rm grnd}^2}{2c^2} = \frac{\omega_\oplus^2 R_\oplus^2}{2c^2} = 1.203\,408\,8 \times 10^{-12}.$$

Difference:

$$\left.\frac{\Delta f}{f}\right|_{\rm vel} = -(8.349\,071\,2 - 0.120\,340\,9)\times 10^{-11} = -8.228\,730\,3 \times 10^{-11}.$$

Per day (`86 400 s`):

$$\left.\Delta t\right|_{\rm vel,\,daily} = -8.2287 \times 10^{-11} \times 86\,400 \;\approx\; -7.110 \times 10^{-6}\ {\rm s} \;=\; \boxed{\;-7.11\ \mu{\rm s/day}\;}.$$

The satellite loses `7.11 μs` per day relative to the ground from velocity dilation, before any gravitational contribution.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; RR=6378137; aGPS=26560000; om=7.2921151467*^-5; day=86400; velRate = -(GMe/aGPS - (om RR)^2)/(2 cc^2); Print["vel rate = ", ScientificForm[N[velRate,8]]]; Print["vel microsec/day = ", N[velRate*day*1*^6,6]]
```

**Result:**
```
vel rate = -8.2287303×10⁻¹¹
vel microsec/day = -7.10962
```
✅ matches §4.

## 6. Comparison with Ashby (2003)

Ashby §3.2 derives the velocity time dilation in Eq. (20). His quoted value is `−7.2 μs/day`, which agrees with ours to within rounding: Ashby uses a slightly different `v_sat ≈ 3.874 km/s` (rounded), giving `v_sat²/(2c²) ≈ 8.34 × 10⁻¹¹`; the `0.1 μs/day` apparent disagreement is the difference between including vs neglecting the ground-station rotation contribution. Our `−7.11 μs/day` includes the equator-station rotation as a credit (it reduces the magnitude of the satellite's deficit); Ashby's `−7.2 μs/day` is the satellite-only contribution.

A common (but incorrect) decomposition treats `v_grnd = 0` (i.e., the ground clock is treated as if at rest in ECI). That gives `−v_sat²/(2c²) × 86 400 = −7.22 μs/day`. The right operational answer includes the ground-station rotation, which contributes `+0.10 μs/day`. The net velocity contribution is `−7.11 μs/day`.

## 7. Verdict

✅ reproduces Ashby (2003) §3.2 Eq. (20) at quoted precision when ground-station rotation is included consistently with how the geoid potential absorbs centrifugal terms in the operational decomposition. The bookkeeping is identical to effect 02's reconciliation: whether you put the centrifugal term in "gravity" (operational geoid) or in "velocity" (this document), the *sum* in effect 04 is the same `+38.5 μs/day`.
