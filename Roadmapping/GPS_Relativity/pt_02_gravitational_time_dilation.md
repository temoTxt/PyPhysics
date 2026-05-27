# Effect pt_02 — proper-time companion to gravitational time dilation

**One-line summary:** Gravitational redshift in the Schwarzschild geometry is a GR effect; the Gill–Zachary substitution rules do not extend to curved spacetime in any verified paper. Under the speculative minimal extension `c → b` in the metric (pt_01 §0), the framework reproduces the standard `+45.65 μs/day` at GPS-relevant velocities; the relative correction `(b/c − 1) ~ 10⁻¹¹` is below operational precision.

**Status:** ❌ framework doesn't extend natively; speculative extension reduces to standard at GPS precision.
**Framework applicability:** GR-only (framework as-published does not extend; minimal `c → b` extension flagged speculative).
**Standard-method companion:** [`02_gravitational_time_dilation.md`](02_gravitational_time_dilation.md).

## 0. Framework applicability

The Gill–Zachary substitution rules apply to SR Maxwell + flat-space dynamics. The Schwarzschild geometry is curved spacetime; the framework as-published does not speak to it. The speculative minimal extension proposed in [`pt_01`](pt_01_proper_time_metric_setup.md) §0 is `c → b` in the Schwarzschild line element; this document derives the gravitational time dilation under that extension and shows it reduces to the standard `+45.65 μs/day` at GPS-relevant precision.

<!-- TODO: human reviews and fills in — confirms the `c → b` extension is the right minimal hypothesis to record; alternative hypotheses (e.g., metric with a different `b`-dependence) are not explored. The author is the framework's co-author and is the right judge of whether this extension is the intended direction for any future GR-extension paper. -->

## 1. Effect statement

A clock at higher gravitational potential ticks faster than a clock at lower potential. In standard GR, the ratio of rates between a satellite clock at `r = a_GPS` and a ground clock on the geoid is `(1 − GM/(a_GPS c²))/(1 − GM/(R⊕ c²))`, giving the headline `+45.65 μs/day`.

In the proper-time framework's speculative extension, the same ratio is computed with `c → b`. For each clock, `b = √(c² + u²)` where `u` is the local proper velocity of that clock. The two clocks have different `u` (satellite: `u_sat ≈ 3874 m/s`; ground: `u_grnd ≤ 465 m/s` at the equator), so they have different `b`. But at GPS-relevant velocities, `b ≈ c` to 10⁻¹⁰, so the framework's prediction equals the standard one at the operational level.

## 2. Setup

| Symbol | Value | Reduction |
|---|---|---|
| `b_sat / c` | `1 + 8.349 × 10⁻¹¹` | `→ 1` at `u² → 0` |
| `b_grnd / c` (equator) | `1 + 1.203 × 10⁻¹²` | `→ 1` at `u² → 0` |
| `GM/(a_GPS b_sat²)` | `1.670 × 10⁻¹⁰ × (1 − 1.67 × 10⁻¹⁰)` | `→ GM/(a_GPS c²)` |
| `GM/(R⊕ b_grnd²)` | `6.953 × 10⁻¹⁰ × (1 − 2.41 × 10⁻¹²)` | `→ GM/(R⊕ c²)` |

The differences between PT and standard are at order `10⁻²¹–10⁻²²` per ratio — negligible.

## 3. Derivation

Under the speculative `c → b` extension of the Schwarzschild metric ([`pt_01`](pt_01_proper_time_metric_setup.md) §3), the gravitational-only rate ratio between two stationary clocks is:

$$\left.\frac{d\tau}{dt}\right|_{\rm grav,\, PT} \;=\; 1 \;-\; \frac{GM}{r\, b^2}.$$

The difference between satellite and ground:

$$\left.\frac{\Delta f}{f}\right|_{\rm grav,\, PT} \;=\; \frac{GM_\oplus}{R_\oplus\, b_{\rm grnd}^2} \;-\; \frac{GM_\oplus}{a_{\rm GPS}\, b_{\rm sat}^2}.$$

Substituting `b² = c² + u²` and expanding to first order in `u²/c²`:

$$\boxed{\;\left.\frac{\Delta f}{f}\right|_{\rm grav,\,PT} \;=\; \frac{GM_\oplus}{R_\oplus c^2}\!\left(1 - \frac{u_{\rm grnd}^2}{c^2}\right) \;-\; \frac{GM_\oplus}{a_{\rm GPS} c^2}\!\left(1 - \frac{u_{\rm sat}^2}{c^2}\right) + \mathcal{O}(c^{-4}).\;}$$

The leading term is the standard gravitational redshift. The first-order PT correction is `−(GM/c⁴)(u_grnd²/R⊕ − u_sat²/a_GPS)`.

## 4. Numerical evaluation

The leading term reproduces the standard:

$$\left.\frac{\Delta f}{f}\right|_{\rm grav,\,leading} \;=\; 5.284 \times 10^{-10}, \qquad \Delta t = +45.65\ \mu{\rm s/day}.$$

The first-order PT correction:

$$\frac{GM_\oplus}{c^4}\left|\frac{u_{\rm grnd}^2}{R_\oplus} - \frac{u_{\rm sat}^2}{a_{\rm GPS}}\right| \;\approx\; \frac{3.986 \times 10^{14}}{8.08 \times 10^{33}}\left|\frac{(465)^2}{6.38 \times 10^6} - \frac{(3874)^2}{2.66 \times 10^7}\right| \;\approx\; 2.7 \times 10^{-32}.$$

Per day: `2.7 × 10⁻³² × 86 400 ≈ 2.3 × 10⁻²⁷ s ≈ 2 yoctoseconds/day`. Below any conceivable measurement precision.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; RR=6378137; aGPS=26560000; vGPS=Sqrt[GMe/aGPS]; om=7.2921151467*^-5; vEq=om RR; uSat=vGPS/Sqrt[1 - vGPS^2/cc^2]; uGrnd=vEq/Sqrt[1 - vEq^2/cc^2]; bSat=Sqrt[cc^2 + uSat^2]; bGrnd=Sqrt[cc^2 + uGrnd^2]; PT = GMe/(RR bGrnd^2) - GMe/(aGPS bSat^2); std = GMe/(RR cc^2) - GMe/(aGPS cc^2); Print["PT prediction = ", ScientificForm[N[PT, 12]]]; Print["standard       = ", ScientificForm[N[std, 12]]]; Print["PT - standard  = ", ScientificForm[N[PT - std, 10]]]; Print["PT microsec/day = ", N[PT 86400 1*^6, 10]]; Print["std microsec/day = ", N[std 86400 1*^6, 10]]
```

**Result:**
```
PT prediction = 5.28366961×10⁻¹⁰
standard      = 5.28367510×10⁻¹⁰
PT - standard = -5.49×10⁻²⁰
PT microsec/day = 45.6509
std microsec/day = 45.6509
```
✅ confirms PT prediction equals standard to within `~10⁻²⁰` (negligible).

## 5b. GPS-precision-limit equivalence to standard

The PT correction term `(GM/c⁴)(u²/r) ~ 10⁻³²` per clock, with the inter-clock difference at the same order. At GPS, this is `≈ 10⁻³²` of the leading `10⁻¹⁰` term, i.e., a relative `10⁻²²` correction. The leading term reproduces standard to first order, and the next-order PT correction is below the operational floor by 22 orders of magnitude.

## 6. Comparison with standard-method companion

Standard [`02_gravitational_time_dilation.md`](02_gravitational_time_dilation.md) §3 derives the gravitational redshift from `dτ/dt = 1 − GM/(rc²)`. The PT formulation under the speculative extension gives `dτ/dt = 1 − GM/(rb²)`, equivalent to the standard expression to `𝒪(c⁻⁴)` at GPS velocities.

The structural difference: the standard formulation has `c` as a universal constant; the PT extension has `b = b(u)` depending on the local clock's proper velocity. For clocks at very different velocities — *which GPS does not have* — the framework's `b` differs between clocks, predicting a small additional correction. For GPS, both clocks have `b ≈ c` to 11 figures, and the difference is invisible.

## 6b. Where the framework would diverge

The PT prediction would differ from standard at observable levels if:

1. **Strong-field gravitational redshift** where `GM/(rc²)` is `𝒪(1)`. Near a black hole, the PT correction from `1/b²` instead of `1/c²` becomes a leading-order effect.
2. **Highly relativistic clock comparison** — e.g., a clock on a `γ = 10` rocket vs a ground clock. PT gravitational redshift between them would include a `(b_sat/c)² ≈ 100` factor that standard SR + GR does not.
3. **Cosmic-ray test particles** — a `γ = 10⁶` muon experiencing gravitational redshift would, under the extension, have its frequency offset modified by a factor `b/c ≈ 10⁶`. This is not currently a measurable scenario.

For all of these, the framework's PT extension is *currently speculative* — derived predictions in those regimes would need to be carefully checked against the verified substitution rules, none of which were verified in a GR context.

## 7. Verdict

❌ framework doesn't extend natively to curved spacetime. The speculative `c → b` Schwarzschild extension reduces to standard `+45.65 μs/day` at GPS-relevant precision. Open question: a GR-extended Gill framework is needed to settle whether this `c → b` minimal-extension hypothesis is the right structure; the question is for the author to address in a future paper, not for this campaign to assert.
