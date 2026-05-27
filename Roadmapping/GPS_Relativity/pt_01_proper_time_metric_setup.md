# Effect pt_01 — proper-time companion to Schwarzschild metric setup

**One-line summary:** The Gill–Zachary substitution rules give an SR-equivalent flat-space master equation `dτ/dt = c/b`, exactly equal to `1/γ`. Extending to Schwarzschild requires a *speculative* substitution `c → b` in the metric, which is not in any verified Gill paper. The minimal extension gives `dτ/dt ≈ 1 − GM/(rb²) − u²/(2b²)`, reducing to the standard `1 − GM/(rc²) − v²/(2c²)` at `b → c` for GPS-relevant velocities.

**Status:** ⚠ partial — kinematic piece is direct framework application; gravitational piece is a speculative minimal extension.
**Framework applicability:** mixed (SR ✅; GR ❌ — extension flagged speculative).
**Standard-method companion:** [`01_schwarzschild_metric_setup.md`](01_schwarzschild_metric_setup.md).

## 0. Framework applicability

The Gill–Zachary substitution rules ([`_proper_time_GPS_cheatsheet.md`](_proper_time_GPS_cheatsheet.md) §1) apply to **Maxwell's equations** and **flat-space relativistic dynamics**. They have **not** been extended to curved spacetime in any verified paper. For the Schwarzschild geometry that the standard derivation uses, the framework as-published does not directly speak — the master equation `dτ/dt ≈ 1 − GM/(rc²) − v²/(2c²)` is a GR result, not an SR + flat-space-EM result.

To produce a proper-time companion to effect 01, this document adopts a *speculative minimal extension*: in the Schwarzschild line element, replace `c` by the local `b = √(c² + u²)` of the test particle whose proper time is being computed. At `u² → 0` this reduces to standard Schwarzschild. The extension is a hypothesis about how the framework might extend to GR, not a derived consequence. It is included so that the side-by-side comparison with the standard derivation is complete.

<!-- TODO: human reviews and fills in — confirms the speculative `c → b` substitution in the Schwarzschild metric is the right minimal-extension hypothesis to record. Alternative extensions (e.g., redefining the metric to involve `b` only in the time-time component, or applying `c → b` only to the matter sector) are not explored. -->

## 1. Effect statement

The proper-time formulation uses the particle's clock `τ` as the natural parameter and `u = dx/dτ` as the natural velocity. For a satellite in orbit, the satellite's atomic clock provides `τ` directly. The question "how does the satellite's `τ` relate to a ground-clock's `τ_grnd`" is the operational question; the answer is to first transform both clocks to a common coordinate time `t` (TCG) and then compose.

The flat-space SR piece is direct: `dt/dτ = b/c`, so `dτ/dt = c/b`. With `b² = c² + u²` and `u = γv`, this equals `c/√(c² + γ²v²) = 1/√(1 + γ²v²/c²) = 1/γ` exactly. So the framework gives the same kinematic time dilation as standard SR, by construction.

The gravitational piece requires the speculative extension. Under `c → b` in the Schwarzschild metric, the time-time component `g_{00} = −(1 − 2GM/rc²)` becomes `g_{00} = −(1 − 2GM/rb²)`. The proper-time line element evaluated on the orbital trajectory then gives the master expression boxed in §3.

## 2. Setup

| Symbol | Value at GPS | Reduction to standard |
|---|---|---|
| `v` (satellite, observer frame) | `3873.96 m/s` | (input) |
| `γ = 1/√(1 − v²/c²)` | `1 + 8.349 × 10⁻¹¹` | (SR identity) |
| `u = γv` (satellite proper velocity) | `3873.96 m/s` (`u − v ≈ 3.2 × 10⁻⁷` m/s) | `u → v` to leading order |
| `b = √(c² + u²)` | `c (1 + 8.349 × 10⁻¹¹)` | `b → c` as `u² → 0` |
| `c/b` | `1 − 8.349 × 10⁻¹¹` | exactly `1/γ` |
| `GM/(rb²)` | `1.670 × 10⁻¹⁰ × (1 − 1.67 × 10⁻¹⁰)` | `GM/(rb²) → GM/(rc²)` |

The relative correction `(b² − c²)/c² = u²/c² ≈ 1.67 × 10⁻¹⁰` is of the same order as the gravitational redshift itself.

## 3. Derivation

**Flat-space SR piece (direct framework application).** From rule (Velocity duality, [cheatsheet §1](_proper_time_GPS_cheatsheet.md)): `w/c = u/b`. The observer-frame proper time satisfies `dt/dτ = γ`. Identifying `γ = b/c` (which holds because `1/√(1 − v²/c²) = 1/√(1 − u²/b²) = b/c` after the velocity-duality substitution):

$$\boxed{\;\frac{d\tau}{dt}\bigg|_{\rm SR} \;=\; \frac{c}{b} \;=\; \frac{c}{\sqrt{c^2 + u^2}} \;=\; \frac{1}{\gamma}.\;}$$

**Schwarzschild piece (speculative minimal extension).** Under `c → b` in the Schwarzschild metric:

$$ds^2 \;=\; -\!\left(1 - \frac{2GM}{r b^2}\right) b^2 dt^2 + \left(1 - \frac{2GM}{r b^2}\right)^{-1} dr^2 + r^2 \, d\Omega^2.$$

For a test particle following a worldline with `dx/dτ = u`, the proper time satisfies (analogous to the standard PN expansion):

$$\boxed{\;\frac{d\tau}{dt} \;\approx\; 1 \;-\; \frac{GM}{r b^2} \;-\; \frac{u^2}{2 b^2} \;+\; \mathcal{O}(b^{-4}).\;}$$

The two pieces combine into a master equation that *structurally* looks like the standard one but with `c → b` everywhere — and that *numerically* reduces to the standard one at GPS-relevant velocities (`b/c − 1 ~ 10⁻¹¹`, all `1/b²` factors equal `1/c²` to that precision).

## 4. Numerical evaluation

| Term | Standard (with `c`) | PT (with `b`) | Difference |
|---|---|---|---|
| `GM/(rc²)` at geoid | `6.953 × 10⁻¹⁰` | `6.953 × 10⁻¹⁰ × (1 + 𝒪(10⁻¹¹))` | `~ 10⁻²⁰` |
| `GM/(rc²)` at GPS | `1.670 × 10⁻¹⁰` | `1.670 × 10⁻¹⁰ × (1 + 𝒪(10⁻¹¹))` | `~ 10⁻²⁰` |
| `v²/(2c²) = u²/(2b²)` at GPS | `8.349 × 10⁻¹¹` | `8.349 × 10⁻¹¹` exactly | `0` (algebraic identity) |
| `c/b` | `1 − 8.349 × 10⁻¹¹` | `1 − 8.349 × 10⁻¹¹` | `0` (exact) |

The `u²/(2b²) = v²/(2c²)` equality is exact (algebraic identity from the velocity duality `w/c = u/b`). The `GM/(rb²) = GM/(rc²)` equality is only true to leading order — the correction is `~ 10⁻²⁰` (the difference `u²/(c² b²)`), well below operational precision.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; RR=6378137; aGPS=26560000; vGPS=Sqrt[GMe/aGPS]; gam=1/Sqrt[1 - vGPS^2/cc^2]; uGPS=gam vGPS; bGPS=Sqrt[cc^2 + uGPS^2]; Print["c/b = ", N[cc/bGPS, 14]]; Print["1/gamma = ", N[1/gam, 14]]; Print["c/b - 1/gamma = ", N[cc/bGPS - 1/gam, 18]]; Print["u^2/b^2 = ", ScientificForm[N[uGPS^2/bGPS^2, 12]]]; Print["v^2/c^2 = ", ScientificForm[N[vGPS^2/cc^2, 12]]]; Print["GM/(R b^2) at geoid = ", ScientificForm[N[GMe/(RR bGPS^2), 12]]]; Print["GM/(R c^2) at geoid = ", ScientificForm[N[GMe/(RR cc^2), 12]]]; Print["relative difference = ", ScientificForm[N[(GMe/(RR bGPS^2) - GMe/(RR cc^2))/(GMe/(RR cc^2)), 8]]]
```

**Result:**
```
c/b = 0.9999999999166
1/gamma = 0.9999999999166
c/b - 1/gamma = 0
u^2/b^2 = 1.66981424×10⁻¹⁰
v^2/c^2 = 1.66981424×10⁻¹⁰
GM/(R b^2) at geoid = 6.95348928×10⁻¹⁰
GM/(R c^2) at geoid = 6.95349044×10⁻¹⁰
relative difference = -1.66981×10⁻¹⁰
```
✅ confirms `c/b = 1/γ` exactly and `GM/(rb²) = GM/(rc²)` to within `~10⁻¹⁰` (negligible at GPS precision).

## 5b. GPS-precision-limit equivalence to standard

Expand `b² = c²(1 + u²/c²)` and `1/b² = (1/c²)(1 − u²/c² + 𝒪(u⁴/c⁴))`. The proper-time master equation becomes:

$$\frac{d\tau}{dt} \;=\; 1 \;-\; \frac{GM}{r c^2}\left(1 - \frac{u^2}{c^2}\right) \;-\; \frac{u^2}{2 c^2}\left(1 - \frac{u^2}{c^2}\right) \;+\; \mathcal{O}(c^{-4}).$$

The cross-terms `GM u²/(r c⁴)` and `u⁴/(2c⁴)` are at order `1/c⁴` and are dropped (effect 01 truncates at `1/c²`). Using `u = γv ≈ v + 𝒪(v³/c²)` to leading order, the equation reduces to:

$$\frac{d\tau}{dt} \;\approx\; 1 \;-\; \frac{GM}{r c^2} \;-\; \frac{v^2}{2c^2}.$$

This is the standard master equation. The reduction error is `𝒪(GM u²/(rc⁴)) ~ 10⁻¹⁹` for GPS — below the operational precision floor.

## 6. Comparison with standard-method companion

Standard [`01_schwarzschild_metric_setup.md`](01_schwarzschild_metric_setup.md) §3 gives the master equation `dτ/dt ≈ 1 − GM/(rc²) − v²/(2c²)`. The PT formulation gives `dτ/dt ≈ 1 − GM/(rb²) − u²/(2b²)`. After applying `b² = c² + u²` and `u = γv`, the two forms are equivalent to `𝒪(c⁻⁴)`, with structural differences only at order `10⁻¹⁹` — well below GPS precision. The two formulations describe the same physics in different bookkeeping.

The structural difference: in the standard formulation, `c` is fundamental and `γ` is a derived correction; in the proper-time formulation, `b` is the fundamental "speed of light in the moving frame" and `γ = b/c` is derived. The bookkeeping privileges the satellite's local clock as the natural parameter, which is operationally what GPS does.

## 6b. Where the framework would diverge

The PT and standard master equations would give numerically distinguishable predictions:

1. **At ultrarelativistic velocities.** When `u² ~ c²`, `b/c` is `𝒪(1)` and the algebraic difference between `1/b²` and `1/c²` is no longer suppressed. This is the regime of cosmic-ray muons (`γ ~ 10`–`10⁶`), where `b/c ~ 10`–`10⁶`. For GPS this is irrelevant.
2. **In strong gravitational fields.** The speculative extension's `1/b²` factor in the Schwarzschild `g_{00}` would matter when `GM/(rc²)` is `𝒪(1)`, e.g., near black-hole horizons. For Earth's field at GPS altitude (`GM/(rc²) ~ 10⁻¹⁰`) the framework's prediction is indistinguishable from standard.
3. **For non-geodesic motion** (acceleration `a` non-zero in the local frame). For GPS the satellite is in free fall on a geodesic; the local-frame acceleration is zero. For a thrusting spacecraft, the framework's third-term contribution to dynamics could appear; this is the regime of [issue #43](https://github.com/temoTxt/PyPhysics/issues/43) (radiation-reaction experiments), not GPS.

## 7. Verdict

⚠ matches standard at GPS precision (algebraic structure differs; numerical reduction is identical to `𝒪(c⁻⁴)`). The kinematic piece (`c/b = 1/γ`) is exact framework application; the gravitational piece requires a speculative `c → b` Schwarzschild extension that is not in any verified Gill paper.
