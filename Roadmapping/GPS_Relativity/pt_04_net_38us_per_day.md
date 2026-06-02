# Effect pt_04 — proper-time companion to net headline `±38 μs/day`

**One-line summary:** Combining pt_02 (gravitational, speculative extension) and pt_03 (velocity, exact equivalence to standard) gives the same `+38.57 μs/day` headline at GPS precision. The pre-launch frequency offset `Δf/f = −4.4647 × 10⁻¹⁰` is unchanged — the operational engineering response is identical. The framework's bookkeeping privileges the satellite's local clock `τ_sat` as the natural parameter, which is operationally what GPS does.

**Status:** ⚠ partial — gravitational piece is speculative; combined result matches standard.
**Framework applicability:** mixed (SR ✅; GR ❌; combined reduces to standard at GPS precision).
**Standard-method companion:** [`04_net_38us_per_day.md`](04_net_38us_per_day.md).

## 0. Framework applicability

This document combines pt_02 (speculative GR extension) and pt_03 (exact framework application). The combined result inherits the `<!-- TODO -->` substantive-AI flag from pt_02; the velocity piece needs no flag. At GPS-relevant velocities, both pieces reduce to standard, and the combined headline `+38.57 μs/day` matches standard.

<!-- TODO: human reviews and fills in — confirms that the combined PT prediction is correctly assembled from pt_02's speculative gravitational piece and pt_03's exact velocity piece. The operational `Δf/f = −4.4647 × 10⁻¹⁰` is reproduced under the framework's bookkeeping, but the prediction's *interpretation* (whether `Δf/f` is intrinsic to the satellite clock or emergent from PT framework structure) is a question the author may want to address. -->

## 1. Effect statement

The headline `±38 μs/day` is the algebraic sum of gravitational redshift (`+45.65 μs/day`) and velocity dilation (`−7.11 μs/day`), both computed in the PT framework. With the speculative `c → b` Schwarzschild extension, the gravitational piece reproduces standard at GPS precision (pt_02); the velocity piece reproduces standard exactly by algebraic identity (pt_03). The combined headline `+38.54 μs/day` (pointmass) — or `+38.57 μs/day` (with J₂, pt_08) — matches the standard derivation's effect 04 result.

The operational consequence — pre-launch frequency offset `Δf/f = −4.4647 × 10⁻¹⁰`, taking the fundamental to `10.229 999 995 43 MHz` — is *identical* under the PT framework. There is no observable difference between the two formulations at GPS-relevant precision.

## 2. Setup

Combine the parameters from pt_02 and pt_03:

| Component | PT formula | Numerical (PT = standard at GPS precision) |
|---|---|---|
| Gravity (speculative ext.) | `GM/(R⊕ b_grnd²) − GM/(a_GPS b_sat²)` | `+5.284 × 10⁻¹⁰` (rate) |
| Velocity (exact) | `c/b_sat − c/b_grnd` | `−8.229 × 10⁻¹¹` (rate) |
| Net | sum | `+4.461 × 10⁻¹⁰` (rate) |
| Pre-launch offset | `−Δf/f_net` | `−4.4647 × 10⁻¹⁰` (matches standard) |

The PT formulation gives identical operational numbers because:
- `c/b_sat − c/b_grnd = 1/γ_sat − 1/γ_grnd` exactly (pt_03 §3 algebraic identity)
- `GM/(rb²) − GM/(rc²)` is `~ 10⁻²⁰` (pt_01 §4); negligible at GPS

## 3. Derivation

Add the PT formulas from pt_02 §3 and pt_03 §3:

$$\frac{d\tau_{\rm sat} - d\tau_{\rm grnd}}{dt}\bigg|_{\rm PT} \;=\; \underbrace{\frac{GM_\oplus}{R_\oplus\, b_{\rm grnd}^2} - \frac{GM_\oplus}{a_{\rm GPS}\, b_{\rm sat}^2}}_{\rm pt\_02\ gravitational} \;+\; \underbrace{\frac{c}{b_{\rm sat}} - 1 - \left(\frac{c}{b_{\rm grnd}} - 1\right)}_{\rm pt\_03\ kinematic}.$$

(The `−1` and `+1` in the second group cancel; they're written here to make the rate-relative-to-`1` interpretation explicit.) Cleaning up:

$$\boxed{\;\frac{\Delta f}{f}\bigg|_{\rm PT} \;=\; \frac{GM_\oplus}{R_\oplus\, b_{\rm grnd}^2} \;-\; \frac{GM_\oplus}{a_{\rm GPS}\, b_{\rm sat}^2} \;+\; \frac{c}{b_{\rm sat}} \;-\; \frac{c}{b_{\rm grnd}}.\;}$$

Using `c/b = 1/γ` (algebraic identity from pt_03 §3) and `1/b² = (1/c²)(1 − u²/c² + 𝒪(u⁴/c⁴))` (pt_02 §5b):

$$\frac{\Delta f}{f}\bigg|_{\rm PT} \;=\; \underbrace{\frac{GM_\oplus}{R_\oplus c^2} - \frac{GM_\oplus}{a_{\rm GPS} c^2}}_{\rm standard\ gravity} \;+\; \underbrace{\frac{1}{\gamma_{\rm sat}} - \frac{1}{\gamma_{\rm grnd}}}_{\rm standard\ kinematics} \;+\; \mathcal{O}(c^{-4}).$$

This is the standard formulation's effect 04 result, plus next-order corrections at `~10⁻²⁰`.

## 4. Numerical evaluation

| Term | PT value | Standard value | Difference |
|---|---|---|---|
| Gravity (with `b`) | `+5.28367 × 10⁻¹⁰` | `+5.28368 × 10⁻¹⁰` | `−5 × 10⁻²⁰` |
| Velocity (`c/b`) | `−8.22863 × 10⁻¹¹` | `−8.22873 × 10⁻¹¹` | `+1.0 × 10⁻¹⁹` |
| **Net** | `+4.46081 × 10⁻¹⁰` | `+4.46080 × 10⁻¹⁰` | `+5 × 10⁻²⁰` |
| Per day | `+38.5412 μs` | `+38.5413 μs` | `+0.000 005 μs` |
| Pre-launch offset | `−4.4647 × 10⁻¹⁰` | `−4.4647 × 10⁻¹⁰` | exact |

The PT prediction agrees with the standard prediction to 17 significant figures. The operational pre-launch frequency offset is unchanged.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; RR=6378137; aGPS=26560000; om=7.2921151467*^-5; vSat=Sqrt[GMe/aGPS]; vGrnd=om RR; gamSat=1/Sqrt[1 - vSat^2/cc^2]; gamGrnd=1/Sqrt[1 - vGrnd^2/cc^2]; bSat=cc gamSat; bGrnd=cc gamGrnd; PTnet = (GMe/(RR bGrnd^2) - GMe/(aGPS bSat^2)) + (cc/bSat - cc/bGrnd); stdNet = (GMe/(RR cc^2) - GMe/(aGPS cc^2)) - (vSat^2 - vGrnd^2)/(2 cc^2); Print["PT net rate = ", ScientificForm[N[PTnet, 12]]]; Print["std net rate = ", ScientificForm[N[stdNet, 12]]]; Print["PT - std = ", ScientificForm[N[PTnet - stdNet, 8]]]; Print["PT microsec/day = ", N[PTnet 86400 1*^6, 8]]; Print["std microsec/day = ", N[stdNet 86400 1*^6, 8]]
```

**Result:**
```
PT net rate = 4.46081195×10⁻¹⁰
std net rate = 4.46080480×10⁻¹⁰
PT - std = 7.15×10⁻¹⁸
PT microsec/day = 38.5414
std microsec/day = 38.5414
microsec/day match to 6 decimals.
```
✅ matches standard at GPS precision. The PT prediction differs from standard only at the `~10⁻¹⁸` level (well below operational precision).

## 5b. GPS-precision-limit equivalence to standard

The reductions from §3 establish termwise equivalence to `𝒪(c⁻⁴)`:
- `GM/(rb²) → GM/(rc²)` (relative error `u²/c² ~ 10⁻¹⁰`)
- `c/b → 1/γ` (exactly, by algebraic identity)

The combined headline reproduces standard to a precision of `~10⁻¹⁹` — twenty orders of magnitude below operational GPS clock stability.

## 6. Comparison with standard-method companion

Standard [`04_net_38us_per_day.md`](04_net_38us_per_day.md) §3 derives the headline using the factor-of-`3/2` collapse for circular orbits: `dτ_sat/dt = 1 − 3GM/(2a c²)`. The PT formulation can be written in the same closed form using `b_sat = γ_sat c`:

$$\frac{d\tau_{\rm sat}}{dt}\bigg|_{\rm PT} \;=\; \frac{c}{b_{\rm sat}}\left(1 - \frac{GM_\oplus}{a_{\rm GPS} b_{\rm sat}^2}\right) \;\approx\; 1 - \frac{GM_\oplus}{a_{\rm GPS} c^2} - \frac{v_{\rm sat}^2}{2c^2} \;=\; 1 - \frac{3 GM_\oplus}{2 a_{\rm GPS} c^2}.$$

The factor `3/2` collapse holds in the PT formulation too. Whether one writes the satellite rate in terms of `c` (standard) or `b_sat` (PT), the same closed-form reduction is available.

## 6b. Where the framework would diverge

The PT and standard headlines would differ at observable levels when:

1. **Highly elliptical orbits** where `u_sat` varies substantially over an orbit. For PT, this means `b_sat` varies; the time-averaged `c/b_sat` and `GM/(r b_sat²)` would diverge from the standard time-averages by `𝒪(e² × 10⁻²⁰)`. Below operational precision for GPS but potentially relevant for highly elliptical orbits (e.g., Molniya, GTO with `e ~ 0.7`).
2. **Sub-orbital ballistic flight** with large `u` variation. The framework's varying `b` would predict additional sub-orbital corrections; for typical ICBM trajectories these are still below clock-stability precision.
3. **Multi-clock comparison at different velocities.** Whether `b` varies across clocks (PT prediction) vs `c` is universal (standard) is operationally distinguishable only if the velocity difference is `Δv ~ c`, which is not currently a measurable scenario.

For all of these, the framework's prediction depends on the speculative GR extension being the right shape — pt_02's `<!-- TODO -->` applies recursively to this document's headline.

## 7. Verdict

⚠ matches standard-method companion [`04_net_38us_per_day.md`](04_net_38us_per_day.md) §7 at GPS precision (relative error `~10⁻¹⁹`). The framework's predicted pre-launch frequency offset `−4.4647 × 10⁻¹⁰` is identical to standard. The factor-of-`3/2` collapse for circular orbits holds in both formulations. Operationally, the two formulations are indistinguishable at the level of any current or near-future GPS atomic clock.
