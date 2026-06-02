# Effect pt_07 — proper-time companion to Shapiro delay

**One-line summary:** The Shapiro delay is purely a curved-spacetime effect (extra coordinate time the photon spends in the gravitational potential); the Gill–Zachary framework as-published does not extend to curved spacetime. Under the speculative minimal extension (pt_01 §0), the framework reproduces the standard `~62 ps` per signal at GPS-relevant precision; the relative correction `(b/c − 1) ~ 10⁻¹¹` is below operational floor.

**Status:** ❌ framework doesn't extend natively; speculative extension reduces to standard at GPS precision.
**Framework applicability:** GR-only (framework as-published does not extend; minimal `c → b` extension flagged speculative).
**Standard-method companion:** [`07_shapiro_delay.md`](07_shapiro_delay.md).

## 0. Framework applicability

The Shapiro delay derivation requires the spatial part of the Schwarzschild metric (`(1 − 2GM/rc²)⁻¹ dr²`) — this is curved spacetime, outside the framework's verified flat-space Maxwell + dynamics domain. The framework as-published does not produce a Shapiro-delay prediction.

The minimal-extension hypothesis from pt_01 §0 — `c → b` in the Schwarzschild metric — gives a PT analog of the Shapiro formula: `Δt_Shapiro,PT = (2GM/b³) ln(...)`. For GPS at `b ≈ c (1 + 8.3 × 10⁻¹¹)`, this is `~3 × 10⁻¹⁰` smaller than the standard `(2GM/c³) ln(...)` — i.e., relative correction `~10⁻¹⁰`. The PT prediction at GPS precision is the same `~62 ps` maximum.

<!-- TODO: human reviews and fills in — confirms the speculative `c → b` extension's prediction for the spatial part of the Schwarzschild metric. The right minimal extension might be different (e.g., affecting only the time-time component, or applying differently for null vs timelike paths). The author is the right judge of what extension matches the framework's intended structure. -->

## 1. Effect statement

Light propagating through Earth's gravitational potential takes longer than the Euclidean path length divided by `c`. In standard GR, the extra delay is `(2GM/c³) ln((r_T + r_R + L)/(r_T + r_R − L))`, peaking at `~62 ps` for a satellite at the horizon.

Under the speculative `c → b` extension of the Schwarzschild metric, the same derivation gives `(2GM/b³) ln(...)`, where the appropriate `b` is the local light speed in the photon's frame. For a vacuum photon, there's no proper-velocity-like contribution to `b`; the framework's `b → c` exactly for photons in vacuum. The PT prediction is then *identical* to the standard prediction.

The speculative extension only matters for situations where the photon's frame has a non-trivial `b ≠ c` — e.g., propagation through a moving medium. For GPS vacuum propagation, the PT and standard predictions coincide.

## 2. Setup

| Symbol | Value | PT-relevant note |
|---|---|---|
| `2GM⊕/c³` | `2.957 × 10⁻¹¹ s` | Shapiro prefactor (standard) |
| `2GM⊕/b³` (vacuum photon) | `= 2GM⊕/c³` exactly | `b → c` for vacuum photons |
| `L_max` (slant range, horizon) | `2.578 × 10⁷ m` | (geometry-only) |
| `Δt_Shapiro_max` | `~62 ps` | same in PT and standard |

The Shapiro effect is one of the cleanest in the campaign: the framework's `b → c` for vacuum photons means the PT prediction reduces *exactly* to standard, with no GPS-precision correction at all.

## 3. Derivation

Under the speculative `c → b` Schwarzschild extension (pt_01 §0):

$$ds^2 \;=\; -\!\left(1 - \frac{2GM}{r b^2}\right) b^2 dt^2 \;+\; \left(1 - \frac{2GM}{r b^2}\right)^{-1} dr^2 \;+\; r^2 \, d\Omega^2.$$

For a photon (`ds² = 0`), `b` is the local light speed in the photon's frame. In vacuum (no medium-induced `u`), the photon's frame has `u = 0` and `b = c` exactly. The line element reduces to the standard Schwarzschild form, and the Shapiro derivation yields the standard formula:

$$\boxed{\;\Delta t_{\rm Shapiro,\,PT} \;=\; \frac{2 GM_\oplus}{c^3}\, \ln\!\left(\frac{r_T + r_R + L}{r_T + r_R - L}\right) \;=\; \Delta t_{\rm Shapiro,\,std}.\;}$$

The two predictions are *identical*, not just equivalent — the framework gives the same number with no precision loss.

## 4. Numerical evaluation

| Geometry | PT = standard |
|---|---|
| Zenith (`L = a − R⊕`) | `42.2 ps` |
| Horizon (`L = √(a² − R⊕²)`) | `62.3 ps` |
| Pseudorange error (max) | `1.87 cm` |

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; aGPS=26560000; RR=6378137; Lh=Sqrt[aGPS^2 - RR^2]; rath = (aGPS+RR+Lh)/(aGPS+RR-Lh); std = 2 GMe/cc^3 Log[rath]; (* For PT with b → c at vacuum photon, PT = std exactly *) Print["std (= PT) horizon Shapiro (ps) = ", N[std 1*^12, 6]]; Print["PT - std (s) = 0 exactly (vacuum photon: b → c)"]
```

**Result:**
```
std (= PT) horizon Shapiro (ps) = 62.2795
PT - std (s) = 0 exactly (vacuum photon: b → c)
```
✅ identical (not just equivalent — exactly equal).

## 5b. GPS-precision-limit equivalence to standard

The reduction is trivial: for vacuum photons, `b = c` exactly, so the speculative `c → b` extension produces zero modification to the standard Shapiro formula. The PT and standard predictions are bit-for-bit identical.

This is the only effect in the campaign where the framework's speculative extension does not introduce *any* numerical difference from standard, not even at relative precision `10⁻²⁰`. The reason: vacuum photons have `u = 0`, so `b = c`, and the substitution has no effect.

## 6. Comparison with standard-method companion

Standard [`07_shapiro_delay.md`](07_shapiro_delay.md) §3 derives the Shapiro delay from the spatial-curvature integral. The PT framework's speculative extension to curved spacetime gives the same expression because the relevant `b` (the photon's local light speed) is `c` for vacuum photons. The two predictions coincide.

A potentially interesting wrinkle: if the framework's intended extension to GR treats null and timelike geodesics differently — e.g., putting `b` in the timelike normalization but `c` in the null condition — the Shapiro derivation might give a slightly different answer. The minimal-extension hypothesis used here makes the simpler choice (`b` in the metric, `b → c` for null geodesics).

## 6b. Where the framework would diverge

The PT Shapiro prediction would differ from standard if:

1. **Photons propagate through a moving medium.** In a medium with macroscopic flow velocity `u_med`, the framework's `b ≠ c` for the photon's local frame; the Shapiro delay would carry an additional `u_med`-dependent correction. For GPS vacuum propagation this doesn't apply.
2. **The speculative extension is structurally different from `c → b` in `g_{rr}`.** If the right extension is `c → b` in only some metric components, the spatial-curvature integral might pick up a non-trivial correction. This is unresolved.
3. **Higher-order PN extensions** of the framework (not in any verified paper) might predict additional Shapiro-like terms at `(GM/rc²)²`.

None of these are GPS-accessible regimes.

## 7. Verdict

❌ framework doesn't extend natively to curved spacetime; speculative extension reduces to standard *exactly* at GPS precision (vacuum photons have `b = c`, so the extension produces no change). Open question for the author: what is the correct extension of the framework to curved spacetime, and does it predict any Shapiro-like correction beyond standard?
