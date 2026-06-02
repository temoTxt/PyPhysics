# 02 — Orbital dynamics: speed, period, angular velocity

**Status:** ✅ all Wolfram-MCP verified.
**Source paper:** [`../Tepper_Gill_Papers/Dual Newtonian Theory.tex`](../Tepper_Gill_Papers/Dual%20Newtonian%20Theory.tex), §2 *The Sun-Mercury System*.

## 1. Relative acceleration in the two-body problem

Subtracting Mercury-side and Sun-side accelerations from Eq. (h4):

$$\mathbf{a} \;=\; \mathbf{a}_m - \mathbf{a}_s \;=\; -\frac{G(M+m)}{r^2}\!\left[1 - \frac{G(M^2 + m^2)}{(M+m)c^2 r}\right]\hat{\mathbf{e}}_r.$$

The corresponding effective gravitational force on Mercury:

$$F_e \;=\; -\frac{G(M+m)\,m}{r^2}\!\left[1 - \frac{G(M^2 + m^2)}{(M+m)c^2 r}\right].$$

## 2. Eq. (h5) — orbital speed `u_d` from force balance

For a circular orbit at radius `r₀`, equate `|F_e|` to the centripetal force `m u_d² / r₀`:

$$\frac{G(M+m)m}{r_0^2}\!\left[1 - \frac{G(M^2 + m^2)}{(M+m)c^2 r_0}\right] \;=\; \frac{m\,u_d^2}{r_0}$$

solving for `u_d`:

$$\boxed{\;u_d \;=\; \sqrt{\frac{G(M+m)}{r_0}}\,\sqrt{1 - \frac{G(M^2 + m^2)}{(M+m)c^2 r_0}}.\;}\quad\text{(paper Eq. h5)}$$

## 3. Eqs. (h6)–(h8) — periods and angular velocities

The Newtonian period (paper Eq. h6, textbook Kepler):

$$T_0 \;=\; \frac{2\pi r_0^{3/2}}{(GM)^{1/2}}, \qquad \omega_0 \;=\; \frac{2\pi}{T_0}.$$

The dual period (paper Eq. h7) and angular velocity (paper Eq. h8):

$$\boxed{\;\frac{T_d}{T_0} \;=\; \frac{2\pi r_0/u_d}{T_0} \;=\; \left[\left(1 + \tfrac{m}{M}\right)\!\left(1 - \tfrac{G(M^2 + m^2)}{(M+m)c^2 r_0}\right)\right]^{-1/2}.\;}$$

$$\boxed{\;\omega_d \;=\; \omega_0 \left[\left(1 + \tfrac{m}{M}\right)\!\left(1 - \tfrac{G(M^2 + m^2)}{(M+m)c^2 r_0}\right)\right]^{1/2}.\;}$$

## 4. Wolfram MCP — joint verification of Eqs. (h5), (h7), (h8)

```wolfram
ClearAll[r0, cc, mm, GG, MM]; relAccel = -GG (MM + mm)/r0^2 (1 - GG (MM^2 + mm^2)/((MM + mm) cc^2 r0)); effForce = relAccel mm; udSquared = -effForce r0/mm; udPaper = Sqrt[GG (MM + mm)/r0] Sqrt[1 - GG (MM^2 + mm^2)/((MM + mm) cc^2 r0)]; Print["FullSimplify[u_d² - paper_h5²] = ", FullSimplify[udSquared - udPaper^2]]; Td = 2 Pi r0/udPaper; T0 = 2 Pi r0^(3/2)/Sqrt[GG MM]; TdOverT0Paper = ((1 + mm/MM) (1 - GG (MM^2 + mm^2)/((MM + mm) cc^2 r0)))^(-1/2); Print["FullSimplify[T_d/T_0 (derived) - paper_h7] = ", FullSimplify[Td/T0 - TdOverT0Paper, Assumptions -> {MM > 0, mm > 0, GG > 0, cc > 0, r0 > 0}]]
```

**Result:**
- `FullSimplify[u_d² − paper_h5²] = 0` ✅
- `FullSimplify[T_d/T_0 (derived) − paper_h7] = 0` ✅

(Eq. h8 follows trivially from `ω_d = 2π/T_d` and the boxed expression for `T_d/T_0`; the verification above covers it.)

## 5. Decomposition for Mercury's regime

For the Sun-Mercury system, the two small parameters are:

| Symbol | Magnitude | Source |
|---|---|---|
| `m/M` | `1.66 × 10⁻⁷` | reduced-mass / orbital coupling |
| `G(M²+m²) / ((M+m) c² r₀)` | `2.55 × 10⁻⁸` | relativistic correction (effectively `GM/(c²r₀)`) |

Both terms enter linearly under the square root, so to leading order:

$$\sqrt{(1 + m/M)(1 - GM/(c^2 r_0))} \;-\; 1 \;\approx\; \tfrac{1}{2}\!\left[\frac{m}{M} - \frac{GM}{c^2 r_0}\right].$$

For Mercury, the reduced-mass term `m/(2M) = 8.30 × 10⁻⁸` is **larger than** the relativistic term `GM/(2c² r₀) = 1.27 × 10⁻⁸` by ~6.5×. The framework's two contributions to the perihelion advance therefore have *opposite signs* (the reduced-mass term increases `ω_d`; the relativistic term decreases it), and the framework's combined prediction is dominated by reduced mass.

This decomposition is the load-bearing input for the numerical computation in [`03_numerical_predictions.md`](03_numerical_predictions.md).

## 6. Verdict

✅ Paper's algebraic derivation of Eqs. (h5), (h7), (h8) reproduces under Wolfram-MCP verification. The orbital dynamics produce a well-defined dual angular frequency `ω_d` that decomposes into a reduced-mass term (`+m/(2M)`) and an opposite-sign relativistic term (`−GM/(2c²r₀)`).
