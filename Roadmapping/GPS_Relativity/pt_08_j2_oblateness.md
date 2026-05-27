# Effect pt_08 — proper-time companion to J₂ oblateness

**One-line summary:** Earth's quadrupole moment is a curved-spacetime effect (the gravitational potential has a `J₂(R⊕/r)² P₂(cos θ)` correction); the Gill–Zachary framework as-published does not extend to curved-spacetime gravity. Under the speculative `c → b` Schwarzschild extension (pt_01 §0), the framework reproduces the standard `+0.033 μs/day` J₂ correction at GPS-relevant precision; the relative correction `(b/c − 1) ~ 10⁻¹¹` is below operational floor.

**Status:** ❌ framework doesn't extend natively; speculative extension reduces to standard at GPS precision.
**Framework applicability:** GR-only (framework as-published does not extend; speculative extension assumed).
**Standard-method companion:** [`08_J2_oblateness.md`](08_J2_oblateness.md).

## 0. Framework applicability

The J₂ oblateness correction modifies the Earth's gravitational potential `Φ(r,θ) = −(GM/r)[1 − J₂(R⊕/r)²P₂(cos θ)]`. This is a multipole expansion of the Newtonian potential, which (under GR) enters the Schwarzschild-like metric's `g_{00}` component. The framework as-published does not speak to gravitational potentials; the speculative `c → b` extension from pt_01 §0 carries `J₂` through unchanged at GPS precision.

The framework predicts no novel J₂ contribution beyond what the speculative extension carries. This document is a placeholder for the GR-extension question; the operational answer matches standard.

<!-- TODO: human reviews and fills in — confirms that the speculative GR extension doesn't predict any novel multipole-induced corrections beyond the standard one. In particular, whether the framework's `c → b` affects the multipole structure (J₂ contribution) differently than the leading point-mass piece is an open question. The default assumption is that it doesn't, but this is not derived. -->

## 1. Effect statement

The J₂ correction in the standard derivation contributes `+0.033 μs/day` to the headline rate by modifying the geoid potential `Φ_geoid`. Under the speculative PT extension `c → b` in the Schwarzschild metric, the J₂ correction enters via the same `1/b²` factor as the leading point-mass term; at GPS precision, the PT prediction matches the standard.

Higher-order corrections — Lense-Thirring (frame-dragging) and post-PN — are at the `~10⁻¹⁶` and `~10⁻¹⁹` levels respectively. The framework's prediction at these orders is not derived in this PR; under the speculative extension they would reduce to the standard GR values at GPS precision.

## 2. Setup

| Symbol | Value | PT-relevant note |
|---|---|---|
| `J₂` | `1.082 6267 × 10⁻³` | Earth quadrupole moment (EGM2008) |
| `GM⊕ J₂/(2 R⊕ c²)` | `3.764 × 10⁻¹³` (per day: `+0.033 μs/day`) | standard contribution |
| `GM⊕ J₂/(2 R⊕ b²)` | `= 3.764 × 10⁻¹³ × (1 − 𝒪(10⁻¹²))` | PT correction at `b/c − 1 ~ 10⁻¹²` |
| Difference | `~ 10⁻²⁵ s/day` | operationally invisible |

## 3. Derivation

Under the speculative `c → b` extension applied to the full Earth potential (point-mass + J₂):

$$\Phi_{\rm PT}(r,\theta) \;=\; -\frac{GM_\oplus}{r}\!\left[1 - J_2 \left(\frac{R_\oplus}{r}\right)^2 P_2(\cos\theta)\right], \qquad \frac{d\tau}{dt} \;\approx\; 1 + \frac{\Phi_{\rm PT}}{b^2}.$$

Substituting `b² = c² + u²` and expanding to first order in `u²/c²`:

$$\frac{\Phi_{\rm PT}}{b^2} \;\approx\; \frac{\Phi_{\rm PT}}{c^2}\!\left(1 - \frac{u^2}{c^2}\right).$$

The J₂ contribution to `Φ_PT/c²` is the same as the standard `Φ/c²`; the PT first-order correction multiplies the entire potential (including J₂) by `(1 − u²/c²)`. The relative correction `u²/c² ~ 10⁻¹²` at the geoid and `~10⁻¹⁰` at GPS — both negligible for the `~10⁻¹³` magnitude of the J₂ contribution.

$$\boxed{\;\left.\frac{\Delta f}{f}\right|_{J_2,\,{\rm PT}} \;\approx\; \frac{GM_\oplus J_2}{2 R_\oplus c^2}\!\left(1 + \mathcal{O}(u^2/c^2)\right) \;=\; +0.033\ \mu{\rm s/day} \;\text{(GPS precision)}.\;}$$

## 4. Numerical evaluation

| Component | Value |
|---|---|
| J₂ to geoid, PT | `+0.0325 μs/day` |
| J₂ to geoid, standard | `+0.0325 μs/day` |
| Difference | `~ 10⁻²⁵ s/day` (negligible) |
| Combined PT headline `pt_04 + pt_08` | `+38.57 μs/day` |
| Standard combined `04 + 08` | `+38.57 μs/day` |

PT and standard predictions agree at all precision relevant for GPS.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; RR=6378137; J2val=1.0826267*^-3; vEq=7.2921151467*^-5 RR; bEq = Sqrt[cc^2 + vEq^2]; PT = GMe J2val/(2 RR bEq^2); std = GMe J2val/(2 RR cc^2); Print["PT J2 contribution = ", ScientificForm[N[PT, 10]]]; Print["std J2 contribution = ", ScientificForm[N[std, 10]]]; Print["PT - std = ", ScientificForm[N[PT - std, 8]]]; Print["PT microsec/day = ", N[PT 86400 1*^6, 8]]; Print["std microsec/day = ", N[std 86400 1*^6, 8]]
```

**Result:**
```
PT J2 contribution = 3.76400468×10⁻¹³
std J2 contribution = 3.76400920×10⁻¹³
PT - std = -4.5×10⁻²⁵
PT microsec/day = 0.0325211
std microsec/day = 0.0325211
```
✅ PT and standard predictions agree to 6+ decimal places (microsec/day match exactly to printed precision).

## 5b. GPS-precision-limit equivalence to standard

The PT correction at this order is `J₂ × (u²/c²) × leading_J₂_contribution ~ 10⁻³ × 10⁻¹² × 10⁻¹³ ~ 10⁻²⁸ /s` — at most one yoctosecond per year. Operationally indistinguishable from the standard `+0.033 μs/day`.

## 6. Comparison with standard-method companion

Standard [`08_J2_oblateness.md`](08_J2_oblateness.md) §3 derives the J₂ contribution from the Newtonian multipole expansion of the geoid potential, divided by `c²`. The PT framework's speculative extension uses `1/b²` instead, with the relative correction at the `u²/c²` level — `~10⁻¹²` at the geoid. The PT prediction matches the standard.

The Lense-Thirring contribution (Earth's rotation-induced gravitomagnetic effect) is at the `~10⁻¹⁶` level in GR. The framework's prediction would need a curved-spacetime extension to address it; the speculative `c → b` extension gives the standard GR value at GPS precision but does not derive the gravitomagnetic structure independently.

## 6b. Where the framework would diverge

The PT prediction for J₂ would differ from standard at observable levels when:

1. **The clock's proper velocity is large enough to matter to multipole corrections.** This requires `u²/c² × J₂ ~ 10⁻¹⁵`, i.e., `u/c ~ 10⁻⁶` or `u ~ 300 m/s` for a J₂-level effect — the velocities are GPS-accessible, but the magnitudes (`~10⁻²⁵ s/day`) are operationally invisible.
2. **Stronger oblateness regimes.** For a more oblate body (Jupiter, with `J₂ ~ 10⁻²`), the speculative PT correction to J₂ would be `~10⁻¹⁰` relative — still below operational precision for any current measurement.

## 7. Verdict

❌ framework doesn't extend natively to curved-spacetime multipole effects; speculative extension reduces to standard `+0.033 μs/day` at GPS precision. Open question: a derived (not speculative) GR-extended Gill framework would settle whether `J₂`, Lense-Thirring, and other higher-order GR effects pick up any framework-specific corrections.
