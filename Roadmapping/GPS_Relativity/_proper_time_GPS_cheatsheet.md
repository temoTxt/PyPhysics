# Proper-time substitution rules — GPS applicability cheat sheet

This is the GPS-specific cheat sheet for the proper-time companion campaign. It restates the [Gill–Zachary substitution rules](../Electromagnetism/_proper_time_cheatsheet.md) verified in [`Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) and adds the GPS-specific applicability decisions: which rules apply to which effect, and where curved-spacetime extension is needed.

## 1. The substitution rules (recap)

For a particle with observer-frame velocity `w = dx/dt` and proper-frame velocity `u = dx/dτ`, the framework defines

$$b^{2} = c^{2} + u^{2},$$

so `b ≥ c` always, and `b = c` for `u = 0`. The verified rules are:

| Rule | Form | Origin |
|---|---|---|
| Velocity duality | `w/c = u/b` | Maxwell-paper Eq. (1) |
| Time-derivative duality | `(1/c) ∂_t = (1/b) ∂_τ` | Eq. (2) |
| Maxwell equations (PT form) | `∇·E = 4πρ`, `∇·B = 0`, `∇×E = −(1/b)∂_τ B`, `∇×B = (1/b)∂_τ E + (4π/b)ρu` | Eq. (3′) |
| Chain rule | `∂_τ(1/b) = −(u·a)/b³` | derived |
| Dissipative coefficient | `−(u·a)/b⁴ · ∂_τ` | Eq. (4) |
| Modified Lorentz force | `F = e[E + (u/b) × B] + ...` | Eq. (18) |

The rules are SR + flat-space Maxwell + flat-space dynamics. They do **not** include a curved-spacetime extension; the Schwarzschild geometry of effects 02, 04, 07, 08 lies outside the verified corpus.

## 2. GPS-specific reductions

At GPS satellite kinematics:

| Quantity | Symbol | Value | Reduction to standard |
|---|---|---|---|
| Satellite orbital speed (observer) | `v` | `3 874 m/s` | (input) |
| Lorentz factor | `γ` | `1 + 8.349 × 10⁻¹¹` | `γ − 1 = v²/(2c²) + 𝒪(v⁴/c⁴)` |
| Satellite proper velocity | `u = γv` | `≈ v + 𝒪(v³/c²)` | `u ≈ v` to first order in `v²/c²` |
| Effective speed of light | `b = √(c² + u²)` | `c · (1 + 8.349 × 10⁻¹¹)` | `b → c` as `u² → 0` |
| `b/c − 1` | | `≈ v²/(2c²) ≈ 8.35 × 10⁻¹¹` | matches the velocity-time-dilation magnitude |
| `b² − c² = u²` | `≈ v²` | `1.50 × 10⁷ m²/s²` | (the SR correction itself) |

The numerical consequence: every quantity in the proper-time formulation that is naively "`c` for standard, `b` for PT" differs by exactly the magnitude of the velocity-time-dilation correction. **This is why the PT formulation reproduces standard SR predictions for GPS at the operational level** — the same `1 − v²/(2c²) = c/b` factor appears in both, just attributed to different bookkeeping.

## 3. GPS-effect applicability matrix

| Effect | Standard derivation cites | Framework applies? | PT prediction at GPS |
|---|---|---|---|
| **01** Schwarzschild metric setup | Schwarzschild line element + PN expansion | ⚠ partial (kinematic ✅; gravitational ❌ — minimal extension `c → b` in metric flagged speculative) | `dτ/dt = 1 − GM/(rc²) − u²/(2b²)` reduces to standard at `b ≈ c` |
| **02** Gravitational time dilation | Schwarzschild redshift only | ❌ GR-only | `+45.65 μs/day` (reduces from speculative extension) |
| **03** Velocity time dilation | SR `1/γ` | ✅ yes (direct PT application) | `dτ/dt = c/b` — exact equivalent to `1/γ` |
| **04** Net headline `±38 μs/day` | combined Schwarzschild + SR | ⚠ partial | `+38.57 μs/day` (bookkeeping reorganised; same number) |
| **05** Eccentricity periodic | Kepler orbit + SR | ⚠ partial (gravity ❌, orbital dynamics ✅, third term contributes) | matches standard; third term `(u·a)/b⁴` vanishes for circular orbit, ~`10⁻³⁴` for `e=0.02` |
| **06** Sagnac effect | Flat-space rotating-frame Maxwell | ✅ yes (direct PT application) | `Δt = −2ω·A/b²`; reduces to standard at `b ≈ c` |
| **07** Shapiro delay | Schwarzschild null geodesic | ❌ GR-only | `~62 ps` (reduces from speculative extension) |
| **08** J₂ oblateness | Schwarzschild + multipole | ❌ GR-only | `+0.033 μs/day` (reduces from speculative extension) |
| **09** Full clock equation | (synthesis) | (synthesis) | full PT operational equation |

## 4. The two open framework questions for this campaign

These are flagged in [`_proper_time_cheatsheet.md` §4a](../Electromagnetism/_proper_time_cheatsheet.md) as load-bearing open questions for the framework. Both surface in this campaign:

1. **Covariance of `u · a`.** The third-term contribution `(u · a)/b⁴` in effect pt_05 depends on `u · a` having a frame-independent operational meaning. The covariance of this 3-vector dot product under the proper-time boost (Maxwell-paper Eq. (11)) has not been derived in the campaign. **For GPS**, `u · a = 0` for circular orbits (radial acceleration ⊥ tangential velocity), so the third term vanishes regardless of how it transforms. **For elliptical GPS orbits** (`e ≤ 0.02`), `u · a ≠ 0` but the contribution is `~ 10⁻³⁴`, well below operational precision.

2. **GR extension of the framework.** The framework as-published is SR + flat-space. The Schwarzschild-based effects (02, 04, 07, 08) require an extension that has not been published. The minimal-extension hypothesis used in pt_02, pt_04, pt_07, pt_08 is: *the Schwarzschild metric, written in any local frame, has its `c` replaced by the local `b = √(c² + u²)`*. This is a substantive AI proposal (not a verified Gill–Zachary result) — every doc that invokes it carries the `<!-- TODO -->` Crocco substantive-AI block.

## 5. The GPS-precision-limit equivalence theorem (claimed; not formally proven)

The proper-time formulation reduces to standard SR + flat-space Maxwell predictions whenever:

- The particle's `u² ≪ c²` (so `b ≈ c` to leading order in `u²/c²`).
- The particle's acceleration `a` satisfies `(u · a) τ_signal / b⁴ ≪ 1` (so the third term integrates to below the noise floor over the signal time `τ_signal`).

For GPS:
- `u²/c² ≈ 1.67 × 10⁻¹⁰` ✅
- `(u · a) τ_sig / b⁴ ~ (10⁻³⁴) × 0.07 s ~ 10⁻³⁵` ✅ (the satellite's acceleration is gravitational, so `u · a = u · g`; for circular orbit, `u ⊥ g` and the product vanishes identically)

Hence at GPS-relevant precision, the framework reproduces standard predictions. The campaign records the explicit per-effect reductions in §5b of each pt_*.md document.
