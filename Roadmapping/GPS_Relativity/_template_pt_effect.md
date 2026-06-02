# Per-effect proper-time document template

This template defines the structure of every per-effect proper-time document under `Roadmapping/GPS_Relativity/pt_0N_*.md`. It mirrors the standard-method [`_template_effect.md`](_template_effect.md) with three additions: §0 (Framework applicability), §5b (GPS-precision-limit equivalence), and §6b (Where the framework would diverge).

## Per-effect template (the body of every `pt_0N_*.md`)

````markdown
# Effect NN — proper-time companion to <standard-method effect title>

**One-line summary:** {Plain-English statement of what the proper-time framework says about this effect. Always cite the standard-method companion file.}

**Status:** ✅ matches standard at GPS precision / ⚠ structural difference / ❌ framework doesn't extend natively.
**Framework applicability:** {SR + flat-space EM / GR-only / mixed}.
**Standard-method companion:** [`0N_<title>.md`](0N_<title>.md).

## 0. Framework applicability

{One paragraph: which Gill–Zachary substitution rules apply to this effect, and which don't. If the effect is GR-only (Schwarzschild geometry, gravitational potential), state explicitly that the framework as-published does not extend to curved spacetime, and that the §3 derivation carries a *speculative minimal extension* — typically `c → b` substituted into the Schwarzschild metric, with `b → c` at low velocity. This paragraph carries `<!-- TODO: human reviews and fills in -->` per [Crocco §1](../Tooling/CROCCO_COMPLIANCE.md) when the doc is in the substantive (GR-extension) category.}

## 1. Effect statement

{Mirror the standard-method companion's §1. State what's being computed and in what frame. The proper-time framework uses the particle's clock `τ` as the natural parameter and `u = dx/dτ` as the natural velocity (proper velocity); restate the problem in those terms.}

## 2. Setup

{Mirror the standard-method §2. Add a small table of the substitution-rule consequences: typically `b² = c² + u²`, and for satellite circular orbit `u ≈ γv ≈ v + 𝒪(v³/c²)`, so `b ≈ c(1 + v²/(2c²))`.}

## 3. Derivation

{Apply the substitution rules (cited from [`_proper_time_GPS_cheatsheet.md`](_proper_time_GPS_cheatsheet.md)) to the standard-method derivation. Box the proper-time result in its as-derived form (typically with `b`, `u`, `τ` not yet reduced).}

## 4. Numerical evaluation

{Substitute parameters (`b = √(c² + u²)`, `u = √(GM/a)`, etc.) and reduce to a number. Report to the same precision as the standard-method companion.}

## 5. Wolfram MCP check

```wolfram
{Single-line Wolfram Language expression with the same gotcha rules. The check should reproduce the same number as the standard-method companion (within the rounding-precision of `b ≈ c` reduction).}
```

**Result:** `{numerical output}` ✅ matches standard at GPS precision.

## 5b. GPS-precision-limit equivalence to standard

{Explicit algebraic reduction showing that the proper-time prediction matches the standard prediction when `u² ≪ c²` is applied. Typical form: expand `b² = c² + u²` as `b² = c²(1 + u²/c²)`, take `b/c ≈ 1 + u²/(2c²)`, and show that the proper-time expression reduces to the standard one term-by-term. State the order of the reduction error (typically `𝒪(u⁴/c⁴) ≈ 10⁻²⁰` for GPS).}

## 6. Comparison with standard-method companion

{One paragraph: does the PT derivation match the standard formulation step-by-step, or take a structurally different path? If different, are the algebraic forms convertible via the substitution rules? Cite the standard-method companion's verdict line for the operational number.}

## 6b. Where the framework would diverge

{One paragraph: identify a regime where the framework's prediction is *operationally distinguishable* from standard. Typical answers: (a) at speeds where `u² ~ c²` (relativistic particle, not GPS), (b) in strong gravitational fields where the GR extension matters, (c) at high accelerations where the third-term `(u·a)/b⁴` contribution becomes observable. Cite the relevant cross-experiment investigation (e.g., issue [#43](https://github.com/temoTxt/PyPhysics/issues/43)) if applicable.}

## 7. Verdict

✅ matches standard-method companion [`0N_<title>.md`](0N_<title>.md) §7 at GPS precision.
or ⚠ matches numerically at GPS precision but algebraic structure differs — see §6.
or ❌ framework doesn't extend natively; speculative extension reduces to standard at GPS precision (see §0 and §3).
````

## Voice compliance

Prose in §§0, 1, 6, 6b, 7 follows the voice guide at [`Roadmapping/Tooling/VOICE_MATCH_GILL.md`](../Tooling/VOICE_MATCH_GILL.md). The substantive-AI disclosure in §0 (for GR-extension docs) is a Crocco compliance requirement, not a voice-style choice — the `<!-- TODO -->` block is mandatory.

## A note on the template's stability

This template is provisional. Refinements identified during the nine-doc draft are folded back here.
