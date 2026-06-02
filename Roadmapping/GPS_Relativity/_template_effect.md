# Per-effect document template

This template defines the structure of every per-effect document under `Roadmapping/GPS_Relativity/0N_*.md`. Each effect is presented in this fixed format so that a reader can locate the same kind of content in the same place across the campaign. Refinements identified during the campaign's draft are folded back into this file.

## Per-effect template (the body of every `0N_*.md`)

````markdown
# Effect NN — short title

**One-line summary:** {plain-English what this effect is and what its magnitude is for GPS}.

**Status:** ✅ drafted / ⚠ flagged / ❌ blocked.
**Ashby (2003) cross-reference:** §X.Y, Eq. (NN).

## 1. Effect statement

{One paragraph stating the physical effect, the frame in which it is naturally computed, and the parameters it depends on. Avoid algebra here — this is the orientation paragraph.}

## 2. Setup

{Geometry, frame choice, parameters held fixed, what counts as the "ground reference". Include a small table of any parameters specific to this effect not already in the README's parameter table.}

## 3. Derivation

{Step-by-step algebra. Each intermediate form on its own display-math line. Use `\boxed{...}` for the final result. Cite Ashby (2003) or MTW or Weinberg at each non-trivial step.}

## 4. Numerical evaluation

{Substitute the parameter values from the README's table. Compute the magnitude. Report to the precision Ashby reports (typically 4–5 sig figs).}

## 5. Wolfram MCP check

```wolfram
{Single-line Wolfram Language expression reproducing the numerical answer from scratch. Follows the three gotcha rules from CLAUDE.md: (a) single-line code joined with `;`; (b) avoid `V` as a variable name (Vanadium); (c) avoid `e` as a variable name (Euler).}
```

**Result:** `{numerical output from Wolfram MCP}` ✅ matches §4.

## 6. Comparison with Ashby (2003)

{One paragraph: which section / equation in Ashby is the canonical reference? Does our derivation match step-by-step or does it take a shortcut? Are there any sign conventions that differ?}

## 7. Verdict

✅ reproduces Ashby (2003) §X.Y Eq. (NN) at quoted precision.
or ⚠ minor numerical disagreement at the {N}-th sig fig — likely traceable to {convention}.
or ❌ structural disagreement — see [`FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md).
````

## Voice compliance

Prose in §§1, 6, 7 follows the voice guide at [`Roadmapping/Tooling/VOICE_MATCH_GILL.md`](../Tooling/VOICE_MATCH_GILL.md). The derivation steps in §3 are operational and need not match Gill's voice. The numerical evaluation in §4 and the Wolfram check in §5 are mechanical and need not match.

## A note on the template's stability

This template is provisional. The nine per-effect documents in the first PR will exercise it on the full range of derivation types — closed-form algebraic (effects 02–04), orbit-averaged (effect 05), path-dependent line integrals (effects 06–07), and orbit-secular (effect 08). Any refinement identified during drafting is folded back into this file before the proper-time companion campaign begins.
