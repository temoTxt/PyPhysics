# Implementation plan: GPS relativity — standard-method derivation

**Tracks:** issue [#57](https://github.com/temoTxt/PyPhysics/issues/57).
**Status:** plan; per-effect documents to follow in this same PR.
**Successor:** a proper-time companion campaign will repeat the derivation under the Gill–Zachary substitution rules ([[two_mathematically_equivalent_versions_of_maxwells_equations]]); that work is **not** in scope here.

This is the *standard* relativistic derivation — Schwarzschild geometry, post-Newtonian expansion to order `1/c²`, and the textbook decomposition into nine separately-named effects. The campaign's value lies in covering the full Ashby (2003) treatment in one place so the proper-time companion has a complete point-for-point reference to confront.

---

## 1. Goal and scope

Reproduce the full standard-method derivation of the relativistic clock correction for a GPS satellite, from the Schwarzschild line element through to the operational clock equation programmed into the GPS ground segment. The campaign covers nine separately-named effects (the headline ±38 μs/day decomposition plus six smaller corrections), each in its own per-effect document with a Wolfram-MCP numerical check, a comparison against the canonical Ashby (2003) treatment, and a verdict.

The deliverable is a complete reference document for the *standard* answer. The companion proper-time campaign (deferred to a future PR) will work the same nine effects under the Gill–Zachary substitution rules and record where the predictions agree numerically and where the proper-time formulation surfaces extra structure.

### Why this earns its own thread

- **Self-contained reference.** The repository has no consolidated GPS-relativity derivation. The verification campaign covers the Gill papers' equations; the Electromagnetism and Quantum Mechanics campaigns work textbook problems. A standalone clock-correction derivation belongs in neither.
- **Setup for the proper-time companion.** The proper-time predictions are mostly the same numbers (the dominant effects are independent of the third-term structure), but the *framework's* clock equation will look structurally different. We cannot do that comparison without first having the standard derivation written down.
- **Pedagogical value.** The standard derivation is scattered across Ashby's *Living Reviews* article, Misner-Thorne-Wheeler §40, and Weinberg's *Gravitation and Cosmology* §9; collecting it in one document is itself useful.

### Explicit non-goals

- Not a re-derivation of the Schwarzschild metric (cited, not re-derived).
- Not the proper-time companion — that is a separate planned PR.
- Not an empirical-comparison sub-investigation (no GPS performance data analysis).
- Not a critique of Ashby (2003) — the campaign reproduces his treatment with full agreement; any disagreement would be flagged for re-derivation, not for publication.

---

## 2. Repository structure

New campaign folder, sibling to `Electromagnetism/` / `Quantum_Mechanics/` / `Equation_Verification/`:

```
Roadmapping/GPS_Relativity/
├── README.md                              # status table + scope + conventions
├── _template_effect.md                    # per-effect document template
├── 01_schwarzschild_metric_setup.md       # PN-expanded line element
├── 02_gravitational_time_dilation.md      # static-field clock rate
├── 03_velocity_time_dilation.md           # kinematic-frame clock rate
├── 04_net_38us_per_day.md                 # headline combined offset
├── 05_eccentricity_correction.md          # e·sin(E) periodic term
├── 06_sagnac_effect.md                    # rotating-frame propagation
├── 07_shapiro_delay.md                    # gravitational propagation delay
├── 08_J2_oblateness.md                    # quadrupole + higher PN
└── 09_full_clock_equation.md              # synthesis
```

Bibliography note scaffolded under `Roadmapping/History/Bibliography/Retrospective/`:

- `ashby2003_relativity_in_gps.md` — Ashby, N. (2003). *Relativity in the Global Positioning System*. Living Reviews in Relativity, 6(1).

---

## 3. Per-effect document template

Every per-effect document follows the structure in [`_template_effect.md`](../../Roadmapping/GPS_Relativity/_template_effect.md):

1. **Effect statement** — one-paragraph plain-English summary of what's being computed.
2. **Setup** — geometry, frame choice, parameters, what's held constant.
3. **Derivation** — step-by-step algebra, intermediate forms boxed.
4. **Numerical evaluation** — substitute GPS parameters; record the number to the precision Ashby reports.
5. **Wolfram MCP check** — one-line `WolframLanguageEvaluator` invocation reproducing the number from scratch, with the gotcha rules from [CLAUDE.md](../../CLAUDE.md) baked in (single-line code; avoid `V`/`e` as variable names).
6. **Comparison with Ashby (2003)** — section/equation reference, agreement statement.
7. **Verdict** — ✅ reproduces / ⚠ minor numerical disagreement / ❌ structural disagreement.

The Crocco compliance disclosure for this campaign: this is a *pragmatic* AI use under [§5 of CROCCO_COMPLIANCE.md](../../Roadmapping/Tooling/CROCCO_COMPLIANCE.md). The mathematics is textbook physics being reproduced; no novel theoretical claims are made. Per-paragraph `<!-- TODO -->` blocks are *not* required for the derivations themselves; the *framing* paragraphs at the top of the README and §9 (synthesis) carry the standard substantive `<!-- TODO -->` block.

---

## 4. Decision points (deferred to author)

- **D1.** Do we want a Manim animation tying the nine effects together? Default answer: no, until the proper-time companion exists — the animation's value is in the side-by-side comparison.
- **D2.** Do we adopt a single units convention throughout (SI) or also record a Gaussian conversion table? Default: SI only — the GPS literature is uniformly SI, and the proper-time companion's substitution rules are unit-agnostic.

---

## 5. Honest scoping

This is the *standard* derivation. It does **not** validate the proper-time formulation; it does not deny it either. The proper-time companion will revisit each effect, and any structural difference between the two formulations will be flagged at that time. The standard derivation has been worked out by many authors over many decades (Ashby, Allan, Hodge, Will); this document's value is consolidation, not novelty.
