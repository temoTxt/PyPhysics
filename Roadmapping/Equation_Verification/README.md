# Equation Verification Campaign

Tracks issue [#3](https://github.com/temoTxt/PyPhysics/issues/3) — systematic equation-by-equation verification of the Gill paper corpus under [../Tepper_Gill_Papers/](../Tepper_Gill_Papers/).

## Goals

1. **Verify** every equation using a Wolfram Mathematica MCP tool as the reproducible backend.
2. **Expand** each derivation to a level a first-year physics graduate student can follow end-to-end.
3. **Compare** every Gill-framework equation to the standard textbook form a grad student would already know (Jackson, Sakurai, Peskin & Schroeder, Goldstein, Weinberg).

## Status Table

Legend: ⬜ not started · 🟨 in progress · ✅ complete · ⚠ blocked

| Paper | File | Verified / Total | Status | Notes |
|---|---|---:|---|---|
| Two Mathematically Equivalent Versions of Maxwell's Equations | [Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md](Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) | 23 / 24 | ✅ + 🔴 | Eqs (1)–(23) verified; **Eq (24) typos identified** (missing `c` in `eℏΣ·B` denominator + missing `V²/(2mc²)` term) |
| Dual Relativistic Quantum Mechanics I | [Dual_Relativistic_Quantum_Mechanics_I.md](Dual_Relativistic_Quantum_Mechanics_I.md) | ~28 / ~30 | ✅ + 🔴 | Eqs (I.1)–(III.22) verified; **Eq (III.22) g-factor claim fails** — paper's $r_e = 0.499857150068631\,r_0$ does not reproduce experimental $g_e = -2.00231930436256$ (formula gives $-2.0005714$ instead). Required $r_e \approx 0.4994205\,r_0$. |
| FOUNDATIONS FOR QED I: MATHEMATICAL | — | 0 / ? | ⬜ | |
| FoundationsII-Classical | — | 0 / ? | ⬜ | |
| The Classical Electron Problem | [The_Classical_Electron_Problem.md](The_Classical_Electron_Problem.md) | novel content ✅ | ✅ + ⚠️ | Liénard-Wiechert full chain (3.36, 3.40, 3.41, 3.45-3.46) ✅; Sec 3.3 Larmor analog ✅; **Eq (4.16) sign typo** flagged |
| Analytic Representation of The Dirac Equation | [Analytic_Representation_of_The_Dirac_Equation.md](Analytic_Representation_of_The_Dirac_Equation.md) | core results ✅ | ✅ | Green's-function-based analytical separation of Dirac (Sec II); operator-algebra chain (Sec V) parallels DRQM I Sec III. No new findings. |
| Analytic Representation of The Square-Root Operator | [Analytic_Representation_of_The_Square-Root_Operator.md](Analytic_Representation_of_The_Square-Root_Operator.md) | core integral identities ✅ | ✅ | Yosida + Schulman path-integral construction; Bessel-function identities (9b, 11, 18, 32, 44) all ✅ via Wolfram MCP. No new findings. |
| Relativistic Transformations of Thermodynamics | — | 0 / ? | ⬜ | |
| Mathematical Concepts in Physics | — | 0 / ? | ⬜ | |
| On the physical and mathematical foundations of quantum physics via functional integrals | — | 0 / ? | ⬜ | |
| Constructive Representation Theory for the Feynman Operator Calculus | — | 0 / ? | ⬜ | |
| Foundations for Relativistic Quantum Theory I — Feynman's Operator Calculus and the Dyson Conjectures | — | 0 / ? | ⬜ | |
| A sufficiency class for global (in time) solutions to the 3D Navier-Stokes equations II | — | 0 / ? | ⬜ | Pure-math; verify only load-bearing results (see scope below) |
| Global (in Time) Solutions to the 3D-Navier-Stokes Equations on R³ | — | 0 / ? | ⬜ | Pure-math |
| Global solutions to the homogeneous and inhomogeneous Navier-Stokes equations | — | 0 / ? | ⬜ | Pure-math |
| Adjoint Operators on Banach Spaces | — | 0 / ? | ⬜ | Pure-math |
| Adjoint for Operators in Banach Spaces | — | 0 / ? | ⬜ | Pure-math |
| Banach spaces for the Schwartz distributions | — | 0 / ? | ⬜ | Pure-math |
| Constructive Analysis in Infinitely many variables | — | 0 / ? | ⬜ | Pure-math |
| Note on the Spectral Theorem | — | 0 / ? | ⬜ | Pure-math |
| Some Banach spaces are almost Hilbert | — | 0 / ? | ⬜ | Pure-math |
| The Jones Strong Distribution Banach Spaces | — | 0 / ? | ⬜ | Pure-math |
| The S-basis and M-basis Problems for Separable Banach Spaces | — | 0 / ? | ⬜ | Pure-math |

## Methodology

### Per-equation template

````markdown
### Eq. (X.Y) — {short name}

**As printed:**
{LaTeX exactly as in the paper}

**Context / claim:**
{1–3 sentences on what the paper is asserting}

**Mathematica check:**
```mathematica
(* commands that reproduce or refute the claim *)
```
Result: {Simplify output / True / counterexample}

**Expanded derivation (grad-student level):**
- Step 1. …
- Step 2. …
- …

**Standard-equation comparison:**
In the limit {…} this reduces to {reference, e.g. Jackson Eq. 11.140}:
  {standard equation}
Identification: {variable mapping}.

**Verdict:** ✅ / ⚠ / ❌
````

### Wolfram Mathematica MCP

**Status:** license being procured by repo owner. Until the MCP server is wired in, every `Mathematica check` block is committed in source-form (the commands that *will* be run) marked with `(* PENDING: Wolfram MCP not yet configured *)`. Once the MCP is online, the campaign back-fills each result.

**Candidate MCP servers** (decision pending):

- Official Wolfram MCP (preferred if released)
- A `wolframscript`-backed shim (community projects exist; cf. `mathematica-mcp`-style wrappers)

**Companion notebooks:** each paper's checks will also be committed as a `.wl` (or `.nb`) under [../Mathematica_Notebooks/](../Mathematica_Notebooks/) so a human can re-run verification independent of the MCP.

### Scope cut on pure-math papers

The Banach-space and Navier–Stokes papers contain heavy functional analysis that doesn't all feed downstream. Default policy:

- **Physics papers** (Maxwell, Dirac/QED, dual relativistic QM, thermodynamics, classical electron, square-root operator, functional-integral foundations): **full coverage** — every numbered equation.
- **Pure-math papers**: **load-bearing results only** — lemmas/theorems whose conclusions are imported by a physics paper. Verify the remainder on demand.

### Audience note

Derivations are written for a **first-year physics graduate student**, not for the repo owner (who is a co-author of *Dual Relativistic Quantum Mechanics I*). Expect explicit algebra and references to standard texts.

## Known issues in existing equation-error catalog

[`../Equation_Errors_Dual_Theory_of_Relativity_and_Quantum_Mechanics.md`](../Equation_Errors_Dual_Theory_of_Relativity_and_Quantum_Mechanics.md) flags several "errors" that are actually **intentional dual-theory modifications** — specifically the appearance of the collaborative speed `b` (or `b/c` factors) where standard Maxwell has `c`. Replacing `b → c` in those equations does not "correct" them; it un-does the dual formulation that is the entire subject of Gill & Zachary (2011). The verification campaign will explicitly classify these as ✅ (correct in the dual framework) rather than ❌.

Likely affected entries to revisit:
- Error 4 — "J' transformation (2.13)" : the `(b/c)ρv` factor is intentional.
- Error 5 — "Four-vector source (3.12)" : `Jα = (Jx, Jy, Jz, ibρ)` is intentional under the dual theory's `c → b` substitution.
- Error 6 — "Vector potential wave equation (4.1)" : a notation/clarity issue at most, not an error.
- Error 7 — "Gradient operator (4.13)" : likely a typesetting artifact from PDF → Markdown, not a paper error.
