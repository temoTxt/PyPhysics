# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A **physics research repository** centered on the Gill *dual theory of relativity and quantum mechanics*. The bulk of the repo is **Markdown verification docs** and **Manim animation scenes**, not application code. There are no tests and no lint configuration — output quality is measured by whether (a) the algebra reproduces in Wolfram Mathematica and (b) the rendered animation looks right.

The repo owner (`temoTxt` / Trey Morris) is a co-author of *Dual Relativistic Quantum Mechanics I* (2021). When findings concern that paper, frame them as "for author review" rather than corrections.

## Two Python environments — both `uv`-managed, different Python versions

Always invoke Python via `uv run …` — never bare `python`/`pip`/`python3`.

| Project | Path | Python | Purpose |
|---|---|---|---|
| `pyphysics` (root) | `/` | **3.14** (`requires-python = ">=3.14"`) | Paper-download (`arxiv`) and PDF→Markdown conversion (`marker-pdf`/`torch`) |
| `pyphysics-animations` | `Roadmapping/Animations/` | **3.13** (`>=3.13,<3.14`) | Manim CE scenes. **Intentionally NOT a workspace member** of the root project — Manim's PyAV transitive lacks Python-3.14 wheels |

When working in `Roadmapping/Animations/`, always pass `--project Roadmapping/Animations` to `uv` (or `cd` first). Running `uv init` inside that subdir from earlier versions of `uv` has been observed to overwrite the *root* `pyproject.toml`; create sub-project pyprojects by hand if needed.

## Common commands

Paper pipeline (root project):
```bash
uv run python Roadmapping/download_paper.py        # arxiv → Roadmapping/Tepper_Gill_Papers/*.pdf
uv run python Roadmapping/parse_papers.py          # PDFs → Roadmapping/Converted_Markdown/<paper>/<paper>.md (+ images/)
```
`parse_papers.py` uses `marker-pdf` in CPU mode (set `TORCH_DEVICE=cpu` if needed) and skips already-converted papers. It can take hours on a full corpus.

Manim animations (sub-project):
```bash
uv sync --project Roadmapping/Animations           # installs Manim 0.20+, numpy, into Roadmapping/Animations/.venv/
cd Roadmapping/Animations
uv run manim -ql --media_dir rendered manim_scenes/<file>.py <SceneClass>   # 480p preview (~10s render)
uv run manim -qh --media_dir rendered manim_scenes/<file>.py <SceneClass>   # 1080p — required to satisfy issue #5 acceptance
uv run manim -pqh ...                                                       # add -p to auto-play after render
```

System requirements for Manim: `ffmpeg`, `texlive-latex-base`+`texlive-latex-extra`+`texlive-pictures`, `libcairo2-dev`, `libpango1.0-dev`. Hardcoded as system deps (not in `pyproject.toml`).

## Repository structure

```
Roadmapping/
├── download_paper.py / parse_papers.py     # paper-corpus pipeline (root pyproject)
├── Tepper_Gill_Papers/                     # 22 PDFs from arxiv (gitignored if large; check)
├── Converted_Markdown/<paper>/<paper>.md   # marker-pdf output; line-numbers in these files
│                                           # are cited by the verification docs
├── Equation_Verification/                  # ⭐ THE HEART OF THE REPO
│   ├── README.md                           # status table + per-equation template
│   ├── FINDINGS_for_author_review.md       # ⚠ 3 substantive findings flagged for authors
│   └── <paper>.md                          # one verification doc per physics paper
├── Mathematica_Notebooks/                  # .wl files runnable independent of MCP
├── Animations/                             # Manim sub-project (see above)
│   ├── manim_scenes/<paper>_<eq>_<topic>.py
│   └── rendered/videos/<scene>/480p15/     # 480p previews committed (~400-500 KB each)
│                              /1080p60/    # gitignored (rebuild via -qh)
└── Research Roadmap …, Equation_Errors_…, Experimental_Design_Plan_… (.md notes)

Hydrogen_Visualizer/                        # standalone scratch script (plotly); not part of pipeline
main.py                                     # placeholder ("Hello from pyphysics!")
```

## How equation verification works

Each paper has a verification doc under `Roadmapping/Equation_Verification/` that follows the per-equation template (full template in `Roadmapping/Equation_Verification/README.md`):

```
### Eq. (X.Y) — short name
**As printed:**   (LaTeX exactly as in the paper, with line ref into Converted_Markdown)
**Mathematica check:** (single-line Wolfram Language; run via Wolfram MCP; record "Result: 0 ✅" inline)
**Expanded derivation (grad-student level):** step-by-step
**Standard-equation comparison:** (Jackson / Sakurai / Peskin / Goldstein / Weinberg)
**Verdict:** ✅ / ⚠ / ❌
```

Verification is performed via a **Wolfram MCP server** (`mcp-remote https://services.wolfram.com/api/mcp`). The Mathematica blocks in the docs are the canonical record of what was actually run.

**Three Wolfram MCP gotchas that have bitten previous runs** — bake into any new Mathematica check:
1. **Multi-line code is silently broken at line boundaries** by the MCP transport. Always put each `ClearAll`/`Print`/`FullSimplify` chain on a single line, joined with `;`.
2. **`V` is a reserved symbol** (Vanadium element) — uppercase `V` causes silent term-dropping. Use `potV` or similar.
3. **`e` may resolve to Euler's number** in some contexts — use `ee` for the electron charge.
4. Mathematica's `Dot` is **non-commutative** for symbolic 3-vectors. Either substitute scalar surrogates (e.g. `udota` for `u·a`) or apply rewrite rules like `a[t] . u[t] -> u[t] . a[t]`.

## Scope conventions

- **Physics papers:** full per-equation coverage (every numbered equation).
- **Pure-math papers** (Banach spaces, Navier–Stokes — the bottom of the status table in `Equation_Verification/README.md`): only load-bearing results referenced by a physics paper. The full math content is intentionally deferred.

Some entries in `Roadmapping/Equation_Errors_Dual_Theory_of_Relativity_and_Quantum_Mechanics.md` flag the `b`/`b/c` factors as "errors". They are **not** errors — they are the intentional dual-theory modifications. Don't propose those as corrections.

## Manim scene conventions

When adding a new scene under `Roadmapping/Animations/manim_scenes/`:
- One file per equation cluster, named `<paper>_eq<N>_<short>.py`.
- Docstring **must** include the markdown anchor link to the relevant section of the verification doc and the `uv run manim ...` command.
- Walk through the derivation step-by-step using `Transform(old, new)` between intermediate `MathTex` forms — this preserves the "this becomes that" narrative.
- Box the final result with `\boxed{...}` colored `YELLOW`, and cite the Wolfram MCP confirmation in a closing card.
- Render at `-ql` first to catch syntax errors, then `-qh` to confirm the 1080p acceptance criterion.
- Wrap LaTeX in raw strings (`r"..."`).

## Three open findings flagged for author review

These live in `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` and are *reproducible from the recorded MCP inputs*. If a future task touches related papers, treat these as known:

1. **Maxwell paper Eq. (24)** — missing factor of `c` in `eℏΣ·B/(2m)` (should be `/(2mc)`) and missing entire `+V²/(2mc²)` term.
2. **DRQM I Eq. (III.22)** — paper's stated `r_e = 0.499857150068631 r_0` gives `g = -2.0005714`, not the claimed experimental `-2.00231930436256`. Required `r_e ≈ 0.499420510 r_0`.
3. **TCEP Eq. (4.16)** — paper writes `v_g = v_g' - v`; algebra and paper's own commentary give `v_g = v_g' + v` (sign typo).

Open Q3 from the original Maxwell paper open-questions list (the `r_0` critical-point claim) was *resolved* during the campaign by cross-reference to `FoundationsII-Classical.md` Sec 2.2.1.
