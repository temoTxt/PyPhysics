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
├── Equation_Verification/                  # ⭐ THE HEART OF THE REPO (issue #3)
│   ├── README.md                           # status table + per-equation template
│   ├── FINDINGS_for_author_review.md       # ⚠ 3 substantive findings flagged for authors
│   └── <paper>.md                          # one verification doc per physics paper
├── Mathematica_Notebooks/                  # .wl files runnable independent of MCP
├── Animations/                             # Manim sub-project (see above)
│   ├── manim_scenes/<paper>_<eq>_<topic>.py
│   └── rendered/videos/<scene>/480p15/     # 480p previews committed (~400-500 KB each)
│                              /1080p60/    # gitignored (rebuild via -qh)
├── Historical_Papers/{Primary,Retrospective}/        # PDFs for History campaign (issue #7)
│                                                     # gitignored by default; public-domain
│                                                     # files allow-listed via `git add -f`
├── Historical_Converted_Markdown/{Primary,Retrospective}/   # marker-pdf output; committed
└── History/                                # History-of-physics multi-chapter project (issue #7)
    ├── README.md, PLAN.md                  # methodology + 12-PR roadmap
    ├── _template_chapter.md                # historical-chapter scaffold
    ├── 01_…md through 06_…md               # historical + synthesis chapters
    ├── Forward/                            # forward-looking chapters (Ch 7-9)
    │   ├── _template_forward_open_questions.md     # Variant A (speculative)
    │   ├── _template_forward_derivational.md       # Variant B (derivational)
    │   └── 07_…md, 08_…md, 09_…md
    ├── Bibliography/                       # YAML-frontmatter cite cards
    │   ├── README.md                       # cite-key + tag conventions
    │   ├── {Primary,Retrospective}/<cite_key>.md
    │   ├── build_bibtex.py                 # YAML → bibliography.bib (gitignored)
    │   ├── scaffold_bib_note.py            # generate stub from --cite-key (+optional --from-doi)
    │   └── update_acquisition_tracker.py   # regenerate Historical_Papers/Acquisition_Tracker.md
    ├── Podcast/                            # 3-voice dialogue scripts per chapter
    │   ├── README.md                       # persona cast (Historian/Physicist/Experimentalist)
    │   └── episode_NN_*.md
    └── _tools/                             # per-chapter helpers (added in PR C)

Hydrogen_Visualizer/                        # standalone scratch script (plotly); not part of pipeline
main.py                                     # placeholder ("Hello from pyphysics!")
```

### Obsidian vault conventions

The repo doubles as an **Obsidian vault** (`.obsidian/` is checked in; `workspace*.json`/`cache/` are ignored). Cross-references between verification docs, bibliography notes, and chapters use Obsidian wikilinks (`[[name]]`, optionally `[[name#anchor]]`) rather than relative URLs. Cite-keys follow `firstauthor + year + slug` snake-case (e.g., `[[maxwell1865_dynamical_theory]]`). Full convention spec: [`Roadmapping/History/Bibliography/README.md`](Roadmapping/History/Bibliography/README.md).

### PDF acquisition + storage policy (History campaign)

Two-tier policy keyed on each bibliography note's YAML `pdf_status`:
- `out_of_copyright_public` — PDF committed to `Historical_Papers/<Primary|Retrospective>/` via `git add -f`.
- `acquired` — in-copyright PDF kept local only; only the marker-pdf-converted markdown is committed (fair-use academic quotation), under `Historical_Converted_Markdown/`.
- `pending` / `unavailable` — not yet sourced.

`Historical_Papers/` is gitignored by default (see [.gitignore](.gitignore)); each public-domain PDF must be force-added. Run `uv run python Roadmapping/History/Bibliography/update_acquisition_tracker.py` after any batch of bib-note edits to refresh `Historical_Papers/Acquisition_Tracker.md`.

### Tooling locations

- **Root project tools** (Python 3.14): `Roadmapping/download_paper.py`, `Roadmapping/parse_papers.py` (now accepts `--input` / `--output` / `--skip-existing`).
- **Bibliography tools**: `Roadmapping/History/Bibliography/{build_bibtex,scaffold_bib_note,update_acquisition_tracker}.py`.
- **Per-chapter helpers** (added in PR C): `Roadmapping/History/_tools/{fetch_pdf,validate_wikilinks,qa_converted_markdown,chapter_status}.py`.
- **Podcast tools** (introduced piecewise): `Roadmapping/History/Podcast/lint_episode.py` (PR C), optional `build_audio.py` / `build_episode_video.py` (PR L).
- **Manim scenes** (Python 3.13 sub-project): `Roadmapping/Animations/manim_scenes/*.py`.

All scripts follow `uv run python <path> [--dry-run] …` and ship with a top-of-file docstring + argparse. Don't add top-level dependencies without noting it in the PR.

## Crocco et al. (2026) AI-use compliance (load-bearing)

The repo's AI usage is governed by [Crocco, Rasdi, Garavan (2026), *Responsible AI in Non-Empirical Research*](https://doi.org/10.1177/15344843261445531). Full mapping: [`Roadmapping/Tooling/CROCCO_COMPLIANCE.md`](Roadmapping/Tooling/CROCCO_COMPLIANCE.md). Short rules every Claude Code session should follow:

1. **Pragmatic vs substantive.** Tag every AI use as one or the other. *Pragmatic* = grammar, formatting, translation, reference management, validator runs — doesn't shape what the manuscript argues. *Substantive* = generating themes, proposing frameworks, drafting prose that conveys core arguments. Substantive uses require fuller disclosure (the prompt + the model + a human-acceptance section in any synth report).
2. **Never invent citations.** Every `cite_key` in any committed chapter must correspond to an existing file under `Roadmapping/History/Bibliography/{Primary,Retrospective}/`. If you reference a paper that isn't there yet, scaffold the bib stub first via `scaffold_bib_note.py`. AI-fabricated citations are the most common Crocco-flagged failure mode.
3. **`human_reviewed` is binary and honest.** The YAML field flips to `true` only after a human has *read the primary source*. Don't infer "true" from a hand-written abstract; the body summary must reflect a reading, not an AI paraphrase.
4. **Substantive prompts live in `Roadmapping/Tooling/_prompts/`.** When invoking AI in a substantive role (Phase 4 synth tools, future literature-review work), use a committed prompt from there — don't compose prompts inline in synth tools. The prompt-of-record must be version-controlled.
5. **Synth reports leave a "human acceptance" section blank for the human.** When Claude produces a cluster, cross-ref suggestion, or claim-comparison report under `Roadmapping/Tooling/synth_reports/`, it includes that section as `<!-- TODO: human reviews and fills in -->`. The synth report does not commit chapter changes; only the human's subsequent PR does.
6. **AI is never an author.** Use `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>` trailers only. Don't list Claude in `Author:` fields, PR `--author` flags, or chapter byline attributions.

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
