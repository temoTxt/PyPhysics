# History of Physics 1800–1965 — Multi-Step Plan

**Tracks:** GitHub issue #7
**Status:** plan-only (this document). No execution PRs opened yet.

A multi-step plan to recreate the history of physics from 1800 to 1965 — through the birth of particle physics — comparing the **standard development** with **Tepper Gill's proper-time / dual-theory framework**. Five chronological chapters covering three parallel threads (electromagnetism, quantum mechanics, solid-state / superconductivity / transistors), each ending in 1965 at the cusp of the standard model and Moore's law.

## Scoping decisions (already settled)

| Question | Decision |
|---|---|
| Primary output format | Markdown chapters + companion Manim scenes (mirrors verification campaign #3 + animations campaign #5). |
| What when Gill has no published paper covering a topic | Extrapolate cautiously; mark every extrapolation `#inferred`. |
| Solid-state / transistor track | Include with speculative proper-time commentary, also flagged `#inferred`. |
| Existing verification + animations docs | Retroactively converted to Obsidian wikilinks in a small preparatory PR (PR A below). |
| BibTeX bridge | Built in Phase 0 (PR B), so chapters can cite via wikilinks from day one. |
| PDF storage | `Historical_Papers/` gitignored by default. Allow-list small public-domain PDFs via `git add -f`. In-copyright PDFs are never committed; only the converted markdown is, under fair-use academic quotation. |

## Output directory layout

```
.obsidian/                                       # vault config (PR A)
Roadmapping/
├── Tepper_Gill_Papers/                          # existing
├── Converted_Markdown/                          # existing
├── Equation_Verification/                       # existing — wikilinked in PR A
├── Animations/                                  # existing — wikilinked in PR A
├── Mathematica_Notebooks/                       # existing
├── Historical_Papers/                           # NEW (PR B) — gitignored allow-list
│   ├── Primary/
│   │   ├── maxwell1865_dynamical_theory.pdf     # public-domain → committed
│   │   ├── michelson1887_relative_motion.pdf    # public-domain → committed
│   │   └── ...
│   ├── Retrospective/
│   │   ├── pais1982_subtle_is_the_lord.pdf      # in-copyright → NOT committed
│   │   └── ...
│   └── Acquisition_Tracker.md                   # NEW — master status table
├── Historical_Converted_Markdown/               # NEW (PR B) — output of parse_papers.py
│   ├── Primary/
│   │   └── maxwell1865_dynamical_theory/
│   │       ├── maxwell1865_dynamical_theory.md
│   │       └── images/
│   └── Retrospective/...
└── History/                                     # NEW
    ├── README.md                                # methodology + Obsidian conventions
    ├── PLAN.md                                  # this file
    ├── _template_chapter.md                     # copyable scaffold
    ├── 01_early_electromagnetism_1800_1860.md
    ├── 02_classical_synthesis_1860_1900.md
    ├── 03_old_quantum_theory_1900_1925.md
    ├── 04_quantum_mechanics_1925_1948.md
    ├── 05_QED_renormalization_solid_state_1948_1965.md
    ├── 06_synthesis_divergence_map.md
    └── Bibliography/
        ├── README.md                            # cite-key conventions, YAML schema, tags
        ├── build_bibtex.py                      # YAML → bibliography.bib generator
        ├── bibliography.bib                     # generated; gitignored
        ├── Primary/
        │   └── *.md                             # citation cards
        └── Retrospective/
            └── *.md
```

## Five eras + divergence axis

| Ch | Era | EM / QM thread | Solid-state thread | Where proper-time diverges |
|---|---|---|---|---|
| 1 | **1800–1860** | Volta → Ørsted → Ampère → Faraday | (none) | `#inferred`: quasi-static sources, `b ≈ c` throughout. **Null comparison.** |
| 2 | **1860–1900** | Maxwell (1865) → Hertz (1888) → Michelson-Morley (1887) → Thomson electron (1897) | Hall effect (1879), photoconductivity (1873), Braun rectification (1874) | `#verified`: Eqs (3)→(3′) verified in PR #4. **Michelson-Morley** is the natural launch point for the dual-theory contrast. |
| 3 | **1900–1925** | Planck → photoelectric → Bohr → Compton → de Broglie → matrix/wave mechanics | Onnes He liquefaction (1908) → superconductivity (1911) | `#inferred`: proper-time Bohr levels (K replaces H); Compton on free electron unchanged (μ_photon=0 for inertial sources). |
| 4 | **1925–1948** | Schrödinger → Klein-Gordon → Dirac → Anderson positron | Bloch (1928), band theory (Wilson 1931), Meissner (1933) | `#verified`: DRQM I verified in PR #4. **Klein-Gordon's negative-probability problem dissolved** in proper-time. Antimatter via Santilli isoduals. |
| 5 | **1948–1965** | Lamb shift → Schwinger (1947) → Feynman/Tomonaga (1948–50) → Dyson → parity violation (Wu 1957) | Transistor (1947), BCS (1957), Josephson (1962) | `#verified`: DRQM I **Finding 2** — the **headline payoff** of the whole project. `#inferred`: Dyson conjecture dissolved by Gill's time-ordered operator calculus. |
| 6 | **Synthesis** | Divergence map across all eras; reuses + extends [`synthesis_tour.py`](../Animations/manim_scenes/synthesis_tour.py). |

## Obsidian + citation conventions

### Cite-key format

`firstauthor + year + slug`, lowercase, snake-case:
- `[[maxwell1865_dynamical_theory]]`
- `[[lamb1947_fine_structure]]`
- `[[bcs1957_superconductivity]]` (BCS as conventional alias)
- `[[schweber1994_qed_and_men]]`

### Bibliography note YAML schema

```yaml
---
title: "A Dynamical Theory of the Electromagnetic Field"
authors: ["James Clerk Maxwell"]
year: 1865
type: primary                    # primary | retrospective
era: "1860-1900"
tags: [electromagnetism, maxwell-equations, foundational]
journal: "Philosophical Transactions of the Royal Society"
volume: 155
pages: "459-512"
doi: "10.1098/rstl.1865.0008"
pdf_status: out_of_copyright_public   # acquired | pending | unavailable | out_of_copyright_public
pdf_path: "../../../Historical_Papers/Primary/maxwell1865_dynamical_theory.pdf"
converted_md: "../../../Historical_Converted_Markdown/Primary/maxwell1865_dynamical_theory/maxwell1865_dynamical_theory.md"
gill_corpus_overlap: ["Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations"]
---
```

### Tag taxonomy

- **By framework claim type**: `#verified` (cross-ref to a verification doc), `#inferred` (extrapolation), `#gill-silent` (standard physics only).
- **By thread**: `#thread/electromagnetism`, `#thread/quantum`, `#thread/solid-state`.
- **By era**: `#era/1800-1860`, `#era/1860-1900`, `#era/1900-1925`, `#era/1925-1948`, `#era/1948-1965`.
- **By topic**: free-form, e.g. `#hyperfine-splitting`, `#renormalization`, `#superconductivity`.

### Chapter template structure

Each chapter has these sections (full template in `_template_chapter.md`):
1. **Overview** — 1-paragraph framing of the era.
2. **Historical narrative** — standard physics, key experiments highlighted; wikilinks to primary sources inline.
3. **Proper-time commentary** — section per thread (EM / QM / solid-state); each claim tagged `#verified`, `#inferred`, or `#gill-silent`.
4. **Key derivations worth animating** — list of which animations exist (cross-ref) or are proposed for this chapter.
5. **Primary sources cited** — wikilink list to `Bibliography/Primary/*.md`.
6. **Retrospective reviews drawn on** — wikilink list to `Bibliography/Retrospective/*.md`.

## Eight-PR execution plan

### PR A — Obsidian retrofit *(prep, small)*

- Add `.obsidian/` vault config; enable plugins: **Tag Pane**, **Backlinks**, **Outgoing Links**, **Dataview**, **Templates**.
- `.gitignore` for `.obsidian/`: commit `app.json`, `core-plugins.json`, `community-plugins.json`, `hotkeys.json`, `templates.json`; ignore `workspace.json`, `workspace-mobile.json`, `cache/`.
- Convert relative-URL cross-references in:
  - `Roadmapping/Equation_Verification/*.md` (11 files)
  - `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md`
  - `Roadmapping/Animations/README.md`
- Manual spot-check of the graph view (Trey verifies that the existing 12 docs form an interconnected cluster).
- **Out of scope:** rewriting any content; just link syntax.

### PR B — Phase 0 infrastructure

- `Roadmapping/History/README.md` — methodology + Obsidian conventions + comparison template.
- `Roadmapping/History/_template_chapter.md` — copyable scaffold.
- `Roadmapping/History/Bibliography/README.md` — cite-key conventions, YAML schema, tag taxonomy.
- `Roadmapping/History/Bibliography/build_bibtex.py` — small `uv run python` script that reads YAML frontmatter from every `Bibliography/**/*.md` and emits `bibliography.bib`. Generated file gitignored.
- Empty chapter files (`01_…md` through `06_…md`) with template + section headings.
- **Refactor `Roadmapping/parse_papers.py`** to accept `--input` / `--output` / `--skip-existing` CLI args. Existing behavior preserved when called with defaults.
- Create `Roadmapping/Historical_Papers/Primary/` + `Retrospective/` + `Roadmapping/Historical_Converted_Markdown/Primary/` + `Retrospective/` (with `.gitkeep`).
- Add `Roadmapping/Historical_Papers/Acquisition_Tracker.md` skeleton listing all ~95 cite-keys with `pdf_status: pending`.
- Add `Historical_Papers/` to `.gitignore` with allow-list pattern for `*.pdf` matching public-domain entries (revisited per-chapter).
- Update root `CLAUDE.md` with one paragraph on Obsidian conventions + PDF-acquisition workflow.

### PR C–G — Phases 1–5 (one chapter per PR)

Each chapter PR follows the same **five-step** pattern:

1. **Bibliography stubs** — YAML frontmatter only, for ~10–24 era sources.
2. **PDF acquisition pass** — fetch what's freely available (Royal Society archive, Wikisource, arxiv, archive.org); update `Acquisition_Tracker.md`; commit only public-domain PDFs (allow-listed).
3. **Conversion pass** — `uv run python Roadmapping/parse_papers.py --input Roadmapping/Historical_Papers/Primary --output Roadmapping/Historical_Converted_Markdown/Primary`. Commit converted markdown + images.
4. **Narrative pass** — write the chapter. Inline wikilinks resolve to the now-filled bibliography notes; quote directly from the converted markdown where useful.
5. **QA + retrospective summaries + animation** — eyeball converted equations for OCR errors and fix the load-bearing ones; fill in 2–3-paragraph summaries on retrospective notes; render any new Manim scenes.

Per-chapter sizing:

| PR | Phase | Chapter | New bib notes | New Manim scenes | Est. lines (excl. converted PDFs) |
|---|---|---|---|---|---|
| C | 1 | 1: 1800–1860 | ~13 | 1 | ~600 |
| D | 2 | 2: 1860–1900 (load-bearing) | ~16 | 2 | ~1100 |
| E | 3 | 3: 1900–1925 | ~20 | 2 | ~1000 |
| F | 4 | 4: 1925–1948 | ~17 | 2 | ~1100 |
| G | 5 | 5: 1948–1965 (headline) | ~24 | 2 | ~1300 |
| H | 6 | 6: Synthesis + indexes | 0 | 1 (extends `synthesis_tour`) | ~500 |

**Totals:** ~95 bibliography notes, ~10 new Manim scenes, 6 chapters, ~5,600 lines of narrative + bibliography (PDFs and converted markdown not counted).

### PR H — Phase 6 synthesis

- `06_synthesis_divergence_map.md` — closing chapter: where standard and proper-time agree, where they diverge sharply.
- Dataview query pages: `_index_by_year.md` (all bibliography by year), `_index_by_tag.md` (all by tag), `_index_inferred_claims.md` (every `#inferred` extrapolation across the project, for periodic review).
- Extend `Roadmapping/Animations/manim_scenes/synthesis_tour.py` to span the full 165-year arc, or add `hist_full_synthesis.py` as a complement.
- Optional: a small Python script that emits a citation-network graph PNG from the wikilinks (for embedding in the synthesis chapter).

## Per-chapter scene proposals

| Chapter | Proposed new scenes |
|---|---|
| 1 | `hist_faraday_induction.py` — field lines through a moving loop; induced EMF. |
| 2 | `hist_michelson_morley.py` — predicted fringe shift vs measured null. `hist_maxwell_synthesis.py` — Ampère + Faraday → wave equation → `c` emerges. |
| 3 | `hist_bohr_proper_time.py` — Bohr levels rederived with `K = H²/(2mc²) + mc²/2`. `hist_compton_null.py` — Gill's framework recovers standard Compton for inertial electrons. |
| 4 | `hist_klein_gordon_vs_dual.py` — negative-probability problem and how the dual K dissolves it. `hist_positron_isodual.py` — Anderson + Santilli isodual reinterpretation. |
| 5 | `hist_lamb_shift_contrast.py`. Reuses [`drqm_eq22_g_factor_finding.py`](../Animations/manim_scenes/drqm_eq22_g_factor_finding.py) as the chapter payoff. |
| 6 | Extends [`synthesis_tour.py`](../Animations/manim_scenes/synthesis_tour.py). |

## OCR-quality + copyright caveats

1. **Modern PDFs (post-1980, arxiv-style)** convert cleanly via `marker-pdf`.
2. **Mid-century journal scans (1947 Lamb-Retherford, 1957 BCS)** convert OK with some equation cleanup needed. Plan a per-chapter QA pass.
3. **Pre-1900 scans (Maxwell 1865, Hertz 1893)** require significant manual cleanup of equations and figures. Budget extra QA time.
4. **Copyright reality** — two-tier policy:
   - `pdf_status: out_of_copyright_public` → commit the PDF.
   - `pdf_status: acquired` (in copyright) → DO NOT commit the PDF. Convert locally, commit *only* the markdown if licensing allows (likely fair-use for academic quotation), and note source in YAML.
   - For books we cite but can't redistribute, convert and commit only the chapters relevant to our citations, with an explicit `pages_quoted` field in YAML.

## Cross-referencing the verification campaign

The 11 papers verified in PR #4 are the "anchor texts" for proper-time commentary. Each chapter cross-references them via wikilinks (post-PR-A):

| Chapter | Verification docs cross-referenced |
|---|---|
| 2 | `[[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]` |
| 3 | `[[FoundationsII-Classical]]`, `[[Analytic_Representation_of_The_Square-Root_Operator]]` |
| 4 | `[[Dual_Relativistic_Quantum_Mechanics_I]]`, `[[Analytic_Representation_of_The_Dirac_Equation]]`, `[[FoundationsII-Classical]]` |
| 5 | `[[Dual_Relativistic_Quantum_Mechanics_I]]`, `[[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]]`, `[[Feynman_Operator_Calculus_Papers]]`, `[[Relativistic_Transformations_of_Thermodynamics]]`, `[[FINDINGS_for_author_review]]` |

## Open items left for execution time

- Final per-chapter bibliography lists (specific cite-keys with paper titles) — drafted during the bibliography-stub step of each chapter PR.
- Specific PDF source URLs for each cite-key — populated into `Acquisition_Tracker.md` as PDFs are sourced.
- The exact YAML schema for the `gill_corpus_overlap` field — finalized in PR B once the Obsidian retrofit (PR A) sets the canonical filenames.
- Whether to LaTeX-typeset a final bound document (e.g. via pandoc) at the end — deferred; not in this plan.
