# History of Physics 1800–1965

Working campaign tracking the history-of-physics multi-chapter project. Full plan: [[PLAN]]. Execution tracker: GitHub epic [#7](https://github.com/temoTxt/PyPhysics/issues/7) with sub-issues #9–#20.

## Framing principle (load-bearing)

> **We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics.**

Every chapter narrative, every Manim scene, every podcast episode opens (or closes) on an explicit version of this disclaimer. The dual-theory framework reproduces all experimentally confirmed predictions of standard SR + QM + QED within current measurement precision; the contrast is at the level of mathematical conventions, not predictions. Findings flagged for author review (see [[FINDINGS_for_author_review]]) are *internal-consistency* failures of cited papers — algebra that doesn't reproduce — not "Gill is wrong about physics".

## Layout

```
History/
├── README.md                                # this file
├── PLAN.md                                  # 12-PR roadmap
├── _template_chapter.md                     # historical-chapter scaffold (Ch 1–6)
├── 01_…md through 06_…md                    # historical chapters
├── Forward/
│   ├── _template_forward_open_questions.md  # Variant A (speculative)
│   ├── _template_forward_derivational.md    # Variant B (derivational)
│   ├── 07_PNT_GPS_SLR_QKD.md                # Variant B
│   ├── 08_quantum_computing_open_questions.md     # Variant A
│   └── 09_fusion_open_questions.md          # Variant A
├── Podcast/
│   ├── README.md                            # 3-voice cast conventions
│   ├── _template_episode.md
│   └── episode_NN_*.md                      # dialogue scripts per chapter
├── Bibliography/
│   ├── README.md                            # cite-key, YAML schema, tag taxonomy
│   ├── Primary/                             # contemporary primary sources
│   └── Retrospective/                       # secondary literature
└── _tools/                                  # per-chapter helper scripts (PR C+)
```

The PDF + converted-markdown tree lives one level up:

```
Roadmapping/
├── Historical_Papers/{Primary,Retrospective}/   # PDFs; gitignored-by-default with allow-list
└── Historical_Converted_Markdown/{Primary,Retrospective}/    # marker-pdf output; committed
```

## Methodology

### Per-chapter 6-step pattern

1. **Bibliography stubs** — populate Bibliography/Primary + Retrospective for the era (~13–24 sources). Use `uv run python Roadmapping/History/Bibliography/scaffold_bib_note.py --cite-key … --from-doi …`.
2. **PDF acquisition** — pull what's freely available. Public-domain PDFs are committed to `Historical_Papers/` via `git add -f`; in-copyright PDFs stay local.
3. **Conversion** — `uv run python Roadmapping/parse_papers.py --input Roadmapping/Historical_Papers/Primary --output Roadmapping/Historical_Converted_Markdown/Primary`.
4. **Narrative** — write the chapter `.md` per the template, opening with the framing-principle disclaimer.
5. **QA + retrospectives + animations** — eyeball converted equations for OCR errors; fill 2–3-paragraph summaries on retrospective notes; render the new Manim scenes.
6. **Podcast episode script** — write the dialogue version of the chapter as `Podcast/episode_NN_*.md`.

### Confidence tiers (tags)

- `#verified` — cross-references a Wolfram-confirmed verification doc.
- `#inferred` — direct/mechanical extrapolation of a published Gill equation.
- `#speculative` — beyond mechanical extrapolation; framework *might* bear, derivation deferred.
- `#gill-silent` — no claim about proper-time relevance.

### Cross-references

Use Obsidian `[[wikilinks]]` for in-vault citations. Bibliography entries follow the `firstauthor + year + slug` snake-case convention. The verification campaign's 11 Gill-corpus papers (PR #4) are the load-bearing anchors:

- [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]
- [[Dual_Relativistic_Quantum_Mechanics_I]]
- [[The_Classical_Electron_Problem]]
- [[Analytic_Representation_of_The_Dirac_Equation]]
- [[Analytic_Representation_of_The_Square-Root_Operator]]
- [[Relativistic_Transformations_of_Thermodynamics]]
- [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]]
- [[FoundationsII-Classical]]
- [[Feynman_Operator_Calculus_Papers]]
- [[FINDINGS_for_author_review]]

## Tooling

- `Roadmapping/parse_papers.py` — PDF → markdown via marker-pdf (CLI: `--input` / `--output` / `--skip-existing`).
- `Roadmapping/History/Bibliography/build_bibtex.py` — YAML frontmatter → `.bib`.
- `Roadmapping/History/Bibliography/scaffold_bib_note.py` — generate a YAML stub from `--cite-key`+ optional `--doi`.
- `Roadmapping/History/Bibliography/update_acquisition_tracker.py` — regenerate `Historical_Papers/Acquisition_Tracker.md` from bibliography frontmatter.

Per-chapter helpers land in PR C: `_tools/fetch_pdf.py`, `_tools/validate_wikilinks.py`, `_tools/qa_converted_markdown.py`, `_tools/chapter_status.py`, `Podcast/lint_episode.py`. Synthesis tools land in PR H. Optional audio pipeline lands in PR L (Phase 10).
