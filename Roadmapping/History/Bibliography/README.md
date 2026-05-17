# Bibliography

Cite-source-of-truth for the History campaign. Every paper cited in any chapter, podcast episode, or animation traces back to a single YAML-frontmatter markdown file in this tree.

## Cite-key convention

`firstauthor + year + slug`, lowercase, snake-case, no diacritics. The slug is 1–3 words from the title.

Examples:
- `maxwell1865_dynamical_theory`
- `lamb1947_fine_structure`
- `bcs1957_superconductivity` (BCS as a conventional alias; bibtex `bcs1957`)
- `schweber1994_qed_and_men`
- `verrier1859_perihelie_mercure`

Cite-keys must be globally unique within the `Bibliography/` tree. Wikilinks use the bare cite-key: `[[maxwell1865_dynamical_theory]]`.

## YAML frontmatter schema

Every bibliography note opens with this block. Fields not yet known stay as the empty string or `pending`; never delete a field, just leave it blank.

```yaml
---
cite_key: maxwell1865_dynamical_theory
title: "A Dynamical Theory of the Electromagnetic Field"
authors: ["James Clerk Maxwell"]
year: 1865
type: primary                                # primary | retrospective
era: "1860-1900"                             # one of: 1800-1860, 1860-1900, 1900-1925, 1925-1948, 1948-1965, forward
tags: [electromagnetism, maxwell-equations, foundational]
journal: "Philosophical Transactions of the Royal Society"
volume: 155
issue: ""
pages: "459-512"
doi: "10.1098/rstl.1865.0008"
url: ""
arxiv_id: ""
pdf_status: out_of_copyright_public          # acquired | pending | unavailable | out_of_copyright_public
pdf_path: "../../Historical_Papers/Primary/maxwell1865_dynamical_theory.pdf"
converted_md: "../../Historical_Converted_Markdown/Primary/maxwell1865_dynamical_theory/maxwell1865_dynamical_theory.md"
pages_quoted: ""                             # populated only for books we excerpt under fair use
gill_corpus_overlap: ["Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations"]
chapters_citing: []                          # filled by chapters that wikilink this note
---
```

After the frontmatter, the body is freeform — typically a 2–3-paragraph summary written during the chapter's "QA + retrospectives" step, plus any notes the author wants to keep.

## Tag taxonomy

Three orthogonal axes; tags are namespace-scoped so an Obsidian Tag Pane query can filter by any combination.

**By framework claim type** (confidence tiers, strongest → weakest):

- `#verified` — cross-references a Wolfram-confirmed verification doc under [../Equation_Verification/](../../Equation_Verification/).
- `#inferred` — direct/mechanical extrapolation of a published Gill equation to a topic he didn't cover.
- `#speculative` — beyond mechanical extrapolation; framework *might* bear, derivation deferred.
- `#gill-silent` — no claim about proper-time relevance.

**By thread**: `#thread/electromagnetism`, `#thread/quantum`, `#thread/solid-state`.

**By era**: `#era/1800-1860`, `#era/1860-1900`, `#era/1900-1925`, `#era/1925-1948`, `#era/1948-1965`, `#era/forward`.

**By topic** (free-form): `#hyperfine-splitting`, `#renormalization`, `#superconductivity`, etc.

## Layout

```
Bibliography/
├── README.md                       # this file
├── Primary/                        # contemporary sources from the era under study
│   └── <cite_key>.md
└── Retrospective/                  # later historical / pedagogical reviews
    └── <cite_key>.md
```

## Tooling

All scripts run from the repo root via `uv run python …`.

- `build_bibtex.py` — walks every `Bibliography/**/*.md`; emits `bibliography.bib` (gitignored — YAML is canonical). Idempotent.
- `scaffold_bib_note.py --cite-key <key>` — generates a skeleton bib note with the full YAML schema. Optional `--from-doi <doi>` auto-fills `title`, `authors`, `year`, `journal`, `volume`, `pages` via the Crossref REST API.
- `update_acquisition_tracker.py` — reads every bib note's `pdf_status` + `pdf_path` and regenerates `../../Historical_Papers/Acquisition_Tracker.md`.

Per-chapter helpers (introduced in PR C) live in [`../_tools/`](../_tools/): `fetch_pdf.py`, `validate_wikilinks.py`, `qa_converted_markdown.py`, `chapter_status.py`.

## Workflow for adding a bibliography entry

1. `uv run python Roadmapping/History/Bibliography/scaffold_bib_note.py --cite-key foo1880_slug --from-doi 10.…` — generates the stub.
2. Hand-edit any fields the Crossref lookup missed (tags, gill_corpus_overlap, era).
3. If the PDF is freely available: `uv run python Roadmapping/History/_tools/fetch_pdf.py --cite-key foo1880_slug --url <…>` (PR C tool). Otherwise leave `pdf_status: pending` and acquire manually.
4. `uv run python Roadmapping/History/Bibliography/update_acquisition_tracker.py` — refresh the tracker.
5. Write the 2–3-paragraph body during the chapter's QA pass (step 5 of the per-chapter pattern).
