# Author-review packet — Bethe–Salpeter campaign + r_e/r_0 findings

A PDF for **Tepper Gill** (DRQM I senior author) to fine-tooth-comb review. Implements issue [#93](https://github.com/temoTxt/PyPhysics/issues/93).

This is an **author-review packet, not a journal submission.** Per [CLAUDE.md](../../../CLAUDE.md), findings concerning DRQM I are "for author review" rather than corrections. The Branch A vs Branch C tension surfaced in the packet is content for Tepper to adjudicate, not a position the packet takes.

## What it contains

- The 28-result Bethe–Salpeter campaign inventory, with the campaign-closing honest-scope statement preserved verbatim from [`10_CrossComparison.md`](../../Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md) §1.
- The Z=1 triangulated `r_e/r_0 = 0.4994205099128317` and the §III.D closed-form characterisation per [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md).
- The Z-scan refinement and the Branch A vs Branch C tension (from open PRs [#84](https://github.com/temoTxt/PyPhysics/pull/84), [#85](https://github.com/temoTxt/PyPhysics/pull/85), [#86](https://github.com/temoTxt/PyPhysics/pull/86), [#87](https://github.com/temoTxt/PyPhysics/pull/87)).
- Nine Wolfram-MCP `.wl` notebooks translated to per-cell LaTeX under `equations/`, each with its source line and the Wolfram-MCP confirmation preserved.
- The mis-provenance audit log from PR bodies #84/#85/#86.
- A Crocco §5 AI-use disclosure section tied to the prompt-of-record under `prompts/`.
- A numbered closing question list for Tepper.

## Build

```bash
# From the repo root.
cd Roadmapping/Manuscripts/author_review_packet
make
# -> paper.pdf
```

Requires: `pdflatex`, `bibtex` (TeX Live; system packages already required for the Manim animations pipeline per CLAUDE.md).

The `references.bib` is built from the YAML stubs under `Roadmapping/History/Bibliography/{Primary,Retrospective}/` via `build_bibtex.py`; the Makefile invokes it automatically.

## Layout

```
author_review_packet/
├── README.md          # this file
├── Makefile           # pdflatex + bibtex (latexmk not required)
├── paper.tex          # main document
├── paper.pdf          # built artifact (committed)
├── references.bib     # generated from YAML stubs
├── sections/          # 11 input files, one per packet section
├── equations/         # 9 LaTeX renderings of the load-bearing .wl notebooks
├── figures/           # cross-comparison table, Z-scan curve, optional Manim stills
└── prompts/           # Crocco-disclosed substantive AI prompts
```

## Writing standard

This packet is governed by the same writing standard that governs the
Author Reports, codified in
[`Roadmapping/Tooling/VOICE_MATCH_GILL.md`](../../Tooling/VOICE_MATCH_GILL.md):

- **§3.bis** — every load-bearing equation is a numbered display equation.
  LaTeX `equation` / `align` environments give numbering by default.
- **§3.ter** — short complete sentences; one claim per sentence; the
  mathematics does most of the talking.
- **§3.quater** — numbered footnotes plus a bibliography section; no inline
  links scattered through paragraphs.
- **§3.quinquies** — no emojis or pseudo-emoji symbols. Use bold text
  labels (**Pass** / **Marginal** / **Blocked** / **Refuted**) instead.

Status notes (verdict markers in tables, derivation-chain confirmation
markers) use bold text labels rather than verdict emoji. Findings recorded
elsewhere in the repository may still use the emoji vocabulary; the
packet's transcriptions preserve content but normalise the rendering.

Further prose tightening against §3.ter (cutting paragraphs whose content
the equations already carry) is tracked as follow-up work and is
deliberately left out of this packet revision, which focuses on the
emoji-removal pass and the explicit standard reference.

## Honest framing — do not soften

- The campaign-closing statement *"Zero of 28 results constitute independent corroboration of the dual theory's content distinct from textbook QM"* is preserved verbatim. Any edit that softens this leaves the honest scope.
- The Z-scan empirical finding (PR #87) and the Branch A "by construction" reading (PRs #84/#85/#86) are presented in parallel without endorsement; the packet asks Tepper to adjudicate.
- Finding 2's verdict trajectory `Refuted → Marginal-characterised → Marginal/Pass-conditional` (per `FINDINGS_for_author_review.md:209`) is transcribed verbatim, with [issue #75](https://github.com/temoTxt/PyPhysics/issues/75) referenced as the unconditional-Pass hypothesis-(i) path.
