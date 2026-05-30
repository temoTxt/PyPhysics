# Task 93: Author-review PDF for Tepper Gill: BS campaign equations + r_e/r_0 findings (incl. Z-scan / Branch A vs C tension)

## Objective

Build a single PDF (`Roadmapping/Manuscripts/author_review_packet/paper.pdf`) that Tepper Gill — DRQM I senior author — can fine-tooth-comb review, because **he reads PDFs of equation derivations, not Mathematica code**. The packet translates the Wolfram-MCP-verified `.wl` notebooks into human-readable LaTeX derivations, consolidates the 28-result Bethe–Salpeter campaign, presents the Z=1 triangulated `r_e/r_0 = 0.4994205099128317` and the Z-scan refinement honestly, and closes with an explicit numbered question list for the author. Per [CLAUDE.md](../../CLAUDE.md) Crocco-compliance rules, the packet is **for author review**, not "corrections," and the Branch A vs Branch C tension is *content the PDF surfaces for Tepper to adjudicate*, not a position the PDF takes.

## Background

The repo already carries the substantive material; this task is a *transcription and consolidation* effort, not a new physics derivation. The honest-scope and verdict state of the underlying campaign is **richer than the issue body conveys** and must be preserved faithfully.

**Campaign source-of-truth on `main`:**

- The Bethe–Salpeter campaign closes in [Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md). Its §1 inventory (28 results across PRs A–I) and the campaign-closing statement *"Zero of 28 results constitute independent corroboration of the dual theory's content distinct from textbook QM"* are load-bearing and must be transcribed verbatim. The §2 six-observable table now carries a `NIST MCP cross-check` column (per [PR #91](https://github.com/temoTxt/PyPhysics/pull/91), already merged) that the PDF should reproduce as the per-observable provenance audit. The campaign README at [Bethe_Salpeter/README.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/README.md) anchors the proper-time/dual-theory framing.
- The per-PR results live in [Bethe_Salpeter/{01..09}_*.md](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/) (`01_NonRelHydrogen.md` through `09_HeliumExcited.md`), per the `find` listing performed in this session — each section corresponds to a numbered result block the PDF will consume.

**Finding 2 verdict state on `main` (the load-bearing nuance the issue body underplays):**

[Roadmapping/Equation_Verification/FINDINGS_for_author_review.md:209](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md#L209) currently records Finding 2 at verdict **⚠/✅ at hypothesis-(ii) framework precision**, with the trajectory `🔴 → ⚠ characterised → ⚠/✅-conditional`. The closed-form algebraic inversion `r_e/r_0 = (2 - a_e)/(2(2 + a_e))` (documented in [Dual_Relativistic_Quantum_Mechanics_I.md:58](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md#L58) §III.D-extension) reproduces the triangulated value to `3.45e-13` (within precision floor `2.5e-13`) when CODATA-full `a_e` is supplied. **Unconditional ✅** requires the hypothesis-(i) proper-time photon-propagator re-derivation tracked in **[issue #75](https://github.com/temoTxt/PyPhysics/issues/75)**. The PDF must present this `(i)/(ii)` characterisation precisely as worded on `main`, not summarised away.

**The Z-scan tension (additional layer the PDF adds on top of the above):**

The four open BS PR branches carry `.md` and `.wl` material that has *not* been merged to `main` — confirmed by `gh pr view --json files` and `git ls-tree origin/<branch>` listings performed in this session:

| Branch | New section doc | New Wolfram notebook |
|---|---|---|
| `78-bethe-salpeter-z-extension-li2plus` (PR #84) | `Bethe_Salpeter/11_Li2plus_HydrogenicIon.md` | `Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_joint_fit.wl` |
| `78-li2plus-spectroscopy` (PR #85) | `Bethe_Salpeter/12_Li2plus_Spectroscopy.md` | `Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_spectroscopy.wl` |
| `78-li2plus-hyperfine` (PR #86) | `Bethe_Salpeter/13_Li2plus_Hyperfine.md` | `Mathematica_Notebooks/Quantum_Mechanics/r_e_Li2plus_hyperfine.wl` |
| `82-hydrogenic-z-scan-g-factor` (PR #87) | `Bethe_Salpeter/14_HydrogenicIon_Zscan.md` | `Mathematica_Notebooks/Quantum_Mechanics/r_e_Zscan_fit.wl` |

PRs #84/#85/#86 reach "Branch A by construction" but in the *absence* of valid cross-Z residuals (3 of 4 Li²⁺ measurements were mis-provenanced — see PR bodies). PR #87 supplies the Z-scan residual table at Z = 1, 2, 6, 8, 14 and finds Branch A ruled out at `~10⁶–10⁷σ`, recommending Branch C (per-Z QED inheritance of `a_e(Zα)`). The packet must present both readings without taking sides.

**Existing `.wl` notebooks already on `main` that must also be translated:**

The `ls` of [Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/) confirms `r_e_triangulation.wl`, `r_e_derivation_self_energy.wl`, `r_e_derivation_variational.wl`, `r_e_schwinger_residual_test.wl`, and `BetheSalpeter_S3.wl` — these are the existing PR #62/#64/#70 outputs the packet must also render in LaTeX. Plus the Griffiths notebooks (not in scope for this packet).

**Tooling already in place this task reuses:**

- Bibliography pipeline: [Roadmapping/History/Bibliography/scaffold_bib_note.py](../../Roadmapping/History/Bibliography/scaffold_bib_note.py), [build_bibtex.py](../../Roadmapping/History/Bibliography/build_bibtex.py), [audit_human_reviewed.py](../../Roadmapping/History/Bibliography/audit_human_reviewed.py) — all confirmed present by `ls`.
- Crocco-compliance framework: [Roadmapping/Tooling/CROCCO_COMPLIANCE.md](../../Roadmapping/Tooling/CROCCO_COMPLIANCE.md) (Read in this session) and the substantive prompts directory [Roadmapping/Tooling/_prompts/](../../Roadmapping/Tooling/_prompts/) (`README.md`, `chapter_qa_review.md`, `synth_cluster_claims.md`, `synth_suggest_cross_refs.md`, `triage_summary.md` — all listed).
- TeX Live system deps (`texlive-latex-base`, `texlive-latex-extra`, `texlive-pictures`) — already required by the Manim animations pipeline per CLAUDE.md.
- Animation stills available under [Roadmapping/Animations/rendered/videos/](../../Roadmapping/Animations/rendered/videos/) — `drqm_eq18_g_factor_derivation`, `drqm_eq22_g_factor_finding`, and `hist_lamb_shift_contrast` are confirmed present and would aid the argument.

**Manuscript directory does NOT exist yet** (verified by `ls Roadmapping/Manuscripts/` returning absent). The plan creates it.

## Implementation Plan

1. **Scaffold the packet directory** under `Roadmapping/Manuscripts/author_review_packet/`:
   ```
   author_review_packet/
   ├── README.md          # build instructions, packet purpose, "for author review" framing
   ├── Makefile           # latexmk targets: paper.pdf, clean
   ├── paper.tex          # main document; \input{sections/...}; loads bibliography
   ├── references.bib     # generated by build_bibtex.py (do NOT hand-edit)
   ├── paper.pdf          # built artifact, committed
   ├── sections/
   │   ├── 00_for_author_review.tex
   │   ├── 01_abstract.tex
   │   ├── 02_campaign_overview.tex
   │   ├── 03_z1_triangulation.tex
   │   ├── 04_zscan_refinement.tex
   │   ├── 05_branch_a_vs_branch_c.tex
   │   ├── 06_mis_provenance_audit.tex
   │   ├── 07_findings_per_observable.tex
   │   ├── 08_open_theoretical_question.tex
   │   ├── 09_ai_use_disclosure_crocco.tex
   │   └── 10_questions_for_tepper.tex
   ├── equations/         # one .tex per load-bearing .wl derivation
   │   ├── BetheSalpeter_S3.tex
   │   ├── r_e_triangulation.tex
   │   ├── r_e_derivation_self_energy.tex
   │   ├── r_e_derivation_variational.tex
   │   ├── r_e_schwinger_residual.tex
   │   ├── r_e_Li2plus_joint_fit.tex
   │   ├── r_e_Li2plus_spectroscopy.tex
   │   ├── r_e_Li2plus_hyperfine.tex
   │   └── r_e_Zscan_fit.tex
   ├── figures/
   │   ├── cross_comparison_table.tex          # reproduces 10_CrossComparison.md §2 (with NIST-MCP cross-check column)
   │   ├── zscan_r_e_curve.{pdf,tex}           # reproduced from r_e_Zscan_fit.wl
   │   └── stills/                             # optional Manim stills (drqm_eq22_g_factor_finding, hist_lamb_shift_contrast)
   └── prompts/
       └── README.md      # Crocco-disclosed substantive prompts used during drafting (if any)
   ```

2. **Consume the BS PR branches without merging.** Per the issue body's reframed answer, draft from the four open branches directly. Mechanically:
   ```bash
   git fetch origin 78-bethe-salpeter-z-extension-li2plus 78-li2plus-spectroscopy \
                   78-li2plus-hyperfine 82-hydrogenic-z-scan-g-factor
   for branch in 78-bethe-salpeter-z-extension-li2plus 78-li2plus-spectroscopy \
                 78-li2plus-hyperfine 82-hydrogenic-z-scan-g-factor; do
     git show origin/$branch:Roadmapping/Quantum_Mechanics/Bethe_Salpeter/...   # read the new section doc
     git show origin/$branch:Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/... # read the new .wl
   done
   ```
   Record per-PR cell→equation cross-references in each `equations/*.tex` header so Tepper can trace any equation back to its Wolfram-MCP symbolic check.

3. **Translate each `.wl` cell to a LaTeX block.** Per-cell pattern, applied uniformly across the nine `equations/*.tex` files:
   ```latex
   % --- Cell N (source: r_e_Li2plus_joint_fit.wl line LLL) ---
   % Wolfram-MCP confirmation: Result: 0  (verified YYYY-MM-DD)
   \begin{equation} ... \end{equation}
   ```
   Honour the three Wolfram MCP gotchas from [CLAUDE.md](../../CLAUDE.md) when re-deriving anything to verify the transcription: single-line `ClearAll;Print;FullSimplify` chains; use `potV` not `V`; use `ee` for electron charge; non-commutative `Dot` for symbolic 3-vectors.

4. **Write paper sections** as individual `.tex` files inputted from `paper.tex`. Each section has a fixed scope:
   - `00_for_author_review.tex` — short cover note: "for author review (not corrections)," names Tepper as senior author of DRQM I, points to the closing question list.
   - `01_abstract.tex` — single paragraph: campaign + Z=1 triangulation + Z-scan refinement; transcribe the honest-scope statement from `10_CrossComparison.md` §1.
   - `02_campaign_overview.tex` — the 28-result inventory table reproduced from `10_CrossComparison.md` §1; honest-scope statement verbatim.
   - `03_z1_triangulation.tex` — `r_e/r_0 = 0.4994205099128317` framed as a Z=1 joint-fit calibration; cite PR #62; reproduce the §III.D closed-form `r_e/r_0 = (2 - a_e)/(2(2 + a_e))` from `Dual_Relativistic_Quantum_Mechanics_I.md:58` and the 3.45e-13 agreement; record the verdict trajectory `🔴 → ⚠ characterised → ⚠/✅-conditional` verbatim from `FINDINGS_for_author_review.md:209`; flag the hypothesis-(i)/(ii) characterisation and reference issue #75 for unconditional ✅.
   - `04_zscan_refinement.tex` — the Z-scan dataset table at Z = 1, 2, 6, 8, 14 from `r_e_Zscan_fit.wl` (consumed from branch `82-hydrogenic-z-scan-g-factor`); back-fit `r_e^(Z)/r_0` curve along QED `-(Zα)²/3`; cite the four measurements (Schneider 2022 *Nature* 606,878; Sturm 2011 *PRL* 107,023002; Verdú 2004 *PRL* 92,093002; Sturm 2013 *PRL* 110,263002).
   - `05_branch_a_vs_branch_c.tex` — present **both** readings without endorsement: Branch A (Z-universal cutoff, from PRs #84/#85/#86, blocked-residual context) and Branch C (per-Z QED inheritance, from PR #87, ~10⁶σ Branch-A refutation). Quote the PR #87 verdict line verbatim. State the tension explicitly. Defer adjudication to the closing question list.
   - `06_mis_provenance_audit.tex` — the audit log already documented in PR bodies #84/#85/#86: `Sturm 2014 Nature` cited for #1 actually measured ¹²C⁵⁺ (Z=6), not Li²⁺ (Z=3); the brief's `7,367 MHz` for #3 is helium-like Li⁺ 2P FS, not hydrogenic Li²⁺; `Beckmann 1974` cited for #4 is a nuclear-moment paper, not a Li²⁺ HFS measurement. Ask Tepper to confirm or correct.
   - `07_findings_per_observable.tex` — the §2 six-observable cross-comparison table reproduced from `10_CrossComparison.md:61-68`, with the `NIST MCP cross-check` column from [PR #91](https://github.com/temoTxt/PyPhysics/pull/91) preserved (CODATA g-factor, ASD-theoretical H 2P, curated list_dirac_targets entries — provenance audit, not corroboration).
   - `08_open_theoretical_question.tex` — first-principles `r_e` derivation per issue [#54](https://github.com/temoTxt/PyPhysics/issues/54) and the hypothesis-(i) proper-time photon propagator per issue [#75](https://github.com/temoTxt/PyPhysics/issues/75). State both readings of what the derivation would need to reproduce (Z-invariant constant vs Z-curve), per the Branch A vs C split.
   - `09_ai_use_disclosure_crocco.tex` — Crocco §5 disclosure (per [Roadmapping/Tooling/CROCCO_COMPLIANCE.md](../../Roadmapping/Tooling/CROCCO_COMPLIANCE.md)): Wolfram MCP and CODATA/PDG/NIST-ASD lookups are pragmatic; any prose-generation steps are substantive and must list the prompt-of-record from `prompts/`.
   - `10_questions_for_tepper.tex` — the six numbered questions from the issue body (Z-scan validity; Eq. (III.23) freedom; Branch C as "borrowed QED a_e(Zα)" acceptability; verdict downgrade question; valid Li²⁺ measurements; triangulated value as Z=1 calibration). Add a seventh question if hypothesis-(i) re-derivation (#75) overlaps with the requested verdict.

5. **Build `references.bib` via the existing pipeline.**
   ```bash
   uv run python Roadmapping/History/Bibliography/build_bibtex.py \
     --output Roadmapping/Manuscripts/author_review_packet/references.bib
   ```
   For any `\cite{...}` in the packet whose key has no YAML stub under `Roadmapping/History/Bibliography/{Primary,Retrospective}/`, scaffold it via `scaffold_bib_note.py` **before** the cite appears in the packet. Per CLAUDE.md rule "Never invent citations," no `\cite{...}` may resolve to a missing stub at PDF build time.

6. **Build figures.**
   - `cross_comparison_table.tex`: LaTeX `tabular` rendering of `10_CrossComparison.md:61-68` plus the NIST-MCP cross-check column.
   - `zscan_r_e_curve.pdf`: produced from `r_e_Zscan_fit.wl` (consumed from branch `82-hydrogenic-z-scan-g-factor`) — back-fit `r_e^(Z)/r_0` vs Z curve overlaid with QED `-(Zα)²/3` reference.
   - Optional Manim stills under `figures/stills/`: candidates from `Roadmapping/Animations/rendered/videos/` are `drqm_eq22_g_factor_finding/`, `drqm_eq18_g_factor_derivation/`, and `hist_lamb_shift_contrast/` — confirmed present by `ls`. Use only stills that improve the argument at print resolution; skip others.

7. **Crocco-compliance enforcement.** For each substantive AI use during drafting, commit the prompt-of-record to `Roadmapping/Manuscripts/author_review_packet/prompts/` and reference it from §09. Confirm against the pragmatic vs substantive split documented in [Roadmapping/Tooling/CROCCO_COMPLIANCE.md](../../Roadmapping/Tooling/CROCCO_COMPLIANCE.md) §2.

8. **Build and commit.**
   ```bash
   cd Roadmapping/Manuscripts/author_review_packet/
   make paper.pdf       # latexmk via Makefile
   git add paper.pdf references.bib equations/ sections/ figures/ prompts/ paper.tex Makefile README.md
   ```
   The committed `paper.pdf` is the artifact Tepper opens; reviewers don't have to rebuild it.

## Files to Modify

| File | Change |
|---|---|
| `Roadmapping/Manuscripts/author_review_packet/README.md` | New. Build/purpose/for-author-review overview. |
| `Roadmapping/Manuscripts/author_review_packet/Makefile` | New. `paper.pdf`, `clean` targets via latexmk. |
| `Roadmapping/Manuscripts/author_review_packet/paper.tex` | New. Main document; `\input{sections/*}`; loads `references.bib`. |
| `Roadmapping/Manuscripts/author_review_packet/references.bib` | New. Generated by `build_bibtex.py`. |
| `Roadmapping/Manuscripts/author_review_packet/paper.pdf` | New. Committed build artifact (the file Tepper opens). |
| `Roadmapping/Manuscripts/author_review_packet/sections/00_for_author_review.tex` | New. Cover note. |
| `Roadmapping/Manuscripts/author_review_packet/sections/01_abstract.tex` | New. Abstract with honest-scope statement verbatim. |
| `Roadmapping/Manuscripts/author_review_packet/sections/02_campaign_overview.tex` | New. 28-result inventory + closing statement. |
| `Roadmapping/Manuscripts/author_review_packet/sections/03_z1_triangulation.tex` | New. Triangulated value + §III.D closed form + verdict trajectory. |
| `Roadmapping/Manuscripts/author_review_packet/sections/04_zscan_refinement.tex` | New. Z-scan dataset + back-fit curve from PR #87. |
| `Roadmapping/Manuscripts/author_review_packet/sections/05_branch_a_vs_branch_c.tex` | New. Both readings, no adjudication. |
| `Roadmapping/Manuscripts/author_review_packet/sections/06_mis_provenance_audit.tex` | New. Per-citation audit log from PRs #84/#85/#86 bodies. |
| `Roadmapping/Manuscripts/author_review_packet/sections/07_findings_per_observable.tex` | New. §2 cross-comparison table + NIST-MCP cross-check column. |
| `Roadmapping/Manuscripts/author_review_packet/sections/08_open_theoretical_question.tex` | New. Issues #54 + #75 framing under both readings. |
| `Roadmapping/Manuscripts/author_review_packet/sections/09_ai_use_disclosure_crocco.tex` | New. Crocco §5 disclosure tied to `prompts/`. |
| `Roadmapping/Manuscripts/author_review_packet/sections/10_questions_for_tepper.tex` | New. Numbered closing-question list. |
| `Roadmapping/Manuscripts/author_review_packet/equations/BetheSalpeter_S3.tex` | New. LaTeX translation of `BetheSalpeter_S3.wl`. |
| `Roadmapping/Manuscripts/author_review_packet/equations/r_e_triangulation.tex` | New. LaTeX translation of `r_e_triangulation.wl`. |
| `Roadmapping/Manuscripts/author_review_packet/equations/r_e_derivation_self_energy.tex` | New. LaTeX translation of `r_e_derivation_self_energy.wl`. |
| `Roadmapping/Manuscripts/author_review_packet/equations/r_e_derivation_variational.tex` | New. LaTeX translation of `r_e_derivation_variational.wl`. |
| `Roadmapping/Manuscripts/author_review_packet/equations/r_e_schwinger_residual.tex` | New. LaTeX translation of `r_e_schwinger_residual_test.wl`. |
| `Roadmapping/Manuscripts/author_review_packet/equations/r_e_Li2plus_joint_fit.tex` | New. LaTeX translation of `r_e_Li2plus_joint_fit.wl` (consumed from branch `78-bethe-salpeter-z-extension-li2plus`). |
| `Roadmapping/Manuscripts/author_review_packet/equations/r_e_Li2plus_spectroscopy.tex` | New. LaTeX translation of `r_e_Li2plus_spectroscopy.wl` (from branch `78-li2plus-spectroscopy`). |
| `Roadmapping/Manuscripts/author_review_packet/equations/r_e_Li2plus_hyperfine.tex` | New. LaTeX translation of `r_e_Li2plus_hyperfine.wl` (from branch `78-li2plus-hyperfine`). |
| `Roadmapping/Manuscripts/author_review_packet/equations/r_e_Zscan_fit.tex` | New. LaTeX translation of `r_e_Zscan_fit.wl` (from branch `82-hydrogenic-z-scan-g-factor`). |
| `Roadmapping/Manuscripts/author_review_packet/figures/cross_comparison_table.tex` | New. LaTeX `tabular` of `10_CrossComparison.md:61-68` + NIST-MCP column. |
| `Roadmapping/Manuscripts/author_review_packet/figures/zscan_r_e_curve.pdf` | New. Back-fit `r_e^(Z)/r_0` curve from `r_e_Zscan_fit.wl`. |
| `Roadmapping/Manuscripts/author_review_packet/figures/zscan_r_e_curve.tex` | New. TikZ/pgfplots source for the curve. |
| `Roadmapping/Manuscripts/author_review_packet/figures/stills/` | New (directory). Optional Manim stills (candidates: `drqm_eq22_g_factor_finding`, `hist_lamb_shift_contrast`); include only those that aid the argument at print resolution. |
| `Roadmapping/Manuscripts/author_review_packet/prompts/README.md` | New. Index of Crocco-disclosed substantive prompts (if any). |
| `Roadmapping/History/Bibliography/{Primary,Retrospective}/<missing_cite_key>.md` | Potentially new. Scaffold via `scaffold_bib_note.py` any `\cite{...}` whose stub does not exist; zero-or-more new files. |

## Dependencies

- **System packages** (already required by the Manim animations pipeline per [CLAUDE.md](../../CLAUDE.md)): `texlive-latex-base`, `texlive-latex-extra`, `texlive-pictures`, plus `latexmk` for the Makefile build.
- **Existing repo scripts** (no new dependencies): `Roadmapping/History/Bibliography/build_bibtex.py`, `scaffold_bib_note.py`, `audit_human_reviewed.py` — invoked via `uv run python ...` per CLAUDE.md.
- **No new Python deps.** No new MCP servers. No changes to root or sub-project `pyproject.toml` files.
- **External access:** none required at build time (the bibliography YAML and the four BS PR branches are local once `git fetch origin` runs).

## Acceptance Criteria

- [ ] `Roadmapping/Manuscripts/author_review_packet/paper.pdf` exists and renders cleanly with no LaTeX errors.
- [ ] Every load-bearing `.wl` derivation listed in the `Files to Modify` table appears as a per-cell LaTeX block under `equations/`, each with a `% Cell N (source: <file> line LLL) — Wolfram-MCP confirmation: Result: 0 (verified <date>)` header.
- [ ] The §2 six-observable cross-comparison table (with the NIST-MCP cross-check column from [PR #91](https://github.com/temoTxt/PyPhysics/pull/91)) is reproduced in `figures/cross_comparison_table.tex` and rendered in the PDF.
- [ ] The Z-scan back-fit `r_e^(Z)/r_0` curve from `r_e_Zscan_fit.wl` appears as `figures/zscan_r_e_curve.pdf`.
- [ ] The honest-scope statement *"Zero of 28 results constitute independent corroboration of the dual theory's content distinct from textbook QM"* appears verbatim in the abstract.
- [ ] Both Branch A and Branch C readings are presented in `sections/05_branch_a_vs_branch_c.tex` without endorsement; the PDF does not adjudicate.
- [ ] Finding 2's verdict trajectory `🔴 → ⚠ characterised → ⚠/✅-conditional` is preserved verbatim from `FINDINGS_for_author_review.md:209`, and issue [#75](https://github.com/temoTxt/PyPhysics/issues/75) (hypothesis-(i)) is referenced as the unconditional-✅ path.
- [ ] The mis-provenance audit (Sturm 2014 / `7,367 MHz` / Beckmann 1974) is in `sections/06_mis_provenance_audit.tex` exactly as documented in PR bodies #84/#85/#86.
- [ ] `sections/10_questions_for_tepper.tex` contains the six numbered questions from the issue body, plus any seventh question added to disambiguate hypothesis-(i) (#75) vs Z-scan (PR #87) verdicts.
- [ ] A Crocco §5 AI-use disclosure section (`sections/09_ai_use_disclosure_crocco.tex`) is present and references any prompts committed under `prompts/`.
- [ ] Every `\cite{...}` in the packet resolves to an existing YAML stub under `Roadmapping/History/Bibliography/{Primary,Retrospective}/`; running `build_bibtex.py` produces `references.bib` with no missing-key errors.
- [ ] No edits to files outside `Roadmapping/Manuscripts/author_review_packet/` and the conditional bibliography stubs row in the `Files to Modify` table.

## Testing

Reviewer commands (run from the repo root):

```bash
# 1. Build the PDF (no network needed; consumes BS PR branches via git fetch).
git fetch origin 78-bethe-salpeter-z-extension-li2plus 78-li2plus-spectroscopy \
                 78-li2plus-hyperfine 82-hydrogenic-z-scan-g-factor
(cd Roadmapping/Manuscripts/author_review_packet/ && make paper.pdf)
test -s Roadmapping/Manuscripts/author_review_packet/paper.pdf || echo "BUILD FAILED"

# 2. Confirm every \cite{...} resolves to a YAML stub before bibtex builds.
grep -rohE '\\cite\{[^}]+\}' Roadmapping/Manuscripts/author_review_packet/ \
  | sed -E 's/\\cite\{([^}]+)\}/\1/' | tr ',' '\n' | sort -u \
  | while read k; do
      test -f Roadmapping/History/Bibliography/Primary/${k}.md \
        || test -f Roadmapping/History/Bibliography/Retrospective/${k}.md \
        || echo "MISSING STUB: $k"
    done

# 3. Rebuild references.bib and confirm no missing-key errors.
uv run python Roadmapping/History/Bibliography/build_bibtex.py \
  --output Roadmapping/Manuscripts/author_review_packet/references.bib

# 4. Honest-scope statement transcription check.
grep -q "Zero of 28 results constitute independent corroboration" \
  Roadmapping/Manuscripts/author_review_packet/sections/01_abstract.tex \
  || echo "honest-scope statement missing or paraphrased"

# 5. Verdict-trajectory check (Finding 2).
grep -q "characterised.*conditional\|hypothesis-(ii)" \
  Roadmapping/Manuscripts/author_review_packet/sections/03_z1_triangulation.tex \
  || echo "Finding 2 verdict trajectory not preserved"

# 6. Both Branch A and Branch C presented without adjudication.
grep -q "Branch A" Roadmapping/Manuscripts/author_review_packet/sections/05_branch_a_vs_branch_c.tex \
  && grep -q "Branch C" Roadmapping/Manuscripts/author_review_packet/sections/05_branch_a_vs_branch_c.tex \
  || echo "Branch A vs C section incomplete"

# 7. Question list present.
grep -cE '^\\(item|question)' Roadmapping/Manuscripts/author_review_packet/sections/10_questions_for_tepper.tex
# expect >= 6

# 8. Crocco AI-use disclosure present.
test -s Roadmapping/Manuscripts/author_review_packet/sections/09_ai_use_disclosure_crocco.tex \
  || echo "Crocco disclosure missing"
```

Tests added:

- A make-target smoke check (`make paper.pdf` builds without LaTeX errors) is the primary acceptance gate.
- The `grep` checks in steps 2-8 above are the structural / content-fidelity gates. They are intended to be run by the reviewer; they are not a permanent test suite.
- No Python unit tests are added (the repo has no pytest suite at the root; per CLAUDE.md "There are no tests").

Reviewer additionally opens `paper.pdf` and confirms (a) the 28-result inventory matches `10_CrossComparison.md` §1, (b) the Z-scan curve is legible at print resolution, (c) the six numbered questions appear as the closing section, (d) "For author review" framing is present in the cover note, and (e) no `\cite{...}` is unresolved in the rendered bibliography.
