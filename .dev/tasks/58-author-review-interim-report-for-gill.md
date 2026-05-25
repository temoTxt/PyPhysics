# Implementation plan: Interim report for Tepper Gill (author review)

**Tracks:** issue [#58](https://github.com/temoTxt/PyPhysics/issues/58).
**Status:** plan; report drafting + LaTeX/PDF build pipeline to follow.
**Dependencies:** [#42 PR #51](https://github.com/temoTxt/PyPhysics/pull/51), [#49 PR #52](https://github.com/temoTxt/PyPhysics/pull/52), [#50 PR #53](https://github.com/temoTxt/PyPhysics/pull/53) — all three campaign PRs delivered as of 2026-05-25; this work consolidates them into an author-engageable summary.

This is **not** a multi-PR campaign. It is a single report-authoring task with three deliverables (markdown source + LaTeX + PDF) plus a small amount of supporting tooling (build script + `Author_Reports/` folder convention). The plan inherits the [§13 honest-framing discipline of #42](42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix) verbatim — the report's value to Tepper is that he can trust its scoping.

---

## 1. Goal and scope

Produce an interim report for Tepper Gill summarising the repository's verification + campaign work since 2026-05-11, with explicit author-side questions and decision points. The report's audience is a single known collaborator (a co-author of *DRQM I* and the foundational dual-theory papers); its purpose is to surface honest scoping, flagged findings, and open questions in a form Tepper can act on in one sitting.

**Honest framing of the request.** Tepper has not asked for this report. The repository's verification campaign has surfaced three errata-level findings in his published work plus six downstream operational signatures of Finding 2; the campaign documents (#51, #52, #53) record these in a per-PR form, but no consolidated author-engageable summary exists. The report is the author's (Trey's) initiative to surface this material to Tepper in a form that does not require him to read the per-PR campaign documents. It is unsolicited, and the report's cover note in §1 of the deliverable acknowledges this explicitly. Phase 4's human-acceptance pass is the place to decide whether the report is ready to send.

### Why this earns its own thread

- **Author-engageable format.** The verification campaign's [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) covers three substantive errata findings; the per-PR campaign documents (50+ files across #42 / #49 / #50) cover the load-bearing structural results. Neither is an author-engageable summary.
- **PDF deliverable.** Tepper receives a PDF, not a GitHub URL trail. The LaTeX/PDF build pipeline is the bridge between the markdown source-of-truth (which keeps `<!-- TODO -->` blocks for human-acceptance) and the deliverable (which does not).
- **Folder structure for any successor reports.** `Roadmapping/Author_Reports/` becomes the location *if* further author reports are written; the README documents naming conventions so the folder remains navigable. *Not* a claim that this report's build pipeline is the standard convention for hypothetical future reports — that generalisation is deferred until a second use-case actually appears (see §5 scope note).

### Explicit non-goals

- Not a journal manuscript. Working document for a known collaborator.
- Not a re-summary of any per-PR campaign document — cross-reference, do not duplicate.
- Not a pre-judgement of Tepper's answers. The report frames the questions; it does not propose Tepper's preferred answer for him.
- Not a new-campaign plan. Suggested next steps in the report's §7 are tactical (resolve specific questions), not strategic (open new campaigns).

---

## 2. Repository structure

New top-level folder, sibling to `Equation_Verification/` / `Electromagnetism/` / `Quantum_Mechanics/` / `History/`:

```
Roadmapping/Author_Reports/
├── README.md                              # convention for future reports
├── build_report.sh (or .py)               # single-command markdown → LaTeX → PDF
├── 2026-05_interim_for_gill.md            # markdown source-of-truth
├── 2026-05_interim_for_gill.tex           # LaTeX source (derived)
└── 2026-05_interim_for_gill.pdf           # built PDF (committed)
```

The folder's `README.md` documents the convention so subsequent author reports follow the same naming + build pipeline.

---

## 3. Report structure (the markdown source)

**Length budget: 3–7 pages in the final PDF, soft target 5.** Earlier drafts of this plan set a hard `[3, 5]` cap as a load-bearing constraint; that cap is relaxed because the campaign's honest-scoping caveats (back-fit caveat on branch (c), K-identity-is-by-definition caveat, Lamb-shift-is-reproduction-not-endorsement caveat, asymmetric absorbing-correction caveat) are *paragraph-scale* and the steel-man edits to PRs #51/#53 confirmed that compressing them costs fidelity. The honest-scoping discipline is load-bearing for the report; the page count is not. If the natural drafted length comes in at 6–7 pages because the caveats demand the space, the report is the right length and the caveats stay.

Soft section targets (at the 5-page nominal): §§1, 2, 7 ≈ ½ page each; §§3, 4, 5 ≈ ¾ page each; §6 the load-bearing Q1–Q8 section ≈ 1 page. If the report goes long, the page budget yields *before* the caveats do.

The build-step defensive check (§5 below) fails loudly **above 7** (overshoot is real) or **below 3** (undershoot suggests missing substance, most likely §6 Q1–Q8 needs more context per question or §4 conditional predictions are under-cited). Inside `[3, 7]` the check passes; the human-acceptance pass owns the final judgement on length. **Honest scoping never gets sacrificed to hit a length target.**

`Roadmapping/Author_Reports/2026-05_interim_for_gill.md` has the following sections, in order, per [issue #58 acceptance criterion 2](https://github.com/temoTxt/PyPhysics/issues/58):

1. **§1 Cover note.** One paragraph — what the report is, why now, what response we'd value.
2. **§2 Scope of the campaigns to date.** The four threads (verification + #42 + #49 + #50): what each set out to do, what it actually delivered. Cross-references to the PRs (#51, #52, #53) rather than re-summary.
3. **§3 What we have actually shown.** Load-bearing structural results, with explicit honest scoping:
   - Algebraic identities of `K` (trivial by construction — the [`BetheSalpeter_S3.wl` honest-scope note](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S3.wl) is the model: the identity `K = mc² + p²/(2m)` is true by *definition*, not a physical claim about relativistic content).
   - Dual-Dirac FW reduction at leading `(Zα)²` (non-trivial; DRQM I §II.D dependent).
   - Bethe (1947) self-energy estimate reproduced at the Bethe-estimate precision floor ([#50 PR E](../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/05_LambShift.md)).
   - Non-relativistic correspondence: 20 of 28 #50 results ✅ unconditionally.
4. **§4 Conditional predictions.** What the campaigns predict *if* the framework's load-bearing claims hold:
   - The 6 `r_e`-discriminator observables and their branched verdicts (b/c).
   - The Lamb-shift at full one-loop precision ([#55](https://github.com/temoTxt/PyPhysics/issues/55); apparatus-limited).
   - The radiation-reaction dissolution claim of J3e-P16.1 ([#42 PR D](../../Roadmapping/Electromagnetism/Jackson/Ch16_Radiation_Damping.md)) — conditional on author intent.
5. **§5 The three flagged findings.** Verbatim or near-verbatim recap from [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md), plus #50's six operational signatures of Finding 2.
6. **§6 Open questions for the author.** Load-bearing section. Eight specific decision points (Q1–Q8) as enumerated in [issue #58 acceptance criterion 2 §6](https://github.com/temoTxt/PyPhysics/issues/58):
   - Q1 — `r_e` resolution (branch (b) vs (c)) per [#54](https://github.com/temoTxt/PyPhysics/issues/54)
   - Q2 — `r_μ` and `r_p` numerical values (DRQM I Eq. III.23)
   - Q3 — Maxwell paper Eq. (24) erratum confirmation
   - Q4 — TCEP Eq. (4.16) sign-typo confirmation
   - Q5 — Operator-ordering choice in DRQM I §II.D FW reduction
   - Q6 — Radiation-reaction dissolution interpretation (J3e-P16.1)
   - Q7 — *Roadmap question (not a commitment ask):* is proper-time one-loop dual-QED on the framework's medium-term roadmap, or honestly out of scope? ([#55](https://github.com/temoTxt/PyPhysics/issues/55)). The earlier framing — "framework's intent to develop" — was downgraded to a roadmap-yes/no/maybe question per §10 honest framing, because asking a senior co-author to commit to a multi-year programme via a PDF checkbox is too heavy a request.
   - Q8 — Dual-Dirac propagator sufficiency for one-loop apparatus
7. **§7 Suggested next steps.** Short ordered list, ordered by what would unblock the most downstream work. Tactical, not strategic.

---

## 4. Per-section authoring discipline

The report is **substantive AI** end-to-end. Per the Crocco compliance protocol, every paragraph carries a per-paragraph `<!-- TODO: human reviews and fills in -->` block in the markdown source. The committed `.md` retains these for the human-acceptance pass; the built PDF does not (per the build-script's TODO-stripping step in §5 below).

Voice follows [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) **as a default**: Gill's published-paper voice, "we plural", connective openers, narrative attribution. The Phase 4 human-acceptance pass should **reconsider** this choice for the report specifically. Two countervailing concerns: (a) `VOICE_MATCH_GILL.md` is Claude+Trey's reading of Gill's published voice, not Gill's confirmation of his own voice; (b) addressing a senior co-author in their own writing style can read as parody or pandering more easily than as deference. The default is the inherited campaign voice; the override is a more neutral co-author-to-co-author register if the human reviewer prefers. Flag explicitly rather than assume.

Honest scoping in every section: algebraic identities labelled as such; conditional predictions labelled as such; out-of-scope claims labelled as such. The model is the [`BetheSalpeter_S3.wl` honest-scope linter note](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S3.wl) — verbatim discipline.

---

## 5. Build pipeline

A single-command build from the markdown source to the deliverable PDF.

**Scope note (revised).** Earlier drafts justified the six-step pipeline as "convention-setting for future reports" (cf. §1.2). That justification is downgraded: there is no second author report planned, no second recipient, and no second use case. Per the CLAUDE.md *don't design for hypothetical future requirements* discipline, the build script is scoped to *this* report only. If a second report appears later, the script will be generalised then with the benefit of an actual second use-case. The pipeline is still six steps because each one solves a concrete problem (TODO-strip, math passthrough, page-count guard, etc.), not because it is forward-engineered for hypothetical reports.

Pipeline:

1. **Markdown → LaTeX (pandoc)**. `pandoc 2026-05_interim_for_gill.md -o 2026-05_interim_for_gill.tex` with appropriate template + filter flags. Math passes through unchanged (`$...$`, `$$...$$`). Tables convert via pandoc's reliable markdown-table → LaTeX-table conversion. Code blocks via `listings` or `verbatim`.
2. **TODO-comment stripping**. Markdown HTML comments `<!-- TODO: human reviews and fills in -->` must not leak into the PDF. Either:
   - Use a pandoc filter (Lua filter) that drops HTML comments at conversion time, or
   - Pre-process the `.md` to strip `<!-- ... -->` lines before pandoc, or
   - Use pandoc's `--strip-comments` if available.
   The build script must **fail loudly** if any `<!-- TODO -->` is present at PDF-build time after the strip step (defensive check).
3. **LaTeX → PDF**. `pdflatex` or `xelatex`, run twice (for cross-references). Document class `article` with `geometry`, `hyperref`, `amsmath`, `booktabs`. PDF metadata set via pandoc `--metadata` or LaTeX `\hypersetup`: Title, Author = "Trey Morris with Claude Opus 4.7", Date, Subject = "Interim report for author review".
4. **Cross-reference handling**. GitHub-relative `../blob/main/...` links cannot be followed from a PDF. **Every cross-reference in the markdown source must appear in the PDF**, either as a `\footnote{See \texttt{path} in the repository.}` (default, for inline citations) or in a closing references-table appendix (for high-density sections like §5 where multiple findings cite the same files). Earlier drafts of this plan allowed "less-cited cross-references to be dropped" from the PDF; that policy is withdrawn — if a cross-reference is in the report, the report's claim depends on it, and dropping it leaves the PDF making claims without supporting references. The page-count budget accommodates the additional footnote/appendix space (per §3's relaxed `[3, 7]` window).
5. **Bibliography integration** (if cited). Reuse the existing [`Roadmapping/History/Bibliography/build_bibtex.py`](../../Roadmapping/History/Bibliography/build_bibtex.py) pipeline to produce a `bibliography.bib` and reference it from the LaTeX via `\bibliography{}`. Cite only existing bib stubs; do not fabricate citations.
6. **Page-count defensive check.** After PDF build completes, the script runs `pdfinfo` (or `pdftk … dump_data`, or `qpdf --show-npages`) to extract the page count and fails loudly if outside `[3, 7]` (per the relaxed §3 budget). Sample check:

   ```bash
   pages=$(pdfinfo "$PDF" | awk '/^Pages:/ {print $2}')
   if [ "$pages" -lt 3 ] || [ "$pages" -gt 7 ]; then
     echo "ERROR: PDF page count $pages is outside the [3, 7] budget per plan §3." >&2
     exit 1
   fi
   ```

   This check is symmetric with the TODO-comment defensive check: both fail the build rather than silently shipping a non-conforming PDF. The window is `[3, 7]` rather than the tighter `[3, 5]` an earlier draft proposed; honest scoping is load-bearing and is allowed to push the page count into the 6–7 range when needed (see §3).

System dependencies (reused from Manim sub-project per [`CLAUDE.md`](../../CLAUDE.md)): `pandoc`, `texlive-latex-base`, `texlive-latex-extra`, `texlive-pictures`, plus `poppler-utils` (for `pdfinfo`). Document install command in the script's docstring.

---

## 6. PR sequencing

This task is delivered as a **single PR** (not a multi-PR campaign). The phase ordering is **deliberate**: no PDF is built and committed before the human-acceptance pass completes, so Tepper never sees a finished-looking PDF that is actually a pre-review AI draft (which would materially violate Crocco §5 — see §10 honest framing).

| Phase | Scope | Status |
|---|---|---|
| Phase 1 | `Roadmapping/Author_Reports/` scaffold + `README.md` convention + build-script skeleton (script tested on a stub `.md`, but PDF *not* committed yet) | pending |
| Phase 2 | `2026-05_interim_for_gill.md` draft — §1 through §7 with per-paragraph TODO blocks | pending |
| Phase 3 | Build pipeline *dry-runs* against the draft `.md` (verify pandoc + LaTeX + page-count + TODO-strip mechanics on the actual content) but does **not** commit a PDF yet | pending |
| Phase 3.5 | **Devil's-advocate self-review** of the Phase 2 draft against the campaign's §13 honest-framing discipline + **steel-man patches** applied to the `.md`. Analogous to the workflow that pulled PRs #51 and #53's per-chapter docs back into alignment with §13 after they drifted above it. The AI-drafted report is the most likely place this drift recurs; the gate is explicit. The self-review's output (the critique + the steel-man diff) is committed alongside the draft so the human reviewer in Phase 4 can see what was already caught | pending |
| Phase 4 | Human-acceptance pass: review remaining TODO blocks, edit prose, confirm honest scoping after Phase 3.5 patches, **build and commit final PDF** | requires human review |
| Phase 5 | PR opened (closes #58) | pending |

Phases 1–3 are AI-authored work that does not ship a PDF; Phase 3.5 is the AI-authored devil's-advocate / steel-man gate that catches honest-framing drift before human review; Phase 4 is human-authored (the load-bearing review pass — the committed `.md` is a draft until this pass completes, and the PDF is built and committed only at the end of Phase 4); Phase 5 closes the loop.

---

## 7. Acceptance criteria

Per [issue #58](https://github.com/temoTxt/PyPhysics/issues/58):

1. Markdown source at `Roadmapping/Author_Reports/2026-05_interim_for_gill.md` (criterion 1 of issue).
2. Report structure per §3 above (criterion 2 of issue): §§1–7 with the eight Q1–Q8 open questions explicit.
3. Honest scoping discipline throughout (criterion 3): per-paragraph TODO blocks, no overselling, scoping labels on every identity / prediction / out-of-scope claim.
4. Voice match per [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) (criterion 4).
5. Crocco compliance (criterion 5): substantive AI end-to-end, per-paragraph TODO blocks, no AI-fabricated citations.
6. Cross-references present (criterion 6): the seven specific cross-references named in #58's body.
7. LaTeX/PDF deliverable (criterion 7): `.tex`, `.pdf`, and `build_report.sh` (or `.py`) all committed; build is single-command and rebuildable (see §13 reproducibility scope; *not* byte-identical determinism); PDF carries correct metadata; no `<!-- TODO -->` leaks into the PDF.
8. **PDF page count is in `[3, 7]`** per §3 length budget (relaxed from an earlier `[3, 5]` cap to protect honest-scoping caveats); the build's `pdfinfo`-based defensive check (§5) passes.
9. **The PDF is committed only at the end of Phase 4 (human-acceptance pass)**, not during Phases 1–3. The build script's dry-run capability is exercised in Phase 3 to verify the pipeline mechanics, but the deliverable PDF — the artifact Tepper would actually receive — does not enter the repository until a human has reviewed the AI-drafted `.md` (and the Phase 3.5 devil's-advocate self-review's steel-man patches) and accepted them.
10. **Phase 3.5 devil's-advocate self-review artifact** is committed alongside the draft `.md`. The artifact is a short markdown file (one or two pages) recording the critiques and the corresponding steel-man patches applied; it lets the Phase 4 human reviewer see what was already caught vs what still needs attention. Analogous to the workflow that produced the `913802d` and earlier steel-man commits on PRs #51/#53.

The PR is mergeable when all ten hold and a human-acceptance pass has reviewed (and signed off on) the report's prose.

---

## 8. Files to modify / create

| File | Change |
|---|---|
| `Roadmapping/Author_Reports/README.md` | Create — convention for future author reports (naming, build pipeline, expectations) |
| `Roadmapping/Author_Reports/build_report.sh` (or `.py`) | Create — single-command markdown → LaTeX → PDF build script with TODO-stripping + defensive check |
| `Roadmapping/Author_Reports/2026-05_interim_for_gill.md` | Create — report markdown source, §§1–7 per §3 of this plan |
| `Roadmapping/Author_Reports/2026-05_interim_for_gill.tex` | Create (derived) — LaTeX from pandoc, human-edited as needed; *built and committed at end of Phase 4, after human-acceptance* |
| `Roadmapping/Author_Reports/2026-05_interim_for_gill.pdf` | Create (derived) — built PDF; *committed at end of Phase 4, after human-acceptance* (NOT during Phases 1–3, per §6 ordering + §7 acceptance criterion 9) |
| `Roadmapping/Author_Reports/2026-05_interim_for_gill_self_review.md` | Create — Phase 3.5 devil's-advocate self-review artifact + steel-man patch summary (one to two pages), committed alongside the draft `.md` per §7 acceptance criterion 10 |

No source-code or campaign-document changes. The report only *consolidates* existing campaign documents; it does not modify them.

---

## 9. Out-of-scope

- **Re-authoring or correcting any campaign document.** The report references the campaign documents; it does not edit them. If a reviewer identifies an error in a per-PR document while drafting the report, that error is filed as a separate issue, not fixed inline.
- **Building infrastructure for arbitrary report types.** The build script is scoped to *this* report — academic article style, with the specific TODO-stripping discipline. Earlier drafts of this plan justified the script as "convention-setting for future reports"; that framing is downgraded (see §5 scope note) because no second report is planned. A more general report-pipeline (with options for different document classes, journal templates, etc.) is a future thread that will be designed when a second use-case actually appears.
- **Author-side response to the report.** Tepper's answers to Q1–Q8 are recorded in follow-on threads (#54 already opened for Q1; #55 for Q7; further as needed); they are out of scope for this PR.
- **Translation of the campaign verdicts.** The report records what the campaigns said; it does not re-derive their verdicts or change them. Disagreements with the recorded verdicts are filed as separate issues against the relevant campaign.

---

## 10. Honest framing

Per [§13 of #42](42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix), inherited verbatim. The report is the campaign's most direct test of the honest-framing discipline: it is the document Tepper will read first, and if it oversells the campaigns' results, the rest of the work is worth less to him.

**The report is unsolicited by its recipient.** This is the load-bearing framing the deliverable must respect. Tepper has not asked for an interim report; the campaign has produced one because the verification work has surfaced material (three errata findings, six downstream signatures, three load-bearing structural results) that benefits from an author-engageable summary. The report's §1 cover note acknowledges this — *"we have not been asked for this, but the work has reached a point where surfacing it to you is useful; here is what we have, here is what we are honestly uncertain about, here is what would help us most to hear from you."* That framing is the antithesis of a demand letter, a status report, or a deliverable; it is a peer note from a co-author to a senior author, acknowledging the asymmetry honestly.

**Q7 specifically.** Earlier drafts of this plan listed Q7 as *"framework's intent to develop proper-time one-loop dual-QED."* That phrasing reads as asking Tepper to commit to a multi-year research programme via a checkbox response. The §6 listing now downgrades Q7 to *"is proper-time one-loop dual-QED on the framework's medium-term roadmap, or out of scope?"* — a yes/no/maybe-with-context question rather than a commitment ask.

**Crocco §5 substantive-AI compliance, on the report-as-deliverable.** The report is substantive AI end-to-end. The PDF must not be built and committed before the human-acceptance pass (Phase 4) is complete (see §6's revised phase ordering). The earlier draft of this plan ordered Phase 3 (PDF build) before Phase 4 (human review); that ordering would have stripped per-paragraph `<!-- TODO -->` markers from the deliverable before a human reviewed the draft, materially violating Crocco §5's *"the human-acceptance section is explicitly blank for the human to fill in"* rule. The revised ordering — dry-run the build pipeline in Phase 3, run the devil's-advocate self-review + steel-man in Phase 3.5, then human review + final PDF in Phase 4 — keeps the AI markers visible to the human reviewer and ships the PDF only after acceptance.

The three honest caveats explicit in #50 carry over to the report:

1. **The Lamb-shift route is the Bethe (1947) estimate, not a one-loop QED calculation.** The report's §4 records this explicitly when summarising the Lamb-shift verdict; the §6 Q7 question asks Tepper whether the framework intends to develop the apparatus needed to push beyond Bethe-estimate precision.
2. **`r_e`-dependent observables are conditional on the `r_e` finding's resolution.** The report's §4 records both branches (b) and (c); §6 Q1 asks Tepper to resolve which is intended.
3. **The campaign's "passes the precision-spectroscopy test" claim is conditional on branch (c).** The report does not assert that the framework agrees with measurement; it asserts that *if* branch (c) is the intended `r_e`, the framework agrees at the Bethe-estimate precision floor across six independent observables.

These caveats are load-bearing in the report. They are not buried in a closing-disclaimer paragraph; they are recorded at the point in §3 / §4 / §5 where the relevant claim appears.

---

## 11. Decision points

Inherited from [§13.5 of #42](42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24) and [§8 of #49](49-quantum-mechanics-griffiths-proper-time.md#8-decision-points):

| Inherited | Override for #58 | Status |
|---|---|---|
| Per-paragraph TODO for substantive content | same — applies end-to-end (the report is substantive end-to-end) | confirmed |
| Voice = Gill's published-paper voice | same — and load-bearing because the audience IS Gill | confirmed |
| Honest scoping discipline | same — every identity/prediction/out-of-scope claim labelled | confirmed |
| Mechanical deferral on flagged findings | does not engage — the report does not modify or re-derive findings | confirmed |
| AI is never an author | same — Tepper's "author response" recipient is human; report's "Author = Trey Morris with Claude Opus 4.7" is honest pragmatic-AI tooling attribution, not journal authorship | confirmed |

One #58-specific decision: **how many cross-references in §2 to cite by URL vs by relative path.** Relative paths assume the reader has the repository checked out; URL paths assume PDF-only reading. Default: footnote-style URL references for the dozen most-cited paths (so the PDF stands alone), relative paths in the markdown source-of-truth (so internal navigation works on GitHub).

---

## 12. Linked PRs / branches

- Parent issue: [#58](https://github.com/temoTxt/PyPhysics/issues/58)
- Branch: `58-author-review-interim-report-for-tepper-gill-goals-progress-and-open-questions-across-the-verification-campaign-work` (current)
- PR (to be opened after Phase 3 build is operational + Phase 4 human-acceptance pass): closes #58

The branch base is `origin/main`. Implementation requires the campaign content (Bethe-Salpeter, Griffiths, Jackson docs) to be present on this branch for the report to reference local files. By the time this PR is being implemented, [#51](https://github.com/temoTxt/PyPhysics/pull/51) / [#52](https://github.com/temoTxt/PyPhysics/pull/52) / [#53](https://github.com/temoTxt/PyPhysics/pull/53) should have merged to main; if not, the implementer rebases this branch onto the most-recent campaign branch (stacked, as #53 stacks on #52).

---

## 13. Working notes

- The build script is scoped to *this* report's format; generalising it to handle multiple document classes / journal templates is a future thread.
- The `Author_Reports/` folder may host multiple reports over time (`2026-05_interim_for_gill.md`, `2026-XX_followup_for_gill.md`, etc.); the `README.md` documents naming convention (`{YYYY}-{MM}_{report_type}_for_{recipient}.{ext}`) so the folder remains navigable.
- The committed PDF is rebuilt from the markdown source-of-truth on every meaningful change to the `.md`. The build script's **reproducibility scope** (not the same as byte-identical determinism — pandoc + pdflatex output is rarely byte-identical across pandoc/texlive versions, font availability, or `/tmp` paths embedded by some LaTeX packages) is: *rebuildable from the same `.md` on a clean checkout without manual intervention*. The build script pins the PDF metadata date via `--metadata date=YYYY-MM-DD` (no timestamp leakage), documents the pandoc + texlive versions it was developed against in its docstring, and accepts that byte-identical output across machines requires a pinned container that is out of scope for this script. Earlier drafts overclaimed byte-identical determinism; the relaxed requirement is the honest one.

---

## 14. Disclosure: this plan is substantive AI

Per Crocco §5 substantive-AI compliance, this plan document — making architectural decisions about the report's structure, build pipeline, phase ordering, honest-framing discipline, and acceptance criteria — is **substantive AI** (Trey Morris with Claude Opus 4.7, 1M context), not pragmatic AI. The decisions encoded here (Phase 3.5 self-review gate, relaxed `[3, 7]` page budget, PDF-after-human-acceptance ordering, voice-match-as-flag-not-assertion, Q7 downgrade, byte-identical-determinism withdrawal, cross-reference no-drops policy, build-pipeline scoped to this report only) are AI-proposed and require a human-acceptance pass on the plan itself before Phase 1 work begins. The plan's substantive-AI status is symmetric with the report's: both require human acceptance before the deliverable ships.

This disclosure was added after a devil's-advocate self-review (analogous to the workflow that produced commits `393efb4` on #51 and `913802d` on #53). The earlier draft of this plan carried no substantive-AI marker; this section corrects that omission. The committed self-review artifact for the *plan* (the critiques + the steel-man patches) is the commit that introduces this disclosure plus the §3 / §5.4 / §5.6 / §6 / §7 / §10 / §13 edits — they are the load-bearing patches that pulled the plan back into alignment with the §13 honest-framing discipline of #42 that the plan claims (and now actually) inherits verbatim.
