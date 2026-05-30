# Substantive AI prompt-of-record — author-review packet implementation

**Date:** 2026-05-30.
**Issue:** [#93](https://github.com/temoTxt/PyPhysics/issues/93).
**Plan:** [`.dev/tasks/93_Tepper_Review_PDF_BS_Equations_r_e_Findings.md`](../../../../.dev/tasks/93_Tepper_Review_PDF_BS_Equations_r_e_Findings.md).
**Model:** Claude Opus 4.7 (1M context).

## What this file is

Per the [Crocco-compliance policy](../../../Tooling/CROCCO_COMPLIANCE.md): substantive AI use (synthesis, framing, prose drafting) requires a committed prompt-of-record. This file captures the human-authored instructions and the AI work scope for the implementation of issue #93's author-review packet.

## User-authored instruction sequence (verbatim, in temporal order)

1. *(after triage of #93)*: "now implement the plan"

That's the entire user prompt for this implementation. The "plan" referenced is the triaged plan committed at [`.dev/tasks/93_Tepper_Review_PDF_BS_Equations_r_e_Findings.md`](../../../../.dev/tasks/93_Tepper_Review_PDF_BS_Equations_r_e_Findings.md), which Claude generated under the `/triage-ticket 93` skill earlier in the same session, against issue #93's body.

The plan itself reads as a structured instruction with eight numbered implementation steps. The implementation followed it.

## Substantive AI scope in the resulting packet

1. **Selection of which `.wl` cells to surface in `equations/*.tex`.** Not every cell of every notebook is rendered — the cells judged most relevant to the senior author's verdict questions are. Specifically: the cell with the load-bearing result equation, the cell with the Wolfram-MCP confirmation note, and any cell whose interpretation bears on the Branch A vs Branch C question. Other cells (intermediate simplifications, sanity checks) are summarised in prose rather than rendered.

2. **Branch A vs Branch C synthesis** in [`sections/05_branch_a_vs_branch_c.tex`](../sections/05_branch_a_vs_branch_c.tex). The non-contradiction reading — "the algebra of (III.22) is Z-independent (Branch A correct as algebra); the empirical residual at $\chi^2 \sim 10^{16}$ rejects (Z-i) (Branch C correct as empirics); the synthesis is that the framework's apparatus supplies (III.22) and leaves the cutoff free per (III.23), so the per-Z back-fit inherits QED's $a_e^{\rm bound}(Z\alpha)$ structure without the framework deriving the binding coefficient" — is a substantive AI interpretive move. It is explicitly flagged inline in §05 and surfaced in the AI-use disclosure ([§09](../sections/09_ai_use_disclosure_crocco.tex)).

3. **Outcome-matrix mapping** in [`sections/08_open_theoretical_question.tex`](../sections/08_open_theoretical_question.tex): A excluded by Z-scan; B excluded structurally; C the consistent residual; D unlikely without revising (III.22). This is substantive synthesis across PRs #84/#85/#86/#87, the master #67 outcome-matrix issue, and the §III.D-extension reading. Flagged inline.

4. **Selection of the seven closing questions** in [`sections/10_questions_for_tepper.tex`](../sections/10_questions_for_tepper.tex). The questions are drawn from the BS PR bodies (#84/#85/#86/#87), the FINDINGS Finding 2 verdict trajectory, the open theoretical issues (#54, #75), and the original issue #94 (now closed as duplicate of #93). The question selection and ordering are substantive AI.

5. **Prose drafting throughout the sections.** Section 00 (cover note), abstract (§01), and the connective prose in §§02–08 are AI-drafted under the user's explicit instruction. The substance — the 28-result inventory, the honest-scope statement, Finding 2's verdict trajectory, the Z-scan numerics, the mis-provenance audit — is transcribed verbatim from the source documents (`10_CrossComparison.md`, `FINDINGS_for_author_review.md`, `Dual_Relativistic_Quantum_Mechanics_I.md`, PR bodies, `.wl` notebook headers). Where AI-drafted prose makes an interpretive move beyond transcription, it is flagged with `\authorreview{...}` in the LaTeX source.

## Substantive AI NOT present in the packet

- No new physics derivations are introduced.
- No verdict on Branch A vs Branch C is announced; the packet asks the senior author to adjudicate.
- No abstract claim that the packet itself adjudicates the $r_e/r_0$ story.

## Constraints honoured

- **"For author review" framing** preserved per [CLAUDE.md](../../../../CLAUDE.md) — packet asks for the senior author's verdict, does not announce one.
- **Honest-scope statement preserved verbatim** ("Zero of 28 results constitute independent corroboration of the dual theory's content distinct from textbook QM," `10_CrossComparison.md` §1).
- **Finding 2 verdict trajectory preserved verbatim** from `FINDINGS_for_author_review.md:209`.
- **Mathematica → LaTeX faithful transcription**: each load-bearing `.wl` cell is rendered with source line + Wolfram-MCP confirmation marker. No equation content is invented.
- **No fabricated citations.** The packet uses one `\cite{...}` (`nist_mcp_pr91`) which resolves to the merged PR. Any subsequent citation must be scaffolded via `scaffold_bib_note.py` before it appears.

## Human review section

`<!-- TODO: human reviews and fills in — confirms (a) the AI-substantive scope above accurately describes the moves Claude made, (b) the substantive synthesis in §§05 and 08 is correctly framed as AI interpretation (not human position), (c) the closing-question selection in §10 surfaces the right verdict-blockers for senior-author adjudication, and (d) the packet's overall framing as "an author-review vehicle, not a journal submission" matches the user's intent for the implementation. -->`
