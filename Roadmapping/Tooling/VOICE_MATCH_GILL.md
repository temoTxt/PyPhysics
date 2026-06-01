# Voice and style: matching Gill's published prose

This document records the writing-voice constraint for substantive physics prose in this repository. It applies to the Electromagnetism / Jackson campaign (issue [#42](https://github.com/temoTxt/PyPhysics/issues/42)), to any future per-problem documents continuing that thread, and to the Physicist-voice lines in podcast scripts when they convey dual-theory derivations.

Source for the voice characterisation: direct sampling of Tepper Gill's published papers in [`Roadmapping/Converted_Markdown/`](../Converted_Markdown/), specifically [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], [[Dual_Relativistic_Quantum_Mechanics_I]], [[FoundationsII-Classical]], and [[The_Classical_Electron_Problem]].

---

## 1. Why this constraint exists

Gill is the framework's co-author. The Electromagnetism campaign documents are continuous with his corpus — derivational work *in* the framework, not commentary *about* it. Voice continuity signals that framing. Drift into a third-party voice (Claude's default register, a textbook-summary register, or a journalistic register) undermines the framing and is also one of the easiest tells of AI-generated text even when the underlying derivations are correct. Crocco-compliance-adjacent: see [`CROCCO_COMPLIANCE.md`](CROCCO_COMPLIANCE.md) §6 on the substantive/pragmatic distinction; voice-drift in substantive prose is one of the failure modes that disclosure cannot retroactively fix.

## 2. Scope — where this applies

| Document type | Match Gill's voice? |
|---|---|
| Per-problem documents under `Roadmapping/Electromagnetism/Jackson/Ch*.md` | Yes |
| `Roadmapping/Electromagnetism/_proper_time_cheatsheet.md` | Yes |
| Companion `.wl` notebook header comments | Yes |
| Podcast-script Physicist-voice lines that convey dual-theory derivations | Yes |
| `Roadmapping/Electromagnetism/README.md` status tables and headers | Short-form, but adopt Gill's diction (no hyperbole) |
| Equation Verification per-equation documents | Already in this voice; maintain it |
| **`Roadmapping/Author_Reports/*.md`, `.tex`, `.pdf`** (per [#115](https://github.com/temoTxt/PyPhysics/issues/115)) | **Yes — and §3.bis + §3.ter below are load-bearing here** |
| **`Roadmapping/Quantum_Mechanics/Bethe_Salpeter/*.md`** result documents | **Yes — §3.bis + §3.ter apply** |
| GitHub issue bodies, PR descriptions | No — operational tone |
| `.dev/tasks/` planning documents | No |
| Devil's-advocate / honest-limitations sections | No — plain operational tone preserves their sharpness |
| TODOs and internal notes | No |

## 3. Voice markers to use

Drawn from sampling Gill's prose:

- **First person.** "We" as plural of modesty throughout. Never "I". Never address the reader as "you"; use "one can construct", "it is easy to show", "one observes".
- **Connective sentence openers.** "Thus,", "Indeed,", "Furthermore,", "However,", "On the other hand,", "It is clear that…", "It is easy to show that…", "We observe that…", "It follows that…", "In order to see…".
- **Sentence rhythm.** Long subordinate-clause sentences linked by semicolons; not fragmentary punchy prose.
- **Attribution.** Narrative, by name and year, embedded in the prose ("Wheeler and Feynman [13] showed that…", "Liénard, in 1898, computed…"). Footnote-style citations supplement but do not replace.
- **Equation-to-interpretation hand-off.** *Separate sentence* after the equation, never inline. Name the physical mechanism (e.g., "radiation reaction") with attribution to prior work.
- **Reformative framing of classical theory.** "Reformulation", "alternative", "extension", "parallel image", "recovers lost physical insight". Acknowledge the success of the classical formulation before contrasting.
- **Hedges.** "To the extent that", "it is not obvious that", "under the assumption that", "if one accepts", "it is not clear how", "one should be able to show".
- **Conceptual anchor phrase.** "*Mathematically equivalent but not physically equivalent*" is a Gill signature; reuse it verbatim whenever the contrast between classical and proper-time formulations is being introduced.
- **Parenthetical asides** for clarification: "(in c.g.s. units)", "(what we mean by) physically equivalent". Used sparingly but characteristically.

## 3.bis. Equation discipline (added per [#115](https://github.com/temoTxt/PyPhysics/issues/115))

Every mathematical expression that conveys a result or carries a label is set as a **numbered display equation**, not inline. The numbering scheme follows the existing per-chapter / per-section convention (e.g. `Eq. (5.7)` in verification docs, `(\theequation)` in LaTeX manuscripts).

Inline math is permitted only for **single-symbol references** to a quantity already established elsewhere, for example "the cutoff `r_e`" or "the factor `(1 + V/mc²)`". Multi-term expressions, derivatives, integrals, and any equation a reader might want to cite by number must be display-set and numbered.

| Forbidden | Required |
|---|---|
| "…the proper-time Hamiltonian `K = H²/(2mc²) + mc²/2`…" | "…the proper-time Hamiltonian is\\n\\n    K = H²/(2mc²) + mc²/2,    (5.7)\\n\\nwhich, by direct algebra, …" |
| "…differentiating gives `dp/dt = -∇V(1 + V/mc²)`…" | "…differentiating gives\\n\\n    dp/dt = -∇V(1 + V/mc²),    (2.11)\\n\\nwhich vanishes at `r = r_0`." |

The point of the discipline is that a reader skimming for results can find them. An equation buried in a sentence is invisible. A numbered display equation can be cross-linked, cited, and returned to.

## 3.ter. Sentence-length discipline (added per [#115](https://github.com/temoTxt/PyPhysics/issues/115))

Every sentence in substantive prose must be **short and complete**.

- **Short**: one or two clauses, rarely more. A sentence carrying three or more subordinate clauses is almost always doing the work of two sentences and should be split.
- **Complete**: subject, verb, and the load-bearing claim. No sentence fragments. No telegraphic shorthand.
- **One claim per sentence**: if a sentence makes claim A and claim B, split into two sentences. A numbered list is the right structure only when the enumeration is genuinely list-shaped, not when it is a paragraph in disguise.

The discipline is not about choppy fragments. It is about making every sentence end with the reader knowing what claim was just made. Gill's long subordinate-clause rhythm is preserved when each clause is itself short and the chain points at a single conclusion.

**The mathematics does most of the talking.** In a substantive prose document, the load-bearing technical claims live in the numbered display equations of §3.bis. The role of the surrounding sentences is to *illustrate* the equation — to name the symbols, to flag the limit, to point at the next step. Sentences that paraphrase what an equation already says, or that try to *replace* an equation with English, are doing the wrong job. A reader who skims the equations alone should still come away with the document's technical narrative; the prose layer adds context and connective tissue, not redundant restatement.

A prose sentence earns its place if removing it would leave the equations less interpretable. A prose sentence that could be cut without losing meaning should be cut.

| Forbidden ("AI dribble") | Required |
|---|---|
| "We observe that, given the proper-time Hamiltonian formulation, which differs from the standard Hamiltonian formulation by the algebraic rearrangement of the kinetic-energy term, the resulting equations of motion, when expanded to the appropriate order in `(u/c)²`, reproduce the textbook Sommerfeld–Dirac fine-structure result up to a multiplicative factor that depends on the chosen value of the cutoff `r_e/r_0`." | "We observe that the proper-time Hamiltonian formulation differs from the standard one by an algebraic rearrangement of the kinetic-energy term. Expanding the resulting equations of motion to order `(u/c)²` reproduces the Sommerfeld–Dirac fine-structure result. The reproduction agrees up to a multiplicative factor depending on the cutoff `r_e/r_0`." |

The forbidden example is one 65-word sentence carrying four distinct claims. The required example is three sentences, each carrying one claim, each readable in a single breath.

## 3.quater. Citation discipline (added per [#117](https://github.com/temoTxt/PyPhysics/issues/117))

Every cross-reference is a **numbered footnote** in the prose and is resolved in a **bibliography section at the end of the document**. This applies to references to repository files, to GitHub issues or pull requests, to external papers, and to verification documents. No inline links. No `#NN` tags scattered through paragraphs. No verbose file paths embedded in sentences.

The reasoning: when citations are footnoted, the prose carries the argument on its own. The reader who skims for the claim sees the claim; the reader who wants the source consults the references list. Gill's published papers follow the same pattern.

| Forbidden | Required |
|---|---|
| "The candidates note ([`Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md`](../Author_Reports/2026-05_re_derivation_candidates_for_gill.md)) found, in [issue #65](https://github.com/temoTxt/PyPhysics/issues/65) under master [#67](https://github.com/temoTxt/PyPhysics/issues/67), that the published §III.D NR-expansion is structurally inadequate." | "The candidates note[^1] found, in the §III.D investigation[^2], that the published NR-expansion is structurally inadequate.\\n\\n## References\\n[^1]: \`Roadmapping/Author_Reports/2026-05_re_derivation_candidates_for_gill.md\`.\\n[^2]: Issue #65 under master #67 (\`github.com/temoTxt/PyPhysics/issues/65\`)." |

### Bibliography section conventions

The bibliography section at the end of each document:

- Is titled `## References` or `## Bibliography`. Pick one and use it consistently within a document.
- Lists footnoted items in the order they appear, numbered to match the footnote markers.
- For repository files: bare path (`Roadmapping/...`), and an optional one-sentence description.
- For GitHub issues or pull requests: bare number plus the bare URL (`Issue #65, github.com/temoTxt/PyPhysics/issues/65`).
- For primary-source papers: the existing `cite_key` from `Roadmapping/History/Bibliography/`, optionally with the full citation expanded.
- For external URLs: bare URL, no rich-text linkification.

Pandoc renders `[^1]` as a real LaTeX footnote, so the discipline survives the `.md` → `.tex` → `.pdf` build pipeline without modification.

## 3.quinquies. No emojis or pseudo-emoji symbols (added per [#120](https://github.com/temoTxt/PyPhysics/issues/120))

Emojis must not appear in any committed `.md`, `.tex`, or `.pdf` release governed by §2. Verdict markers, bullet decorators, and any Unicode glyph in the Emoji block of Unicode are all out.

The reason is concrete. pdflatex cannot render Unicode emoji glyphs in its default T1 encoding. The build pipeline substitutes them to ASCII tags (`[OK]`, `[warn]`, `[X]`). Both renderings read as visual noise. A reader sees a tag, not a claim. The ASCII fallbacks are therefore also out — the rule is on emojis *and* their substitutes.

This rule does **not** restrict:

- Mathematical symbols (≈, ×, −, →, ∂, ∇, ℏ, …). These are the language of physics and pandoc handles them natively in math mode.
- Greek letters in math mode (α, β, γ, …). Same.

| Forbidden | Required |
|---|---|
| "Acceptance verdicts: ✅ predicts residual; ⚠ matches at floor only; 🔴 framework refuted." | "Acceptance verdicts: **Pass** — the framework predicts the residual within the campaign's precision. **Marginal** — the framework matches at precision floor only (consistent but not discriminating). **Refuted** — the framework's predicted signature is absent at the predicted sensitivity." |

The recommended substitution is a short bold label, in line with §3 voice markers and §3.ter sentence-length discipline. The label is a recommendation, not a binding template; authors may choose other text labels that fit the document.

## 4. Anti-patterns

Phrases and rhetorical moves that **do not appear** in Gill's prose and must be removed when content is elevated into per-problem documents:

- **Hyperbole.** "Headline payoff", "smoking gun", "dissolves a 120-year-old pathology", "the framework's superiority", "campaign's emotional arc".
- **Casual or idiomatic English.** "let's", "the cool thing about", "kind of", "pretty much", "gotcha".
- **Markdown-heavy bullet structures** for substantive content; Gill writes paragraphs.
- **Bold / italic for rhetorical force.** Reserve emphasis for genuine technical terms.
- **"Claude voice" tells.** "I'll", "let me", "we should consider", "interestingly", "it's worth noting".

## 5. Compliance check

Before any per-problem document is committed, the author (or Claude, when authoring) reads each prose paragraph aloud. If a paragraph could plausibly appear in a Gill paper — same register, same hedging, same connective rhythm — it passes. If it sounds like a third-party summary of a Gill paper, it is revised.

When in doubt, find the closest passage in Gill's corpus that does the same rhetorical work (transition, equation-to-interpretation hand-off, attribution) and pattern-match against it.

## 6. Representative voice fingerprints

Direct quotes from Gill's corpus that exemplify the target voice:

1. *"It is known that Maxwell's electrodynamics as usually understood at the present time — when applied to moving bodies, leads to asymmetries which do not appear to be inherent in the phenomena."* — opening with received wisdom and historical framing.
2. *"However, we must be careful that we don't use the description of what these waves are composed of (i.e., solutions of the wave equation) as an interpretation of how they travel."* — measured caution, parenthetical clarification.
3. *"We observe that the representation `dτ = (Mc²/H)dt` does not depend on the number of particles in the system."* — collective observation, mathematical certainty.
4. *"Thus, the collaborative use of the observer's coordinate system and the local clock of the observed system provides intrinsic information about the field dynamics."* — affirmative clause with causal weight.
5. *"In order to see their impact on Maxwell's equations, write them (in c.g.s. units)…"* — pedagogical framing, directive without imperative.
6. *"This result is not surprising given the close relation between the two groups."* — reassurance of logical coherence.
7. *"It follows that no consideration of the action of a particle on itself or the problematic Lorentz–Dirac equation is required."* — logical closure, semantic distance from emotion.
8. *"Although the form of the Dirac equation serves to make space and time appear on an equal footing mathematically, it is clear that they are still not on an equal footing from a physical point of view."* — balanced concession, then assertion.
9. *"These observations imply that the Dirac Hamiltonian and the square-root Hamiltonian are mathematically, but not (what we mean by) physically equivalent."* — parenthetical definition, technical precision.

When authoring new prose, the test is: would a passage of this length and rhetorical purpose, written by Gill, look like one of the above? If yes, pass.

## 7. Status

- **2026-05-23.** Original voice characterisation written by sampling four of Gill's papers (see source list above). Initially stored as Claude memory; relocated into the repo on the same date so the constraint is version-controlled and reviewable.
- **Open:** at PR A of issue #42, this document should be reviewed against the prose actually written for Ch. 6 problems. If the prose passes the §5 read-aloud test but the markers in §3 are not the load-bearing ones, this document is updated accordingly.

---

Related:
- [`CROCCO_COMPLIANCE.md`](CROCCO_COMPLIANCE.md) — the broader AI-use compliance framework.
- [`.dev/tasks/42-electromagnetism-jackson-proper-time.md`](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md) §14 — campaign-plan section that points back here.
- [`Roadmapping/History/Podcast/README.md`](../History/Podcast/README.md) — podcast persona style guidance (Historian, Experimentalist voices have their own conventions; only Physicist voice is bound by this document).
