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
