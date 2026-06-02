# Task 74: Pre-show conversation topic menu with Tepper Gill — proper time vs standard time

**Tracks:** issue [#74](https://github.com/temoTxt/PyPhysics/issues/74).
**Status:** plan; topic-menu drafting + bib-note scaffolding to follow.
**Dependencies:** verification docs [`The_Classical_Electron_Problem.md`](../../Roadmapping/Equation_Verification/The_Classical_Electron_Problem.md) and [`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) (both on `main`); podcast conventions [`Roadmapping/History/Podcast/README.md`](../../Roadmapping/History/Podcast/README.md); bib-note scaffolding tool [`scaffold_bib_note.py`](../../Roadmapping/History/Bibliography/scaffold_bib_note.py).

This is **not** a multi-PR campaign and **not** a full podcast episode. It is a single topic-menu authoring task that produces one new markdown file under `Roadmapping/History/Podcast/`, plus any History-bucket primary-source bib stubs that don't yet exist. The deliverable is a *prep document* for a future 2-hour recorded conversation with Tepper Gill — the topic menu that the eventual episode(s) will draw from.

---

## 1. Goal and scope

Produce a structured topic menu under `Roadmapping/History/Podcast/preshow_tepper_proper_time.md` that supports a 2-hour pre-show recorded conversation with Tepper Gill, aimed at an amateur-physicist / first-year-grad-student audience, structured around the four buckets in [issue #74](https://github.com/temoTxt/PyPhysics/issues/74): **History** (Newton → Lagrange → Hamilton → ... → Gill), **Society/Culture** (historical figures + Tepper's IAS / Howard experience), **Experiments**, **Thought experiments**.

The format keeps the standing 3-voice Podcast cast (Historian / Physicist / Experimentalist) per [`Roadmapping/History/Podcast/README.md`](../../Roadmapping/History/Podcast/README.md), with **Tepper Gill as guest fourth voice and centerpiece interlocutor**.

### Why this earns its own thread

- **Different doc type from existing episodes.** The nine `episode_NN_*.md` scripts are full ~5,000–7,000-word dialogue scripts targeting one chapter each. This is a *pre-show prep menu* — opening prompts + opener-voice tags + citations, not a script. New convention within the Podcast directory.
- **Source weighting is non-obvious.** The issue explicitly directs the topic mine to favor [`The_Classical_Electron_Problem.md`](../../Roadmapping/Equation_Verification/The_Classical_Electron_Problem.md) and [`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) over the Dirac-equation paper. That weighting needs to live in the deliverable, not be re-derived each future drafting session.
- **History bucket spans 1687 → present.** Newton → Lagrange → Hamilton → Jacobi → Maxwell → Einstein/Minkowski → Dirac → Schwinger → Stueckelberg/Feynman → Gill. Each waypoint needs a primary-source citation; bib stubs that don't already exist must be scaffolded before the menu commits.

### Explicit non-goals

- Not a podcast episode script. The eventual episode(s) drawn from this menu are a separate downstream task.
- Not recording / audio production / TTS pipeline. The `audio/` subdirectory and `lint_episode.py` flow remain untouched.
- Not a topic mine of the Dirac-equation verification doc. That paper is in-scope only as a History waypoint (Dirac 1938) for primary-source citation, not as a topic-content source.
- Not a re-derivation of any verification-doc finding. The menu cites; it does not modify.
- Not a commitment of biographical content about Tepper (IAS / Howard threads). The menu provides *open-ended prompts* with named gaps for Tepper to fill; it does not pre-answer them.

---

## 2. Deliverable structure

`Roadmapping/History/Podcast/preshow_tepper_proper_time.md` — one new file, target ~600–1,000 lines. Structure:

1. **YAML frontmatter** — a new variant on the episode frontmatter shape, with `doc_type: preshow_prep` and `guest: Tepper Gill` fields added; `speakers` retains `[Historian, Physicist, Experimentalist]` plus a new `guest` slot. The frontmatter is documented in the `Roadmapping/History/Podcast/README.md` update (see §5 files-to-modify).
2. **Opening note** (~½ page) — what the doc is, who it's for, how to use it during the live conversation. Acknowledges that the menu is a *prompt bank*, not a script; the regulars are expected to follow Tepper's threads rather than walk the menu top-to-bottom.
3. **Running order** (~¼ page) — proposed sequencing: History (40 min, expanded from the nominal 30) → Society/Culture (30 min, with two interlocking threads) → Experiments (25 min) → Thought experiments (25 min). Conceptual hooks (Maxwell Eq. 1 / Eq. 2 / Eq. 6 / "b is not a speed") woven in across all four; not their own segment.
4. **Bucket 1 — History** (~1½ pages, the longest bucket). Nine waypoints per the issue's spine: Newton (1687) → Maupertuis/Euler/Lagrange (1740s–1788) → Hamilton/Jacobi (1830s–1840s) → Maxwell (1861) → Einstein/Minkowski (1905/1908) → Dirac (1932, 1938) → Schwinger (1951) → Stueckelberg/Feynman → Gill (1990s–present). Each waypoint:
   - 1-paragraph framing
   - Primary-source citation as `[[cite_key]]` (scaffolded if absent — see Phase 1)
   - Verification-doc cross-reference where applicable
   - 2–4 candidate opener questions, each tagged with voice (`[H]` / `[P]` / `[E]`)
   - Plus the two bucket-closers: the "free Hamiltonian looks Newtonian" payoff (Maxwell Eq. 16 / Classical Electron Eq. 3.4) and the Dresden-renormalization word-origin via Maxwell Eq. (17). The Newtonian-Hamiltonian payoff is marked as the **bucket's conceptual climax** — narratively load-bearing.
5. **Bucket 2 — Society/Culture, two threads** (~1¼ pages). Thread A (historical figures, H-led) and Thread B (Tepper's IAS + Howard experience, Tepper-led / regulars cross-examine) printed as separate prompt lists in the same bucket, clearly labeled. Thread B prompts are deliberately **open-ended and biographical** — names of mentors, cafeteria-culture moments, the "what kept you in the program" arc. No pre-answers; gaps left for Tepper.
6. **Bucket 3 — Experiments** (~¾ page). Four prompt clusters per the issue: GPS clock corrections (Maxwell Eq. 9), muon storage-ring / g−2 (Classical Electron Eq. 3.51), classical electron radius as critical point (Maxwell Eq. 18 + cross-reference to active campaign [#54](https://github.com/temoTxt/PyPhysics/issues/54) / [#61](https://github.com/temoTxt/PyPhysics/issues/61) / [#64](https://github.com/temoTxt/PyPhysics/issues/64) / [#65](https://github.com/temoTxt/PyPhysics/issues/65) / [#67](https://github.com/temoTxt/PyPhysics/issues/67)), synchrotron / cyclotron angle-dependent correction (Classical Electron Eq. 3.54 vs. 3.55).
7. **Bucket 4 — Thought experiments** (~¾ page). Five prompt clusters: accelerated radiating twin (Classical Electron Eq. 4.5), preacceleration / Liénard–Wiechert third term (Classical Electron Eq. 3.45–3.46), "soft sphere" electron (Maxwell Eq. 18), group-velocity-adds vs. subtracts (Classical Electron Eq. 4.16 — flag the known sign-typo from [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md)), Newtonian-looking free Hamiltonian (Maxwell Eq. 16 — explicit cross-ref to the History bucket payoff).
8. **Cross-cutting conceptual hooks** (~¼ page). Four hooks per the issue: `u/b = w/c`, `dt/dτ = b/c`, "b is not a speed," dynamical photon mass (Maxwell Eq. 6). Each with a "when to invoke" sentence — these hooks are not their own segment; they get sprinkled when a regular's question naturally surfaces them.
9. **FINDINGS-touching items index** (~¼ page). Explicit list of menu items that touch open entries in [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) — currently the Maxwell Eq. (24) erratum, the DRQM I §III.D `r_e` triangulation, and the TCEP Eq. (4.16) sign typo — so Tepper can react to them in his own words during the conversation rather than discovering them in post-production.
10. **Honest-scoping note** (~¼ page). Per the campaign-wide discipline: identify which menu items are *unconditional* (algebraic identities, established history) vs *conditional* (predictions that depend on the `r_e` resolution from issues [#54](https://github.com/temoTxt/PyPhysics/issues/54) / [#61](https://github.com/temoTxt/PyPhysics/issues/61)). The doc itself does not assert framework-experiment agreement; it presents prompts.

---

## 3. Implementation plan

| Phase | Scope | Status |
|---|---|---|
| Phase 1 | **Bib-note audit + scaffolding.** Walk the nine History-bucket waypoints; for each, check whether the primary-source bib note exists under `Roadmapping/History/Bibliography/{Primary,Retrospective}/`; scaffold any missing stubs via `uv run python Roadmapping/History/Bibliography/scaffold_bib_note.py --cite-key <key>` (with `--from-doi` where DOIs are available — Schwinger 1951, Dirac 1938). No body content is written into the stubs in this phase; the `human_reviewed: false` field stays false. | pending |
| Phase 2 | **Draft `preshow_tepper_proper_time.md`** — §§1–10 per §2 above. Per-paragraph TODO blocks for substantive content per Crocco §5. Topic mine sourced from the [Agent-Explore harvest in the issue-creation conversation](https://github.com/temoTxt/PyPhysics/issues/74) (24 items, weighted toward the Maxwell paper). | pending |
| Phase 3 | **Update [`Roadmapping/History/Podcast/README.md`](../../Roadmapping/History/Podcast/README.md)** with a short "Preshow prep documents" subsection documenting the new `doc_type: preshow_prep` frontmatter variant and the `preshow_*.md` naming convention. One paragraph + the frontmatter example block. | pending |
| Phase 3.5 | **Devil's-advocate self-review** against the campaign's §13 honest-framing discipline. Specifically: verify no biographical content about Tepper is asserted (only prompted), no menu item overstates a framework-experiment agreement, every primary-source citation has a corresponding scaffolded bib note. Critique + steel-man patches written to a *transient* file `preshow_tepper_proper_time_self_review.md` for the Phase 4 human reviewer; deleted before PR open per the [#58 plan](58-author-review-interim-report-for-gill.md#7-acceptance-criteria) precedent. | pending |
| Phase 4 | **Human-acceptance pass.** Trey reviews the menu + the transient self-review, edits prompts toward `VOICE_MATCH_GILL.md` register where appropriate, confirms the History primary-source citations are honest (no Dirac-1938-page-number guesses, no Schwinger-1951-equation-number guesses), deletes the transient self-review file. | requires human review |
| Phase 5 | **PR opened** (closes #74). | pending |

The phasing mirrors the [#58 interim-report plan](58-author-review-interim-report-for-gill.md#6-pr-sequencing) — substantive AI authoring with a devil's-advocate gate before human review — because the same Crocco §5 substantive-AI compliance posture applies (the menu's selection of *which* prompts shape what gets discussed is substantive, not pragmatic).

---

## 4. Files to modify / create

| File | Change |
|---|---|
| `Roadmapping/History/Podcast/preshow_tepper_proper_time.md` | Create — topic menu per §2 above |
| `Roadmapping/History/Podcast/README.md` | Append "Preshow prep documents" subsection documenting the new `doc_type: preshow_prep` frontmatter variant + `preshow_*.md` naming convention |
| `Roadmapping/History/Bibliography/{Primary,Retrospective}/*.md` | Create — scaffold stubs for any History-bucket primary sources that don't already have a bib note (audit happens in Phase 1; expected candidates include Newton's *Principia*, Lagrange's *Mécanique Analytique*, Hamilton 1834/1835, Jacobi 1842, Dirac 1938 proper-time paper, Schwinger 1951 vacuum-polarization paper, Stueckelberg parameterized-paths papers). Stubs are scaffold-only with `human_reviewed: false`; no AI-generated body content. |
| `Roadmapping/History/Podcast/preshow_tepper_proper_time_self_review.md` | Create (transient) — Phase 3.5 devil's-advocate artifact; **deleted before PR open** per Phase 4. Content survives in the Phase 3.5 commit message. |

No verification-doc changes. No source-code changes. No `lint_episode.py` changes (the linter targets `episode_NN_*.md`; the new doc type is out of its scope by design — documented in the README update).

---

## 5. Acceptance criteria

Per [issue #74](https://github.com/temoTxt/PyPhysics/issues/74)'s "Deliverable" section:

1. [ ] `Roadmapping/History/Podcast/preshow_tepper_proper_time.md` exists with the §2 structure and YAML frontmatter (`doc_type: preshow_prep`, `guest: Tepper Gill`, `speakers: [Historian, Physicist, Experimentalist]`).
2. [ ] Every menu item carries the underlying verification-doc citation (file + line range or section heading) **and**, for History-bucket items, the primary-source `[[cite_key]]` wikilink resolving to an existing bib note.
3. [ ] Every menu item carries 2–4 candidate opener questions tagged with voice (`[H]` / `[P]` / `[E]`).
4. [ ] The Society/Culture bucket prints **Thread A (historical figures)** and **Thread B (Tepper's IAS + Howard experience)** as separately-labeled prompt lists; Thread B prompts are open-ended and biographical with named gaps for Tepper (no pre-answered names, dates, or anecdotes).
5. [ ] The running-order section explicitly marks the **Newton → Lagrange → Hamilton → Jacobi → Maxwell → Einstein/Minkowski → Dirac → Schwinger → Stueckelberg/Feynman → Gill** arc as the spine that makes the rest of the conversation land.
6. [ ] The Newtonian-looking-Hamiltonian payoff (Maxwell Eq. 16 / Classical Electron Eq. 3.4) is explicitly cross-referenced between the History bucket and Bucket 4 — the payoff is the conceptual climax of the History bucket and the natural hand-off into the rest of the show.
7. [ ] The "FINDINGS-touching items" index lists every menu item that touches an open entry in [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) so Tepper can react in his own words during the conversation.
8. [ ] [`Roadmapping/History/Podcast/README.md`](../../Roadmapping/History/Podcast/README.md) carries a "Preshow prep documents" subsection documenting `doc_type: preshow_prep` and `preshow_*.md` naming.
9. [ ] All primary-source bib stubs cited as `[[cite_key]]` in the menu exist under `Roadmapping/History/Bibliography/`. No `[[cite_key]]` is fabricated (Crocco §2).
10. [ ] The menu's substantive paragraphs carry per-paragraph `<!-- TODO: human reviews and fills in -->` blocks per Crocco §5; the Phase 3.5 self-review file is created and **deleted before PR open**.
11. [ ] Honest-scoping section labels each menu item as *unconditional* (algebraic identity, established history) vs *conditional* (depends on the `r_e` resolution from [#54](https://github.com/temoTxt/PyPhysics/issues/54) / [#61](https://github.com/temoTxt/PyPhysics/issues/61)).

---

## 6. Out-of-scope

- **The pre-show recording itself.** Audio capture, scheduling with Tepper, post-production are a separate downstream task.
- **The eventual podcast episode(s) drawn from the conversation.** Those follow the existing `episode_NN_*.md` convention and are a separate downstream task.
- **Mining the Dirac-equation verification doc for additional topics.** Per the issue's "Out of scope" — that paper is in-scope only as a History waypoint (Dirac 1938) for primary-source citation.
- **AI-generated biographical content about Tepper.** Per Crocco §2 and the issue's triage notes, names / dates / anecdotes from Tepper's IAS / Howard career are *primary-source* and must not be AI-fabricated. The menu prompts; Tepper answers.
- **Bib-note body content.** Phase 1 scaffolds stubs only; `human_reviewed: false`. Body summaries follow when Trey has read the primary source.
- **`lint_episode.py` changes.** The linter targets `episode_NN_*.md`; the new preshow-prep doc type is out of its scope by design.

---

## 7. Honest framing

Inherits the §13 discipline from [#42](42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix) verbatim. Specifics:

- **The menu prompts, it does not assert.** Every menu item is a prompt for one of the four voices (H / P / E / Tepper); the doc itself does not claim that the framework agrees with experiment, does not pick a branch on the `r_e` question, does not pre-answer Tepper's biographical threads. Where a prompt invokes a conditional prediction (e.g., the GPS-clock-correction prompt, the muon-storage-ring prompt), the conditionality is stated inline rather than tucked into a closing paragraph.
- **Primary-source citations are scaffolded honestly.** Phase 1's bib-stub scaffolding uses `--from-doi` where DOIs exist (Dirac 1938 Proc. Roy. Soc. A 167, Schwinger 1951 Phys. Rev. 82); for older works without DOIs (Newton 1687, Lagrange 1788), the scaffold relies on standard bibliographic data. No equation numbers, no page numbers, no paragraph quotes are written into the bib bodies until Trey has read the primary source — `human_reviewed: false` is honest, not a placeholder.
- **The "Newtonian-Hamiltonian payoff" framing is load-bearing.** The conceptual climax of the History bucket — that the proper-time free Hamiltonian `K = mu²/2 + mc²` (Maxwell Eq. 16 / Classical Electron Eq. 3.4) recovers a Newtonian form — is presented as Tepper's pedagogical hook to underline, not as the menu's editorial claim about what the dual theory "really is." The prompt sets him up; he chooses how strongly to lean.
- **The biographical thread acknowledges asymmetry.** Tepper's lived experience at IAS and Howard is primary-source material that the menu cannot generate. The prompts name the rooms and roles (cafeteria culture, math-vs-physics-department dynamics, conferences that welcomed vs sidelined the work) but leave the specifics blank — the menu's job is to give the regulars opening questions sharp enough that Tepper's answers carry the segment.

---

## 8. Decision points

Inherited from earlier campaign plans:

| Inherited | Override for #74 | Status |
|---|---|---|
| Per-paragraph TODO for substantive content | same — applies end-to-end (menu authoring is substantive) | confirmed |
| Voice = `VOICE_MATCH_GILL.md` register | **not applicable** — the deliverable is a topic menu, not prose addressed to Tepper. The regulars' voices follow the existing Podcast cast conventions (Historian / Physicist / Experimentalist personas in [`Roadmapping/History/Podcast/README.md`](../../Roadmapping/History/Podcast/README.md)). | overridden |
| Honest scoping — every claim labelled unconditional vs conditional | same — applies per menu item, recorded in §11 honest-scoping note | confirmed |
| AI is never an author | same — Tepper is the guest, the regulars are personas, Claude is tooling | confirmed |
| Mechanical deferral on flagged findings | same — menu cites FINDINGS, does not re-derive or modify | confirmed |

One #74-specific decision: **whether to commit the menu before the recording happens, or hold it until post-recording for revision.** Default: **commit pre-recording**, treat the menu as a prep artifact whose value is the pre-conversation alignment. Post-recording, the eventual episode(s) drafted from the conversation are a separate task that can cite the menu but does not modify it. Rationale: the menu's value is destroyed if it is rewritten after the fact to match what was actually said — that turns the menu into a transcript, which is not what it is.

---

## 9. Linked PRs / branches

- Parent issue: [#74](https://github.com/temoTxt/PyPhysics/issues/74)
- Branch: `74-podcast-pre-show-conversation-topics-with-tepper-gill-2-hour-conceptual-tour-of-proper-time-vs-standard-time-cep-maxwell-weighted` (current)
- PR (to be opened after Phase 4 human-acceptance pass): closes #74

The branch base is `main`. No sequencing dependencies on other open PRs — the two source verification docs are already on `main`. The active `r_e` triangulation campaign ([#54](https://github.com/temoTxt/PyPhysics/issues/54) / [#61](https://github.com/temoTxt/PyPhysics/issues/61) / [#64](https://github.com/temoTxt/PyPhysics/issues/64) / [#65](https://github.com/temoTxt/PyPhysics/issues/65) / [#67](https://github.com/temoTxt/PyPhysics/issues/67)) is *referenced* by the Experiments bucket's classical-electron-radius prompt but does not block this work — the prompt is "Tepper's intuition for why r₀ is special," which is a conversation-opener whether or not the campaign has converged.

---

## 10. Disclosure: this plan is substantive AI

Per Crocco §5, this plan document — making architectural decisions about which topics make the menu, how the History bucket spine is structured (Newton → ... → Gill), how the Society/Culture bucket splits into Thread A / Thread B, what the phase ordering looks like, and what the acceptance criteria are — is **substantive AI** (Trey Morris with Claude Opus 4.7, 1M context), not pragmatic AI. The decisions encoded here require a human-acceptance pass on the plan itself before Phase 1 work begins. The plan's substantive-AI status is symmetric with the menu's: both require human acceptance before the deliverable ships.

The plan inherits the substantive-AI disclosure precedent from [§14 of the #58 plan](58-author-review-interim-report-for-gill.md#14-disclosure-this-plan-is-substantive-ai).
