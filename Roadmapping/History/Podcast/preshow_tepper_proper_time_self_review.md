# Devil's-advocate self-review — `preshow_tepper_proper_time.md`

**Phase 3.5 artifact** of [`.dev/tasks/74-podcast-preshow-topics-tepper-proper-time.md`](../../../.dev/tasks/74-podcast-preshow-topics-tepper-proper-time.md). Transient file — *to be deleted before PR open*. Content preserved in the Phase 3.5 commit message body so git history retains the record without cluttering `Podcast/`.

The review pass tested the draft menu against the §13 honest-framing discipline of [#42](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix) (inherited via the task plan's §7), the Crocco §2 no-fabricated-citations rule, and the issue's acceptance criteria. Two patches were applied; six concerns are surfaced for Phase 4 human consideration.

---

## Patches applied during Phase 3.5

### Patch 1 — Finding 1 attribution corrected

**Where:** §8 FINDINGS-touching items index.

**Before:** "The Maxwell Eq. (24) erratum (**Finding 1**) is not on the menu — it lives in the Dirac-equation verification doc and the issue explicitly defers Dirac-paper topic mining."

**Problem:** Wrong attribution. Finding 1 is in the Maxwell paper (Eq. 24 of `Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`), not in the Dirac-equation verification doc. The "defers Dirac-paper topic mining" justification doesn't apply — the Maxwell paper is one of the two source-weighted docs the issue *centers*.

**After:** Patched to attribute Finding 1 correctly to the Maxwell paper, with a pointer to §3.7 Schwinger / Lamb-shift cluster as the natural slot where it would come up (V²/(2mc²) is the proper-time Lamb-shift contribution Schwinger's regulator would handle in standard QED). The menu still does not pre-empt Tepper's framing with a dedicated prompt — that judgment call is left for Phase 4.

### Patch 2 — Ashby 2003 hedge removed

**Where:** §5.1 GPS clock corrections.

**Before:** "Cross-reference: [`ashby2003_gps_relativity`](../Bibliography/Retrospective/) (if exists; this is Neil Ashby's standard GR-corrections review)."

**Problem:** The "(if exists; ...)" hedge was authored before bib-stub audit. `ashby2003_gps_relativity.md` does exist in `Bibliography/Primary/` (note: Primary, not Retrospective — Ashby's paper is contemporary research, not retrospective historical review).

**After:** Patched to clean `[[ashby2003_gps_relativity]]` wikilink with no hedge.

---

## Concerns surfaced for Phase 4 human consideration

### Concern 1 — Should Finding 1 (Maxwell Eq. 24 typo) get a dedicated prompt?

Patch 1 fixes the attribution but leaves the underlying question open: does Tepper need a direct prompt about the missing `c` in `eℏΣ·B/(2m)` and the missing `V²/(2mc²)` term, or is leaving it to organic surfacing in §3.7 Schwinger / Lamb-shift the right call?

**Argument for adding a prompt:** Findings 2 and 3 are flagged as FINDINGS-touching prompts (§5.3, §6.4). Symmetry suggests Finding 1 should too — otherwise Tepper might not realize the menu is *aware* of the typo when he's narrating proper-time Lamb-shift physics in §3.7.

**Argument against:** Finding 1 is a transcription erratum, not a structural physics result. Asking Tepper "did you mean `/(2mc)` or `/(2m)`?" is a typo-confirmation question more suited to the formal author-review report ([Roadmapping/Author_Reports/2026-05_interim_for_gill.md](../../Author_Reports/2026-05_interim_for_gill.md)) than to a pre-show prep menu about *concepts*.

**Recommendation for Phase 4:** Leave as-is. Trey can add a one-line "Q: 2m vs 2mc, was that a proofs-stage typo?" parenthetical in §3.7 during the human pass if the recording session feels like it has the airtime.

### Concern 2 — Bucket 4 §6.4 prompt double-asks the FINDINGS sign question

§6.4 (group velocity adds, not subtracts) has two opener prompts. The `[P]` prompt asks for the conceptual / pedagogical walkthrough; the `[H]` prompt asks for the proofs-stage typo confirmation. The second prompt is a sub-instance of Concern 1 — it duplicates the "ask about a printed typo" question that Concern 1 advises against having dedicated prompts for elsewhere.

**Asymmetry:** Concern 1 advises *not* making Finding 1 a dedicated prompt. §6.4 *does* make Finding 3 a dedicated prompt. Is the asymmetry honest?

**Defense:** Finding 3 (sign typo) bears on a *thought-experiment-level* claim about causality (vg = vg′ + v not breaking causality). The signage matters to the *story* in a way Eq. 24's prefactor does not — a reader following the thought experiment along can't act on it without knowing which sign is intended. Finding 1's prefactor doesn't affect any conceptual narrative — it's just a transcription correction.

**Recommendation for Phase 4:** Asymmetry is defensible. Leave §6.4 as-is. If Trey disagrees, the patch is one-line: collapse §6.4's two prompts into one and move the typo-confirmation question to a unified "errata pass" section.

### Concern 3 — §3.9 Gill section has no `[[cite_key]]`

The History bucket's spine has 9 waypoints. Waypoints §3.1–§3.8 each carry one or more `[[cite_key]]` wikilinks. §3.9 (Gill 1990s–present) does not — it points only at the verification docs, with the in-body note that "the framework papers themselves are in the Gill verification corpus, not in the History bibliography."

**Asymmetry honest?** Yes — the History `Bibliography/` tree is for the *historical-canon* corpus (1687 → 1948), not for the Gill framework papers. The Gill papers live in their own verification docs because they are the *subject* of the verification campaign, not historical context. The bib stubs for `gill1990s_*` would be circular — citing them as historical context for the verification of *themselves*.

**Defense:** Acceptance criterion 9 says "no `[[cite_key]]` is fabricated." Not having a Gill cite-key is honest — fabricating one would violate the criterion. The pointer to the verification docs is the right level of indirection.

**Recommendation for Phase 4:** Leave as-is. Confirm with Trey that this asymmetry is intentional — the convention is "Bibliography is for historical context, verification docs are for framework primary sources" — and consider documenting that convention in the [Bibliography README](../Bibliography/README.md) if it isn't already.

### Concern 4 — §4.2 Thread B has named gaps but no explicit "if you don't remember, that's fine" framing

The biographical thread asks for "the load-bearing collaborators — the ones whose names show up in the acknowledgments most often" and "the Howard student who walked into your office and pushed the dual-theory program in a direction you hadn't planned." These are good open-ended prompts, but they don't acknowledge the asymmetry of expecting Tepper to remember specific names on the spot.

**Risk:** A prompt that asks "name the load-bearing collaborators" without an off-ramp can put a guest on the spot in a way that produces silence rather than a story.

**Recommendation for Phase 4:** Add a one-line preamble at the top of §4.2 like *"These prompts are openers; Tepper, if a question doesn't land or you'd rather circle back to it, just say so — the regulars will pick it back up later in the segment."* This is standard interviewer-asymmetry hygiene and isn't currently in the menu.

### Concern 5 — Conceptual hook §7.4 (dynamical photon mass) is borderline-speculative

The hook is grounded in Maxwell Eq. (6) and the verification doc cites a specific algebraic identity for µ². But the framing "Is that a real physical claim, or a bookkeeping artifact of the proper-time rewriting?" invites Tepper toward a *speculative* answer rather than an *unconditional* algebraic one. The hook is marked unconditional in §9 honest-scoping, but the prompt phrasing nudges conditional.

**Tension:** Algebraic identity (unconditional) vs. physical interpretation (conditional / speculative).

**Recommendation for Phase 4:** Either (a) reword §7.4's prompt to ask Tepper to explain the algebra first, *then* the physical reading, or (b) move §7.4 from "unconditional" to "borderline / conditional" in §9. (a) is the more honest fix.

### Concern 6 — Voice tagging undercounts the Experimentalist in History

The History bucket has 9 waypoints; the Experimentalist gets `[E]` opener prompts in 6 of them. The Society/Culture and Thought-experiments buckets have far fewer `[E]` openers. This is partly correct (history is naturally `[H]`-led and thought experiments are naturally `[P]`-led), but the Experimentalist's role in keeping the conversation grounded was framed as load-bearing in the running-order (§2) and the bucket-1 lede.

**Risk:** If the Experimentalist gets crowded out of bucket-3, the apparatus-grounding role might land thinner than the menu intends.

**Recommendation for Phase 4:** Audit the §5 Experiments bucket — the `[E]` openers there are sometimes paired with `[P]` openers in a way that lets the Physicist drive the segment. Consider re-tagging some §5 prompts as `[E]`-primary, `[P]`-secondary to honor the Experimentalist-leads convention in the buckets where they should lead.

---

## Acceptance-criteria checklist (against task plan §5)

| # | Criterion | Status |
|---|---|---|
| 1 | `preshow_tepper_proper_time.md` exists with §2 structure + YAML frontmatter | ✓ |
| 2 | Every menu item carries verification-doc citation + History primary-source `[[cite_key]]` | ✓ (§3.9 carries verification-doc pointer only — see Concern 3) |
| 3 | 2–4 candidate opener questions per item tagged with voice (`[H]`/`[P]`/`[E]`) | ✓ |
| 4 | Society/Culture prints Thread A + Thread B with biographical gaps left for Tepper | ✓ |
| 5 | Running-order marks Newton → ... → Gill arc as the spine | ✓ |
| 6 | Newtonian-Hamiltonian payoff cross-referenced between §3.10 Closer 1 and §6.5 | ✓ |
| 7 | FINDINGS-touching items indexed in §8 | ✓ (Finding 1 attribution patched during Phase 3.5) |
| 8 | Podcast README has "Preshow prep documents" subsection | ✓ |
| 9 | All `[[cite_key]]` resolve to existing bib notes — none fabricated | ✓ (12 stubs scaffolded in Phase 1; Maxwell/Einstein/Minkowski pre-existed) |
| 10 | Per-paragraph TODO blocks; Phase 3.5 self-review file created (this file) and to-be-deleted-before-PR | ✓ |
| 11 | Honest-scoping §9 labels unconditional vs conditional vs FINDINGS-touching | ✓ |

All 11 criteria meet at Phase 3.5 close. Concerns 1–6 are surfaced for Phase 4 human-acceptance judgment but do not block Phase 4 entry.

---

## Recommended Phase 4 actions

1. Read the menu end-to-end as if you were the Physicist trying to drive a recorded conversation from it. Where do you stop because a prompt isn't quite right?
2. Address Concern 4 (interviewer-asymmetry hygiene in §4.2) — recommended as a one-line preamble add.
3. Address Concern 5 (§7.4 prompt phrasing) — recommended as a (a)-route rewording.
4. Audit Concern 6 (`[E]` distribution in §5) and re-tag if the Experimentalist is being crowded out.
5. Decide on Concerns 1, 2, 3 — current recommendations are "leave as-is" but Trey owns the call.
6. Confirm the History primary-source bib stubs have honest enough metadata for Phase 4 use (the no-DOI stubs — Newton, Maupertuis, Euler, Lagrange, Jacobi, Stueckelberg — are bare scaffolds with `human_reviewed: false`; that's correct per the plan but it means the wikilinks render as titles only).
7. Delete this file (`preshow_tepper_proper_time_self_review.md`) before opening the PR.
