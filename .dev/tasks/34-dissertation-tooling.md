# Implementation plan: continuous-crawl → human-triage → continuous-synthesis tool suite

**Tracks:** issue [#34](https://github.com/temoTxt/PyPhysics/issues/34) (dissertation tooling survey, now merged).
**Status:** plan only. No code committed yet. Execution should be staged across multiple focused PRs, not one monolithic PR.

## 1. Goal and scope

Build an end-to-end pipeline so the repo's content (chapters, bibliography, animations) **stays current automatically** as new physics literature appears, while keeping the human in the loop on triage and synthesis decisions. Three pipelines feeding each other:

```
┌─────────────────────┐    ┌─────────────────────┐    ┌──────────────────────┐
│  Pipeline 1: CRAWL  │ →  │ Pipeline 2: TRIAGE  │ →  │ Pipeline 3: SYNTHESIZE│
│  (autonomous)       │    │ (human-in-the-loop) │    │  (Claude-assisted)   │
│                     │    │                     │    │                      │
│  Daily: S2 + ArXiv  │    │  Weekly digest of   │    │  Bib stubs, chapter  │
│  + Crossref + Zotero│    │  candidates with    │    │  cross-refs, claim   │
│  → candidate queue  │    │  why-surfaced score │    │  clustering, mss diff│
└─────────────────────┘    └─────────────────────┘    └──────────────────────┘
         ↑                          ↑                            ↑
         │                          │                            │
    cron / /loop              human reviews              /loop /synthesize
```

The end state: Trey opens Claude Code on Monday morning and finds a fresh triage queue waiting; flicks through it in 20 minutes; commits the kept entries to Zotero; Claude then takes 10 minutes to scaffold bib stubs, propose cross-refs, and identify which chapter sections need revision. Compounded weekly, the dissertation/book stays current without anyone scheduling "literature review days."

**Explicit non-goals:**
- Auto-publishing or auto-committing chapter content. Synthesis suggestions land as PRs Trey reviews; nothing merges without human approval.
- Replacing the existing YAML bibliography or chapter workflow. Zotero is *added* as a capture-and-store layer; YAML stays canonical.
- A web UI. Everything is CLI + markdown files + Claude Code + Obsidian.
- Real-time / sub-hour latency. Daily crawl + weekly triage is the design target.

## 2. Architecture overview

### Three pipelines

| Pipeline | Trigger | Frequency | Output |
|---|---|---|---|
| **Crawl** | `/loop 1d /crawl-papers` (or cron) | daily | `Roadmapping/Tooling/triage/queue.md` — candidates with provenance |
| **Triage** | `/loop 1w /triage-papers` + human | weekly | Updated `queue.md` with `decision: keep|skip|later`; "keep" entries pushed to Zotero |
| **Synthesize** | `/loop 1w /synthesize-claims` | weekly | Per-chapter "what changed" markdown + new bib stubs (PR draft) |

### Component inventory

**Existing (no work needed; just leveraged):**

- YAML bibliography pipeline (`Roadmapping/History/Bibliography/*.py`, 115 stubs).
- Chapter-level cross-reference structure (Obsidian wikilinks, validated by `validate_wikilinks.py`).
- Acquisition tracker + Dataview indexes (auto-regenerated).
- `scaffold_bib_note.py --from-doi` (Crossref lookup).
- `parse_papers.py` (PDF → markdown via marker-pdf).
- 9 historical / forward chapters with `era` / `tags` / `verification_anchors` frontmatter.

**New tools to write (~1500 lines Python total, split across 5 PRs):**

| Component | Path | Purpose | Est. lines |
|---|---|---|---|
| `crawl/from_s2.py` | `Roadmapping/Tooling/crawl/` | Semantic Scholar recommendations from existing library | 200 |
| `crawl/from_arxiv.py` | same | Daily ArXiv listings filtered by physics + tags | 180 |
| `crawl/from_crossref.py` | same | New papers citing entries already in our bib | 150 |
| `crawl/from_zotero.py` | same | Newly-added Zotero entries that aren't yet in YAML | 120 |
| `triage/build_queue.py` | `Roadmapping/Tooling/triage/` | Merge crawl outputs; score relevance; produce digest | 250 |
| `triage/digest.py` | same | Render the human-facing triage digest | 100 |
| `synth/cluster_claims.py` | `Roadmapping/Tooling/synth/` | Group related `#inferred`/`#speculative` claims across chapters | 180 |
| `synth/suggest_cross_refs.py` | same | Find unlinked references between chapters | 150 |
| `synth/manuscript_diff.py` | same | Per-chapter "what changed since last build" | 120 |
| `pyphysics-mcp` | `mcp_servers/pyphysics_mcp/` | Custom MCP exposing the repo pipeline as Claude tools | 400 |
| Sync script | `Roadmapping/History/Bibliography/sync_from_zotero.py` | Promised in earlier Q&A | 150 |

**3rd-party apps to install (manual; document in install guide):**

- Zotero 7 + Better BibTeX + Zotero Connector (browser).
- Obsidian + Citations plugin (already in use).
- Optional: Elicit account (AI literature search), Quarto (manuscript output).

**MCP servers to wire up:**

- `pyphysics-mcp` (new, custom — described below).
- `zotero-mcp` (community PyPI package, `--local` mode).
- Optional: `semantic-scholar-mcp` (community) if the crawl/from_s2.py needs interactive use during synthesis.

**Slash commands (`.claude/commands/`):**

- `/crawl-papers` → invokes `crawl/build_queue.py`.
- `/triage-papers` → opens the current `queue.md` for human review; on `--apply`, pushes keep-decisions into Zotero.
- `/synthesize-claims` → runs the three `synth/*.py` tools + composes a PR draft.
- `/update-tracker` → runs `update_acquisition_tracker.py` + `build_dataview_indexes.py` + `validate_wikilinks.py`.
- `/chapter-status` → thin wrapper around the existing `chapter_status.py`.

**Scheduled execution:**

- `/loop` for in-session recurring tasks (when you're actively at the keyboard).
- `CronCreate` for autonomous runs when you're not (daily crawl, weekly digest).

## 3. Detailed component designs

### 3.1 Crawl pipeline

Four independent crawlers, each producing a normalised candidate format. All write to the same intermediate format (one JSON line per candidate) so the merge step is trivial.

```jsonc
// Roadmapping/Tooling/triage/_inbox/2026-05-21-from_arxiv.jsonl
{"source": "arxiv", "doi": "...", "arxiv_id": "2505.12345", "title": "...",
 "authors": ["..."], "year": 2026, "abstract": "...",
 "why_surfaced": "matches tag #superconductivity; cites bcs1957",
 "score": 0.78, "candidate_cite_key": "smith2026_qubit_coherence"}
```

**`crawl/from_s2.py`** — Semantic Scholar recommendations. Input: list of cite-keys from `Roadmapping/History/Bibliography/{Primary,Retrospective}/*.md` that have DOIs. For each, hit `https://api.semanticscholar.org/recommendations/v1/papers/forpaper/DOI:<doi>` (free, no API key required). Score = citation overlap with rest of library. Output: top-N candidates per seed paper.

**`crawl/from_arxiv.py`** — Daily ArXiv listings. Query the ArXiv RSS endpoint for `physics.gen-ph`, `physics.hist-ph`, `quant-ph`, `cond-mat.supr-con` (matching the repo's threads). Filter by abstract keywords drawn from the chapter `tags:` frontmatter. Score by keyword density + author overlap with existing bib.

**`crawl/from_crossref.py`** — Crossref event data. The Crossref Event Data API surfaces newly-published papers that cite a given DOI. Iterate over our 115 cite-keys with DOIs; collect "new papers citing one of these." Score = number of our entries cited (a paper citing 5 of ours is more relevant than one citing 1).

**`crawl/from_zotero.py`** — Newly-added Zotero entries. Read the local Zotero API (`http://localhost:23119/api/users/0/items?since=<last_run>`). Surfaces papers Trey added to Zotero (e.g., via Zotero Connector browser capture) that don't yet have a YAML stub in the repo. These are the *highest-priority* candidates because Trey already curated them.

**Run mode.** `crawl/from_*.py` are individually runnable for debugging; the main entry point is `triage/build_queue.py` which invokes all four and merges.

### 3.2 Triage pipeline

**`triage/build_queue.py`** — merges the four crawler outputs, deduplicates by DOI/ArXiv-id, applies a relevance threshold (default: keep top 30 per week), and writes the human-facing digest.

**`triage/digest.py`** — produces `Roadmapping/Tooling/triage/queue.md`. **Writes raw candidate metadata only** (no Claude-generated summaries — decision E). When `/triage-papers` runs inside Claude Code, the slash command instructs Claude to generate inline summaries on-the-fly for the top-N candidates; summaries are ephemeral, not written back to the file. Format:

```markdown
# Triage queue — 2026-05-21 to 2026-05-28

**29 candidates**, sorted by relevance score.

---

## 1. [smith2026_qubit_coherence] Coherence limits in transmon qubit arrays
- Authors: Smith, Liu, Brown (Nature Physics, May 2026)
- DOI: 10.1038/s41567-026-12345
- Sources: 🟢 surfaced by ArXiv crawl, 🟢 also surfaced by Crossref (cites bcs1957)
- Score: **0.78**
- Why surfaced: matches `#superconductivity` tag; cites `bcs1957_superconductivity`; 4 of the author's previous papers are already in our bib.
- Abstract: [first 300 chars]…
- Suggested chapter: [[08_quantum_computing_open_questions]] §4.1 (decoherence)

### Decision
<!-- Set one of: keep / skip / later -->
decision:

### Notes
<!-- Optional Trey notes -->

---

## 2. …
```

Trey edits the file in Obsidian, fills in `decision:`, optionally adds notes, saves. Re-running `/triage-papers --apply` reads the decisions, pushes `keep`s to Zotero, archives `skip`s, leaves `later`s in the queue.

**Promotion to Zotero.** On `keep`, the script uses the Zotero local API to create the item with the appropriate collection tag (`era/forward`, `thread/superconductivity`, etc., pulled from the suggested-chapter's frontmatter). Better BibTeX auto-exports to `zotero.bib`; the next synthesis run scaffolds the YAML stub.

### 3.3 Synthesis pipeline

**`synth/cluster_claims.py`** — reads every chapter's `#inferred` and `#speculative` claims (already extracted by `_tools/build_dataview_indexes.py`). Clusters by tag co-occurrence. Output: a markdown report grouping related claims so Trey can see, e.g., "all four chapters touching `#fine-structure` make compatible claims, but Ch 3 §3.2 and Ch 5 §3.2 use different precision regimes — reconcile?"

**`synth/suggest_cross_refs.py`** — for each chapter, find references *that exist in the bibliography but are not yet wikilinked from the chapter*, with a relevance score based on tag overlap. Output: a per-chapter "should you cite this?" list. Trey reviews; updates wikilinks; runs `validate_wikilinks.py`.

**`synth/manuscript_diff.py`** — runs after a `git fetch`. Identifies which chapters have unmerged changes since the last manuscript build; for each, summarises the diff in plain English ("Ch 5 added 3 new bib entries and 1 new derivation in §3.2"). Useful for the weekly synthesis digest and eventually for a Quarto/Pandoc auto-rebuild.

### 3.4 Custom MCP server: `pyphysics-mcp`

A small MCP server (Python, FastMCP framework) that exposes the repo's pipeline as Claude tools. Lives in `mcp_servers/pyphysics_mcp/`. Useful so Claude can act on the bibliography without shelling out to scripts each time.

**Tools exposed:**

| Tool | Args | Returns |
|---|---|---|
| `search_bib` | `query: str`, `tags?: list[str]`, `era?: str`, `limit?: int=10` | matching YAML notes |
| `get_bib_note` | `cite_key: str` | full YAML frontmatter + body |
| `scaffold_bib_note` | `cite_key`, `doi?`, `type`, `era` | invokes existing `scaffold_bib_note.py` |
| `list_chapters` | (none) | chapter metadata |
| `get_chapter` | `chapter: str` | full chapter markdown |
| `search_claims` | `tag: str` (e.g., `#inferred`), `era?: str` | locations + section context |
| `validate_wikilinks` | `paths?: list[str]` | broken-link report |
| `regenerate_indexes` | (none) | runs `build_dataview_indexes.py` |
| `regenerate_tracker` | (none) | runs `update_acquisition_tracker.py` |
| `render_scene` | `scene_name`, `quality: ql|qh` | shells out to manim; returns video path |

**Why custom MCP rather than just letting Claude shell out?** Three reasons: (a) repeatable tool surface across sessions (Claude knows these tools exist without having to be told); (b) clean JSON-typed returns instead of raw stdout parsing; (c) lets the user grant once-and-done permission to the MCP tools instead of per-Bash-call confirms.

**Implementation.** ~400 lines using FastMCP. Each tool method wraps an existing script or directly reads/writes files. Distribution: `uv tool install pyphysics-mcp` from this repo (locally; not published to PyPI).

### 3.5 Sync layer: `sync_from_zotero.py`

Bridges Zotero's BBT-exported `.bib` to the YAML pipeline. Specced in detail in the earlier conversation. Recap:

1. Read `Roadmapping/History/Bibliography/zotero.bib`.
2. For each entry whose cite-key has no corresponding `Primary/<cite_key>.md` or `Retrospective/<cite_key>.md`, scaffold a stub by translating BibTeX fields → YAML schema.
3. Detect *modified* entries (different title, new authors) — flag for manual reconciliation rather than auto-overwriting.
4. Optional: copy PDFs from Zotero storage into `Historical_Papers/{Primary,Retrospective}/<cite_key>.pdf` when `pdf_status` flips to `acquired` or `out_of_copyright_public`.

### 3.6 Slash commands

Each is a thin markdown file in `.claude/commands/` invoking the underlying scripts. Format example:

```markdown
---
description: Run a full crawl across S2, ArXiv, Crossref, and Zotero; merge into the triage queue.
---

Run the crawl pipeline:

1. `uv run python Roadmapping/Tooling/crawl/from_s2.py`
2. `uv run python Roadmapping/Tooling/crawl/from_arxiv.py`
3. `uv run python Roadmapping/Tooling/crawl/from_crossref.py`
4. `uv run python Roadmapping/Tooling/crawl/from_zotero.py`
5. `uv run python Roadmapping/Tooling/triage/build_queue.py`
6. `uv run python Roadmapping/Tooling/triage/digest.py`

Then summarise the digest in 5 bullets, highlighting any candidate with score > 0.85.
```

Slash commands invoked from `/loop` produce a running stream of summary outputs while the pipelines work in the background.

### 3.7 `/loop` and `CronCreate` usage

Two execution modes:

**Interactive `/loop`** — for when Trey is at the keyboard and wants the agent to keep working on something for a session.

- `/loop 6h /crawl-papers` — re-crawl every 6 hours during an active research session.
- `/loop 1d /synthesize-claims` — daily synthesis update inside a long-running session.
- `/loop /triage-papers` (no interval; self-paced) — agent decides when to re-run based on whether the queue has new entries.

**Autonomous `CronCreate`** — for when Trey isn't at the keyboard.

- Cron: `0 7 * * *` → `/crawl-papers` — daily crawl at 7 a.m.
- Cron: `0 8 * * MON` → `/triage-papers` — Monday morning digest waiting.
- Cron: `0 9 * * FRI` → `/synthesize-claims` — Friday end-of-week synthesis.

Autonomous runs produce committed-to-the-repo artefacts (the queue.md, the synthesis reports). Nothing modifies chapter content without human review.

## 4. Implementation phases

Staged across 5 PRs so each one is reviewable and the system is incrementally usable.

### Phase 1 — Install + sync layer (1 PR)

- Install guide: `Roadmapping/Tooling/install_zotero_obsidian.md` — Linux install of Zotero + BBT + Connector, config matching repo cite-key convention, Obsidian Citations plugin setup, BBT auto-export to `zotero.bib`.
- `Bibliography/sync_from_zotero.py` — Zotero → YAML one-way sync.
- `.gitignore` adds `Roadmapping/History/Bibliography/zotero.bib`.
- Smoke test: capture one paper via Zotero Connector → run sync → confirm YAML stub appears with correct cite-key.

**Estimated effort:** ~200 lines + docs. Half a day.

### Phase 2 — Custom MCP server (1 PR)

- `mcp_servers/pyphysics_mcp/` — FastMCP-based server exposing the 10 tools listed in §3.4.
- Installation docs (`uv tool install` from this repo).
- Optional: register with Claude Code via `.claude/mcp.json` so it's available on next session start.
- Smoke test: ask Claude "search the bib for entries tagged #superconductivity" and confirm it goes via the MCP rather than reading files directly.

**Estimated effort:** ~400 lines. One day.

### Phase 3 — Crawl pipeline (1 PR)

- `Roadmapping/Tooling/crawl/{from_s2,from_arxiv,from_crossref,from_zotero}.py`.
- `triage/build_queue.py` + `triage/digest.py`.
- Tests with a small fixture (cached responses from each source) so the pipeline can be re-run deterministically.
- `.gitignore` for `Roadmapping/Tooling/triage/_inbox/` (raw crawler outputs are regeneratable).

**Estimated effort:** ~750 lines + tests. Two days.

### Phase 4 — Synthesis pipeline (1 PR)

- `Roadmapping/Tooling/synth/{cluster_claims,suggest_cross_refs,manuscript_diff}.py`.
- Each runs against the current repo state; outputs go to `Roadmapping/Tooling/synth_reports/<date>-*.md` (gitignored unless explicitly committed).
- The reports are designed to *seed* PR drafts, not auto-commit anything.

**Estimated effort:** ~450 lines. One day.

### Phase 5 — Slash commands + cron + docs (1 PR)

- `.claude/commands/{crawl-papers,triage-papers,synthesize-claims,update-tracker,chapter-status}.md`.
- Documentation for `/loop` usage (interactive) and `CronCreate` setup (autonomous).
- A `Roadmapping/Tooling/AUTOMATION.md` operator-guide covering the weekly rhythm: Monday triage, mid-week ad-hoc work, Friday synthesis.

**Estimated effort:** ~100 lines of slash-command markdown + docs. Half a day.

### Total

~1900 lines of new code + docs across 5 PRs, ~5 days of focused work. Each PR is independently mergeable and incrementally useful: after Phase 1 you have working Zotero sync; after Phase 2 Claude can act on the bibliography via MCP; after Phase 3 you have a triage queue; etc.

## 5. Locked decisions (2026-05-21)

All seven open decisions are now locked. Implementation reflects these choices.

**A. Zotero tags, not collections.** The crawl pipeline maps chapter `tags:` frontmatter → Zotero tags 1:1 (`#era/forward`, `#thread/quantum`, etc.). Collections are unused; Trey can add them later for visual organisation if desired.

**B. Hard cap: 30/day, 10/source.** `triage/build_queue.py` enforces a per-source cap of 10 and a total cap of 30. Surplus candidates drop with a "N more dropped at cap" note in the digest so nothing is silently lost.

**C. `triage_decisions.jsonl` audit log + `queue.md` working UI.** Decisions append-only to `Roadmapping/Tooling/triage/triage_decisions.jsonl` (committed; permanent record). `queue.md` is the live working surface Trey edits weekly.

**D. Synthesis emits summary reports only.** `synth/*.py` writes markdown reports to `Roadmapping/Tooling/synth_reports/<date>-*.md`. No auto-PR generation in v1; Trey applies suggestions by hand. Auto-PR generation is deferred to a possible Phase 6.

**E. No paid API summarisation — Claude Code handles it ephemerally.** `triage/digest.py` writes the queue *without* Claude-generated summaries (raw candidate metadata only — title, authors, abstract, why-surfaced, score, suggested chapter). When `/triage-papers` runs inside Claude Code, the slash command instructs Claude to generate a one-paragraph "why this matters for [suggested chapter]" summary inline in the chat for the top-N candidates by score. Summaries are ephemeral (regenerated each `/triage-papers` invocation, not written back to disk). Billed via Trey's existing Claude Code session, no extra API cost, no `ANTHROPIC_API_KEY` management. Architecture impact: §3.2 digest description and §3.4 MCP tool list updated accordingly.

**F. Only Unpaywall-confirmed open-access PDFs auto-fetched.** The crawl pipeline never attempts paywalled sources. `fetch_pdf.py` already uses Unpaywall, which by construction only returns open-access URLs. If Unpaywall returns nothing for a `keep`-decided candidate, the pipeline leaves `pdf_status: pending` for Trey to manually source (or skip). No paywalled-publisher API integration. No institutional-access scraping.

**G. Both CronCreate and local cron supported.** The plan ships:
- A `CronCreate` invocation snippet Trey pastes into Claude Code once (cloud-side; runs when laptop is off; uses Anthropic credits).
- A systemd-timer + `.service` file template under `Roadmapping/Tooling/cron/` (machine-side; free; only runs when laptop is on).
Same slash commands; choice is per-user-schedule.

## 6. What this gives Trey day-to-day

Once Phase 5 is complete, the routine is:

**Monday morning.** Open Claude Code in the PyPhysics repo. The autonomous Sunday-night crawl has already populated `Roadmapping/Tooling/triage/queue.md` with 20–30 candidates. Run `/triage-papers`. Claude opens the queue, summarises the top-5 with the most-cited papers in the repo's existing bibliography, and waits for Trey to fill in `decision:` fields. Trey takes 20 minutes; saves. Re-runs `/triage-papers --apply`. Keep-decisions push to Zotero; PDFs fetch where available; `sync_from_zotero.py` scaffolds YAML stubs for the new entries.

**Mid-week.** Normal chapter work proceeds. Claude can be asked things like "search the bib for entries tagged `#fine-structure` published since 2020" via the `pyphysics-mcp` tools; "check whether the new entries cite anything in Chapter 4 that should be cross-linked back."

**Friday afternoon.** Run `/synthesize-claims`. Three reports land:
- `cluster_claims_2026-05-22.md` — groupings of related `#inferred`/`#speculative` claims across chapters; flags potential inconsistencies.
- `suggest_cross_refs_2026-05-22.md` — proposed wikilink additions between chapters and new bib entries.
- `manuscript_diff_2026-05-22.md` — per-chapter "what changed this week" digest.

Trey skims, applies the useful suggestions as small commits, opens a PR for the week's chapter updates.

**Once a quarter.** Run `chapter_status.py` + `build_dataview_indexes.py` + `validate_wikilinks.py` to confirm the campaign-wide invariants. (These can also be cronned.)

## 7. Out of scope (explicitly)

- A web UI for the triage queue. Markdown + Obsidian is sufficient.
- Real-time crawl (sub-hour latency). Daily is the design target.
- Auto-merging anything. Every chapter change goes via PR + human review.
- Multi-user collaboration features. This is a single-author dissertation/book pipeline.
- Replacing the YAML bibliography or Obsidian vault. Both stay as the canonical layer.
- Audio podcast generation (PR L from the history campaign). Deferred per the original campaign plan.

## 8. Cross-reference

- Survey: [`Roadmapping/Tooling/dissertation_tooling_survey.md`](../../Roadmapping/Tooling/dissertation_tooling_survey.md).
- Recap of Zotero install + MCP setup: see the conversation log preceding this plan.
- Existing campaign tools: [`Roadmapping/History/_tools/`](../../Roadmapping/History/_tools/), [`Roadmapping/History/Bibliography/*.py`](../../Roadmapping/History/Bibliography/).
- Issue: [#34](https://github.com/temoTxt/PyPhysics/issues/34) (now closed; this plan picks up from the survey's "follow-up issues worth opening" list).

---

## Status

- [ ] Phase 1 — Zotero install guide + sync_from_zotero.py
- [ ] Phase 2 — pyphysics-mcp custom MCP server
- [ ] Phase 3 — Crawl pipeline (S2, ArXiv, Crossref, Zotero)
- [ ] Phase 4 — Synthesis pipeline (claim clustering, cross-ref suggestion, manuscript diff)
- [ ] Phase 5 — Slash commands + /loop + CronCreate docs

Each phase is one PR. Recommend opening as a new GitHub epic issue ("Continuous-crawl synthesis pipeline") with the five phases as sub-issues, the same way the history-of-physics campaign was structured (epic #7 + sub-issues #9–#20).

When Trey is ready to proceed, the recommended order is: open the epic + sub-issues first, then Phase 1 (which is the cheapest and most immediately useful — Zotero capture working end-to-end with the repo).
