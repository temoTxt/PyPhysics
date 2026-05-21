# Zotero + Obsidian + this repo — install guide

Phase 1 of the dissertation-tooling implementation plan ([`.dev/tasks/34-dissertation-tooling.md`](../../.dev/tasks/34-dissertation-tooling.md)). End state after following this guide:

1. Zotero 7 running on Linux with Better BibTeX (BBT) auto-exporting your library to `Roadmapping/History/Bibliography/zotero.bib` whenever an entry changes.
2. Zotero Connector browser extension capturing papers in one click from journal pages, ArXiv, JSTOR, Google Scholar.
3. Zotero's local HTTP API enabled (port 23119) — bridge for the future MCP integration and the existing `sync_from_zotero.py` script.
4. Obsidian's Citations plugin reading the BBT-exported `.bib` directly from the vault, enabling `[[cite:@maxwell1865_dynamical_theory]]`-style live citations.
5. `sync_from_zotero.py` translating new Zotero entries into YAML stubs under `Roadmapping/History/Bibliography/{Primary,Retrospective}/`, keeping the existing pipeline canonical.

Estimated time: 20–30 minutes including verification.

## 1. Install Zotero on Linux

The community **zotero-deb** repository (maintained by Emiliano Heyns, who also maintains BBT) is the de facto standard on Debian/Ubuntu — updates land within hours of upstream.

```bash
wget -qO- https://raw.githubusercontent.com/retorquere/zotero-deb/master/install.sh | sudo bash
sudo apt update
sudo apt install zotero
```

Launch Zotero. Sign in at zotero.org/user/register if you want cross-device sync; the local library works without an account. Either way, **create the local library** before continuing.

## 2. Install Better BibTeX (BBT)

1. Download the latest `.xpi` from [github.com/retorquere/zotero-better-bibtex/releases/latest](https://github.com/retorquere/zotero-better-bibtex/releases/latest).
2. In Zotero: `Tools → Add-ons → ⚙ → Install Add-on From File…` → select the `.xpi`.
3. Restart Zotero.

## 3. Configure cite-keys to match this repo's convention

The repo's convention is `firstauthor + year + slug`, lowercase, snake_case (e.g., `maxwell1865_dynamical_theory`). Configure BBT to produce keys in exactly this format:

`Edit → Settings → Better BibTeX → Citation Keys → Citation key formula:`

```
auth.lower + year + ('_' + shorttitle(3,3).lower) | random
```

`shorttitle(3,3)` takes the first three significant words of the title (after stopwords like "the", "a", "of" are stripped) and joins them. The `lower` and snake-case formatting match the existing 115 stubs exactly.

`Edit → Settings → Better BibTeX → Citation Keys → ✅ Pin on first import`. **Important** — without pinning, BBT will *regenerate* the cite-key any time a title or author field changes, which would silently break wikilinks in the chapters.

## 4. Configure BBT auto-export to the repo vault

`Edit → Settings → Better BibTeX → Automatic export → On change`.

Then in Zotero's main UI:
1. Right-click the top-level "My Library" → `Export Library…`.
2. Format: `Better BibTeX`. Tick `✅ Keep updated`.
3. Save path: `~/PyPhysics/Roadmapping/History/Bibliography/zotero.bib`.

From now on every Zotero edit auto-rewrites that file. The file is gitignored (added in this PR) — it's a mirror of Zotero, regeneratable.

## 5. Install Zotero Connector (browser)

Download the extension for your browser from [zotero.org/download](https://www.zotero.org/download/). Firefox, Chrome, Edge, and Brave are all supported. After install, the toolbar icon shows a context-aware capture button — a single click on most journal pages, ArXiv listings, JSTOR pages, and Google Scholar results captures the metadata + (often) the open-access PDF.

## 6. Enable Zotero's local HTTP API

`Edit → Settings → Advanced → ✅ Allow other applications on this computer to communicate with Zotero`.

Exposes the local API on `http://localhost:23119`. The `pyphysics-mcp` server (Phase 2) and Phase 3's `crawl/from_zotero.py` both use this. The `sync_from_zotero.py` in this PR doesn't need it — it reads the BBT-exported `.bib` file directly — but enabling it now means later phases don't need a config step.

## 7. Tag convention for the integration

To let `sync_from_zotero.py` correctly auto-fill the `type:` and `era:` YAML fields, **tag your Zotero entries with the same tag namespaces the repo uses**:

| Repo YAML field | Zotero tag | Example |
|---|---|---|
| `type: retrospective` | `repo:retrospective` | review articles, biographies, textbooks |
| `type: primary` | (default; no tag needed) | original research papers |
| `era: "1860-1900"` | `era/1860-1900` | the per-era buckets from the campaign |
| `era: "forward"` | `era/forward` | Ch 7–9 forward chapters |
| `tags: [electromagnetism]` | `thread/electromagnetism` | the campaign's three threads |
| `tags: [quantum]` | `thread/quantum` | |
| `tags: [solid-state]` | `thread/solid-state` | |
| (any free-form topic) | the tag verbatim, e.g. `superconductivity` | topic-level tags |

Tags applied in Zotero are exported in the BBT `keywords` field and re-parsed by `sync_from_zotero.py` into the YAML schema. If you don't apply tags, the YAML stubs get scaffolded with empty fields you fill in by hand.

## 8. Install Obsidian Citations plugin

In Obsidian: `Settings → Community plugins → Browse → Citations → Install → Enable`.

Configure:
- `Citation database format`: BibTeX
- `Citation database path`: `Roadmapping/History/Bibliography/zotero.bib`
- `Literature note location`: `Roadmapping/History/Bibliography/Primary` (or leave blank if you don't want auto-created notes)
- `Literature note title template`: `{{citekey}}` (matches the repo convention)

After save, the command palette (`Cmd/Ctrl-P`) gains `Citations: Insert literature note link` and `Citations: Open literature note`. Type `[[` in a note and search by cite-key to insert wikilinks; the Citations plugin auto-completes from the live `zotero.bib`.

## 9. Verify the integration

End-to-end smoke test. Should take under 2 minutes:

1. **Capture a test paper.** Open a paper page in your browser (ArXiv or a journal). Click the Zotero Connector icon. Confirm Zotero opens to the new entry.
2. **Tag it.** In Zotero, add tags: `repo:retrospective` (or skip for primary), `era/forward`, `thread/quantum` (whichever fits).
3. **Confirm BBT auto-export.** Within ~5 seconds, `Roadmapping/History/Bibliography/zotero.bib` should contain the new entry. Verify with:
   ```bash
   tail -20 Roadmapping/History/Bibliography/zotero.bib
   ```
4. **Run the sync.** From the repo root:
   ```bash
   uv run python Roadmapping/History/Bibliography/sync_from_zotero.py --dry-run
   ```
   Expected: a "would scaffold" line for the new cite-key.
5. **Apply the sync.** Drop `--dry-run`:
   ```bash
   uv run python Roadmapping/History/Bibliography/sync_from_zotero.py
   ```
   Expected: a new YAML stub at `Roadmapping/History/Bibliography/{Primary,Retrospective}/<cite_key>.md` with `title`, `authors`, `year`, `journal`, `doi`, `tags`, and `era` auto-filled from Zotero.
6. **Regenerate the tracker.**
   ```bash
   uv run python Roadmapping/History/Bibliography/update_acquisition_tracker.py
   ```
   Confirm the new cite-key appears in `Roadmapping/Historical_Papers/Acquisition_Tracker.md`.
7. **Validate wikilinks.**
   ```bash
   uv run python Roadmapping/History/_tools/validate_wikilinks.py
   ```
   Expected: still green; the new stub doesn't break anything.

If all 7 steps pass, you're done with Phase 1.

## 10. Day-to-day workflow

- **Reading a paper**: capture via Zotero Connector → entry auto-syncs to `zotero.bib` → run `sync_from_zotero.py` to scaffold the YAML stub → edit the stub to add the chapter-specific fields (`gill_corpus_overlap`, `chapters_citing`, 2–3-paragraph body).
- **Citing in a chapter**: type `[[cite_key]]` in Obsidian; the Citations plugin auto-completes from `zotero.bib`. Wikilink validation catches typos.
- **Tracking acquisition status**: after a batch of new captures, `update_acquisition_tracker.py` regenerates the per-cite-key PDF status table.

## 11. What's *not* done in Phase 1

Deferred to later phases:

- **PDF auto-copy from Zotero storage into `Historical_Papers/` tree.** Currently the sync script only generates YAML stubs; PDFs stay in Zotero's storage. Phase 3 (or a small follow-up) will add: if `pdf_status` is `out_of_copyright_public`, copy the PDF from Zotero storage to `Historical_Papers/{Primary,Retrospective}/<cite_key>.pdf` so it commits to git per the existing two-tier policy.
- **`pyphysics-mcp` server.** Phase 2 wires up the custom MCP so Claude can search the bibliography, scaffold stubs, validate wikilinks, etc. via typed MCP tools rather than shelling out.
- **Crawl pipeline.** Phase 3 adds the S2 / ArXiv / Crossref / Zotero-delta crawlers that build the weekly triage queue.

## 12. Cross-reference

- Plan doc: [`.dev/tasks/34-dissertation-tooling.md`](../../.dev/tasks/34-dissertation-tooling.md).
- Tooling survey: [`Roadmapping/Tooling/dissertation_tooling_survey.md`](dissertation_tooling_survey.md).
- Existing bibliography pipeline: [`Roadmapping/History/Bibliography/README.md`](../History/Bibliography/README.md).
- Existing per-chapter helper tools: [`Roadmapping/History/_tools/`](../History/_tools/).
