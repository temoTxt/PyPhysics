# pyphysics-mcp

MCP server exposing the [PyPhysics](https://github.com/temoTxt/PyPhysics) repo's bibliography, chapter, validation, and animation tools as Claude-callable MCP tools.

Phase 2 of the dissertation-tooling implementation plan ([`.dev/tasks/34-dissertation-tooling.md`](../../.dev/tasks/34-dissertation-tooling.md)).

## What it does

Without this MCP, Claude (or Claude Code) interacts with the repo's pipeline by reading files directly and shelling out to scripts via Bash. That's permission-prompt-heavy and produces ad-hoc stdout-parsing in every conversation.

With this MCP, Claude has 10 typed tools:

| Tool | Action |
|---|---|
| `search_bib` | Search YAML bib stubs by query + tags + era |
| `get_bib_note` | Fetch a single bib note's frontmatter + body |
| `scaffold_bib_note` | Create a new YAML stub (wraps `scaffold_bib_note.py`) |
| `list_chapters` | List all 9 chapters with metadata |
| `get_chapter` | Fetch a single chapter's frontmatter + body |
| `search_claims` | Find every `#verified` / `#inferred` / `#speculative` / `#gill-silent` claim across chapters |
| `validate_wikilinks` | Run the wikilink validator + report broken count |
| `regenerate_indexes` | Refresh the three Dataview-like index pages |
| `regenerate_tracker` | Regenerate `Historical_Papers/Acquisition_Tracker.md` from current YAML |
| `render_scene` | Render a Manim scene at quality `ql`/`qm`/`qh` |

Each returns structured JSON; no stdout parsing.

## Install

```bash
# From this repo's root.
uv tool install ./mcp_servers/pyphysics_mcp
```

This installs `pyphysics-mcp` as a console script in `~/.local/share/uv/tools/`. Verify:

```bash
pyphysics-mcp --help 2>&1 | head -3
```

(The `--help` won't print MCP-server help — the server expects stdio from a client — but it should at least not error.)

## Register with Claude Code

Edit `~/.claude/mcp.json` (create if missing):

```json
{
  "mcpServers": {
    "pyphysics": {
      "command": "pyphysics-mcp",
      "env": {
        "PYPHYSICS_REPO_PATH": "/home/tmorris/PyPhysics"
      }
    }
  }
}
```

Replace `/home/tmorris/PyPhysics` with the absolute path to your repo clone. Restart Claude Code.

After restart, the deferred tools list will include `mcp__pyphysics__*` entries — load them via `ToolSearch` when you want to use them, or just ask Claude something like:

> Use the pyphysics MCP to find all bib entries tagged #superconductivity from before 1900.

## How `PYPHYSICS_REPO_PATH` is resolved

In order:

1. `PYPHYSICS_REPO_PATH` env var, if set and points to a directory containing `Roadmapping/History/Bibliography/`.
2. Walk up from `cwd` looking for the same directory. (Lets you run `pyphysics-mcp` directly from the repo without setting the env var.)
3. Error with a clear remediation message.

## Architecture

```
mcp_servers/pyphysics_mcp/
├── pyproject.toml
├── README.md
└── src/pyphysics_mcp/
    ├── __init__.py       # re-exports main
    ├── __main__.py       # python -m pyphysics_mcp
    ├── server.py         # FastMCP server + tool decorators
    ├── repo.py           # path resolution + YAML helpers
    └── tools/
        ├── bibliography.py   # search_bib, get_bib_note, scaffold_bib_note
        ├── chapters.py       # list_chapters, get_chapter, search_claims
        ├── validation.py     # validate_wikilinks, regenerate_indexes, regenerate_tracker
        └── animations.py     # render_scene
```

The server is intentionally thin: tools that wrap existing scripts (`scaffold_bib_note`, the regenerators, `render_scene`) shell out via `uv run python …` and return structured `{status, returncode, stdout, stderr, …}`. Tools that read the repo state directly (`search_bib`, `get_bib_note`, `list_chapters`, `get_chapter`, `search_claims`) parse YAML frontmatter + markdown in-process and return typed dicts.

Adding a new tool: write the implementation in the appropriate `tools/*.py` module, then expose it via `@mcp.tool()` in `server.py`.

## Verifying it works

After registering with Claude Code, in any conversation inside the repo:

```
Use the pyphysics MCP search_bib tool to find entries with the cite_key prefix 'maxwell'.
```

Should return 4 entries (`maxwell1855_lines_of_force`, `maxwell1861_physical_lines`, `maxwell1865_dynamical_theory`, `maxwell1873_treatise`).

```
Use the pyphysics MCP list_chapters tool.
```

Should return 9 chapter entries with full metadata.

```
Use the pyphysics MCP search_claims with tag='inferred' and era='1860-1900'.
```

Should return the `#inferred` claims from Chapter 2.

## Troubleshooting

**"PYPHYSICS_REPO_PATH not set and no Roadmapping/History/Bibliography/ found"** — set the env var in your Claude Code MCP config (above).

**`uv tool install` fails with "module 'pyphysics_mcp' not found"** — the `pyproject.toml` build-system is `hatchling`. If installation fails, try `uv tool install --reinstall ./mcp_servers/pyphysics_mcp`.

**`render_scene` returns non-zero** — the Animations sub-project has its own venv (Python 3.13). Tools that shell out via `uv run` will pick it up automatically when run from inside `Roadmapping/Animations/` (which the tool does via `cwd=animations_dir()`). If first-time setup hasn't happened, run `uv sync --project Roadmapping/Animations` from the repo root once.

## Cross-reference

- Plan: [`.dev/tasks/34-dissertation-tooling.md`](../../.dev/tasks/34-dissertation-tooling.md) §3.4.
- Phase 1: [`Roadmapping/Tooling/install_zotero_obsidian.md`](../../Roadmapping/Tooling/install_zotero_obsidian.md) + `sync_from_zotero.py`.
- Issue: [#34](https://github.com/temoTxt/PyPhysics/issues/34).
