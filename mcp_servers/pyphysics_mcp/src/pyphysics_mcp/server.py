"""FastMCP server registering the 10 PyPhysics tools.

Run as:
    pyphysics-mcp                  # via console_script entry point
    python -m pyphysics_mcp        # equivalent

The server speaks MCP over stdio. Register it with Claude Code by adding
to your MCP config (~/.claude/mcp.json or via the Claude Code UI):

    {
      "mcpServers": {
        "pyphysics": {
          "command": "pyphysics-mcp",
          "env": {"PYPHYSICS_REPO_PATH": "/home/tmorris/PyPhysics"}
        }
      }
    }

Tools registered (10):
  Bibliography:  search_bib, get_bib_note, scaffold_bib_note
  Chapters:      list_chapters, get_chapter, search_claims
  Validation:    validate_wikilinks, regenerate_indexes, regenerate_tracker
  Animations:    render_scene
"""

from mcp.server.fastmcp import FastMCP

from pyphysics_mcp.tools import animations, bibliography, chapters, validation

mcp = FastMCP("pyphysics")


# ───── Bibliography ─────

@mcp.tool()
def search_bib(query: str = "", tags: list[str] | None = None,
                era: str | None = None, limit: int = 10) -> list[dict]:
    """Search YAML bibliography stubs by free-text query + optional tag/era filters.

    Args:
        query: substring matched case-insensitively against title, authors, cite_key.
        tags: required tags (entry must have all listed).
        era: required era string (e.g. "1860-1900" or "forward").
        limit: max results.
    """
    return bibliography.search_bib(query=query, tags=tags, era=era, limit=limit)


@mcp.tool()
def get_bib_note(cite_key: str) -> dict:
    """Return full YAML frontmatter + body for a single bib note."""
    return bibliography.get_bib_note(cite_key)


@mcp.tool()
def scaffold_bib_note(cite_key: str, type: str = "primary",
                       doi: str | None = None, era: str | None = None) -> dict:
    """Create a new YAML bibliography stub via scaffold_bib_note.py.

    Args:
        cite_key: firstauthor+year+slug, snake_case.
        type: 'primary' or 'retrospective'.
        doi: optional; triggers Crossref auto-fill.
        era: optional; written into the YAML after creation.
    """
    return bibliography.scaffold_bib_note(cite_key=cite_key, type=type, doi=doi, era=era)


# ───── Chapters ─────

@mcp.tool()
def list_chapters() -> list[dict]:
    """List all chapters (01-09) with metadata: title, era, status, scenes, verification anchors."""
    return chapters.list_chapters()


@mcp.tool()
def get_chapter(chapter: str) -> dict:
    """Return full frontmatter + body of a chapter.

    Args:
        chapter: chapter number ('01', '07') or filename stem.
    """
    return chapters.get_chapter(chapter)


@mcp.tool()
def search_claims(tag: str, era: str | None = None) -> list[dict]:
    """Find every claim tagged with a confidence-tier tag across chapters.

    Args:
        tag: one of 'verified', 'inferred', 'speculative', 'gill-silent'.
        era: optional era filter.
    """
    return chapters.search_claims(tag, era=era)


# ───── Validation ─────

@mcp.tool()
def validate_wikilinks(paths: list[str] | None = None) -> dict:
    """Run the wikilink validator. Returns broken-link count + raw output."""
    return validation.validate_wikilinks(paths)


@mcp.tool()
def regenerate_indexes() -> dict:
    """Regenerate the three Dataview-like index pages under Roadmapping/History/."""
    return validation.regenerate_indexes()


@mcp.tool()
def regenerate_tracker() -> dict:
    """Regenerate Historical_Papers/Acquisition_Tracker.md from current bib YAML."""
    return validation.regenerate_tracker()


# ───── Animations ─────

@mcp.tool()
def render_scene(scene_name: str, quality: str = "ql") -> dict:
    """Render a Manim scene from Roadmapping/Animations/manim_scenes/.

    Args:
        scene_name: scene filename stem OR ClassName.
        quality: 'ql' (480p, ~10s), 'qm' (720p), 'qh' (1080p, slow).
    """
    return animations.render_scene(scene_name, quality=quality)


# ───── entry point ─────

def main() -> None:
    """stdio MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
