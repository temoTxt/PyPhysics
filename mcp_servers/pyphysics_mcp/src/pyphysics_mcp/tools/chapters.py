"""Chapter tools: list_chapters, get_chapter, search_claims."""

import re
from pathlib import Path

from pyphysics_mcp.repo import (
    forward_dir,
    history_dir,
    parse_frontmatter,
    relpath,
    split_frontmatter,
)

CLAIM_TAG_RE = re.compile(r"`#(verified|inferred|speculative|gill-silent)`")
HEADING_RE = re.compile(r"^(#+)\s+(.+?)\s*$", re.MULTILINE)


def _chapter_files() -> list[Path]:
    """All numbered chapter files: 01-06 under History/, 07-09 under Forward/."""
    out: list[Path] = []
    out.extend(sorted(history_dir().glob("0*.md")))
    out.extend(sorted(forward_dir().glob("0*.md")))
    return out


def list_chapters() -> list[dict]:
    """List all chapters with metadata.

    Returns: list of {chapter, title, era, status, path, scenes, verification_anchors}.
    """
    out: list[dict] = []
    for ch in _chapter_files():
        meta = parse_frontmatter(ch)
        if not meta:
            continue
        out.append({
            "chapter": str(meta.get("chapter", ch.stem)),
            "title": meta.get("title", ""),
            "era": meta.get("era", ""),
            "status": meta.get("status", ""),
            "scenes": meta.get("animations") or [],
            "verification_anchors": meta.get("verification_anchors") or [],
            "path": relpath(ch),
        })
    return out


def get_chapter(chapter: str) -> dict:
    """Return full frontmatter + body for a chapter.

    Args:
        chapter: chapter number ('01', '07') or filename stem.
    """
    target: Path | None = None
    for ch in _chapter_files():
        meta = parse_frontmatter(ch)
        if str(meta.get("chapter", "")).zfill(2) == chapter.zfill(2):
            target = ch
            break
        if ch.stem == chapter:
            target = ch
            break
    if target is None:
        return {"error": f"chapter {chapter!r} not found"}

    meta, body = split_frontmatter(target)
    return {
        "chapter": str(meta.get("chapter", target.stem)),
        "path": relpath(target),
        "frontmatter": meta,
        "body": body,
        "body_bytes": len(body.encode("utf-8")),
    }


def search_claims(tag: str, era: str | None = None) -> list[dict]:
    """Find every claim tagged with a confidence-tier tag across chapters.

    Args:
        tag: one of 'verified', 'inferred', 'speculative', 'gill-silent'.
        era: optional era filter (matches chapter's era frontmatter).

    Returns: list of {chapter, line, section, snippet, path}.
    """
    tag_norm = tag.lstrip("#")
    out: list[dict] = []
    for ch in _chapter_files():
        meta, body = split_frontmatter(ch)
        if era and str(meta.get("era", "")) != era:
            continue

        text = body
        # Map line numbers to nearest preceding heading.
        headings: list[tuple[int, str]] = []
        for m in HEADING_RE.finditer(text):
            line_no = text[: m.start()].count("\n") + 1
            headings.append((line_no, m.group(2)))

        for m in CLAIM_TAG_RE.finditer(text):
            matched_tag = m.group(1)
            if matched_tag != tag_norm:
                continue
            line_no = text[: m.start()].count("\n") + 1
            # Snippet: ±100 chars around the match.
            start = max(0, m.start() - 100)
            end = min(len(text), m.end() + 200)
            snippet = text[start:end].replace("\n", " ").strip()
            section = "(top)"
            for h_line, h_text in headings:
                if h_line <= line_no:
                    section = h_text
                else:
                    break
            out.append({
                "chapter": str(meta.get("chapter", ch.stem)),
                "line": line_no,
                "section": section,
                "snippet": snippet,
                "path": relpath(ch),
            })
    return out
