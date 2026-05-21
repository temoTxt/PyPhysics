"""Bibliography tools: search_bib, get_bib_note, scaffold_bib_note."""

import subprocess
from pathlib import Path

from pyphysics_mcp.repo import bib_dir, parse_frontmatter, relpath, repo_root, split_frontmatter

ALL_TYPES = ("Primary", "Retrospective")


def _all_stubs() -> list[Path]:
    out: list[Path] = []
    for sub in ALL_TYPES:
        d = bib_dir() / sub
        if d.is_dir():
            out.extend(sorted(d.glob("*.md")))
    return out


def search_bib(query: str = "", tags: list[str] | None = None,
                era: str | None = None, limit: int = 10) -> list[dict]:
    """Search bibliography stubs by free-text query + optional tag/era filters.

    Args:
        query: substring matched case-insensitively against title + authors + cite_key.
        tags: required tags (entry must have all listed tags).
        era: required era string (e.g., "1860-1900" or "forward").
        limit: max results (default 10).

    Returns: list of {cite_key, title, year, authors, era, type, tags, path}.
    """
    tags = tags or []
    q = (query or "").lower().strip()
    out: list[dict] = []

    for md in _all_stubs():
        meta = parse_frontmatter(md)
        if not meta:
            continue

        if era and str(meta.get("era", "")) != era:
            continue
        meta_tags = set(meta.get("tags") or [])
        if tags and not all(t in meta_tags for t in tags):
            continue

        if q:
            haystack = " ".join([
                str(meta.get("cite_key", "")),
                str(meta.get("title", "")),
                " ".join(str(a) for a in (meta.get("authors") or [])),
            ]).lower()
            if q not in haystack:
                continue

        out.append({
            "cite_key": meta.get("cite_key", md.stem),
            "title": meta.get("title", ""),
            "year": meta.get("year", ""),
            "authors": meta.get("authors") or [],
            "era": meta.get("era", ""),
            "type": meta.get("type", ""),
            "tags": list(meta_tags),
            "path": relpath(md),
        })
        if len(out) >= limit:
            break

    return out


def get_bib_note(cite_key: str) -> dict:
    """Return full YAML frontmatter + body for a single bib note."""
    for sub in ALL_TYPES:
        candidate = bib_dir() / sub / f"{cite_key}.md"
        if candidate.exists():
            meta, body = split_frontmatter(candidate)
            return {
                "cite_key": cite_key,
                "path": relpath(candidate),
                "frontmatter": meta,
                "body": body,
            }
    return {"error": f"cite_key {cite_key!r} not found under Bibliography/{{Primary,Retrospective}}/"}


def scaffold_bib_note(cite_key: str, type: str = "primary",
                       doi: str | None = None, era: str | None = None) -> dict:
    """Shell out to scaffold_bib_note.py to create a new YAML stub.

    Args:
        cite_key: firstauthor+year+slug, snake_case (e.g. "smith2026_qubit").
        type: 'primary' or 'retrospective'.
        doi: optional; triggers Crossref lookup for auto-filled metadata.
        era: optional; passed through to the script if not derivable from DOI.

    Returns: {status, path, stdout, stderr}.
    """
    script = repo_root() / "Roadmapping" / "History" / "Bibliography" / "scaffold_bib_note.py"
    if not script.exists():
        return {"error": f"scaffold_bib_note.py not found at {script}"}

    cmd = ["uv", "run", "python", str(script), "--cite-key", cite_key, "--type", type]
    if doi:
        cmd.extend(["--from-doi", doi])

    res = subprocess.run(cmd, capture_output=True, text=True, cwd=repo_root())
    target_sub = "Primary" if type == "primary" else "Retrospective"
    target = bib_dir() / target_sub / f"{cite_key}.md"

    out = {
        "status": "ok" if res.returncode == 0 else "error",
        "returncode": res.returncode,
        "path": relpath(target) if target.exists() else None,
        "stdout": res.stdout.strip(),
        "stderr": res.stderr.strip(),
    }

    # Pass era through by editing the YAML if specified — scaffold_bib_note.py doesn't
    # take --era directly, so set it post-creation if the file landed.
    if era and target.exists():
        text = target.read_text(encoding="utf-8")
        text = text.replace("era: ''", f"era: '{era}'", 1)
        target.write_text(text, encoding="utf-8")
        out["era_set"] = era

    return out
