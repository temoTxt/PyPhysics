"""Shared helpers for inspecting the existing YAML bibliography.

Used by every crawler to determine 'which papers do we already have?'
and 'which DOIs/cite-keys should we use as seeds?'
"""

import re
from functools import lru_cache
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
BIB_DIR = REPO_ROOT / "Roadmapping" / "History" / "Bibliography"

FRONTMATTER_RE = re.compile(r"^---\n(.*?)\n---", re.DOTALL)


@lru_cache(maxsize=1)
def all_entries() -> list[dict]:
    """Return every YAML bib note's frontmatter as a list of dicts."""
    out: list[dict] = []
    for sub in ("Primary", "Retrospective"):
        d = BIB_DIR / sub
        if not d.is_dir():
            continue
        for md in sorted(d.glob("*.md")):
            text = md.read_text(encoding="utf-8")
            m = FRONTMATTER_RE.match(text)
            if not m:
                continue
            try:
                meta = yaml.safe_load(m.group(1)) or {}
            except yaml.YAMLError:
                continue
            meta["_path"] = str(md.relative_to(REPO_ROOT))
            meta["_type_dir"] = sub.lower()
            out.append(meta)
    return out


@lru_cache(maxsize=1)
def cite_keys() -> set[str]:
    return {e.get("cite_key", "") for e in all_entries() if e.get("cite_key")}


@lru_cache(maxsize=1)
def dois() -> set[str]:
    return {str(e.get("doi", "")).strip().lower() for e in all_entries() if e.get("doi")}


@lru_cache(maxsize=1)
def arxiv_ids() -> set[str]:
    return {str(e.get("arxiv_id", "")).strip() for e in all_entries() if e.get("arxiv_id")}


@lru_cache(maxsize=1)
def known_authors() -> set[str]:
    """Last names of all authors in our bib, lowercased."""
    out: set[str] = set()
    for e in all_entries():
        for a in e.get("authors") or []:
            # Take the last whitespace-separated token as a rough surname.
            tokens = str(a).strip().split()
            if tokens:
                out.add(tokens[-1].lower())
    return out


def chapter_tags() -> set[str]:
    """All free-form tags (excluding namespaced ones) used across the bib.

    Used by ArXiv crawler to filter on keyword density.
    """
    out: set[str] = set()
    for e in all_entries():
        for t in e.get("tags") or []:
            ts = str(t).strip().lstrip("#")
            if "/" in ts:  # skip era/* and thread/* — they're navigation, not topic
                continue
            out.add(ts.lower())
    return out
