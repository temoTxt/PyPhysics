"""Shared helpers used by every crawler.

- get_inbox_path: today's JSONL output path under triage/_inbox/.
- emit_jsonl: write a list of dicts as JSONL.
- slugify_cite_key: build a 'firstauthor + year + slug' candidate cite-key.
- score_clamp: clip score into [0, 1].
"""

import json
import re
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
INBOX_DIR = REPO_ROOT / "Roadmapping" / "Tooling" / "triage" / "_inbox"

STOPWORDS = {
    "the", "a", "an", "of", "on", "in", "for", "to", "and", "or", "with",
    "from", "by", "is", "are", "as", "at", "be",
}


def get_inbox_path(source: str) -> Path:
    """Return today's per-source JSONL path under triage/_inbox/."""
    INBOX_DIR.mkdir(parents=True, exist_ok=True)
    today = date.today().isoformat()
    return INBOX_DIR / f"{today}-from_{source}.jsonl"


def emit_jsonl(path: Path, records: list[dict]) -> int:
    """Append-write records as JSONL to path. Returns the count written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False, default=str) + "\n")
    return len(records)


def slugify_cite_key(authors: list[str] | None, year: int | str | None, title: str) -> str:
    """Build 'firstauthor + year + slug' candidate matching the repo convention.

    Example: ['Jane Smith', 'Wei Liu'], 2026, "Coherence in Transmon Arrays"
        → "smith2026_coherence_transmon_arrays"
    """
    first_author_surname = ""
    if authors:
        tokens = str(authors[0]).strip().split()
        if tokens:
            first_author_surname = tokens[-1].lower()
    first_author_surname = re.sub(r"[^a-z0-9]", "", first_author_surname)

    year_s = str(year) if year else ""
    year_s = re.sub(r"[^0-9]", "", year_s)[:4]

    # Slug from title: take first 3 significant words (stopwords stripped).
    words = re.findall(r"[A-Za-z]+", title.lower())
    significant = [w for w in words if w not in STOPWORDS][:3]
    slug = "_".join(significant)

    parts = [p for p in (first_author_surname, year_s) if p]
    base = "".join(parts) if parts else "unknown"
    return f"{base}_{slug}" if slug else base


def score_clamp(x: float) -> float:
    return max(0.0, min(1.0, float(x)))
