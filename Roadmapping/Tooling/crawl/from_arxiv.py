"""ArXiv daily-listing crawler.

Pulls the most recent submissions from a set of physics ArXiv categories
matching the campaign's three threads (electromagnetism, quantum,
solid-state), filters by topic-keyword density drawn from chapter
tags:, and scores by author overlap with the existing bibliography.

Emits one JSONL candidate per line to:
    Roadmapping/Tooling/triage/_inbox/<YYYY-MM-DD>-from_arxiv.jsonl

Uses the ArXiv Atom export endpoint (no rate limit; no API key).

Usage:
    uv run python Roadmapping/Tooling/crawl/from_arxiv.py
    uv run python Roadmapping/Tooling/crawl/from_arxiv.py --dry-run
    uv run python Roadmapping/Tooling/crawl/from_arxiv.py --categories quant-ph,cond-mat.supr-con
"""

import argparse
import json
import re
import sys
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bib import all_entries, arxiv_ids, dois, cite_keys, chapter_tags, known_authors  # noqa: E402
from _common import (  # noqa: E402
    REPO_ROOT, get_inbox_path, slugify_cite_key, emit_jsonl, score_clamp,
)
from _http import get_text  # noqa: E402
from _state import set_cursor  # noqa: E402

# ArXiv Atom feed (better-structured than the RSS feed).
ARXIV_QUERY_URL = "http://export.arxiv.org/api/query"

DEFAULT_CATEGORIES = ["physics.hist-ph", "quant-ph", "cond-mat.supr-con", "physics.gen-ph"]
ATOM_NS = {"a": "http://www.w3.org/2005/Atom",
           "arxiv": "http://arxiv.org/schemas/atom"}


def fetch_listing(category: str, max_results: int) -> list[dict]:
    """Fetch recent submissions in a category via the ArXiv API."""
    params = {
        "search_query": f"cat:{category}",
        "sortBy": "submittedDate",
        "sortOrder": "descending",
        "max_results": max_results,
    }
    body = get_text(ARXIV_QUERY_URL, params=params, accept="application/atom+xml")
    return parse_atom(body)


def parse_atom(xml_text: str) -> list[dict]:
    """Parse an ArXiv Atom feed into a list of raw entry dicts."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        print(f"  XML parse error: {e}", file=sys.stderr)
        return []
    entries: list[dict] = []
    for entry in root.findall("a:entry", ATOM_NS):
        e_id = entry.findtext("a:id", default="", namespaces=ATOM_NS)
        # ArXiv IDs come as URLs http://arxiv.org/abs/2505.12345v1
        m = re.search(r"abs/([^/]+?)(v\d+)?$", e_id)
        arxiv_id = m.group(1) if m else ""

        title = (entry.findtext("a:title", default="", namespaces=ATOM_NS) or "").strip()
        title = re.sub(r"\s+", " ", title)

        abstract = (entry.findtext("a:summary", default="", namespaces=ATOM_NS) or "").strip()
        abstract = re.sub(r"\s+", " ", abstract)

        published = entry.findtext("a:published", default="", namespaces=ATOM_NS) or ""
        year = None
        if published:
            try:
                year = int(published[:4])
            except ValueError:
                year = None

        authors = []
        for au in entry.findall("a:author", ATOM_NS):
            name = au.findtext("a:name", default="", namespaces=ATOM_NS) or ""
            if name.strip():
                authors.append(name.strip())

        # DOI sometimes present in arxiv:doi.
        doi = entry.findtext("arxiv:doi", default="", namespaces=ATOM_NS) or ""

        entries.append({
            "arxiv_id": arxiv_id,
            "title": title,
            "abstract": abstract,
            "year": year,
            "authors": authors,
            "doi": doi.lower(),
            "url": e_id,
        })
    return entries


def score_entry(entry: dict, topic_tags: set[str], our_authors: set[str]) -> tuple[float, list[str]]:
    """Compute (score, why_surfaced_tokens) for one entry.

    Score = baseline + keyword density + author-overlap bonus, clipped to [0, 1].
    """
    text = f"{entry['title']} {entry['abstract']}".lower()
    matched_tags = [t for t in topic_tags if t and t in text]
    keyword_density = min(0.4, 0.08 * len(matched_tags))

    matched_authors = []
    for a in entry["authors"]:
        surname = a.split()[-1].lower() if a.split() else ""
        if surname in our_authors:
            matched_authors.append(a)
    author_bonus = min(0.3, 0.1 * len(matched_authors))

    baseline = 0.1
    score = score_clamp(baseline + keyword_density + author_bonus)

    tokens = []
    if matched_tags:
        tokens.append(f"matches tag {', '.join(f'#{t}' for t in matched_tags[:3])}")
    if matched_authors:
        tokens.append(f"author overlap: {', '.join(matched_authors[:3])}")
    if not tokens:
        tokens.append("recent submission in scope")

    return score, tokens


def normalize(entry: dict, score: float, why_tokens: list[str]) -> dict:
    return {
        "source": "arxiv",
        "doi": entry["doi"],
        "arxiv_id": entry["arxiv_id"],
        "title": entry["title"],
        "authors": entry["authors"],
        "year": entry["year"],
        "abstract": entry["abstract"],
        "candidate_cite_key": slugify_cite_key(entry["authors"], entry["year"], entry["title"]),
        "why_surfaced": "; ".join(why_tokens),
        "score": score,
        "seed_paper": None,
        "raw_url": entry["url"],
        "tags_inferred": [],
        "captured_at": datetime.utcnow().isoformat() + "Z",
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--categories", default=",".join(DEFAULT_CATEGORIES),
                   help="comma-separated ArXiv categories (default: %(default)s)")
    p.add_argument("--per-category", type=int, default=40,
                   help="max recent submissions to scan per category (default 40)")
    p.add_argument("--min-score", type=float, default=0.20,
                   help="drop candidates with score below this (default 0.20)")
    p.add_argument("--dry-run", action="store_true",
                   help="print candidates to stdout instead of writing to _inbox")
    args = p.parse_args()

    categories = [c.strip() for c in args.categories.split(",") if c.strip()]
    topic_tags = chapter_tags()
    our_authors = known_authors()
    seen_arxiv = set(arxiv_ids())
    seen_dois = set(dois())
    seen_cite_keys = set(cite_keys())

    print(f"ArXiv crawler: {len(categories)} categories, up to {args.per_category} each",
          file=sys.stderr)
    print(f"  filter tags: {sorted(topic_tags)[:10]}{'...' if len(topic_tags) > 10 else ''}",
          file=sys.stderr)

    candidates: list[dict] = []
    for cat in categories:
        try:
            entries = fetch_listing(cat, args.per_category)
        except Exception as e:
            print(f"  {cat}: error {e}", file=sys.stderr)
            continue
        kept = 0
        for entry in entries:
            if entry["arxiv_id"] in seen_arxiv:
                continue
            if entry["doi"] and entry["doi"] in seen_dois:
                continue
            score, tokens = score_entry(entry, topic_tags, our_authors)
            if score < args.min_score:
                continue
            cand = normalize(entry, score, tokens)
            if cand["candidate_cite_key"] in seen_cite_keys:
                continue
            candidates.append(cand)
            seen_arxiv.add(entry["arxiv_id"])
            if entry["doi"]:
                seen_dois.add(entry["doi"])
            seen_cite_keys.add(cand["candidate_cite_key"])
            kept += 1
        print(f"  {cat}: {len(entries)} fetched, {kept} kept (after score >= {args.min_score})",
              file=sys.stderr)

    print(f"\nArXiv crawler: {len(candidates)} total candidates", file=sys.stderr)

    if args.dry_run:
        for c in candidates:
            print(json.dumps(c, default=str))
    else:
        out_path = get_inbox_path("arxiv")
        emit_jsonl(out_path, candidates)
        print(f"Wrote {out_path.relative_to(REPO_ROOT)}", file=sys.stderr)
        set_cursor("arxiv")

    return 0


if __name__ == "__main__":
    sys.exit(main())
