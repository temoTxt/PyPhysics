"""Zotero local-API delta crawler.

Reads Zotero's local HTTP API (http://localhost:23119) for items added
or modified since the last crawl run. Surfaces entries Trey has
captured via Zotero Connector that aren't yet in the YAML bibliography
— the highest-priority crawl source because the human has already
hand-curated these.

Requires Zotero running on the local machine with:
  Edit → Settings → Advanced → ✅ Allow other applications to communicate

Emits one JSONL candidate per line to:
    Roadmapping/Tooling/triage/_inbox/<YYYY-MM-DD>-from_zotero.jsonl

Score: 1.0 (these are already triaged; the human picked them).

Usage:
    uv run python Roadmapping/Tooling/crawl/from_zotero.py
    uv run python Roadmapping/Tooling/crawl/from_zotero.py --dry-run
    uv run python Roadmapping/Tooling/crawl/from_zotero.py --since 2026-05-01
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from urllib.error import URLError

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bib import cite_keys, dois  # noqa: E402
from _common import (  # noqa: E402
    REPO_ROOT, get_inbox_path, slugify_cite_key, emit_jsonl,
)
from _http import get_json  # noqa: E402
from _state import get_cursor, set_cursor  # noqa: E402

ZOTERO_ITEMS_URL = "http://localhost:23119/api/users/0/items"


def fetch_recent_items(since_date: str | None, limit: int = 100) -> list[dict]:
    """Hit the Zotero local API for recent items.

    The local API mirrors the web API: items are returned newest first;
    we filter by `dateAdded` client-side since the local API's `since`
    is a library version, not a date.
    """
    try:
        items = get_json(
            ZOTERO_ITEMS_URL,
            params={"limit": limit, "format": "json", "sort": "dateAdded", "direction": "desc"},
            timeout=15,
        )
    except URLError as e:
        raise RuntimeError(
            "Could not reach Zotero local API at localhost:23119. "
            "Is Zotero running with 'Allow other applications to communicate' enabled? "
            f"(underlying error: {e})"
        )
    if not isinstance(items, list):
        return []
    if since_date:
        items = [i for i in items if i.get("data", {}).get("dateAdded", "") >= since_date]
    return items


def normalize(item: dict) -> dict | None:
    """Translate a Zotero item to our normalized candidate schema."""
    data = item.get("data", {})
    title = data.get("title") or ""
    if not title:
        return None

    creators = data.get("creators", []) or []
    authors: list[str] = []
    for c in creators:
        if c.get("creatorType") not in ("author", "editor", None):
            continue
        first = c.get("firstName", "")
        last = c.get("lastName", "")
        name = (c.get("name") or f"{first} {last}").strip()
        if name:
            authors.append(name)

    date_str = data.get("date") or ""
    year = None
    if date_str:
        m = next((p for p in date_str.split() if p.isdigit() and len(p) == 4), None)
        try:
            year = int(m) if m else int(date_str[:4])
        except (ValueError, TypeError):
            year = None

    doi = (data.get("DOI") or "").lower()
    journal = data.get("publicationTitle") or data.get("publisher") or ""
    cite_key = data.get("citationKey", "") or slugify_cite_key(authors, year, title)
    abstract = data.get("abstractNote") or ""

    tags = []
    for t in data.get("tags", []) or []:
        tag_name = t.get("tag", "") if isinstance(t, dict) else str(t)
        if tag_name:
            tags.append(tag_name)

    return {
        "source": "zotero",
        "doi": doi,
        "arxiv_id": "",
        "title": title.strip(),
        "authors": authors,
        "year": year,
        "abstract": abstract.strip(),
        "candidate_cite_key": cite_key,
        "why_surfaced": "manually captured in Zotero"
                        + (f" with tags: {', '.join(tags[:5])}" if tags else ""),
        "score": 1.0,  # human-curated → highest priority
        "seed_paper": None,
        "raw_url": data.get("url", "") or (f"https://doi.org/{doi}" if doi else ""),
        "tags_inferred": tags,
        "captured_at": datetime.utcnow().isoformat() + "Z",
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--since", help="ISO date (YYYY-MM-DD) to filter on; defaults to last-run cursor")
    p.add_argument("--limit", type=int, default=100, help="max recent items to scan (default 100)")
    p.add_argument("--dry-run", action="store_true",
                   help="print candidates to stdout instead of writing to _inbox")
    args = p.parse_args()

    since = args.since
    if not since:
        cursor = get_cursor("zotero")
        if cursor:
            since = cursor[:10]
    if since:
        print(f"Zotero crawler: filtering to items dateAdded >= {since}", file=sys.stderr)
    else:
        print(f"Zotero crawler: no since cursor; scanning {args.limit} most-recent items",
              file=sys.stderr)

    try:
        items = fetch_recent_items(since_date=since, limit=args.limit)
    except RuntimeError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    seen_cite_keys = set(cite_keys())
    seen_dois = set(dois())
    candidates: list[dict] = []
    for item in items:
        cand = normalize(item)
        if not cand:
            continue
        if cand["candidate_cite_key"] in seen_cite_keys:
            continue
        if cand["doi"] and cand["doi"] in seen_dois:
            continue
        candidates.append(cand)
        seen_cite_keys.add(cand["candidate_cite_key"])
        if cand["doi"]:
            seen_dois.add(cand["doi"])

    print(f"\nZotero crawler: {len(items)} items scanned, {len(candidates)} new candidates",
          file=sys.stderr)

    if args.dry_run:
        for c in candidates:
            print(json.dumps(c, default=str))
    else:
        out_path = get_inbox_path("zotero")
        emit_jsonl(out_path, candidates)
        print(f"Wrote {out_path.relative_to(REPO_ROOT)}", file=sys.stderr)
        set_cursor("zotero")

    return 0


if __name__ == "__main__":
    sys.exit(main())
