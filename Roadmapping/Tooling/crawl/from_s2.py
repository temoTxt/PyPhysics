"""Semantic Scholar recommendations crawler.

For each DOI in our bibliography, fetches up to N papers Semantic Scholar
recommends as related. Deduplicates against entries we already have.
Emits one JSONL candidate per line to:

    Roadmapping/Tooling/triage/_inbox/<YYYY-MM-DD>-from_s2.jsonl

Score: S2's own similarity score (when present) + a small bump for
authors that overlap with our existing bibliography.

Usage:
    uv run python Roadmapping/Tooling/crawl/from_s2.py
    uv run python Roadmapping/Tooling/crawl/from_s2.py --dry-run --max-seeds 3
    uv run python Roadmapping/Tooling/crawl/from_s2.py --per-seed 5

The free S2 API requires no key. The Recommendations API
(https://api.semanticscholar.org/recommendations/v1/) is generous but
rate-limited to ~1 req/sec; we sleep accordingly.
"""

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

# Sibling-module imports without a package layout — add our dir to sys.path.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bib import all_entries, dois, cite_keys, known_authors  # noqa: E402
from _common import (  # noqa: E402
    REPO_ROOT, get_inbox_path, slugify_cite_key, emit_jsonl, score_clamp,
)
from _http import get_json  # noqa: E402
from _state import set_cursor  # noqa: E402

S2_REC_URL = "https://api.semanticscholar.org/recommendations/v1/papers/forpaper/DOI:{doi}"
S2_FIELDS = "title,authors,abstract,year,externalIds,referenceCount,citationCount"
REQUEST_DELAY = 1.1  # seconds, to respect S2's ~1/sec rate limit


def fetch_recommendations(doi: str, per_seed: int) -> list[dict]:
    """Hit S2 Recommendations API for one seed DOI."""
    url = S2_REC_URL.format(doi=doi)
    payload = get_json(url, params={"limit": per_seed, "fields": S2_FIELDS})
    return payload.get("recommendedPapers", []) if isinstance(payload, dict) else []


def normalize(item: dict, seed_doi: str, seed_cite_key: str) -> dict | None:
    """Translate one S2 paper into our normalized candidate schema."""
    title = item.get("title") or ""
    if not title:
        return None
    external = item.get("externalIds") or {}
    item_doi = (external.get("DOI") or "").lower()
    arxiv_id = external.get("ArXiv") or ""
    authors = [a.get("name", "") for a in (item.get("authors") or []) if a.get("name")]

    # Author-overlap bonus.
    our_authors = known_authors()
    overlap = sum(1 for a in authors if a.split()[-1].lower() in our_authors) if authors else 0
    base_score = 0.6  # S2 recommendations are themselves curated, so floor at 0.6
    bonus = min(0.3, 0.05 * overlap)
    score = score_clamp(base_score + bonus)

    return {
        "source": "s2",
        "doi": item_doi,
        "arxiv_id": arxiv_id,
        "title": title.strip(),
        "authors": authors,
        "year": item.get("year"),
        "abstract": (item.get("abstract") or "").strip(),
        "candidate_cite_key": slugify_cite_key(authors, item.get("year"), title),
        "why_surfaced": f"S2 recommends as related to {seed_cite_key}"
                        + (f" (author overlap: {overlap})" if overlap else ""),
        "score": score,
        "seed_paper": seed_cite_key,
        "raw_url": f"https://www.semanticscholar.org/paper/{item.get('paperId', '')}",
        "tags_inferred": [],
        "captured_at": datetime.utcnow().isoformat() + "Z",
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--per-seed", type=int, default=5,
                   help="max recommendations per seed DOI (default 5)")
    p.add_argument("--max-seeds", type=int, default=0,
                   help="cap on seed DOIs to process (0 = all)")
    p.add_argument("--dry-run", action="store_true",
                   help="print candidates to stdout instead of writing to _inbox")
    args = p.parse_args()

    seeds = [(e["cite_key"], str(e["doi"]).strip()) for e in all_entries()
             if e.get("cite_key") and e.get("doi")]
    if args.max_seeds:
        seeds = seeds[: args.max_seeds]

    if not seeds:
        print("No seed DOIs in bibliography; nothing to query.", file=sys.stderr)
        return 0

    print(f"S2 crawler: {len(seeds)} seed DOIs, up to {args.per_seed} recs each",
          file=sys.stderr)

    candidates: list[dict] = []
    seen_dois = set(dois())
    seen_cite_keys = set(cite_keys())

    for i, (seed_key, seed_doi) in enumerate(seeds, 1):
        try:
            recs = fetch_recommendations(seed_doi, args.per_seed)
        except Exception as e:
            print(f"  [{i}/{len(seeds)}] {seed_key}: error {e}", file=sys.stderr)
            time.sleep(REQUEST_DELAY)
            continue

        added = 0
        for item in recs:
            cand = normalize(item, seed_doi, seed_key)
            if not cand:
                continue
            if cand["doi"] and cand["doi"] in seen_dois:
                continue
            if cand["candidate_cite_key"] in seen_cite_keys:
                continue
            candidates.append(cand)
            seen_dois.add(cand["doi"])
            seen_cite_keys.add(cand["candidate_cite_key"])
            added += 1
        print(f"  [{i}/{len(seeds)}] {seed_key}: {added} new candidates",
              file=sys.stderr)
        time.sleep(REQUEST_DELAY)

    print(f"\nS2 crawler: {len(candidates)} total new candidates", file=sys.stderr)

    if args.dry_run:
        for c in candidates:
            print(json.dumps(c, default=str))
    else:
        out_path = get_inbox_path("s2")
        emit_jsonl(out_path, candidates)
        print(f"Wrote {out_path.relative_to(REPO_ROOT)}", file=sys.stderr)
        set_cursor("s2")

    return 0


if __name__ == "__main__":
    sys.exit(main())
