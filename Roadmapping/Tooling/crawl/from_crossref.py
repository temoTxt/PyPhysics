"""Crossref event-data crawler.

For each DOI in our bibliography, queries Crossref's event-data API for
*new papers that cite that DOI* since the last crawl run. Surfaces
papers actively building on the corpus we already track.

Score scales with how many of our existing entries the candidate cites
(a paper citing 5 of our entries is much more relevant than one citing 1).

Emits one JSONL candidate per line to:
    Roadmapping/Tooling/triage/_inbox/<YYYY-MM-DD>-from_crossref.jsonl

API: https://api.eventdata.crossref.org/v1/events (free, no key).

Usage:
    uv run python Roadmapping/Tooling/crawl/from_crossref.py
    uv run python Roadmapping/Tooling/crawl/from_crossref.py --dry-run
    uv run python Roadmapping/Tooling/crawl/from_crossref.py --max-seeds 5
"""

import argparse
import json
import math
import sys
import time
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bib import all_entries, dois as bib_dois, cite_keys  # noqa: E402
from _common import (  # noqa: E402
    REPO_ROOT, get_inbox_path, slugify_cite_key, emit_jsonl, score_clamp,
)
from _http import get_json  # noqa: E402
from _state import get_cursor, set_cursor  # noqa: E402

CROSSREF_EVENTS_URL = "https://api.eventdata.crossref.org/v1/events"
CROSSREF_WORK_URL = "https://api.crossref.org/works/{doi}"
DEFAULT_LOOKBACK_DAYS = 30
REQUEST_DELAY = 0.5


def fetch_events(obj_doi: str, from_date: str, rows: int = 50) -> list[dict]:
    """Query the Crossref event-data API for events where the cited paper is obj_doi."""
    params = {
        "obj-id": f"https://doi.org/{obj_doi}",
        "from-occurred-date": from_date,
        "rows": rows,
    }
    payload = get_json(CROSSREF_EVENTS_URL, params=params)
    if isinstance(payload, dict):
        return payload.get("message", {}).get("events", [])
    return []


def fetch_work_metadata(doi: str) -> dict | None:
    """Fetch full Crossref metadata for a citing paper."""
    try:
        payload = get_json(CROSSREF_WORK_URL.format(doi=doi))
    except Exception:
        return None
    if isinstance(payload, dict):
        return payload.get("message", {})
    return None


def normalize(work: dict, citing_seeds: set[str]) -> dict | None:
    """Translate a Crossref work object + citation-overlap into our schema."""
    title = ""
    if isinstance(work.get("title"), list) and work["title"]:
        title = work["title"][0]
    if not title:
        return None

    authors = []
    for a in work.get("author") or []:
        given = a.get("given", "").strip()
        family = a.get("family", "").strip()
        if family:
            authors.append(f"{given} {family}".strip())

    year_parts = (
        work.get("issued", {}).get("date-parts")
        or work.get("published-print", {}).get("date-parts")
        or work.get("published-online", {}).get("date-parts")
    )
    year = year_parts[0][0] if year_parts and year_parts[0] else None

    journal = ""
    if isinstance(work.get("container-title"), list) and work["container-title"]:
        journal = work["container-title"][0]

    doi = (work.get("DOI") or "").lower()

    # Score: log(1 + N) so 1→0.69, 3→1.39, 5→1.79; normalize by clipping to [0, 1].
    n_overlap = len(citing_seeds)
    score = score_clamp(0.4 + 0.25 * math.log(1 + n_overlap))

    return {
        "source": "crossref",
        "doi": doi,
        "arxiv_id": "",
        "title": title.strip(),
        "authors": authors,
        "year": year,
        "abstract": "",  # Crossref doesn't reliably surface abstracts.
        "candidate_cite_key": slugify_cite_key(authors, year, title),
        "why_surfaced": f"cites {n_overlap} of our entries: {', '.join(sorted(citing_seeds)[:5])}"
                        + (f" (and {n_overlap - 5} more)" if n_overlap > 5 else "")
                        + (f"; published in {journal}" if journal else ""),
        "score": score,
        "seed_paper": next(iter(citing_seeds)) if citing_seeds else None,
        "raw_url": f"https://doi.org/{doi}" if doi else "",
        "tags_inferred": [],
        "captured_at": datetime.utcnow().isoformat() + "Z",
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lookback-days", type=int, default=DEFAULT_LOOKBACK_DAYS,
                   help=f"how far back to query events (default {DEFAULT_LOOKBACK_DAYS} days)")
    p.add_argument("--max-seeds", type=int, default=0,
                   help="cap on seed DOIs to process (0 = all)")
    p.add_argument("--min-overlap", type=int, default=1,
                   help="drop candidates citing fewer than this many of our entries (default 1)")
    p.add_argument("--dry-run", action="store_true",
                   help="print candidates to stdout instead of writing to _inbox")
    args = p.parse_args()

    # Determine from_date: cursor if set, else lookback_days ago.
    cursor = get_cursor("crossref")
    if cursor:
        from_date = cursor[:10]
    else:
        from_date = (datetime.now(timezone.utc) - timedelta(days=args.lookback_days)).date().isoformat()
    print(f"Crossref crawler: looking back to {from_date}", file=sys.stderr)

    seeds = [(e["cite_key"], str(e["doi"]).strip().lower()) for e in all_entries()
             if e.get("cite_key") and e.get("doi")]
    if args.max_seeds:
        seeds = seeds[: args.max_seeds]
    if not seeds:
        print("No seed DOIs; nothing to query.", file=sys.stderr)
        return 0

    # Step 1: collect (citing_doi → set of seeds it cites) across all seeds.
    citing: dict[str, set[str]] = defaultdict(set)
    for i, (seed_key, seed_doi) in enumerate(seeds, 1):
        try:
            events = fetch_events(seed_doi, from_date)
        except Exception as e:
            print(f"  [{i}/{len(seeds)}] {seed_key}: error {e}", file=sys.stderr)
            time.sleep(REQUEST_DELAY)
            continue
        for ev in events:
            subj = ev.get("subj_id") or ""
            if "/" in subj:
                # subj_id is a doi.org URL; extract the DOI.
                citing_doi = subj.split("doi.org/", 1)[-1].lower()
                if citing_doi and citing_doi not in bib_dois():
                    citing[citing_doi].add(seed_key)
        print(f"  [{i}/{len(seeds)}] {seed_key}: {len(events)} events", file=sys.stderr)
        time.sleep(REQUEST_DELAY)

    # Step 2: filter by min-overlap.
    citing = {d: seeds_set for d, seeds_set in citing.items() if len(seeds_set) >= args.min_overlap}
    print(f"  {len(citing)} unique citing papers passing min-overlap >= {args.min_overlap}",
          file=sys.stderr)

    # Step 3: fetch metadata for each citing paper (sorted by overlap descending, to prioritise
    # the most relevant if we get rate-limited part-way through).
    sorted_citing = sorted(citing.items(), key=lambda kv: -len(kv[1]))
    candidates: list[dict] = []
    seen_cite_keys = set(cite_keys())
    for i, (citing_doi, citing_seeds) in enumerate(sorted_citing, 1):
        work = fetch_work_metadata(citing_doi)
        if not work:
            time.sleep(REQUEST_DELAY)
            continue
        cand = normalize(work, citing_seeds)
        if cand and cand["candidate_cite_key"] not in seen_cite_keys:
            candidates.append(cand)
            seen_cite_keys.add(cand["candidate_cite_key"])
        time.sleep(REQUEST_DELAY)

    print(f"\nCrossref crawler: {len(candidates)} total candidates", file=sys.stderr)

    if args.dry_run:
        for c in candidates:
            print(json.dumps(c, default=str))
    else:
        out_path = get_inbox_path("crossref")
        emit_jsonl(out_path, candidates)
        print(f"Wrote {out_path.relative_to(REPO_ROOT)}", file=sys.stderr)
        set_cursor("crossref")

    return 0


if __name__ == "__main__":
    sys.exit(main())
