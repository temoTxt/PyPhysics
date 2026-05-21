"""Merge all crawler outputs into a single triage queue + render queue.md.

Reads every `_inbox/<date>-from_*.jsonl` (default: today's date), deduplicates
candidates by DOI/ArXiv-id/cite-key, applies the per-source cap (10) and
total cap (30) per Decision B, sorts by score descending, and renders
`queue.md` via `digest.py`.

Outputs:
  Roadmapping/Tooling/triage/queue.md           — human-facing weekly UI
  Roadmapping/Tooling/triage/queue_full.jsonl   — pre-render audit trail (gitignored)

Usage:
  uv run python Roadmapping/Tooling/triage/build_queue.py
  uv run python Roadmapping/Tooling/triage/build_queue.py --date 2026-05-21
  uv run python Roadmapping/Tooling/triage/build_queue.py --no-render
  uv run python Roadmapping/Tooling/triage/build_queue.py --crawl-first
      (also invokes the four crawlers before merging)
"""

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
INBOX_DIR = REPO_ROOT / "Roadmapping" / "Tooling" / "triage" / "_inbox"
QUEUE_PATH = REPO_ROOT / "Roadmapping" / "Tooling" / "triage" / "queue.md"
QUEUE_AUDIT_PATH = REPO_ROOT / "Roadmapping" / "Tooling" / "triage" / "queue_full.jsonl"
CRAWL_DIR = REPO_ROOT / "Roadmapping" / "Tooling" / "crawl"

DEFAULT_PER_SOURCE_CAP = 10
DEFAULT_TOTAL_CAP = 30

CRAWLERS = ["from_zotero", "from_s2", "from_arxiv", "from_crossref"]


def run_crawlers() -> None:
    """Shell out to each crawler in turn. Errors are logged but non-fatal."""
    for name in CRAWLERS:
        script = CRAWL_DIR / f"{name}.py"
        cmd = ["uv", "run", "python", str(script)]
        print(f"$ {' '.join(cmd)}", file=sys.stderr)
        res = subprocess.run(cmd, cwd=REPO_ROOT)
        if res.returncode != 0:
            print(f"  {name} exited {res.returncode} — continuing", file=sys.stderr)


def load_inbox(target_date: str) -> dict[str, list[dict]]:
    """Read every `<date>-from_*.jsonl` for the given date. Returns {source: [candidates]}."""
    by_source: dict[str, list[dict]] = defaultdict(list)
    if not INBOX_DIR.is_dir():
        return by_source
    for path in sorted(INBOX_DIR.glob(f"{target_date}-from_*.jsonl")):
        # Filename: 2026-05-21-from_s2.jsonl → source 's2'
        stem = path.stem  # 2026-05-21-from_s2
        source = stem.rsplit("from_", 1)[-1] if "from_" in stem else "unknown"
        with path.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    by_source[source].append(json.loads(line))
                except json.JSONDecodeError:
                    continue
    return by_source


def dedupe(by_source: dict[str, list[dict]]) -> list[dict]:
    """Merge candidates across sources, keyed by DOI > ArXiv-id > cite-key.

    When the same paper surfaces from multiple sources, merge into one
    candidate with `sources: [..., ...]`, summed score (capped at 1.0),
    and concatenated why_surfaced explanations.
    """
    merged: dict[str, dict] = {}
    for source, cands in by_source.items():
        for c in cands:
            key = c.get("doi") or c.get("arxiv_id") or c.get("candidate_cite_key", "")
            if not key:
                continue
            if key in merged:
                existing = merged[key]
                # Track all sources.
                existing_sources = set(existing.get("sources", [existing.get("source", "")]))
                existing_sources.add(source)
                existing["sources"] = sorted(existing_sources)
                # Cumulative score — papers that show up from multiple sources matter more.
                existing["score"] = min(1.0, existing.get("score", 0.0) + c.get("score", 0.0) * 0.3)
                # Concatenate why_surfaced.
                existing_why = existing.get("why_surfaced", "")
                if c.get("why_surfaced") and c["why_surfaced"] not in existing_why:
                    existing["why_surfaced"] = f"{existing_why}; ALSO via {source}: {c['why_surfaced']}"
                # Prefer entries with more metadata.
                if not existing.get("abstract") and c.get("abstract"):
                    existing["abstract"] = c["abstract"]
                if not existing.get("year") and c.get("year"):
                    existing["year"] = c["year"]
            else:
                c = dict(c)
                c.setdefault("sources", [source])
                merged[key] = c

    return list(merged.values())


def apply_caps(by_source: dict[str, list[dict]], total_cap: int, per_source_cap: int
               ) -> tuple[list[dict], dict[str, int]]:
    """Apply per-source cap, then dedupe, then total cap. Returns (kept, dropped_per_source)."""
    capped_per_source: dict[str, list[dict]] = {}
    dropped_per_source: dict[str, int] = {}
    for source, cands in by_source.items():
        sorted_cands = sorted(cands, key=lambda c: -c.get("score", 0.0))
        if len(sorted_cands) > per_source_cap:
            dropped_per_source[source] = len(sorted_cands) - per_source_cap
            sorted_cands = sorted_cands[:per_source_cap]
        capped_per_source[source] = sorted_cands

    merged = dedupe(capped_per_source)
    merged.sort(key=lambda c: -c.get("score", 0.0))

    if len(merged) > total_cap:
        dropped_per_source["__total_cap__"] = len(merged) - total_cap
        merged = merged[:total_cap]
    return merged, dropped_per_source


def write_queue_audit(candidates: list[dict]) -> None:
    QUEUE_AUDIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with QUEUE_AUDIT_PATH.open("w", encoding="utf-8") as f:
        for c in candidates:
            f.write(json.dumps(c, ensure_ascii=False, default=str) + "\n")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--date", default=date.today().isoformat(),
                   help="inbox date to read (YYYY-MM-DD; default today)")
    p.add_argument("--per-source-cap", type=int, default=DEFAULT_PER_SOURCE_CAP,
                   help=f"max candidates per source before dedup (default {DEFAULT_PER_SOURCE_CAP})")
    p.add_argument("--total-cap", type=int, default=DEFAULT_TOTAL_CAP,
                   help=f"max candidates in the final queue (default {DEFAULT_TOTAL_CAP})")
    p.add_argument("--no-render", action="store_true",
                   help="skip queue.md rendering; only write queue_full.jsonl")
    p.add_argument("--crawl-first", action="store_true",
                   help="run all four crawlers before merging")
    args = p.parse_args()

    if args.crawl_first:
        run_crawlers()

    by_source = load_inbox(args.date)
    total_raw = sum(len(v) for v in by_source.values())
    if total_raw == 0:
        print(f"No crawler outputs found for {args.date} under {INBOX_DIR.relative_to(REPO_ROOT)}.",
              file=sys.stderr)
        print("Run a crawler (or `build_queue.py --crawl-first`) first.", file=sys.stderr)
        return 0

    print(f"Loaded {total_raw} raw candidates across {len(by_source)} sources:", file=sys.stderr)
    for source, cands in sorted(by_source.items()):
        print(f"  {source}: {len(cands)}", file=sys.stderr)

    candidates, dropped = apply_caps(by_source, args.total_cap, args.per_source_cap)

    print(f"\nAfter dedup + caps: {len(candidates)} candidates in queue.", file=sys.stderr)
    for source, n in sorted(dropped.items()):
        label = "total cap" if source == "__total_cap__" else f"{source} per-source cap"
        print(f"  {n} dropped at {label}", file=sys.stderr)

    write_queue_audit(candidates)
    print(f"Audit written: {QUEUE_AUDIT_PATH.relative_to(REPO_ROOT)}", file=sys.stderr)

    if not args.no_render:
        from importlib import util
        digest_path = Path(__file__).resolve().parent / "digest.py"
        spec = util.spec_from_file_location("digest", digest_path)
        digest_mod = util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(digest_mod)
        digest_mod.render_queue(candidates, dropped, target_date=args.date)
        print(f"Wrote {QUEUE_PATH.relative_to(REPO_ROOT)}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
