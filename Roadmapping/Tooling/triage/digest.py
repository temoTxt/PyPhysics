"""Render `queue.md` from a list of triage candidates.

Per Decision E: writes raw candidate metadata only — no Claude-API
summaries are baked into the file. Summaries are generated inline when
`/triage-papers` is invoked in Claude Code (Phase 5).

The file is the working surface for human review. Trey edits the
`decision:` field on each candidate (keep / skip / later) and saves;
re-running `build_queue.py --apply-decisions` then pushes `keep`s into
Zotero and appends to `triage_decisions.jsonl`.

Usage:
  uv run python Roadmapping/Tooling/triage/digest.py --from-jsonl queue_full.jsonl
  (typically invoked indirectly via build_queue.py)
"""

import argparse
import json
import sys
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
QUEUE_PATH = REPO_ROOT / "Roadmapping" / "Tooling" / "triage" / "queue.md"

ABSTRACT_TRUNC = 350  # chars; full abstract available in queue_full.jsonl


def render_queue(candidates: list[dict], dropped: dict[str, int],
                  target_date: str | None = None) -> None:
    """Write the human-facing queue.md from candidates + drop stats."""
    target_date = target_date or date.today().isoformat()
    lines: list[str] = []
    lines.append(f"# Triage queue — generated {target_date}\n")
    lines.append(
        f"**{len(candidates)} candidates**, sorted by relevance score descending.\n"
    )
    lines.append(
        "Edit each entry's `decision:` field below to one of `keep` / `skip` / `later`. "
        "Add notes if useful. Then run:\n"
    )
    lines.append("```bash")
    lines.append("uv run python Roadmapping/Tooling/triage/build_queue.py --apply-decisions")
    lines.append("```\n")
    lines.append(
        "to push `keep`s into Zotero, archive `skip`s, leave `later`s in the queue, "
        "and append the audit log to `triage_decisions.jsonl`.\n"
    )

    if dropped:
        lines.append("## Dropped at cap\n")
        for source, n in sorted(dropped.items()):
            label = "total cap" if source == "__total_cap__" else f"{source} per-source"
            lines.append(f"- {n} candidates dropped at {label}")
        lines.append("")

    lines.append("---\n")

    if not candidates:
        lines.append("*No candidates this run.*\n")
    else:
        for idx, c in enumerate(candidates, 1):
            lines.extend(_render_one(idx, c))
            lines.append("---\n")

    QUEUE_PATH.parent.mkdir(parents=True, exist_ok=True)
    QUEUE_PATH.write_text("\n".join(lines), encoding="utf-8")


def _render_one(idx: int, c: dict) -> list[str]:
    lines: list[str] = []
    title = c.get("title", "(no title)").strip()
    cite_key = c.get("candidate_cite_key", "(no cite-key)")
    score = c.get("score", 0.0)
    sources = c.get("sources", [c.get("source", "?")])

    lines.append(f"## {idx}. [{cite_key}] {title}")
    lines.append("")
    authors = c.get("authors") or []
    year = c.get("year") or "?"
    authors_str = ", ".join(authors[:6]) + (" et al." if len(authors) > 6 else "")
    lines.append(f"- **Authors / year:** {authors_str} ({year})")
    doi = c.get("doi") or ""
    arxiv_id = c.get("arxiv_id") or ""
    raw_url = c.get("raw_url") or ""
    if doi:
        lines.append(f"- **DOI:** [`{doi}`](https://doi.org/{doi})")
    if arxiv_id:
        lines.append(f"- **ArXiv:** [`{arxiv_id}`](https://arxiv.org/abs/{arxiv_id})")
    if raw_url and not (doi or arxiv_id):
        lines.append(f"- **URL:** {raw_url}")
    sources_str = ", ".join(f"`{s}`" for s in sources)
    lines.append(f"- **Sources:** {sources_str}")
    lines.append(f"- **Score:** **{score:.2f}**")
    lines.append(f"- **Why surfaced:** {c.get('why_surfaced', '(n/a)')}")
    if c.get("seed_paper"):
        lines.append(f"- **Seed paper:** [[{c['seed_paper']}]]")
    if c.get("tags_inferred"):
        lines.append(f"- **Inferred tags:** {', '.join(c['tags_inferred'][:8])}")

    abstract = (c.get("abstract") or "").strip()
    if abstract:
        if len(abstract) > ABSTRACT_TRUNC:
            abstract = abstract[:ABSTRACT_TRUNC].rsplit(" ", 1)[0] + "…"
        lines.append("")
        lines.append("**Abstract**")
        lines.append("")
        lines.append("> " + abstract.replace("\n", "\n> "))

    lines.append("")
    lines.append("### Decision")
    lines.append("")
    lines.append("<!-- Set one of: keep / skip / later -->")
    lines.append("decision:")
    lines.append("")
    lines.append("### Notes")
    lines.append("")
    lines.append("<!-- optional notes -->")
    lines.append("")
    return lines


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--from-jsonl",
                   default=str(REPO_ROOT / "Roadmapping/Tooling/triage/queue_full.jsonl"),
                   help="JSONL of candidates to render (default: queue_full.jsonl)")
    p.add_argument("--date", default=date.today().isoformat(), help="header date stamp")
    args = p.parse_args()

    jsonl_path = Path(args.from_jsonl)
    if not jsonl_path.exists():
        print(f"Not found: {jsonl_path}", file=sys.stderr)
        return 1

    candidates: list[dict] = []
    with jsonl_path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                candidates.append(json.loads(line))

    render_queue(candidates, dropped={}, target_date=args.date)
    print(f"Wrote {QUEUE_PATH.relative_to(REPO_ROOT)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
