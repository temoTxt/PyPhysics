"""Regenerate Historical_Papers/Acquisition_Tracker.md from bibliography YAML.

Reads every `Bibliography/**/*.md`, pulls `cite_key`, `year`, `era`, `type`,
`pdf_status`, `pdf_path` from the YAML frontmatter, and emits a sorted
status table to `Historical_Papers/Acquisition_Tracker.md`.

Idempotent: re-run after any batch of new bib stubs. Existing prose above and
below the auto-generated table is preserved (delimited by HTML comment markers).

Usage:
    uv run python Roadmapping/History/Bibliography/update_acquisition_tracker.py
    uv run python Roadmapping/History/Bibliography/update_acquisition_tracker.py --dry-run
"""

import argparse
import re
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BIB = REPO_ROOT / "Roadmapping" / "History" / "Bibliography"
DEFAULT_TRACKER = REPO_ROOT / "Roadmapping" / "Historical_Papers" / "Acquisition_Tracker.md"

FRONTMATTER_RE = re.compile(r"^---\n(.*?)\n---", re.DOTALL)
TABLE_BEGIN = "<!-- AUTO-TABLE-BEGIN -->"
TABLE_END = "<!-- AUTO-TABLE-END -->"

STATUS_BADGE = {
    "out_of_copyright_public": "🟢 public-PD (committed)",
    "acquired": "🟡 acquired (local only)",
    "pending": "⚪ pending",
    "unavailable": "🔴 unavailable",
}


def parse_frontmatter(md_path: Path) -> dict | None:
    text = md_path.read_text(encoding="utf-8")
    m = FRONTMATTER_RE.match(text)
    if not m:
        return None
    try:
        return yaml.safe_load(m.group(1)) or {}
    except yaml.YAMLError as e:
        print(f"YAML error in {md_path}: {e}", file=sys.stderr)
        return None


def collect(bib_root: Path) -> list[dict]:
    rows = []
    for sub in ("Primary", "Retrospective"):
        d = bib_root / sub
        if not d.is_dir():
            continue
        for md in sorted(d.glob("*.md")):
            meta = parse_frontmatter(md)
            if not meta:
                continue
            rows.append({
                "cite_key": meta.get("cite_key", md.stem),
                "year": meta.get("year", ""),
                "era": meta.get("era", ""),
                "type": meta.get("type", sub.lower()),
                "pdf_status": meta.get("pdf_status", "pending"),
                "pdf_path": meta.get("pdf_path", ""),
                "title": meta.get("title", ""),
            })
    rows.sort(key=lambda r: (str(r["era"]), str(r["year"]), r["cite_key"]))
    return rows


def render_table(rows: list[dict]) -> str:
    if not rows:
        return "*No bibliography notes yet.*\n"

    header = (
        "| Cite-key | Year | Era | Type | PDF status | Path |\n"
        "|---|---:|---|---|---|---|\n"
    )
    body = []
    for r in rows:
        badge = STATUS_BADGE.get(r["pdf_status"], r["pdf_status"])
        path = f"`{r['pdf_path']}`" if r["pdf_path"] else "—"
        body.append(
            f"| [[{r['cite_key']}]] | {r['year'] or '—'} | {r['era'] or '—'} | "
            f"{r['type']} | {badge} | {path} |"
        )

    counts = {
        "total": len(rows),
        "public_pd": sum(1 for r in rows if r["pdf_status"] == "out_of_copyright_public"),
        "acquired": sum(1 for r in rows if r["pdf_status"] == "acquired"),
        "pending": sum(1 for r in rows if r["pdf_status"] == "pending"),
        "unavailable": sum(1 for r in rows if r["pdf_status"] == "unavailable"),
    }
    summary = (
        f"**Summary**: {counts['total']} entries · "
        f"🟢 {counts['public_pd']} public-PD · "
        f"🟡 {counts['acquired']} acquired · "
        f"⚪ {counts['pending']} pending · "
        f"🔴 {counts['unavailable']} unavailable\n\n"
    )
    return summary + header + "\n".join(body) + "\n"


def update_tracker(tracker_path: Path, new_table: str) -> str:
    """Return the new tracker file contents, preserving prose outside the markers."""
    if tracker_path.exists():
        existing = tracker_path.read_text(encoding="utf-8")
        if TABLE_BEGIN in existing and TABLE_END in existing:
            before, _, rest = existing.partition(TABLE_BEGIN)
            _, _, after = rest.partition(TABLE_END)
            return f"{before}{TABLE_BEGIN}\n{new_table}{TABLE_END}{after}"
        # No markers — append.
        return f"{existing.rstrip()}\n\n{TABLE_BEGIN}\n{new_table}{TABLE_END}\n"
    # New file — bootstrap with framing.
    return (
        "# Historical Papers — Acquisition Tracker\n\n"
        "Auto-regenerated from `Roadmapping/History/Bibliography/**/*.md` via "
        "`update_acquisition_tracker.py`. Edit prose above/below the auto-table markers; "
        "the table itself is overwritten on each run.\n\n"
        "**Two-tier PDF policy:**\n"
        "- `pdf_status: out_of_copyright_public` → PDF committed to `Historical_Papers/<Primary|Retrospective>/` via `git add -f`.\n"
        "- `pdf_status: acquired` → PDF kept local only (in-copyright). Markdown conversion committed under "
        "`Historical_Converted_Markdown/` per fair-use academic quotation.\n"
        "- `pdf_status: pending` → not yet sourced.\n"
        "- `pdf_status: unavailable` → no accessible copy found.\n\n"
        f"{TABLE_BEGIN}\n{new_table}{TABLE_END}\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bib-root", default=str(DEFAULT_BIB),
                        help=f"bibliography root (default: {DEFAULT_BIB})")
    parser.add_argument("--tracker", default=str(DEFAULT_TRACKER),
                        help=f"output tracker markdown (default: {DEFAULT_TRACKER})")
    parser.add_argument("--dry-run", action="store_true", help="print to stdout instead of writing")
    args = parser.parse_args()

    rows = collect(Path(args.bib_root).resolve())
    table = render_table(rows)
    new_contents = update_tracker(Path(args.tracker).resolve(), table)

    if args.dry_run:
        sys.stdout.write(new_contents)
    else:
        tracker_path = Path(args.tracker).resolve()
        tracker_path.parent.mkdir(parents=True, exist_ok=True)
        tracker_path.write_text(new_contents, encoding="utf-8")
        print(f"Wrote {len(rows)} entries to {tracker_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
