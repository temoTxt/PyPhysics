"""Compact campaign dashboard — per-chapter, bib + PDFs + scenes + podcast + links.

Reads every chapter file, every bibliography YAML, the Acquisition_Tracker, the
Animations scene index, and every podcast episode; prints a per-chapter table:

    Ch | Title | Bib stubs | PDFs acquired | Scenes | Podcast | Wikilinks
    01 | Early EM 1800-1860      | 13/13 | 11/13 🟢 | 1/1 | draft | ✅

So a future agent (or Trey, mid-campaign) can see at a glance what's outstanding
without spelunking. Read-only — never writes.

Usage:
    uv run python Roadmapping/History/_tools/chapter_status.py
    uv run python Roadmapping/History/_tools/chapter_status.py --chapter 01
    uv run python Roadmapping/History/_tools/chapter_status.py --json
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
HISTORY_ROOT = REPO_ROOT / "Roadmapping" / "History"
BIB_ROOT = HISTORY_ROOT / "Bibliography"
PODCAST_ROOT = HISTORY_ROOT / "Podcast"
ANIMATIONS_ROOT = REPO_ROOT / "Roadmapping" / "Animations" / "manim_scenes"
FORWARD_ROOT = HISTORY_ROOT / "Forward"

FRONTMATTER_RE = re.compile(r"^---\n(.*?)\n---", re.DOTALL)

# era → expected chapter file glob; also captures the forward chapters.
ERA_CHAPTERS = {
    "1800-1860": "01_*.md",
    "1860-1900": "02_*.md",
    "1900-1925": "03_*.md",
    "1925-1948": "04_*.md",
    "1948-1965": "05_*.md",
}
FORWARD_CHAPTERS = ["07_*.md", "08_*.md", "09_*.md"]


def parse_frontmatter(md: Path) -> dict:
    text = md.read_text(encoding="utf-8")
    m = FRONTMATTER_RE.match(text)
    if not m:
        return {}
    try:
        return yaml.safe_load(m.group(1)) or {}
    except yaml.YAMLError:
        return {}


def bib_for_era(era: str) -> list[dict]:
    out = []
    for sub in ("Primary", "Retrospective"):
        d = BIB_ROOT / sub
        if not d.is_dir():
            continue
        for md in sorted(d.glob("*.md")):
            meta = parse_frontmatter(md)
            if str(meta.get("era", "")).strip() == era:
                out.append(meta)
    return out


def find_chapter_file(pattern: str) -> Path | None:
    for root in (HISTORY_ROOT, FORWARD_ROOT):
        matches = sorted(root.glob(pattern))
        if matches:
            return matches[0]
    return None


def find_podcast_episode(chapter_num: int) -> Path | None:
    matches = sorted(PODCAST_ROOT.glob(f"episode_{chapter_num:02d}_*.md"))
    return matches[0] if matches else None


def scenes_for_chapter(chapter_meta: dict) -> list[tuple[str, bool]]:
    """Returns [(scene_name, rendered)] from the chapter's `animations:` field."""
    out = []
    for name in chapter_meta.get("animations", []) or []:
        scene_file = ANIMATIONS_ROOT / f"{name}.py"
        out.append((name, scene_file.exists()))
    return out


def validate_links_count() -> int | None:
    """Returns count of broken wikilinks; None if validator unavailable."""
    validator = HISTORY_ROOT / "_tools" / "validate_wikilinks.py"
    if not validator.exists():
        return None
    res = subprocess.run(
        ["uv", "run", "python", str(validator), "--quiet"],
        capture_output=True, text=True,
    )
    if res.returncode == 0:
        return 0
    return res.stderr.count("[[")


def chapter_row(num: int, pattern: str, era: str | None) -> dict:
    chapter_file = find_chapter_file(pattern)
    chapter_meta = parse_frontmatter(chapter_file) if chapter_file else {}
    bib = bib_for_era(era) if era else []
    pdfs_acquired = sum(
        1 for m in bib if m.get("pdf_status") in ("out_of_copyright_public", "acquired")
    )
    scenes = scenes_for_chapter(chapter_meta) if chapter_meta else []
    scenes_rendered = sum(1 for _, ok in scenes if ok)
    podcast = find_podcast_episode(num)
    podcast_meta = parse_frontmatter(podcast) if podcast else {}
    return {
        "chapter": f"{num:02d}",
        "title": (chapter_meta.get("title") or (chapter_file.stem if chapter_file else "—")),
        "status": chapter_meta.get("status", "—"),
        "bib_total": len(bib),
        "pdfs_acquired": pdfs_acquired,
        "scenes_total": len(scenes),
        "scenes_rendered": scenes_rendered,
        "podcast_status": podcast_meta.get("status", "—") if podcast else "—",
    }


def collect(only_chapter: str | None) -> list[dict]:
    rows = []
    for i, (era, pattern) in enumerate(ERA_CHAPTERS.items(), 1):
        if only_chapter and only_chapter != f"{i:02d}":
            continue
        rows.append(chapter_row(i, pattern, era))
    # Chapter 6 (synthesis) has no era; 7-9 forward chapters.
    for i, pattern in [(6, "06_*.md")] + list(zip([7, 8, 9], FORWARD_CHAPTERS)):
        if only_chapter and only_chapter != f"{i:02d}":
            continue
        rows.append(chapter_row(i, pattern, None))
    return rows


def render(rows: list[dict], broken_links: int | None) -> str:
    lines = ["Ch | Title                                          | Status | Bib   | PDFs  | Scenes | Podcast",
             "---|------------------------------------------------|--------|-------|-------|--------|--------"]
    for r in rows:
        title = r["title"][:46]
        bib = f"{r['bib_total']:>3}"
        pdfs = f"{r['pdfs_acquired']:>2}/{r['bib_total']:<2}" if r["bib_total"] else "—   "
        scenes = f"{r['scenes_rendered']}/{r['scenes_total']}" if r["scenes_total"] else "—  "
        lines.append(f"{r['chapter']} | {title:<46} | {r['status']:<6} | {bib:<5} | {pdfs:<5} | {scenes:<6} | {r['podcast_status']}")
    lines.append("")
    if broken_links is None:
        lines.append("Wikilink validation: skipped (validate_wikilinks.py not found)")
    elif broken_links == 0:
        lines.append("Wikilink validation: ✅ all resolve")
    else:
        lines.append(f"Wikilink validation: ⚠ {broken_links} broken link(s) — run validate_wikilinks.py for details")
    return "\n".join(lines)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chapter", help="restrict to one chapter (e.g. 01)")
    p.add_argument("--json", action="store_true", help="emit JSON instead of human table")
    args = p.parse_args()

    rows = collect(args.chapter)
    if args.json:
        print(json.dumps({"chapters": rows, "broken_links": validate_links_count()}, indent=2))
    else:
        print(render(rows, validate_links_count()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
