"""Validate a podcast-episode script (YAML schema + wikilinks + scene cues + word count).

Checks:
  1. YAML frontmatter present and parses; required fields: episode, title, era,
     chapter, speakers, target_runtime_min, status.
  2. Every speaker line (`**X:**`) uses one of the canonical personas:
     Historian, Physicist, Experimentalist.
  3. Every `animations_cued:` entry resolves to a real Manim scene file under
     `Roadmapping/Animations/manim_scenes/<name>.py`.
  4. Every `[[wikilink]]` in the body resolves to a real markdown file
     anywhere in the repo (uses the same logic as validate_wikilinks.py).
  5. Rough word-count → runtime cross-check: ~150 wpm spoken; flag if
     `target_runtime_min` is off by more than 30%.

Usage:
    uv run python Roadmapping/History/Podcast/lint_episode.py
    uv run python Roadmapping/History/Podcast/lint_episode.py --episode episode_01_early_electromagnetism.md
"""

import argparse
import re
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
PODCAST_ROOT = REPO_ROOT / "Roadmapping" / "History" / "Podcast"
ANIMATIONS_ROOT = REPO_ROOT / "Roadmapping" / "Animations" / "manim_scenes"

CANONICAL_SPEAKERS = {"Historian", "Physicist", "Experimentalist"}
REQUIRED_FIELDS = ("episode", "title", "era", "chapter", "speakers", "target_runtime_min", "status")
SPEAKER_LINE_RE = re.compile(r"^\*\*([A-Z][a-zA-Z]+)\:\*\*", re.MULTILINE)
WIKILINK_RE = re.compile(r"\[\[([^\]\[\|]+?)(?:\|[^\]]+?)?\]\]")
FRONTMATTER_RE = re.compile(r"^---\n(.*?)\n---\n(.*)$", re.DOTALL)
WORDS_PER_MINUTE = 150


def index_markdown(repo_root: Path) -> set[str]:
    out: set[str] = set()
    for md in repo_root.rglob("*.md"):
        if ".venv" in md.parts:
            continue
        out.add(md.stem)
    return out


def lint(episode: Path, md_index: set[str]) -> list[str]:
    errors: list[str] = []
    text = episode.read_text(encoding="utf-8")

    m = FRONTMATTER_RE.match(text)
    if not m:
        return [f"{episode.name}: no YAML frontmatter"]

    try:
        meta = yaml.safe_load(m.group(1)) or {}
    except yaml.YAMLError as e:
        return [f"{episode.name}: YAML parse error: {e}"]

    body = m.group(2)

    # 1. required fields
    for f in REQUIRED_FIELDS:
        if meta.get(f) in (None, "", []):
            errors.append(f"missing required field: {f}")

    # 2. speakers ∈ canonical cast (frontmatter)
    declared = set(meta.get("speakers") or [])
    invalid_decl = declared - CANONICAL_SPEAKERS
    if invalid_decl:
        errors.append(f"frontmatter speakers contains non-canonical: {sorted(invalid_decl)}")

    # 2b. speaker lines in body
    body_speakers = set(SPEAKER_LINE_RE.findall(body))
    invalid_body = body_speakers - CANONICAL_SPEAKERS
    if invalid_body:
        errors.append(f"body uses non-canonical speakers: {sorted(invalid_body)}")
    unused_decl = declared - body_speakers
    if unused_decl:
        errors.append(f"frontmatter declares speakers not used in body: {sorted(unused_decl)}")

    # 3. animations_cued resolves
    for name in (meta.get("animations_cued") or []):
        scene = ANIMATIONS_ROOT / f"{name}.py"
        if not scene.exists():
            errors.append(f"animations_cued references missing scene: {name}.py")

    # 4. wikilinks — strip fenced + inline code spans first
    body_no_code = re.sub(r"```.*?```", "", body, flags=re.DOTALL)
    body_no_code = re.sub(r"`[^`\n]+`", "", body_no_code)
    for raw in WIKILINK_RE.findall(body_no_code):
        name = raw.split("#", 1)[0].strip().split("/")[-1]
        if name not in md_index:
            errors.append(f"broken wikilink: [[{raw}]]")

    # 5. word-count → runtime sanity
    word_count = len(re.findall(r"\b\w+\b", body))
    expected_min = word_count / WORDS_PER_MINUTE
    target = meta.get("target_runtime_min", 0) or 0
    if target and (expected_min < 0.7 * target or expected_min > 1.3 * target):
        errors.append(
            f"runtime mismatch: ~{word_count} words → ~{expected_min:.1f} min, "
            f"target_runtime_min={target}"
        )

    return errors


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--episode", help="restrict to one episode filename under Podcast/")
    p.add_argument("--quiet", action="store_true", help="only print failing episodes")
    args = p.parse_args()

    episodes = (
        [PODCAST_ROOT / args.episode] if args.episode
        else sorted(PODCAST_ROOT.glob("episode_*.md"))
    )
    if not episodes:
        print("No episodes found to lint.", file=sys.stderr)
        return 0

    md_index = index_markdown(REPO_ROOT)
    failures = 0
    for ep in episodes:
        if not ep.exists():
            print(f"{ep.name}: not found", file=sys.stderr)
            failures += 1
            continue
        errs = lint(ep, md_index)
        if errs:
            failures += 1
            print(f"\n⚠ {ep.name}:")
            for e in errs:
                print(f"    - {e}")
        elif not args.quiet:
            print(f"✅ {ep.name}")

    if failures:
        print(f"\n{failures}/{len(episodes)} episode(s) failed lint.", file=sys.stderr)
        return 1
    if not args.quiet:
        print(f"\nAll {len(episodes)} episode(s) passed lint ✅")
    return 0


if __name__ == "__main__":
    sys.exit(main())
