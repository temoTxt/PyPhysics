"""Walk the repo's markdown files and validate every [[wikilink]].

Checks every Obsidian-style wikilink under:
- `Roadmapping/History/**/*.md`
- `Roadmapping/Equation_Verification/**/*.md`
- `Roadmapping/Animations/README.md`

A wikilink `[[name]]`, `[[name#heading]]`, or `[[name|alias]]` resolves to a
markdown file anywhere in the repo whose stem matches `name`. If a `#heading`
anchor is present, the script also verifies that heading exists in the target.

Exit non-zero if any broken link is found — suitable for a pre-commit hook.

Usage:
    uv run python Roadmapping/History/_tools/validate_wikilinks.py
    uv run python Roadmapping/History/_tools/validate_wikilinks.py --quiet
    uv run python Roadmapping/History/_tools/validate_wikilinks.py --root .
"""

import argparse
import re
import sys
from pathlib import Path

# Inner capture: anything up to the alias separator `|` or `\|` (escaped in tables)
# or the closing `]]`. Strip a trailing `\` left over from an escaped pipe later.
WIKILINK_RE = re.compile(r"\[\[([^\]\[\|]+?)(?:\\?\|[^\]]+?)?\]\]")
HEADING_RE = re.compile(r"^#+\s+(.+?)\s*$", re.MULTILINE)

DEFAULT_SCAN = [
    "Roadmapping/History",
    "Roadmapping/Equation_Verification",
    "Roadmapping/Animations/README.md",
]


def slugify(heading: str) -> str:
    """Match Obsidian's relaxed heading match — case-insensitive, whitespace folded."""
    return re.sub(r"\s+", " ", heading.strip().lower())


def index_files(repo_root: Path) -> dict[str, Path]:
    """Map every .md file's stem to its path. Last-write-wins for collisions."""
    out: dict[str, Path] = {}
    for md in repo_root.rglob("*.md"):
        if ".venv" in md.parts or "node_modules" in md.parts:
            continue
        out[md.stem] = md
    return out


def headings_in(md_path: Path) -> set[str]:
    text = md_path.read_text(encoding="utf-8")
    return {slugify(h) for h in HEADING_RE.findall(text)}


def collect_scan_files(repo_root: Path, scan_paths: list[str]) -> list[Path]:
    files: list[Path] = []
    for rel in scan_paths:
        p = repo_root / rel
        if p.is_file() and p.suffix == ".md":
            files.append(p)
        elif p.is_dir():
            files.extend(sorted(p.rglob("*.md")))
    # Skip template files — their placeholders are intentional, not real wikilinks.
    files = [f for f in files if not f.name.startswith("_template_")]
    return files


def strip_code(text: str) -> str:
    """Replace fenced code blocks and inline code spans with blanks of equal length,
    so wikilinks inside code (illustrative syntax) aren't validated."""
    # Fenced blocks: ```...``` (multiline). Preserve newlines so line numbers stay aligned.
    def fenced(m: re.Match) -> str:
        return re.sub(r"[^\n]", " ", m.group(0))
    text = re.sub(r"```.*?```", fenced, text, flags=re.DOTALL)
    text = re.sub(r"`[^`\n]+`", lambda m: " " * len(m.group(0)), text)
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", default=str(Path(__file__).resolve().parents[3]),
                        help="repo root (default: auto-detected)")
    parser.add_argument("--scan", action="append", default=None,
                        help="paths to scan (relative to root); repeatable; default scans History + Equation_Verification + Animations/README.md")
    parser.add_argument("--quiet", action="store_true", help="only print broken-link summary")
    args = parser.parse_args()

    repo_root = Path(args.root).resolve()
    scan_paths = args.scan or DEFAULT_SCAN
    files = collect_scan_files(repo_root, scan_paths)
    index = index_files(repo_root)

    broken: list[tuple[Path, int, str, str]] = []
    checked = 0
    for md in files:
        text = strip_code(md.read_text(encoding="utf-8"))
        for lineno, line in enumerate(text.splitlines(), 1):
            for m in WIKILINK_RE.finditer(line):
                checked += 1
                # Strip trailing backslash leftover from escaped-pipe table syntax.
                raw = m.group(1).rstrip("\\").strip()
                name, _, anchor = raw.partition("#")
                # Path-style wikilinks: `[[folder/sub/file]]` → resolve by basename.
                name = name.strip().split("/")[-1]
                anchor = anchor.rstrip("\\").strip()
                target = index.get(name)
                if target is None:
                    broken.append((md, lineno, raw, "no such file"))
                    continue
                if anchor:
                    if slugify(anchor) not in headings_in(target):
                        broken.append((md, lineno, raw, f"no heading in {target.name}"))

    if not args.quiet:
        print(f"Scanned {len(files)} files, {checked} wikilinks, "
              f"index covers {len(index)} markdown files.")

    if broken:
        print(f"\n{len(broken)} broken wikilink(s):", file=sys.stderr)
        for src, lineno, raw, reason in broken:
            rel = src.relative_to(repo_root)
            print(f"  {rel}:{lineno}  [[{raw}]]  — {reason}", file=sys.stderr)
        return 1
    if not args.quiet:
        print("All wikilinks resolve ✅")
    return 0


if __name__ == "__main__":
    sys.exit(main())
