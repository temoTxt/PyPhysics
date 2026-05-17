"""Flag marker-pdf OCR failure patterns in Historical_Converted_Markdown/.

Scans every `Historical_Converted_Markdown/**/*.md` for known artifacts of
marker-pdf converting old physics journal scans:

  1. Lone uppercase `V` adjacent to a number or operator — Wolfram-reserved
     (Vanadium element) AND a common marker-pdf misread for `v` or `\nu`.
  2. Long exponents on Greek/math symbols (`c^{22}`, `\pi^{11}`) — usually
     page-number contamination of an inline equation.
  3. Missing `\hbar` adjacent to `/2` in quantum-mechanics contexts (e.g.,
     spin-1/2 systems): a `1/2` near "spin" without an `\hbar` is often
     OCR dropping the bar.
  4. Page-break artifacts — bare 1–4-digit numbers on their own line, bracketed
     by blank lines, in the middle of running text.
  5. Running-header bleed-through — repeated identical 1-line headers showing up
     between pages.

Flags only; does NOT auto-fix. Output is a per-paper checklist the chapter PR
references in step 5 (QA + retrospective summaries).

Usage:
    uv run python Roadmapping/History/_tools/qa_converted_markdown.py
    uv run python Roadmapping/History/_tools/qa_converted_markdown.py --paper maxwell1865_dynamical_theory
    uv run python Roadmapping/History/_tools/qa_converted_markdown.py --json
"""

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ROOT = REPO_ROOT / "Roadmapping" / "Historical_Converted_Markdown"


def find_lone_V(text: str) -> list[tuple[int, str]]:
    """Uppercase V flanked by digits or operators (most marker-pdf misreads of v/nu)."""
    hits = []
    for lineno, line in enumerate(text.splitlines(), 1):
        if re.search(r"(?<![A-Za-z])V(?=[\^_/\d\\\)\]])", line) or \
           re.search(r"(?<=[=+\-/*\\])V(?![A-Za-z])", line):
            hits.append((lineno, line.strip()))
    return hits


def find_runaway_exponents(text: str) -> list[tuple[int, str]]:
    """Exponents > 9 on a single letter or single TeX command — almost always OCR contamination."""
    hits = []
    pat = re.compile(r"(?:[a-zA-Z]|\\[a-zA-Z]+)\s*\^\s*\{?\s*\d{2,}\s*\}?")
    for lineno, line in enumerate(text.splitlines(), 1):
        for m in pat.finditer(line):
            hits.append((lineno, m.group(0)))
    return hits


def find_missing_hbar(text: str) -> list[tuple[int, str]]:
    """Spin context with `1/2` but no `\\hbar` on the same line."""
    hits = []
    spin_pat = re.compile(r"(?i)\bspin\b.{0,50}1\s*/\s*2")
    for lineno, line in enumerate(text.splitlines(), 1):
        if spin_pat.search(line) and r"\hbar" not in line and "ℏ" not in line:
            hits.append((lineno, line.strip()))
    return hits


def find_page_break_artifacts(text: str) -> list[tuple[int, str]]:
    """Bare 1–4-digit numbers on their own line, sandwiched by blank lines, mid-document."""
    lines = text.splitlines()
    hits = []
    pat = re.compile(r"^\s*\d{1,4}\s*$")
    for i, line in enumerate(lines[1:-1], 1):
        if pat.match(line) and not lines[i - 1].strip() and not lines[i + 1].strip():
            hits.append((i + 1, line.strip()))
    return hits


def find_running_headers(text: str) -> list[str]:
    """1-line strings that repeat ≥3 times — usually running-header bleed-through."""
    counts = Counter(line.strip() for line in text.splitlines() if line.strip() and len(line.strip()) < 80)
    return [s for s, c in counts.items() if c >= 3 and not s.startswith("#") and not s.startswith("|")]


def qa(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    return {
        "file": str(path.relative_to(REPO_ROOT)),
        "lone_V": find_lone_V(text),
        "runaway_exponents": find_runaway_exponents(text),
        "missing_hbar_in_spin_context": find_missing_hbar(text),
        "page_break_artifacts": find_page_break_artifacts(text),
        "running_headers_repeating": find_running_headers(text),
    }


def render_human(results: list[dict]) -> str:
    lines = []
    for r in results:
        any_hits = any(r[k] for k in r if k != "file")
        marker = "⚠" if any_hits else "✅"
        lines.append(f"\n{marker} {r['file']}")
        if not any_hits:
            continue
        for key in ("lone_V", "runaway_exponents", "missing_hbar_in_spin_context", "page_break_artifacts"):
            for entry in r[key]:
                if isinstance(entry, tuple):
                    lineno, snippet = entry
                    lines.append(f"    [{key}] L{lineno}: {snippet[:120]}")
        for h in r["running_headers_repeating"]:
            lines.append(f"    [running_header] repeats: {h[:120]}")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", default=str(DEFAULT_ROOT),
                        help=f"converted-markdown root (default: {DEFAULT_ROOT})")
    parser.add_argument("--paper", help="restrict to a single paper's cite-key directory")
    parser.add_argument("--json", action="store_true", help="emit JSON instead of human format")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        print(f"Root not found: {root}", file=sys.stderr)
        return 0  # not an error — just nothing to scan yet

    files: list[Path] = []
    for sub in ("Primary", "Retrospective"):
        d = root / sub
        if not d.is_dir():
            continue
        if args.paper:
            target = d / args.paper / f"{args.paper}.md"
            if target.exists():
                files.append(target)
        else:
            for paper_dir in sorted(d.iterdir()):
                if paper_dir.is_dir():
                    md = paper_dir / f"{paper_dir.name}.md"
                    if md.exists():
                        files.append(md)

    if not files:
        print(f"No converted-markdown files found under {root}", file=sys.stderr)
        return 0

    results = [qa(f) for f in files]
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print(render_human(results))
    flagged = sum(1 for r in results if any(r[k] for k in r if k != "file"))
    print(f"\n{flagged}/{len(results)} paper(s) flagged.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
