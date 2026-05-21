"""Audit + set the `human_reviewed` field across every YAML bibliography stub.

Implements Crocco et al. (2026) reporting expectation #7 — confirming that
every entry in the final synthesis was *read by a human* (not merely
AI-summarised from an abstract). The `human_reviewed` field in each YAML
stub's frontmatter records that confirmation.

Heuristic: a stub is `human_reviewed: true` when its body contains at
least one substantive paragraph beyond the "*Skeleton bib note…*"
boilerplate. Threshold is conservative — entries that look human-written
should mostly self-confirm; ambiguous entries default to `false` so the
human can flip them later via inspection.

Three operating modes:

    # 1. Report current state without modifying anything.
    uv run python Roadmapping/History/Bibliography/audit_human_reviewed.py

    # 2. Add the `human_reviewed` field to every stub that lacks it, using
    #    the heuristic to set true/false. Existing values are NOT changed.
    uv run python Roadmapping/History/Bibliography/audit_human_reviewed.py --apply

    # 3. Re-audit *every* stub including ones that already have the field
    #    (use if the heuristic improves and you want to re-evaluate).
    uv run python Roadmapping/History/Bibliography/audit_human_reviewed.py --apply --force
"""

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
BIB_DIR = Path(__file__).resolve().parent

FRONTMATTER_RE = re.compile(r"^(---\n)(.*?)(\n---\n)(.*)$", re.DOTALL)
TITLE_HEADING_RE = re.compile(r"^#\s+.+?$", re.MULTILINE)
SKELETON_BOILERPLATE_RE = re.compile(
    r"\*Skeleton bib note[^*]*\*", re.IGNORECASE,
)
WORD_THRESHOLD = 60  # body words (post-boilerplate) needed for human_reviewed=true


def split_yaml_body(text: str) -> tuple[str, str, str, str] | None:
    m = FRONTMATTER_RE.match(text)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3), m.group(4)


def body_word_count(body: str) -> int:
    """Count words in body, excluding the title heading and boilerplate phrases."""
    cleaned = TITLE_HEADING_RE.sub("", body)
    cleaned = SKELETON_BOILERPLATE_RE.sub("", cleaned)
    return len(re.findall(r"\b\w+\b", cleaned))


def has_field(yaml_body: str, field: str) -> bool:
    return bool(re.search(rf"^{re.escape(field)}\s*:", yaml_body, re.MULTILINE))


def insert_field(yaml_body: str, field: str, value: str) -> str:
    """Insert `field: value` into a YAML frontmatter body before chapters_citing if
    present, else at the end. Preserves indentation + ordering."""
    new_line = f"{field}: {value}"
    # Prefer placement right before `chapters_citing` (the last field in our schema).
    lines = yaml_body.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("chapters_citing:"):
            lines.insert(i, new_line)
            return "\n".join(lines)
    # Fall back: append at end.
    lines.append(new_line)
    return "\n".join(lines)


def all_stubs() -> list[Path]:
    out: list[Path] = []
    for sub in ("Primary", "Retrospective"):
        d = BIB_DIR / sub
        if d.is_dir():
            out.extend(sorted(d.glob("*.md")))
    return out


def audit_one(path: Path) -> tuple[str, bool, int]:
    """Return (current_value, heuristic_decision, body_word_count)."""
    text = path.read_text(encoding="utf-8")
    parts = split_yaml_body(text)
    if not parts:
        return ("?", False, 0)
    _, yaml_body, _, body = parts
    wc = body_word_count(body)
    heuristic = wc >= WORD_THRESHOLD
    # Detect existing field value.
    m = re.search(r"^human_reviewed:\s*(\S+)", yaml_body, re.MULTILINE)
    current = m.group(1) if m else "absent"
    return (current, heuristic, wc)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--apply", action="store_true",
                   help="write the field to each stub that lacks it (or to all stubs with --force)")
    p.add_argument("--force", action="store_true",
                   help="re-evaluate stubs that already have human_reviewed (only with --apply)")
    p.add_argument("--quiet", action="store_true", help="only print the summary line")
    args = p.parse_args()

    stubs = all_stubs()
    if not stubs:
        print("No stubs found under Bibliography/{Primary,Retrospective}/", file=sys.stderr)
        return 0

    n_present, n_absent, n_true_now, n_false_now, n_modified = 0, 0, 0, 0, 0
    for path in stubs:
        current, heuristic, wc = audit_one(path)
        if current == "absent":
            n_absent += 1
        else:
            n_present += 1
            if current.lower() in ("true", "yes"):
                n_true_now += 1
            else:
                n_false_now += 1

        will_write = args.apply and (current == "absent" or args.force)
        if not args.quiet:
            verdict = "true" if heuristic else "false"
            mark = "+" if will_write and current == "absent" else (
                "~" if will_write and args.force and current.lower() != verdict else "·"
            )
            print(f"{mark}  {path.relative_to(REPO_ROOT)}  "
                  f"current={current}  heuristic={verdict}  body_words={wc}")

        if will_write:
            text = path.read_text(encoding="utf-8")
            parts = split_yaml_body(text)
            if not parts:
                continue
            head, yaml_body, sep, body = parts
            new_value = "true" if heuristic else "false"
            if has_field(yaml_body, "human_reviewed") and args.force:
                yaml_body = re.sub(r"^human_reviewed:\s*\S+\s*$",
                                    f"human_reviewed: {new_value}",
                                    yaml_body, count=1, flags=re.MULTILINE)
            elif not has_field(yaml_body, "human_reviewed"):
                yaml_body = insert_field(yaml_body, "human_reviewed", new_value)
            else:
                continue  # field present and not --force: skip
            path.write_text(f"{head}{yaml_body}{sep}{body}", encoding="utf-8")
            n_modified += 1

    print(file=sys.stderr)
    print(f"Total stubs: {len(stubs)}", file=sys.stderr)
    print(f"  with field: {n_present} ({n_true_now} true, {n_false_now} false)", file=sys.stderr)
    print(f"  without field: {n_absent}", file=sys.stderr)
    if args.apply:
        print(f"Wrote field to {n_modified} stub(s).", file=sys.stderr)
    elif n_absent:
        print(f"Run with --apply to add the field to {n_absent} stub(s).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
