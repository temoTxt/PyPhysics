"""Bridge Zotero's Better-BibTeX export to the YAML bibliography pipeline.

Reads `Roadmapping/History/Bibliography/zotero.bib` (the Better BibTeX
auto-export from the Zotero local library) and, for each entry whose
cite-key has no corresponding `Primary/<cite_key>.md` or
`Retrospective/<cite_key>.md`, scaffolds a YAML stub matching the schema
documented in `Bibliography/README.md`.

For entries that *do* already exist as a YAML stub, the script compares
the BBT-exported metadata against the YAML frontmatter and flags any
*drift* (different title, author, year, DOI) for manual reconciliation
rather than silently overwriting.

The Zotero entry's tags are mapped into the YAML schema:
  - `repo:retrospective`        →  type: retrospective       (default: primary)
  - `era/<YYYY-YYYY|forward>`   →  era
  - `thread/<thread>`           →  tags includes #thread/<thread>
  - any other tag verbatim      →  tags includes that tag

Zotero is treated as a *capture-and-store layer*, not a replacement: the
YAML stubs under `Bibliography/{Primary,Retrospective}/` remain
canonical. This script only *creates* stubs and *flags drift*; it never
modifies existing YAML files.

Usage:
    uv run python Roadmapping/History/Bibliography/sync_from_zotero.py
    uv run python Roadmapping/History/Bibliography/sync_from_zotero.py --dry-run
    uv run python Roadmapping/History/Bibliography/sync_from_zotero.py --bib path/to/zotero.bib
"""

import argparse
import re
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
BIB_DIR = Path(__file__).resolve().parent
DEFAULT_BIB = BIB_DIR / "zotero.bib"

FRONTMATTER_RE = re.compile(r"^---\n(.*?)\n---", re.DOTALL)

# Fields whose drift between Zotero and YAML is worth flagging (the
# canonical metadata that shouldn't change after first capture).
DRIFT_FIELDS = ("title", "year", "doi")


# ───────────────────────── BibTeX parsing ─────────────────────────
# BBT's output is well-structured: one entry per @-block, fields on
# their own lines, values in either {braces} or "quotes". We don't need
# a full BibTeX parser to handle it — a small state machine suffices.

def parse_bibtex(text: str) -> list[dict]:
    """Parse a BBT-format BibTeX file. Returns a list of dicts."""
    entries: list[dict] = []
    # Strip BibTeX comments (% to EOL) and find each @entry{...} block.
    text = re.sub(r"%[^\n]*", "", text)
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,(.*?)\n\}\n", text, re.DOTALL):
        entry_type, cite_key, body = m.group(1), m.group(2), m.group(3)
        meta: dict = {"_type": entry_type.lower(), "cite_key": cite_key.strip()}
        for field_name, field_value in _split_fields(body):
            meta[field_name.lower()] = _unbrace(field_value)
        entries.append(meta)
    return entries


def _split_fields(body: str) -> list[tuple[str, str]]:
    """Split a BibTeX entry body into (field_name, field_value_raw) pairs.

    Handles nested braces (common in titles with {LaTeX} commands) by
    tracking brace depth rather than greedy regex.
    """
    out: list[tuple[str, str]] = []
    pos = 0
    while pos < len(body):
        # Skip whitespace and commas.
        while pos < len(body) and body[pos] in " \t\n\r,":
            pos += 1
        if pos >= len(body):
            break
        # Field name.
        name_match = re.match(r"([A-Za-z][\w-]*)\s*=\s*", body[pos:])
        if not name_match:
            break
        name = name_match.group(1)
        pos += name_match.end()
        # Field value: either {...} (with nested braces tracked) or "..." or bare number.
        if pos < len(body) and body[pos] == "{":
            depth = 1
            pos += 1
            start = pos
            while pos < len(body) and depth > 0:
                if body[pos] == "{":
                    depth += 1
                elif body[pos] == "}":
                    depth -= 1
                pos += 1
            value = body[start: pos - 1]  # exclude closing brace
        elif pos < len(body) and body[pos] == '"':
            pos += 1
            start = pos
            while pos < len(body) and body[pos] != '"':
                pos += 1
            value = body[start: pos]
            pos += 1  # consume closing quote
        else:
            value_match = re.match(r"([^,\n]+)", body[pos:])
            if not value_match:
                break
            value = value_match.group(1).strip()
            pos += value_match.end()
        out.append((name, value))
    return out


def _unbrace(value: str) -> str:
    """Strip outer balanced braces and normalise BibTeX escapes lightly."""
    value = value.strip()
    # Repeatedly strip the outermost matched braces.
    while len(value) > 1 and value[0] == "{" and value[-1] == "}":
        depth = 0
        for i, ch in enumerate(value):
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0 and i < len(value) - 1:
                    return value
        value = value[1:-1].strip()
    # Collapse internal whitespace runs.
    return re.sub(r"\s+", " ", value)


# ───────────────────────── BibTeX → YAML translation ─────────────────────────

def split_authors(field: str) -> list[str]:
    """Split a BibTeX 'author' field on ' and '. Names returned 'First Last'."""
    parts = [p.strip() for p in re.split(r"\s+and\s+", field) if p.strip()]
    out = []
    for p in parts:
        # 'Last, First' → 'First Last'; else pass through.
        if "," in p:
            last, first = p.split(",", 1)
            out.append(f"{first.strip()} {last.strip()}")
        else:
            out.append(p)
    return out


def parse_keywords(field: str) -> list[str]:
    """BBT exports Zotero tags into the 'keywords' field, comma-separated."""
    if not field:
        return []
    return [k.strip() for k in field.split(",") if k.strip()]


def categorise_tags(zotero_tags: list[str]) -> tuple[str, str, list[str]]:
    """Map Zotero tags into (type, era, yaml_tags).

    Per the install guide convention:
      `repo:retrospective`        → type=retrospective (default primary)
      `era/<...>`                 → era field
      `thread/<...>`              → yaml_tags includes `#thread/<...>`
      any other tag (verbatim)    → yaml_tags includes that tag
    """
    type_ = "primary"
    era = ""
    yaml_tags: list[str] = []
    for t in zotero_tags:
        if t == "repo:retrospective":
            type_ = "retrospective"
        elif t.startswith("era/"):
            era = t[len("era/"):]
            yaml_tags.append(f"#{t}")
        elif t.startswith("thread/"):
            yaml_tags.append(f"#{t}")
        elif t.startswith("repo:"):
            continue  # other repo:* tags reserved for future use
        else:
            yaml_tags.append(t)
    return type_, era, yaml_tags


def bib_to_yaml(entry: dict) -> dict:
    """Translate a parsed BBT entry into the YAML schema."""
    zotero_tags = parse_keywords(entry.get("keywords", ""))
    type_, era, yaml_tags = categorise_tags(zotero_tags)

    authors_field = entry.get("author", "")
    authors = split_authors(authors_field) if authors_field else []

    return {
        "cite_key": entry["cite_key"],
        "title": entry.get("title", ""),
        "authors": authors,
        "year": _parse_year(entry.get("year", "") or entry.get("date", "")),
        "type": type_,
        "era": era,
        "tags": yaml_tags,
        "journal": entry.get("journal", "") or entry.get("booktitle", "") or entry.get("publisher", ""),
        "volume": entry.get("volume", ""),
        "issue": entry.get("number", "") or entry.get("issue", ""),
        "pages": entry.get("pages", ""),
        "doi": entry.get("doi", ""),
        "url": entry.get("url", ""),
        "arxiv_id": entry.get("eprint", "") if entry.get("archiveprefix", "").lower() == "arxiv" else "",
        "pdf_status": "pending",
        "pdf_path": "",
        "converted_md": "",
        "pages_quoted": "",
        "gill_corpus_overlap": [],
        # Crocco et al. (2026) reporting expectation #7: human-read-source flag.
        # Zotero capture alone doesn't confirm a human read the source; defaults False.
        "human_reviewed": False,
        "chapters_citing": [],
    }


def _parse_year(field: str) -> int | str:
    """BBT may export year as '1865' or date as '1865-08-01'. Return int."""
    if not field:
        return ""
    m = re.match(r"(\d{4})", str(field))
    return int(m.group(1)) if m else field


# ───────────────────────── existing-stub indexing ─────────────────────────

def index_existing_stubs() -> dict[str, Path]:
    """Map cite_key → path for all existing YAML bibliography stubs."""
    out: dict[str, Path] = {}
    for sub in ("Primary", "Retrospective"):
        d = BIB_DIR / sub
        if d.is_dir():
            for md in sorted(d.glob("*.md")):
                out[md.stem] = md
    return out


def load_stub(md: Path) -> dict:
    text = md.read_text(encoding="utf-8")
    m = FRONTMATTER_RE.match(text)
    if not m:
        return {}
    return yaml.safe_load(m.group(1)) or {}


# ───────────────────────── stub scaffolding ─────────────────────────

def render_stub(meta: dict) -> str:
    """Render a YAML-frontmatter stub file matching the existing schema."""
    yaml_block = yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, default_flow_style=False)
    body_title = meta.get("title") or meta["cite_key"]
    return (
        f"---\n{yaml_block}---\n\n"
        f"# {body_title}\n\n"
        "*Skeleton bib note imported from Zotero via `sync_from_zotero.py`. "
        "Fill in: 2-3-paragraph summary, `gill_corpus_overlap`, `chapters_citing`, "
        "and any Zotero-missing fields (e.g., free-form era if `era/...` tag wasn't applied).*\n"
    )


def write_stub(meta: dict, dry_run: bool) -> Path:
    subdir = "Primary" if meta["type"] == "primary" else "Retrospective"
    out_path = BIB_DIR / subdir / f"{meta['cite_key']}.md"
    if dry_run:
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render_stub(meta), encoding="utf-8")
    return out_path


# ───────────────────────── drift detection ─────────────────────────

def compare_drift(zotero_meta: dict, yaml_meta: dict) -> list[tuple[str, str, str]]:
    """For each DRIFT_FIELDS, return (field, zotero_value, yaml_value) where they differ.

    Empty-vs-set is not drift in one direction (Zotero adds info → fine), only when
    both are populated with different non-empty values.
    """
    drift = []
    for f in DRIFT_FIELDS:
        z = str(zotero_meta.get(f, "")).strip()
        y = str(yaml_meta.get(f, "")).strip()
        if z and y and z != y:
            drift.append((f, z, y))
    return drift


# ───────────────────────── main ─────────────────────────

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bib", default=str(DEFAULT_BIB),
                   help=f"BBT-exported BibTeX file (default: {DEFAULT_BIB})")
    p.add_argument("--dry-run", action="store_true",
                   help="show what would be scaffolded without writing files")
    p.add_argument("--quiet", action="store_true", help="only print summary + drift warnings")
    args = p.parse_args()

    bib_path = Path(args.bib).resolve()
    if not bib_path.exists():
        print(f"Bib file not found: {bib_path}", file=sys.stderr)
        print(f"\nIf this is your first sync, the file is auto-created by Better BibTeX's", file=sys.stderr)
        print(f"'Automatic export → On change' setting. See:", file=sys.stderr)
        print(f"  Roadmapping/Tooling/install_zotero_obsidian.md §4\n", file=sys.stderr)
        return 1

    text = bib_path.read_text(encoding="utf-8")
    entries = parse_bibtex(text)
    if not entries:
        print(f"No entries parsed from {bib_path}.", file=sys.stderr)
        return 0

    existing = index_existing_stubs()
    scaffolded = 0
    skipped_existing = 0
    drift_count = 0

    for entry in entries:
        cite_key = entry["cite_key"]
        zotero_meta = bib_to_yaml(entry)

        if cite_key in existing:
            # Drift check.
            yaml_meta = load_stub(existing[cite_key])
            drift = compare_drift(zotero_meta, yaml_meta)
            if drift:
                drift_count += 1
                print(f"⚠  drift: {cite_key} ({existing[cite_key].relative_to(REPO_ROOT)})", file=sys.stderr)
                for field, z, y in drift:
                    print(f"     {field}: zotero={z!r}  yaml={y!r}", file=sys.stderr)
            else:
                skipped_existing += 1
                if not args.quiet:
                    print(f"·  exists, no drift: {cite_key}")
            continue

        out_path = write_stub(zotero_meta, dry_run=args.dry_run)
        verb = "would scaffold" if args.dry_run else "scaffolded"
        rel = out_path.relative_to(REPO_ROOT)
        print(f"{'+' if not args.dry_run else '?'}  {verb}: {rel}")
        scaffolded += 1

    print(
        f"\nSummary: {len(entries)} entries in {bib_path.name}; "
        f"{scaffolded} scaffolded; "
        f"{skipped_existing} already present (no drift); "
        f"{drift_count} flagged for manual reconciliation.",
        file=sys.stderr,
    )
    if drift_count:
        print(
            "\nDrift means: the entry exists in both Zotero and the YAML stub, with different\n"
            "values in one of the canonical fields (title / year / DOI). Resolve manually;\n"
            "this script never overwrites existing YAML stubs.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
