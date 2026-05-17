"""Download a single primary/retrospective PDF and link it to its bibliography note.

Given a `--cite-key` and a `--url` (or `--doi` for a Crossref unpaywall lookup),
fetches the PDF into `Historical_Papers/<Primary|Retrospective>/<cite_key>.pdf`,
updates the bibliography note's `pdf_status` + `pdf_path` fields, and re-runs
`update_acquisition_tracker.py` so the tracker reflects the new state.

The script does NOT decide whether a PDF is in copyright; that's the operator's
call via `--status out_of_copyright_public` (committed via `git add -f`) vs.
`--status acquired` (kept local only). Default is `acquired`.

Usage:
    uv run python Roadmapping/History/_tools/fetch_pdf.py \\
        --cite-key maxwell1865_dynamical_theory \\
        --url https://royalsocietypublishing.org/doi/pdf/10.1098/rstl.1865.0008 \\
        --status out_of_copyright_public

    uv run python Roadmapping/History/_tools/fetch_pdf.py \\
        --cite-key pais1982_subtle_is_the_lord \\
        --url file:///tmp/pais1982.pdf \\
        --status acquired

    uv run python Roadmapping/History/_tools/fetch_pdf.py \\
        --cite-key foo1900_bar --doi 10.xxxx/yyyy --status acquired --dry-run
"""

import argparse
import re
import shutil
import sys
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
BIB_ROOT = REPO_ROOT / "Roadmapping" / "History" / "Bibliography"
PAPERS_ROOT = REPO_ROOT / "Roadmapping" / "Historical_Papers"
TRACKER_TOOL = BIB_ROOT / "update_acquisition_tracker.py"

USER_AGENT = "PyPhysics-fetch-pdf/0.1 (mailto:morris.trey.j@gmail.com)"
CROSSREF_URL = "https://api.crossref.org/works/{doi}"
UNPAYWALL_URL = "https://api.unpaywall.org/v2/{doi}?email=morris.trey.j@gmail.com"
FRONTMATTER_RE = re.compile(r"^(---\n)(.*?)(\n---\n)(.*)$", re.DOTALL)

VALID_STATUSES = ("out_of_copyright_public", "acquired", "pending", "unavailable")


def find_bib_note(cite_key: str) -> Path:
    for sub in ("Primary", "Retrospective"):
        candidate = BIB_ROOT / sub / f"{cite_key}.md"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"No bibliography note for {cite_key!r} under "
        f"{BIB_ROOT}/{{Primary,Retrospective}}. Run scaffold_bib_note.py first."
    )


def subdir_for(note: Path) -> str:
    return "Primary" if note.parent.name == "Primary" else "Retrospective"


def resolve_url(args: argparse.Namespace) -> str:
    if args.url:
        return args.url
    if args.doi:
        # Try unpaywall for an open-access PDF link.
        req = Request(UNPAYWALL_URL.format(doi=args.doi), headers={"User-Agent": USER_AGENT})
        try:
            import json
            with urlopen(req, timeout=15) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            raise RuntimeError(f"Unpaywall lookup failed for {args.doi}: {e}")
        best = payload.get("best_oa_location") or {}
        url = best.get("url_for_pdf") or best.get("url")
        if not url:
            raise RuntimeError(f"No open-access PDF found via Unpaywall for {args.doi}")
        return url
    raise SystemExit("--url or --doi required")


def download(url: str, dest: Path) -> int:
    parsed = urlparse(url)
    if parsed.scheme == "file":
        src = Path(parsed.path)
        if not src.exists():
            raise FileNotFoundError(src)
        shutil.copy(src, dest)
        return dest.stat().st_size
    req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/pdf,*/*"})
    with urlopen(req, timeout=60) as resp, open(dest, "wb") as out:
        shutil.copyfileobj(resp, out)
    return dest.stat().st_size


def update_frontmatter(note: Path, pdf_path_rel: str, status: str, dry_run: bool) -> None:
    text = note.read_text(encoding="utf-8")
    m = FRONTMATTER_RE.match(text)
    if not m:
        raise RuntimeError(f"{note} has no YAML frontmatter")
    meta = yaml.safe_load(m.group(2)) or {}
    meta["pdf_status"] = status
    meta["pdf_path"] = pdf_path_rel
    new_yaml = yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, default_flow_style=False)
    new_text = f"{m.group(1)}{new_yaml.rstrip()}{m.group(3)}{m.group(4)}"
    if dry_run:
        print(f"--- (dry-run) would update {note}")
        print(new_text[: new_text.find('---\n', 4) + 4])
    else:
        note.write_text(new_text, encoding="utf-8")
        print(f"Updated frontmatter: {note}")


def refresh_tracker(dry_run: bool) -> None:
    import subprocess
    cmd = ["uv", "run", "python", str(TRACKER_TOOL)]
    if dry_run:
        cmd.append("--dry-run")
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cite-key", required=True)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--url", help="direct PDF URL (or file://… for a local copy)")
    src.add_argument("--doi", help="DOI; tries Unpaywall to find an open-access PDF URL")
    p.add_argument("--status", choices=VALID_STATUSES, default="acquired",
                   help="pdf_status to write into the bib note (default: acquired)")
    p.add_argument("--dry-run", action="store_true",
                   help="print actions without downloading or editing")
    p.add_argument("--force", action="store_true",
                   help="overwrite existing PDF if present")
    args = p.parse_args()

    note = find_bib_note(args.cite_key)
    subdir = subdir_for(note)
    dest_dir = PAPERS_ROOT / subdir
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_pdf = dest_dir / f"{args.cite_key}.pdf"

    if dest_pdf.exists() and not args.force and not args.dry_run:
        print(f"{dest_pdf} already exists; pass --force to overwrite.", file=sys.stderr)
        return 1

    url = resolve_url(args)
    print(f"Source: {url}")
    print(f"Destination: {dest_pdf}")

    if args.dry_run:
        print("(dry-run) skipping actual download")
    else:
        size = download(url, dest_pdf)
        print(f"Downloaded {size:,} bytes")

    # Path stored in YAML is repo-root-relative for portability.
    pdf_path_rel = f"../../{dest_pdf.relative_to(REPO_ROOT).as_posix()}"
    update_frontmatter(note, pdf_path_rel, args.status, args.dry_run)
    refresh_tracker(args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
