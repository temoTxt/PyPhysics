#!/usr/bin/env bash
# build_report.sh — markdown → LaTeX → PDF build for Roadmapping/Author_Reports/
#
# Usage:
#   ./build_report.sh <basename>             # full build; writes .tex + .pdf to this folder
#   ./build_report.sh <basename> --dry-run   # builds in a temp dir; verifies mechanics; no commit
#
# Example:
#   ./build_report.sh 2026-05_interim_for_gill
#
# Pipeline:
#   1. Strip <!-- ... --> HTML comments (per-paragraph TODO markers) from a working copy.
#   2. pandoc markdown → LaTeX with a minimal preamble (article class, hyperref, amsmath, booktabs, geometry).
#   3. pdflatex twice (for cross-references).
#   4. Defensive check: fail if any "<!-- TODO" substring survives in the .tex.
#   5. Defensive check: fail if PDF page count is outside [3, 7] per the length-budget discipline.
#
# Reproducibility scope: rebuildable from the same .md on a clean checkout, but not byte-identical-
# deterministic across toolchain versions, font availability, or /tmp path embedding. PDF metadata
# date is pinned via pandoc --metadata date=... to prevent timestamp drift.
#
# Dependencies: pandoc (>= 3.x), texlive-latex-base, texlive-latex-extra, texlive-pictures, poppler-utils.
# Unicode handling: emoji (U+2705 ✅, U+26A0 ⚠, U+1F534 🔴) are substituted to ASCII tags in a
# preprocess pass because pdflatex cannot render them; everything else (≈, ×, −, ½, →, etc.) passes
# through pandoc's utf8 handling unchanged.

set -euo pipefail

# ----------------------------------------------------------------------
# Configuration

PAGE_MIN="${PAGE_MIN:-3}"            # override with PAGE_MIN=N
PAGE_MAX="${PAGE_MAX:-7}"            # override with PAGE_MAX=N (long verification reports)
PDFENGINE="${PDFENGINE:-pdflatex}"   # override with PDFENGINE=lualatex if you have a complete
                                     # lualatex install (needs lualatex-math.sty + texlive
                                     # unicode-math + matching font OTFs). Default pdflatex
                                     # is reliable but requires the per-character Unicode
                                     # substitutions below for Greek + math symbols.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PINNED_DATE="${REPORT_DATE:-$(date -u +%Y-%m-%d)}"   # override with REPORT_DATE=YYYY-MM-DD

# ----------------------------------------------------------------------
# Args

if [ $# -lt 1 ]; then
  echo "usage: $0 <basename> [--dry-run]" >&2
  echo "       <basename> is the filename without extension, e.g. 2026-05_interim_for_gill" >&2
  exit 2
fi

BASENAME="$1"
DRY_RUN=0
if [ "${2:-}" = "--dry-run" ]; then
  DRY_RUN=1
fi

SRC_MD="$SCRIPT_DIR/${BASENAME}.md"

if [ ! -f "$SRC_MD" ]; then
  echo "ERROR: source markdown not found at $SRC_MD" >&2
  exit 1
fi

# ----------------------------------------------------------------------
# Working directory (temp for dry-run, script dir for real build)

if [ "$DRY_RUN" -eq 1 ]; then
  WORK_DIR="$(mktemp -d -t author_report_build.XXXXXX)"
  trap 'rm -rf "$WORK_DIR"' EXIT
  echo "[dry-run] working in $WORK_DIR (will be cleaned up on exit)"
else
  WORK_DIR="$SCRIPT_DIR"
fi

OUT_MD="$WORK_DIR/${BASENAME}.stripped.md"
OUT_TEX="$WORK_DIR/${BASENAME}.tex"
OUT_PDF="$WORK_DIR/${BASENAME}.pdf"

# ----------------------------------------------------------------------
# Step 1: strip HTML comments

echo "[1/5] preprocessing: strip HTML comments + GitHub math fences + Unicode subs..."
# Strip <!-- ... --> blocks; convert GitHub-style ```math ... ``` fences into
# pandoc-native $$...$$ display-math; substitute emoji (no TeX engine renders
# them) + primes (arcsec marks).  Greek + math-symbol substitutions are only
# applied when PDFENGINE is pdflatex; lualatex/xelatex render Unicode natively
# and substituting would produce ugly ASCII names in prose.
perl -CSD -0777 -pe '
  s/<!--.*?-->//gs;
  # GitHub-style ```math ... ``` -> pandoc-style $$...$$ display math
  s/```math\n(.*?)\n```\n/\$\$\n$1\n\$\$\n/sg;
  # Emoji verdicts -> ASCII tags (no TeX engine renders emoji by default)
  s/\x{2705}/[OK]/g;        # U+2705 white heavy check mark
  s/\x{26A0}/[warn]/g;      # U+26A0 warning sign
  s/\x{1F534}/[X]/g;        # U+1F534 large red circle
  s/\x{1F7E1}/[partial]/g;  # U+1F7E1 large yellow circle
  s/\x{274C}/[X]/g;         # U+274C cross mark
  # Primes (arcsec / arcmin marks); \x27 = ASCII apostrophe.
  s/\x{2032}/\x27/g;        # U+2032 PRIME (arcmin)
  s/\x{2033}/\x27\x27/g;    # U+2033 DOUBLE PRIME (arcsec)
  s/\x{2034}/\x27\x27\x27/g;# U+2034 TRIPLE PRIME
' "$SRC_MD" > "$OUT_MD"

# Engine-specific extra substitutions: pdflatex needs Greek + math-symbol -> ASCII
# because its T1 encoding can't render Unicode in prose.  lualatex/xelatex skip
# this and pass Unicode through to fontspec/unicode-math.
if [ "$PDFENGINE" = "pdflatex" ]; then
  perl -CSD -0777 -i -pe '
    # Logical implication / arrows
    s/\x{21D2}/=>/g;          # U+21D2 RIGHTWARDS DOUBLE ARROW (⇒)
    s/\x{21D0}/<=/g;          # U+21D0 LEFTWARDS DOUBLE ARROW (⇐)
    # Superscripts -> ASCII ^N
    s/\x{2070}/^0/g; s/\x{00B9}/^1/g; s/\x{00B2}/^2/g; s/\x{00B3}/^3/g;
    s/\x{2074}/^4/g; s/\x{2075}/^5/g; s/\x{2076}/^6/g; s/\x{2077}/^7/g;
    s/\x{2078}/^8/g; s/\x{2079}/^9/g; s/\x{207B}/^-/g;
    # Subscripts -> ASCII _N
    s/\x{2080}/_0/g; s/\x{2081}/_1/g; s/\x{2082}/_2/g; s/\x{2083}/_3/g;
    s/\x{2084}/_4/g; s/\x{2085}/_5/g; s/\x{2086}/_6/g; s/\x{2087}/_7/g;
    s/\x{2088}/_8/g; s/\x{2089}/_9/g;
    # Math symbols
    s/\x{2248}/~=/g; s/\x{00D7}/x/g; s/\x{2212}/-/g; s/\x{2026}/.../g;
    s/\x{2192}/->/g; s/\x{2190}/<-/g; s/\x{00B7}/./g; s/\x{2022}/*/g;
    s/\x{00B1}/+-/g; s/\x{2264}/<=/g; s/\x{2265}/>=/g; s/\x{2260}/!=/g;
    s/\x{221A}/sqrt/g; s/\x{221E}/inf/g; s/\x{2202}/d/g; s/\x{222B}/int/g;
    s/\x{2207}/grad/g; s/\x{210F}/hbar/g; s/\x{2113}/l/g; s/\x{2211}/sum/g;
    s/\x{220F}/prod/g; s/\x{27E8}/</g; s/\x{27E9}/>/g;
    s/\x{00BD}/1\/2/g; s/\x{00BC}/1\/4/g; s/\x{00BE}/3\/4/g;
    # Greek letters
    s/\x{03B1}/alpha/g; s/\x{03B2}/beta/g; s/\x{03B3}/gamma/g;
    s/\x{03B4}/delta/g; s/\x{03B5}/epsilon/g; s/\x{03B6}/zeta/g;
    s/\x{03B7}/eta/g; s/\x{03B8}/theta/g; s/\x{03B9}/iota/g;
    s/\x{03BA}/kappa/g; s/\x{03BB}/lambda/g; s/\x{03BC}/mu/g;
    s/\x{03BD}/nu/g; s/\x{03BE}/xi/g; s/\x{03C0}/pi/g;
    s/\x{03C1}/rho/g; s/\x{03C3}/sigma/g; s/\x{03C4}/tau/g;
    s/\x{03C5}/upsilon/g; s/\x{03C6}/phi/g; s/\x{03C7}/chi/g;
    s/\x{03C8}/psi/g; s/\x{03C9}/omega/g;
    s/\x{0394}/Delta/g; s/\x{03A3}/Sigma/g; s/\x{03A0}/Pi/g;
  ' "$OUT_MD"
fi

# ----------------------------------------------------------------------
# Step 2: pandoc → LaTeX

echo "[2/5] pandoc markdown -> LaTeX (target: $PDFENGINE)..."
# Font args: lualatex/xelatex need explicit \setmainfont / \setmonofont via fontspec;
# pdflatex uses lmodern from its template and ignores these vars.
FONT_ARGS=()
case "$PDFENGINE" in
  lualatex|xelatex)
    FONT_ARGS=(
      --variable=mainfont:"DejaVu Serif"
      --variable=monofont:"DejaVu Sans Mono"
    )
    ;;
esac

pandoc "$OUT_MD" \
  --from markdown+raw_tex+pipe_tables+grid_tables+yaml_metadata_block \
  --to latex \
  --standalone \
  --pdf-engine="$PDFENGINE" \
  --variable=geometry:margin=1in \
  --variable=fontsize:11pt \
  --variable=documentclass:article \
  --variable=colorlinks:true \
  --variable=linkcolor:blue \
  --variable=urlcolor:blue \
  "${FONT_ARGS[@]}" \
  --metadata author="Trey Morris with Claude Opus 4.7" \
  --metadata date="$PINNED_DATE" \
  --output "$OUT_TEX"
# title and subject are read from the .md's YAML frontmatter (yaml_metadata_block).
# author and date are pinned on the command line so PDF metadata is uniform across
# reports and date is stable across rebuilds (override with REPORT_DATE=YYYY-MM-DD).

# ----------------------------------------------------------------------
# Step 3: defensive check — no <!-- TODO survived

echo "[3/5] defensive check: no <!-- TODO markers in .tex..."
if grep -q '<!-- TODO' "$OUT_TEX"; then
  echo "ERROR: '<!-- TODO' substring found in $OUT_TEX after stripping." >&2
  echo "       The HTML-comment stripper failed to remove all TODO markers." >&2
  grep -n '<!-- TODO' "$OUT_TEX" | head -5 >&2
  exit 1
fi

# Also catch any literal "TODO" labels that might have leaked another way
if grep -qE 'human reviews and fills in' "$OUT_TEX"; then
  echo "ERROR: 'human reviews and fills in' text found in $OUT_TEX after stripping." >&2
  echo "       A TODO marker leaked through; check for non-HTML-comment forms." >&2
  exit 1
fi

# ----------------------------------------------------------------------
# Step 4: $PDFENGINE twice (for cross-references)

echo "[4/5] $PDFENGINE (pass 1)..."
(cd "$WORK_DIR" && "$PDFENGINE" -interaction=nonstopmode -halt-on-error "${BASENAME}.tex" > "${BASENAME}.log1" 2>&1) || {
  echo "ERROR: $PDFENGINE pass 1 failed. Tail of log:" >&2
  tail -30 "$WORK_DIR/${BASENAME}.log1" >&2
  exit 1
}

echo "[4/5] $PDFENGINE (pass 2)..."
(cd "$WORK_DIR" && "$PDFENGINE" -interaction=nonstopmode -halt-on-error "${BASENAME}.tex" > "${BASENAME}.log2" 2>&1) || {
  echo "ERROR: $PDFENGINE pass 2 failed. Tail of log:" >&2
  tail -30 "$WORK_DIR/${BASENAME}.log2" >&2
  exit 1
}

# ----------------------------------------------------------------------
# Step 5: page-count defensive check

echo "[5/5] page-count defensive check [$PAGE_MIN, $PAGE_MAX]..."
if [ ! -f "$OUT_PDF" ]; then
  echo "ERROR: PDF not produced at $OUT_PDF" >&2
  exit 1
fi

PAGES=$(pdfinfo "$OUT_PDF" | awk '/^Pages:/ {print $2}')
if [ -z "$PAGES" ]; then
  echo "ERROR: could not read page count from $OUT_PDF" >&2
  exit 1
fi

if [ "$PAGES" -lt "$PAGE_MIN" ] || [ "$PAGES" -gt "$PAGE_MAX" ]; then
  echo "ERROR: PDF page count $PAGES is outside the [$PAGE_MIN, $PAGE_MAX] budget." >&2
  echo "       See plan §3 for the length-budget discipline. Inspect: $OUT_PDF" >&2
  exit 1
fi

# ----------------------------------------------------------------------
# Done

if [ "$DRY_RUN" -eq 1 ]; then
  echo
  echo "[dry-run SUCCESS] pipeline mechanics verified."
  echo "  source:    $SRC_MD"
  echo "  pages:     $PAGES (within [$PAGE_MIN, $PAGE_MAX])"
  echo "  PDF (temp): $OUT_PDF (cleaned up on exit)"
  echo "  not written to: $SCRIPT_DIR/${BASENAME}.tex / .pdf"
else
  # Clean up pdflatex artifacts that aren't the deliverable
  rm -f "$WORK_DIR/${BASENAME}.aux" \
        "$WORK_DIR/${BASENAME}.log" \
        "$WORK_DIR/${BASENAME}.log1" \
        "$WORK_DIR/${BASENAME}.log2" \
        "$WORK_DIR/${BASENAME}.out" \
        "$WORK_DIR/${BASENAME}.stripped.md"
  echo
  echo "[BUILD SUCCESS]"
  echo "  source:  $SRC_MD"
  echo "  LaTeX:   $OUT_TEX"
  echo "  PDF:     $OUT_PDF"
  echo "  pages:   $PAGES"
fi
