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
#   5. Report PDF page count (no enforcement).
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

echo "[1/5] stripping HTML comments + substituting non-pdflatex unicode..."
# Strip <!-- ... --> blocks; substitute emoji + Unicode superscripts/subscripts that pdflatex
# can't render in its default T1 encoding. Math mode unicode (≈, ×, −, →, ⋯) is handled by
# pandoc's utf8 conversion to LaTeX equivalents.
perl -CSD -0777 -pe '
  s/<!--.*?-->//gs;
  # Emoji verdicts -> ASCII tags
  s/\x{2705}/[OK]/g;        # U+2705 white heavy check mark
  s/\x{26A0}/[warn]/g;      # U+26A0 warning sign
  s/\x{1F534}/[X]/g;        # U+1F534 large red circle
  s/\x{274C}/[X]/g;         # U+274C cross mark
  # Superscripts -> ASCII ^N (renders as literal ^N in \texttt{} contexts, OK)
  s/\x{2070}/^0/g; s/\x{00B9}/^1/g; s/\x{00B2}/^2/g; s/\x{00B3}/^3/g;
  s/\x{2074}/^4/g; s/\x{2075}/^5/g; s/\x{2076}/^6/g; s/\x{2077}/^7/g;
  s/\x{2078}/^8/g; s/\x{2079}/^9/g; s/\x{207B}/^-/g;
  # Subscripts -> ASCII _N
  s/\x{2080}/_0/g; s/\x{2081}/_1/g; s/\x{2082}/_2/g; s/\x{2083}/_3/g;
  s/\x{2084}/_4/g; s/\x{2085}/_5/g; s/\x{2086}/_6/g; s/\x{2087}/_7/g;
  s/\x{2088}/_8/g; s/\x{2089}/_9/g;
  # Math symbols that may appear in verbatim/code spans -> ASCII fallbacks
  # (in math mode, pandoc already converts these to LaTeX commands; this only
  # affects \texttt{} contexts where literal Unicode would break pdflatex.)
  s/\x{2248}/~=/g;          # U+2248 ALMOST EQUAL TO (≈)
  s/\x{00D7}/x/g;           # U+00D7 MULTIPLICATION SIGN (×)
  s/\x{2212}/-/g;           # U+2212 MINUS SIGN (−)
  s/\x{2026}/.../g;         # U+2026 HORIZONTAL ELLIPSIS (…)
  s/\x{2192}/->/g;          # U+2192 RIGHTWARDS ARROW (→)
  s/\x{2190}/<-/g;          # U+2190 LEFTWARDS ARROW (←)
  s/\x{00B7}/./g;           # U+00B7 MIDDLE DOT (·)
  s/\x{2022}/*/g;           # U+2022 BULLET (•)
  s/\x{00B1}/+-/g;          # U+00B1 PLUS-MINUS SIGN (±)
  s/\x{2264}/<=/g;          # U+2264 LESS-THAN OR EQUAL TO (≤)
  s/\x{2265}/>=/g;          # U+2265 GREATER-THAN OR EQUAL TO (≥)
  s/\x{2260}/!=/g;          # U+2260 NOT EQUAL TO (≠)
  s/\x{221A}/sqrt/g;        # U+221A SQUARE ROOT (√)
  s/\x{221E}/inf/g;         # U+221E INFINITY (∞)
  s/\x{2202}/d/g;           # U+2202 PARTIAL DIFFERENTIAL (∂)
  s/\x{222B}/int/g;         # U+222B INTEGRAL (∫)
  s/\x{2207}/grad/g;        # U+2207 NABLA (∇)
  s/\x{210F}/hbar/g;        # U+210F PLANCK CONSTANT OVER TWO PI (ℏ)
  s/\x{2113}/l/g;           # U+2113 SCRIPT SMALL L (ℓ)
  s/\x{2211}/sum/g;         # U+2211 N-ARY SUMMATION (∑)
  s/\x{220F}/prod/g;        # U+220F N-ARY PRODUCT (∏)
  s/\x{27E8}/</g;           # U+27E8 MATHEMATICAL LEFT ANGLE BRACKET (⟨)
  s/\x{27E9}/>/g;           # U+27E9 MATHEMATICAL RIGHT ANGLE BRACKET (⟩)
  s/\x{00BD}/1\/2/g;        # U+00BD VULGAR FRACTION ONE HALF (½)
  s/\x{00BC}/1\/4/g;        # U+00BC VULGAR FRACTION ONE QUARTER (¼)
  s/\x{00BE}/3\/4/g;        # U+00BE VULGAR FRACTION THREE QUARTERS (¾)
  # Greek letters in verbatim contexts -> LaTeX-source notation
  s/\x{03B1}/alpha/g; s/\x{03B2}/beta/g; s/\x{03B3}/gamma/g;
  s/\x{03B4}/delta/g; s/\x{03B5}/epsilon/g; s/\x{03B6}/zeta/g;
  s/\x{03B7}/eta/g; s/\x{03B8}/theta/g; s/\x{03B9}/iota/g;
  s/\x{03BA}/kappa/g; s/\x{03BB}/lambda/g; s/\x{03BC}/mu/g;
  s/\x{03BD}/nu/g; s/\x{03BE}/xi/g; s/\x{03C0}/pi/g;
  s/\x{03C1}/rho/g; s/\x{03C3}/sigma/g; s/\x{03C4}/tau/g;
  s/\x{03C5}/upsilon/g; s/\x{03C6}/phi/g; s/\x{03C7}/chi/g;
  s/\x{03C8}/psi/g; s/\x{03C9}/omega/g;
  s/\x{0394}/Delta/g; s/\x{03A3}/Sigma/g; s/\x{03A0}/Pi/g;
' "$SRC_MD" > "$OUT_MD"

# ----------------------------------------------------------------------
# Step 2: pandoc → LaTeX

echo "[2/5] pandoc markdown -> LaTeX (pdflatex-target, Unicode-native)..."
pandoc "$OUT_MD" \
  --from markdown+raw_tex+pipe_tables+grid_tables+yaml_metadata_block \
  --to latex \
  --standalone \
  --pdf-engine=pdflatex \
  --variable=geometry:left=1.25in \
  --variable=geometry:right=1.25in \
  --variable=geometry:top=1in \
  --variable=geometry:bottom=1in \
  --variable=fontsize:11pt \
  --variable=documentclass:article \
  --variable=colorlinks:true \
  --variable=linkcolor:blue \
  --variable=urlcolor:blue \
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
# Step 4: pdflatex twice

echo "[4/5] pdflatex (pass 1)..."
(cd "$WORK_DIR" && pdflatex -interaction=nonstopmode -halt-on-error "${BASENAME}.tex" > "${BASENAME}.log1" 2>&1) || {
  echo "ERROR: pdflatex pass 1 failed. Tail of log:" >&2
  tail -30 "$WORK_DIR/${BASENAME}.log1" >&2
  exit 1
}

echo "[4/5] pdflatex (pass 2)..."
(cd "$WORK_DIR" && pdflatex -interaction=nonstopmode -halt-on-error "${BASENAME}.tex" > "${BASENAME}.log2" 2>&1) || {
  echo "ERROR: pdflatex pass 2 failed. Tail of log:" >&2
  tail -30 "$WORK_DIR/${BASENAME}.log2" >&2
  exit 1
}

# ----------------------------------------------------------------------
# Step 5: page-count defensive check

echo "[5/5] page count..."
if [ ! -f "$OUT_PDF" ]; then
  echo "ERROR: PDF not produced at $OUT_PDF" >&2
  exit 1
fi

PAGES=$(pdfinfo "$OUT_PDF" | awk '/^Pages:/ {print $2}')
if [ -z "$PAGES" ]; then
  echo "ERROR: could not read page count from $OUT_PDF" >&2
  exit 1
fi

# ----------------------------------------------------------------------
# Done

if [ "$DRY_RUN" -eq 1 ]; then
  echo
  echo "[dry-run SUCCESS] pipeline mechanics verified."
  echo "  source:    $SRC_MD"
  echo "  pages:     $PAGES"
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
