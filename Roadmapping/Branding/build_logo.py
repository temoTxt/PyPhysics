"""Generate the PyPhysics project logo.

Logo concept: the **4-velocity hyperboloid** $b^2 - u^2 = c^2$ -- the algebraic
seed of Gill's dual relativity framework. Visually, the right branch of a
hyperbola sweeping up from the $u$-axis to the $\\tau$-axis, with a single
highlighted point at $u = 0$ where $b = c$ (the "rest" case).

Outputs (committed under `Roadmapping/Branding/`):
- `logo.svg` -- symbol-only, vector
- `logo_horizontal.svg` -- symbol + "PyPhysics" wordmark, vector
- `logo_dark.svg`, `logo_light.svg` -- light/dark variants
- `logo_32.png`, `logo_64.png`, `logo_256.png` -- favicon ladder
- `social_preview_1280x640.png` -- GitHub social preview
- `podcast_cover_3000.png` -- Apple-Podcasts-compatible cover art

PNG generation uses Pillow (already a marker-pdf transitive dependency).
SVG generation is hand-written XML for total control over the geometry.

Usage:
    uv run python Roadmapping/Branding/build_logo.py
    uv run python Roadmapping/Branding/build_logo.py --dry-run
"""

import argparse
import math
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

OUT_DIR = Path(__file__).resolve().parent
YELLOW = "#FFD000"
DARK_BG = "#0B1116"
LIGHT_BG = "#FFFFFF"
DARK_FG = "#0B1116"
LIGHT_FG = "#FFFFFF"


# ───────────────────── geometry ─────────────────────

def hyperbola_points(c: float = 1.0, u_max: float = 2.0, n: int = 200) -> list[tuple[float, float]]:
    """Right branch of b^2 - u^2 = c^2, i.e. b = sqrt(c^2 + u^2), for u in [0, u_max]."""
    pts = []
    for i in range(n + 1):
        u = u_max * i / n
        b = math.sqrt(c * c + u * u)
        pts.append((u, b))
    return pts


# ───────────────────── SVG generation ─────────────────────

def svg_symbol(size: int = 256, bg: str = DARK_BG, fg: str = LIGHT_FG, accent: str = YELLOW,
               include_axes: bool = True, include_labels: bool = True) -> str:
    """Symbol-only SVG. Square viewBox; centered hyperbola."""
    pad = size * 0.08
    # Map (u, b) world coords to SVG pixel coords. World: u in [-2.2, 2.2], b in [-0.6, 2.6].
    # We center horizontally and shift so the hyperbola vertex sits centrally.
    cx = size / 2
    cy = size * 0.62  # vertex slightly below middle for balance
    scale = (size - 2 * pad) / 4.4  # 4.4 = world width

    def to_px(u: float, b: float) -> tuple[float, float]:
        return cx + u * scale, cy - b * scale  # SVG y is inverted

    # Right and left branches.
    pts_right = hyperbola_points(c=1.0, u_max=1.9, n=120)
    pts_left = [(-u, b) for (u, b) in pts_right]

    def path_d(pts: list[tuple[float, float]]) -> str:
        coords = [f"{x:.2f},{y:.2f}" for (u, b) in pts for (x, y) in [to_px(u, b)]]
        return "M " + " L ".join(coords)

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {size} {size}" '
        f'width="{size}" height="{size}" role="img" aria-label="PyPhysics dual-relativity logo">',
        f'  <rect width="100%" height="100%" fill="{bg}"/>',
    ]

    if include_axes:
        # τ axis (vertical) and u axis (horizontal), thin.
        ax_x_y0, _ = to_px(-2.1, 0)
        ax_x_y1, _ = to_px(2.1, 0)
        _, ax_y_x0 = to_px(0, -0.4)
        _, ax_y_x1 = to_px(0, 2.5)
        parts.append(
            f'  <line x1="{ax_x_y0:.2f}" y1="{cy:.2f}" x2="{ax_x_y1:.2f}" y2="{cy:.2f}" '
            f'stroke="{fg}" stroke-width="{max(1.0, size/256):.2f}" opacity="0.45"/>'
        )
        parts.append(
            f'  <line x1="{cx:.2f}" y1="{ax_y_x0:.2f}" x2="{cx:.2f}" y2="{ax_y_x1:.2f}" '
            f'stroke="{fg}" stroke-width="{max(1.0, size/256):.2f}" opacity="0.45"/>'
        )

    # Light cones (asymptotes b = ±u): two dashed diagonals
    asy_w = 2.1
    a1x, a1y = to_px(asy_w, asy_w)
    a2x, a2y = to_px(-asy_w, asy_w)
    parts.append(
        f'  <line x1="{cx:.2f}" y1="{cy:.2f}" x2="{a1x:.2f}" y2="{a1y:.2f}" '
        f'stroke="{fg}" stroke-width="{max(1.0, size/256):.2f}" '
        f'stroke-dasharray="{max(2, size/64):.0f},{max(3, size/64*1.5):.0f}" opacity="0.35"/>'
    )
    parts.append(
        f'  <line x1="{cx:.2f}" y1="{cy:.2f}" x2="{a2x:.2f}" y2="{a2y:.2f}" '
        f'stroke="{fg}" stroke-width="{max(1.0, size/256):.2f}" '
        f'stroke-dasharray="{max(2, size/64):.0f},{max(3, size/64*1.5):.0f}" opacity="0.35"/>'
    )

    # Both hyperbola branches, yellow accent.
    sw = max(2.0, size / 96)
    parts.append(
        f'  <path d="{path_d(pts_right)}" fill="none" stroke="{accent}" '
        f'stroke-width="{sw:.2f}" stroke-linecap="round" stroke-linejoin="round"/>'
    )
    parts.append(
        f'  <path d="{path_d(pts_left)}" fill="none" stroke="{accent}" '
        f'stroke-width="{sw:.2f}" stroke-linecap="round" stroke-linejoin="round"/>'
    )

    # Vertex dot at (u=0, b=c).
    vx, vy = to_px(0, 1.0)
    parts.append(f'  <circle cx="{vx:.2f}" cy="{vy:.2f}" r="{max(3.0, size/56):.2f}" fill="{accent}"/>')

    if include_labels and size >= 96:
        # τ at top of vertical axis; u to the right of horizontal axis.
        font_size = max(11, int(size * 0.075))
        # τ label near top of vertical axis (inside the upper region)
        tau_x, tau_y = to_px(0.18, 2.35)
        parts.append(
            f'  <text x="{tau_x:.2f}" y="{tau_y:.2f}" fill="{fg}" '
            f'font-family="Georgia, serif" font-style="italic" font-size="{font_size}px">τ</text>'
        )
        u_x, u_y = to_px(1.95, -0.15)
        parts.append(
            f'  <text x="{u_x:.2f}" y="{u_y:.2f}" fill="{fg}" '
            f'font-family="Georgia, serif" font-style="italic" font-size="{font_size}px" '
            f'text-anchor="end">u</text>'
        )

    parts.append("</svg>")
    return "\n".join(parts)


def svg_horizontal(width: int = 1024, height: int = 320,
                   bg: str = DARK_BG, fg: str = LIGHT_FG, accent: str = YELLOW) -> str:
    """Horizontal lockup: symbol on the left + 'PyPhysics' wordmark on the right."""
    sym_size = height
    inner = svg_symbol(size=sym_size, bg=bg, fg=fg, accent=accent,
                       include_axes=True, include_labels=True)
    # Strip the outer <svg> wrapper to embed as a group.
    inner_lines = inner.splitlines()
    inner_body = "\n".join(inner_lines[1:-1])  # drop <svg ...> and </svg>

    text_x = sym_size + height * 0.25
    word_size = int(height * 0.34)
    sub_size = int(height * 0.13)
    word_y = height * 0.55
    sub_y = height * 0.78

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="{width}" height="{height}" role="img" aria-label="PyPhysics">\n'
        f'  <rect width="100%" height="100%" fill="{bg}"/>\n'
        f'  <g>{inner_body}</g>\n'
        f'  <text x="{text_x:.2f}" y="{word_y:.2f}" fill="{fg}" '
        f'font-family="Georgia, serif" font-weight="700" font-size="{word_size}px">'
        f'Py<tspan fill="{accent}">Physics</tspan></text>\n'
        f'  <text x="{text_x:.2f}" y="{sub_y:.2f}" fill="{fg}" font-family="Georgia, serif" '
        f'font-style="italic" font-size="{sub_size}px" opacity="0.65">'
        f'dual relativity &amp; quantum mechanics</text>\n'
        f'</svg>\n'
    )


# ───────────────────── PNG generation (Pillow) ─────────────────────

def png_symbol(size: int, bg: str, fg: str, accent: str,
               include_axes: bool = True, include_labels: bool = True,
               include_wordmark: bool = False, canvas_size: tuple[int, int] | None = None) -> Image.Image:
    """Render the same geometry as svg_symbol() using Pillow."""
    if canvas_size is not None:
        canvas_w, canvas_h = canvas_size
    else:
        canvas_w = size * 4 if include_wordmark else size
        canvas_h = size
    img = Image.new("RGB", (canvas_w, canvas_h), color=bg)
    draw = ImageDraw.Draw(img)

    pad = size * 0.08
    # Symbol centered horizontally within its own `size` x size box (which is
    # left-anchored on the canvas). Vertical center aligned to canvas midline.
    cx = size / 2
    cy = canvas_h / 2 + size * 0.12  # slight offset so vertex sits at visual center
    scale = (size - 2 * pad) / 4.4

    def to_px(u: float, b: float) -> tuple[float, float]:
        return cx + u * scale, cy - b * scale

    fg_rgb = _hex_to_rgb(fg)
    accent_rgb = _hex_to_rgb(accent)
    axis_rgb = _alpha_blend(fg_rgb, _hex_to_rgb(bg), 0.45)
    asym_rgb = _alpha_blend(fg_rgb, _hex_to_rgb(bg), 0.35)

    if include_axes:
        # τ axis (vertical), u axis (horizontal)
        axis_w = max(1, int(size / 256))
        ax0x, _ = to_px(-2.1, 0)
        ax1x, _ = to_px(2.1, 0)
        _, ay0y = to_px(0, -0.4)
        _, ay1y = to_px(0, 2.5)
        draw.line([(ax0x, cy), (ax1x, cy)], fill=axis_rgb, width=axis_w)
        draw.line([(cx, ay0y), (cx, ay1y)], fill=axis_rgb, width=axis_w)

    # Asymptotes (dashed) -- only shown alongside axes
    if include_axes:
        _draw_dashed_line(draw, (cx, cy), to_px(2.1, 2.1), asym_rgb,
                          width=max(1, int(size / 256)), dash=max(2, int(size / 64)))
        _draw_dashed_line(draw, (cx, cy), to_px(-2.1, 2.1), asym_rgb,
                          width=max(1, int(size / 256)), dash=max(2, int(size / 64)))

    # Hyperbola branches
    sw = max(2, int(size / 96))
    pts_right = [to_px(u, b) for (u, b) in hyperbola_points(c=1.0, u_max=1.9, n=200)]
    pts_left = [to_px(-u, b) for (u, b) in hyperbola_points(c=1.0, u_max=1.9, n=200)]
    draw.line(pts_right, fill=accent_rgb, width=sw, joint="curve")
    draw.line(pts_left, fill=accent_rgb, width=sw, joint="curve")

    # Vertex dot
    vx, vy = to_px(0, 1.0)
    r = max(3, int(size / 56))
    draw.ellipse([(vx - r, vy - r), (vx + r, vy + r)], fill=accent_rgb)

    # Labels
    if include_labels and size >= 96:
        font_size = max(11, int(size * 0.075))
        font = _load_font(font_size, italic=True)
        tau_x, tau_y = to_px(0.18, 2.35)
        u_x, u_y = to_px(1.95, -0.15)
        draw.text((tau_x, tau_y - font_size * 0.85), "τ", fill=fg_rgb, font=font)
        draw.text((u_x - font_size * 0.6, u_y - font_size * 0.85), "u", fill=fg_rgb, font=font)

    if include_wordmark:
        word_size = int(size * 0.34)
        sub_size = int(size * 0.13)
        word_font = _load_font(word_size, bold=True)
        sub_font = _load_font(sub_size, italic=True)
        text_x = size + size * 0.25
        word_y = canvas_h / 2 - word_size * 0.55  # vertically centered
        sub_y = word_y + word_size * 1.05
        # "Py" in fg, "Physics" in accent — draw separately
        py_w = draw.textlength("Py", font=word_font)
        draw.text((text_x, word_y), "Py", fill=fg_rgb, font=word_font)
        draw.text((text_x + py_w, word_y), "Physics", fill=accent_rgb, font=word_font)
        draw.text((text_x, sub_y), "dual relativity & quantum mechanics",
                  fill=_alpha_blend(fg_rgb, _hex_to_rgb(bg), 0.65), font=sub_font)

    return img


def _hex_to_rgb(h: str) -> tuple[int, int, int]:
    h = h.lstrip("#")
    return tuple(int(h[i: i + 2], 16) for i in (0, 2, 4))


def _alpha_blend(fg: tuple[int, int, int], bg: tuple[int, int, int], alpha: float) -> tuple[int, int, int]:
    return tuple(int(fg[i] * alpha + bg[i] * (1 - alpha)) for i in range(3))


def _draw_dashed_line(draw: ImageDraw.ImageDraw, start: tuple[float, float], end: tuple[float, float],
                       fill, width: int, dash: int):
    x0, y0 = start
    x1, y1 = end
    length = math.hypot(x1 - x0, y1 - y0)
    n_dashes = max(1, int(length / (dash * 2)))
    for i in range(n_dashes):
        t0 = i * 2 * dash / length
        t1 = min(1.0, t0 + dash / length)
        sx, sy = x0 + (x1 - x0) * t0, y0 + (y1 - y0) * t0
        ex, ey = x0 + (x1 - x0) * t1, y0 + (y1 - y0) * t1
        draw.line([(sx, sy), (ex, ey)], fill=fill, width=width)


def _load_font(size: int, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    # Try a sequence of common serif fonts; fall back to default.
    candidates = []
    if italic and bold:
        candidates += ["DejaVuSerif-BoldItalic.ttf"]
    elif italic:
        candidates += ["DejaVuSerif-Italic.ttf", "DejaVuSans-Oblique.ttf"]
    elif bold:
        candidates += ["DejaVuSerif-Bold.ttf"]
    candidates += ["DejaVuSerif.ttf", "DejaVuSans.ttf"]
    for name in candidates:
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


# ───────────────────── orchestration ─────────────────────

def build_all(dry_run: bool = False) -> None:
    targets = [
        # SVG variants
        ("logo_dark.svg", lambda: svg_symbol(size=256, bg=DARK_BG, fg=LIGHT_FG, accent=YELLOW)),
        ("logo_light.svg", lambda: svg_symbol(size=256, bg=LIGHT_BG, fg=DARK_FG, accent=YELLOW)),
        ("logo.svg", lambda: svg_symbol(size=256, bg=DARK_BG, fg=LIGHT_FG, accent=YELLOW)),
        ("logo_horizontal_dark.svg",
         lambda: svg_horizontal(width=1024, height=320, bg=DARK_BG, fg=LIGHT_FG, accent=YELLOW)),
        ("logo_horizontal_light.svg",
         lambda: svg_horizontal(width=1024, height=320, bg=LIGHT_BG, fg=DARK_FG, accent=YELLOW)),
    ]
    for name, gen in targets:
        path = OUT_DIR / name
        if dry_run:
            print(f"(dry-run) would write {path} ({len(gen())} bytes)")
        else:
            path.write_text(gen(), encoding="utf-8")
            print(f"Wrote {path}")

    png_targets = [
        ("logo_32.png", 32, False),
        ("logo_64.png", 64, False),
        ("logo_256.png", 256, False),
        ("social_preview_1280x640.png", 640, True),  # horizontal lockup at h=640 -> w=1920; we crop
        ("podcast_cover_3000.png", 3000, False),
    ]
    for name, size, horiz in png_targets:
        path = OUT_DIR / name
        if dry_run:
            print(f"(dry-run) would render {path} at size {size}, horizontal={horiz}")
            continue
        if name.startswith("social_preview"):
            # Render the horizontal lockup at the exact 1280x640 target.
            target_w, target_h = 1280, 640
            sym_size = int(target_h * 0.55)  # symbol fills 55% of height; leaves room for text
            img = png_symbol(
                size=sym_size, bg=DARK_BG, fg=LIGHT_FG, accent=YELLOW,
                include_wordmark=True, canvas_size=(target_w, target_h),
            )
        else:
            # For small favicons, suppress axes/asymptotes so the curve reads at 32px.
            is_small = size < 96
            img = png_symbol(size=size, bg=DARK_BG, fg=LIGHT_FG, accent=YELLOW,
                              include_wordmark=False,
                              include_axes=not is_small,
                              include_labels=not is_small)
        img.save(path, "PNG", optimize=True)
        print(f"Wrote {path}  ({path.stat().st_size:,} bytes)")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dry-run", action="store_true", help="print actions without writing files")
    args = p.parse_args()
    build_all(dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
