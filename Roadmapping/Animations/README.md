# Animations of equation derivations (Manim)

Manim-CE animations of the load-bearing derivations from the Gill physics corpus. Each scene corresponds to a specific equation cluster in `Roadmapping/Equation_Verification/` and walks through the algebra in animated form.

This sub-project is **isolated from the parent `pyphysics` project**: it has its own `pyproject.toml` and `.python-version` pinned to Python 3.13. The parent project pins Python 3.14, but PyAV (a Manim dependency) does not yet have Python-3.14 wheels.

## Status

**Phase 1 — Maxwell paper** (this PR / branch `5-manim-animations`):

| Scene file | Topic | Verification doc cross-ref |
|---|---|---|
| [`manim_scenes/maxwell_eq01_02_duality.py`](manim_scenes/maxwell_eq01_02_duality.py) | Eqs. (1)–(2): velocity duality $\mathbf{w}/c = \mathbf{u}/b$ and time-derivative duality | [Maxwell verification](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-1--velocity-duality) |
| [`manim_scenes/maxwell_eq03_substitution.py`](manim_scenes/maxwell_eq03_substitution.py) | Eq. (3) → Eq. (3′): standard Maxwell → proper-time Maxwell | [Maxwell verification](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-3-%E2%80%94-proper-time-equivalent-maxwells-equations) |
| [`manim_scenes/maxwell_eq04_dissipative.py`](manim_scenes/maxwell_eq04_dissipative.py) | Eq. (4): emergence of the dissipative term $-(\mathbf{u}\!\cdot\!\mathbf{a})/b^4 \cdot \partial_\tau\mathbf{E}$ via curl-of-curl | [Maxwell verification](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-4--dual-wave-equations-with-dissipative-term) |
| [`manim_scenes/maxwell_eq05_06_photonmass.py`](manim_scenes/maxwell_eq05_06_photonmass.py) | Eqs. (5)–(6): Liouville substitution $\psi = (b/c)^{1/2}\Psi_{\rm new}$ and the effective photon mass $\mu$ | [Maxwell verification](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-6--effective-photon-mass) |
| [`manim_scenes/maxwell_eq10_11_hyperboloid.py`](manim_scenes/maxwell_eq10_11_hyperboloid.py) | Eqs. (10)–(11): boost preserving the 4-velocity hyperboloid $b^2 - \mathbf{u}^2 = c^2$ | [Maxwell verification](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-10--proper-time-boosts-of-position-velocity-acceleration) |

**Phase 2–4** (future work, see [issue #5](https://github.com/temoTxt/PyPhysics/issues/5)): DRQM I g-factor; TCEP Liénard–Wiechert; cross-paper synthesis.

## Build / render

Prerequisites (system-level):
- `ffmpeg` (Manim invokes it as a subprocess for video output).
- TeX Live (`texlive-latex-base`, `texlive-latex-extra`, `texlive-fonts-recommended`, `texlive-pictures`) — Manim renders all `MathTex`/`Tex` via LaTeX.
- Cairo + Pango (`libcairo2-dev`, `libpango1.0-dev`).

Python deps are managed by `uv` from this sub-project's `pyproject.toml`:

```bash
# from the repo root
uv sync --project Roadmapping/Animations
```

This creates `Roadmapping/Animations/.venv/` and installs Manim 0.20.1 + `numpy`.

### Render a single scene

```bash
cd Roadmapping/Animations

# Low-quality preview (480p, ~10 s render):
uv run manim -ql --media_dir rendered manim_scenes/maxwell_eq01_02_duality.py MaxwellEq01_02Duality

# High-quality 1080p (the issue's acceptance criterion):
uv run manim -qh --media_dir rendered manim_scenes/maxwell_eq01_02_duality.py MaxwellEq01_02Duality

# With auto-play after render:
uv run manim -pqh --media_dir rendered manim_scenes/maxwell_eq01_02_duality.py MaxwellEq01_02Duality
```

### Render all Phase 1 scenes

```bash
cd Roadmapping/Animations
for scene_file in manim_scenes/maxwell_eq*.py; do
    scene_name=$(grep -m1 "^class " "$scene_file" | cut -d' ' -f2 | cut -d'(' -f1)
    uv run manim -qh --media_dir rendered "$scene_file" "$scene_name"
done
```

## Layout

```
Roadmapping/Animations/
├── pyproject.toml              # uv project; pinned to Python 3.13
├── .python-version             # 3.13
├── .gitignore                  # excludes .venv, partial_movie_files, 1080p renders
├── uv.lock                     # committed for reproducibility
├── README.md                   # this file
├── manim_scenes/               # one .py per equation cluster
│   ├── maxwell_eq01_02_duality.py
│   ├── maxwell_eq03_substitution.py
│   ├── maxwell_eq04_dissipative.py
│   ├── maxwell_eq05_06_photonmass.py
│   └── maxwell_eq10_11_hyperboloid.py
└── rendered/                   # Manim's output dir
    └── videos/<scene>/480p15/  # 480p preview (committed; ~400KB each)
                  /1080p60/     # high-res (gitignored; rebuild via -qh)
                  /partial_movie_files/  # cache (gitignored)
```

## Adding a new scene

1. Pick an equation cluster from one of the verification docs under `../Equation_Verification/`.
2. Create `manim_scenes/<paper>_eq<N>_<short_name>.py` with a `Scene` subclass.
3. Header docstring must point to the relevant section of the verification doc and include the `uv run manim ...` command.
4. Walk through the derivation step-by-step using `Transform` between intermediate `MathTex` forms — this is the pedagogical core.
5. Cite the Wolfram MCP verification (most scenes end with "Verified by Wolfram MCP: ...").
6. Render at `-ql` to debug, then `-qh` to confirm 1080p success.
7. Update the **Status** table in this README.

## Convention notes

- Use **`MathTex`** for math (LaTeX), **`Tex`** for prose. Always wrap LaTeX in raw strings (`r"..."`) to avoid Python string-escaping the backslashes.
- Keep `font_size` between 22 and 56 — the standard range that reads well at 1080p without crowding.
- Use `Transform(old, new)` rather than `FadeOut` + `FadeIn` for chained equations — this preserves the "this becomes that" narrative.
- Box the final result with `\boxed{...}` and color it `YELLOW` (matching the "headline result" convention from the verification docs).
- Always reference the verification doc in the docstring with a Markdown anchor link.
