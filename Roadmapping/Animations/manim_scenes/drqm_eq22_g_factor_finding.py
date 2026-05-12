"""
DRQM I Eq. (III.22) numerical finding: the paper's r_e value does NOT reproduce the
experimental g_e.

Plots g_r(r/r_0) and shows:
  - the experimental g_e = -2.00231930436256 as a horizontal line
  - the paper's claimed r_e/r_0 = 0.499857150068631 (gives g = -2.0005714, MISS)
  - the correct r_e/r_0 ~= 0.499420510 (gives g_e, HIT)

This is Finding 2 in Roadmapping/Equation_Verification/FINDINGS_for_author_review.md.

Derivation reference:
  Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md#eqs-iii21iii23--the-g-factor-formula--numerical-failure
  Roadmapping/Equation_Verification/FINDINGS_for_author_review.md

Render:
    uv run manim -pql manim_scenes/drqm_eq22_g_factor_finding.py DRQMGFactorFinding
"""

from manim import *
import numpy as np


# Numerical constants
G_EXPERIMENTAL = -2.00231930436256
PAPER_RE_OVER_R0 = 0.499857150068631  # what the paper claims
CORRECT_RE_OVER_R0 = 0.4994205099128318  # what reproduces g_experimental


def g_r(x: float) -> float:
    """Gill's formula:  g_r = 2 [1 - 4 r_0 / (2 r + r_0)],  with x = r/r_0."""
    return 2.0 * (1.0 - 4.0 / (2.0 * x + 1.0))


class DRQMGFactorFinding(Scene):
    """Visualize the numerical finding: paper's r_e doesn't hit g_experimental."""

    def construct(self):
        title = Tex(r"\textbf{Finding 2: the $g$-factor numerical claim doesn't reproduce}", font_size=34)
        subtitle = Tex(
            r"DRQM I Eq.\ (III.22) ----- paper's $r_e = 0.499857150068631\,r_0$ \emph{vs.} the value that fits experiment",
            font_size=20,
        )
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.5)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ---- The formula card ----
        formula_card = MathTex(
            r"g_r(r) \;=\; 2\!\left[1 - \frac{4 r_0}{2r + r_0}\right]",
            font_size=44,
            color=YELLOW,
        ).to_edge(UP)
        self.play(Write(formula_card))
        self.wait(0.5)

        # ---- Axes ----
        # Focus on a small window around r/r_0 = 0.5 where the experimental value sits.
        x_min, x_max = 0.45, 0.55
        y_min, y_max = -2.05, -1.96

        axes = Axes(
            x_range=[x_min, x_max, 0.02],
            y_range=[y_min, y_max, 0.01],
            x_length=10,
            y_length=4.5,
            tips=False,
            axis_config={"font_size": 18, "include_numbers": False},
            x_axis_config={"numbers_to_include": [0.46, 0.48, 0.50, 0.52, 0.54], "decimal_number_config": {"num_decimal_places": 2}},
            y_axis_config={"numbers_to_include": [-2.04, -2.02, -2.00, -1.98], "decimal_number_config": {"num_decimal_places": 2}},
        ).shift(DOWN * 0.5)

        x_label = MathTex(r"r_e / r_0", font_size=28).next_to(axes.x_axis.get_end(), RIGHT, buff=0.2)
        y_label = MathTex(r"g_r", font_size=28).next_to(axes.y_axis.get_end(), UP, buff=0.2)

        self.play(Create(axes), Write(x_label), Write(y_label))
        self.wait(0.3)

        # ---- Plot g_r ----
        # Avoid the singularity at r/r_0 = -1/2 (well outside our window).
        g_curve = axes.plot(g_r, x_range=[x_min, x_max], color=BLUE, stroke_width=4)
        g_curve_label = MathTex(r"g_r(r_e/r_0)", font_size=24, color=BLUE).move_to(axes.c2p(0.535, -1.99))
        self.play(Create(g_curve), Write(g_curve_label))
        self.wait(1.0)

        # ---- Horizontal line: experimental g_e ----
        g_exp_line = axes.plot(
            lambda x: G_EXPERIMENTAL, x_range=[x_min, x_max], color=RED, stroke_width=3
        )
        g_exp_label = Tex(
            rf"$g_e^{{\text{{exp}}}} = {G_EXPERIMENTAL:.11f}$",
            font_size=22,
            color=RED,
        ).move_to(axes.c2p(0.47, -2.005))
        self.play(Create(g_exp_line), Write(g_exp_label))
        self.wait(1.5)

        # ---- The paper's claim: vertical line at r_e = 0.499857... r_0 ----
        paper_x_top = axes.c2p(PAPER_RE_OVER_R0, y_max)
        paper_x_bot = axes.c2p(PAPER_RE_OVER_R0, y_min)
        paper_vline = DashedLine(paper_x_bot, paper_x_top, color=ORANGE, stroke_width=3)

        paper_intersect_y = g_r(PAPER_RE_OVER_R0)  # ~ -2.00057148
        paper_intersect_pt = axes.c2p(PAPER_RE_OVER_R0, paper_intersect_y)
        paper_dot = Dot(paper_intersect_pt, color=ORANGE, radius=0.09)

        paper_label = VGroup(
            Tex(r"\textbf{Paper's $r_e$:}\ $0.499857150068631\,r_0$", font_size=18, color=ORANGE),
            MathTex(rf"\Rightarrow g_r = {paper_intersect_y:.10f}", font_size=18, color=ORANGE),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1).move_to(axes.c2p(0.535, -2.04))

        self.play(Create(paper_vline))
        self.play(FadeIn(paper_dot, scale=2.0))
        self.play(Write(paper_label))
        self.wait(2.0)

        # Highlight the gap (paper's g vs experimental g)
        gap_indicator = DoubleArrow(
            paper_intersect_pt,
            axes.c2p(PAPER_RE_OVER_R0, G_EXPERIMENTAL),
            color=YELLOW,
            stroke_width=3,
            buff=0,
            tip_length=0.15,
        )
        gap_label = Tex(
            rf"$\Delta g \approx {abs(paper_intersect_y - G_EXPERIMENTAL):.6f}$\ \ (\textit{{misses by 7 orders of magnitude}})",
            font_size=18,
            color=YELLOW,
        ).move_to(axes.c2p(0.47, -1.965))
        self.play(GrowArrow(gap_indicator), Write(gap_label))
        self.wait(2.5)

        # ---- The correct value ----
        # Clear out the paper-related callouts to make room
        self.play(
            FadeOut(gap_indicator),
            FadeOut(gap_label),
            FadeOut(paper_label),
        )

        correct_x_top = axes.c2p(CORRECT_RE_OVER_R0, y_max)
        correct_x_bot = axes.c2p(CORRECT_RE_OVER_R0, y_min)
        correct_vline = DashedLine(correct_x_bot, correct_x_top, color=GREEN, stroke_width=3)

        correct_intersect_pt = axes.c2p(CORRECT_RE_OVER_R0, G_EXPERIMENTAL)
        correct_dot = Dot(correct_intersect_pt, color=GREEN, radius=0.11)

        correct_label = VGroup(
            Tex(r"\textbf{Required $r_e$:}\ $\approx 0.499420510\,r_0$", font_size=20, color=GREEN),
            MathTex(rf"\Rightarrow g_r = g_e^{{\text{{exp}}}}", font_size=20, color=GREEN),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1).move_to(axes.c2p(0.47, -2.04))

        self.play(Create(correct_vline))
        self.play(FadeIn(correct_dot, scale=2.5))
        self.play(Write(correct_label))
        self.wait(2.5)

        # ---- Final caption ----
        self.play(FadeOut(formula_card))
        final_card = VGroup(
            Tex(r"\textbf{The formula structure is correct.}\ ", font_size=24),
            Tex(r"\textbf{The published digits of $r_e$ are off in the 4\textsuperscript{th} decimal place.}", font_size=24),
            Tex(r"\textit{Likely a transcription error from a working notebook.}", font_size=20),
        ).arrange(DOWN, buff=0.2).to_edge(UP)
        self.play(Write(final_card))
        self.wait(4.0)
