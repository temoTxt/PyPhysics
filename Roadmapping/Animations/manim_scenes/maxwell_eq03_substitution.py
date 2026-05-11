"""
Maxwell paper Eq. (3) -> Eq. (3'): substitution of duality identities into Maxwell's equations.

Derivation reference:
  Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-3--proper-time-equivalent-maxwells-equations

Standard Maxwell (Gaussian) -> proper-time Maxwell, via:
  (1/c) d/dt  ->  (1/b) d/dtau     (Eq. 2)
  rho w / c   ->  rho u / b        (Eq. 1 with rho)

Render:
    uv run manim -pql manim_scenes/maxwell_eq03_substitution.py MaxwellEq03Substitution
"""

from manim import *


class MaxwellEq03Substitution(Scene):
    """Side-by-side Standard vs Proper-time Maxwell's equations."""

    def construct(self):
        title = Tex(r"\textbf{Maxwell's Equations: standard $\to$ proper-time form}", font_size=40)
        subtitle = Tex(
            r"Gill \& Zachary (2011), Eq.\ (3) $\to$ Eq.\ (3$'$)",
            font_size=24,
        )
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ---- Standard Maxwell (LHS) ----
        std_header = Tex(r"\textbf{Standard (Gaussian)}", font_size=28, color=BLUE)
        std_eqs = MathTex(
            r"\nabla\!\cdot\!\mathbf{B} &= 0 \\",
            r"\nabla\!\cdot\!\mathbf{E} &= 4\pi\rho \\[4pt]",
            r"\nabla\!\times\!\mathbf{E} &= -\tfrac{1}{c}\tfrac{\partial \mathbf{B}}{\partial t} \\[4pt]",
            r"\nabla\!\times\!\mathbf{B} &= \tfrac{1}{c}\!\left[\tfrac{\partial \mathbf{E}}{\partial t} + 4\pi\rho\,\mathbf{w}\right]",
            font_size=28,
        )
        std_group = VGroup(std_header, std_eqs).arrange(DOWN, buff=0.4).to_edge(LEFT, buff=0.8)
        self.play(Write(std_header))
        self.play(Write(std_eqs))
        self.wait(1.5)

        # ---- The substitution rules ----
        rule_box_label = Tex(r"\textbf{Apply Eqs.\ (1) and (2):}", font_size=28).to_edge(UP)
        rule_text = MathTex(
            r"\frac{1}{c}\frac{\partial}{\partial t} \;\to\; \frac{1}{b}\frac{\partial}{\partial \tau}, \qquad \frac{\mathbf{w}}{c} \;\to\; \frac{\mathbf{u}}{b}",
            font_size=36,
            color=YELLOW,
        )
        rule_text.next_to(rule_box_label, DOWN, buff=0.3)
        self.play(Write(rule_box_label))
        self.play(Write(rule_text))
        self.wait(2.0)

        # ---- Proper-time Maxwell (RHS) appears ----
        pt_header = Tex(r"\textbf{Proper-time form}", font_size=28, color=GREEN)
        pt_eqs = MathTex(
            r"\nabla\!\cdot\!\mathbf{B} &= 0 \\",
            r"\nabla\!\cdot\!\mathbf{E} &= 4\pi\rho \\[4pt]",
            r"\nabla\!\times\!\mathbf{E} &= -\tfrac{1}{b}\tfrac{\partial \mathbf{B}}{\partial \tau} \\[4pt]",
            r"\nabla\!\times\!\mathbf{B} &= \tfrac{1}{b}\!\left[\tfrac{\partial \mathbf{E}}{\partial \tau} + 4\pi\rho\,\mathbf{u}\right]",
            font_size=28,
        )
        pt_group = VGroup(pt_header, pt_eqs).arrange(DOWN, buff=0.4).to_edge(RIGHT, buff=0.8)
        self.play(Write(pt_header))
        self.play(Write(pt_eqs))
        self.wait(2.0)

        # ---- Highlight: divergence equations unchanged ----
        self.play(FadeOut(rule_text), FadeOut(rule_box_label))
        unchanged_note = Tex(
            r"\textit{Divergence equations have no $\partial_t$ or $\mathbf{w}$ — unchanged.}",
            font_size=26,
        ).to_edge(UP)
        self.play(Write(unchanged_note))
        self.wait(2.0)

        # ---- Highlight: curl equations transformed ----
        unchanged_note2 = Tex(
            r"\textit{Curl equations: every $1/c$ becomes $1/b$; $\mathbf{w}$ becomes $\mathbf{u}$.}",
            font_size=26,
        ).to_edge(UP)
        self.play(Transform(unchanged_note, unchanged_note2))
        self.wait(2.0)

        # ---- Final caption ----
        final_caption = Tex(
            r"\textit{Algebraically equivalent; physically reinterpreted because $b = b(\tau)$.}",
            font_size=26,
            color=YELLOW,
        ).to_edge(DOWN)
        self.play(Write(final_caption))
        self.wait(3.0)
