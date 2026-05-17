"""
Klein-Gordon equation: derivation, the non-positive-definite conserved current j^0,
and the dual-framework reformulation via positive-definite K = H^2/(2mc^2) + mc^2/2.

Chapter reference:
  Roadmapping/History/04_quantum_mechanics_1925_1948.md#4-key-derivations-worth-animating

The narrative:
  1. E^2 = p^2 c^2 + m^2 c^4 -> Klein-Gordon equation (squaring breaks j^0 sign).
  2. Conserved current j^0 ~ psi* d_t psi - d_t psi* psi: not positive-definite.
  3. Dual reformulation: K = H^2/(2mc^2) + mc^2/2 -- positive-definite by construction.
  4. Probability density from dual framework: psi^dagger psi, manifestly positive.
  5. Final summary: dissolution of the negative-probability problem.

Render (from Roadmapping/Animations):
    uv run manim -pql --media_dir rendered manim_scenes/hist_klein_gordon_vs_dual.py HistKleinGordonVsDual
    uv run manim -qh  --media_dir rendered manim_scenes/hist_klein_gordon_vs_dual.py HistKleinGordonVsDual
"""

from manim import (
    BLUE, DOWN, FadeIn, FadeOut, GREEN, LEFT, MathTex, RED, RIGHT, Scene, Tex,
    Transform, UP, VGroup, WHITE, Write, YELLOW,
)


class HistKleinGordonVsDual(Scene):
    """Klein-Gordon negative-probability problem; dual-framework resolution."""

    def construct(self):
        # ───── Title ─────
        title = Tex(r"\textbf{Klein--Gordon vs. dual: the negative-probability problem}", font_size=36)
        subtitle = Tex(r"Chapter 4, Quantum Mechanics (1925--1948)", font_size=24)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ───── Klein-Gordon derivation ─────
        kg_step1 = MathTex(r"E^2 \;=\; p^2 c^2 + m^2 c^4", font_size=34).to_edge(UP)
        kg_step2 = MathTex(
            r"E \to i\hbar \partial_t, \quad \mathbf{p} \to -i\hbar \nabla",
            font_size=28,
        ).next_to(kg_step1, DOWN, buff=0.5)
        kg_step3 = MathTex(
            r"\left(\Box + \frac{m^2 c^2}{\hbar^2}\right)\psi \;=\; 0",
            color=BLUE, font_size=36,
        ).next_to(kg_step2, DOWN, buff=0.5)
        kg_label = Tex(r"Klein--Gordon equation (Klein 1926, Gordon 1926)",
                       color=BLUE, font_size=22).next_to(kg_step3, DOWN, buff=0.3)

        self.play(Write(kg_step1))
        self.wait(0.6)
        self.play(Write(kg_step2))
        self.wait(0.5)
        self.play(Write(kg_step3), Write(kg_label))
        self.wait(1.5)
        self.play(FadeOut(kg_step1), FadeOut(kg_step2), FadeOut(kg_step3), FadeOut(kg_label))

        # ───── The negative-probability problem ─────
        prob_title = Tex(r"\textbf{The negative-probability problem}", color=RED, font_size=34).to_edge(UP)
        prob_eq = MathTex(
            r"j^0 \;=\; \frac{i\hbar}{2m c^2}\left(\psi^* \partial_t \psi - \psi \, \partial_t \psi^*\right)",
            font_size=30,
        ).next_to(prob_title, DOWN, buff=0.6)
        prob_problem = Tex(
            r"\textit{Not positive-definite.} \\ "
            r"There exist initial $\psi$ for which $j^0 < 0$ in some spacetime region. \\ "
            r"Cannot be interpreted as a probability density.",
            color=RED, font_size=26,
        ).next_to(prob_eq, DOWN, buff=0.6)
        schroedinger_note = Tex(
            r"Schr\"odinger himself encountered this in early 1926 and retreated \\ "
            r"to the non-relativistic limit. Dirac (1928) solved it by going to a \\ "
            r"first-order-in-time equation with a 4-spinor structure.",
            font_size=22, color=WHITE,
        ).next_to(prob_problem, DOWN, buff=0.6)

        self.play(Write(prob_title))
        self.play(Write(prob_eq))
        self.wait(1.0)
        self.play(Write(prob_problem))
        self.wait(1.5)
        self.play(Write(schroedinger_note))
        self.wait(2.5)
        self.play(FadeOut(prob_title), FadeOut(prob_eq), FadeOut(prob_problem), FadeOut(schroedinger_note))

        # ───── Dual-framework reformulation ─────
        dual_title = Tex(r"\textbf{Dual-framework reformulation}", color=GREEN, font_size=34).to_edge(UP)
        dual_step1 = Tex(
            r"Tepper Gill's \textit{Foundations II Classical} construction:",
            font_size=24,
        ).next_to(dual_title, DOWN, buff=0.6)
        dual_K = MathTex(
            r"K \;=\; \frac{H^2}{2 m c^2} + \frac{m c^2}{2}",
            color=YELLOW, font_size=40,
        ).next_to(dual_step1, DOWN, buff=0.5)
        dual_pos = Tex(
            r"\textit{Positive-definite by construction.} Bounded below by $mc^2/2$. \\ "
            r"No negative-energy modes.",
            color=GREEN, font_size=24,
        ).next_to(dual_K, DOWN, buff=0.5)
        dual_density = MathTex(
            r"\rho_{\rm dual} \;=\; \psi^\dagger \psi \;\geq\; 0",
            color=YELLOW, font_size=36,
        ).next_to(dual_pos, DOWN, buff=0.5)

        self.play(Write(dual_title))
        self.play(Write(dual_step1))
        self.wait(0.5)
        self.play(Write(dual_K))
        self.wait(0.8)
        self.play(Write(dual_pos))
        self.wait(0.8)
        self.play(Write(dual_density))
        self.wait(2.0)
        self.play(FadeOut(dual_title), FadeOut(dual_step1), FadeOut(dual_K),
                  FadeOut(dual_pos), FadeOut(dual_density))

        # ───── Closing summary ─────
        summary = Tex(
            r"\textbf{Headline conceptual payoff:} \\ "
            r"\textit{Klein--Gordon's negative-probability problem does not arise} \\ "
            r"\textit{in the dual framework. Positive-definite probability density follows} \\ "
            r"\textit{from the positive-definite Hamiltonian, without spinor-doubling.}",
            color=YELLOW, font_size=28,
        )
        cross_ref = Tex(
            r"See Foundations II Classical for the construction; \\ "
            r"DRQM I Sec III for the relativistic eigenvalue problem.",
            font_size=22, color=WHITE,
        ).next_to(summary, DOWN, buff=0.5)
        self.play(Write(summary))
        self.wait(1.5)
        self.play(Write(cross_ref))
        self.wait(3.0)
        self.play(FadeOut(summary), FadeOut(cross_ref))
