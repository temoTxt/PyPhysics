"""
Phase 4 cross-paper synthesis: the dual-theory framework across Maxwell, DRQM I, TCEP,
plus the three findings flagged for author review.

Closing scene of the Phase 1-4 animation series.  Cross-references:
  - Roadmapping/Equation_Verification/README.md
  - Roadmapping/Equation_Verification/FINDINGS_for_author_review.md

Render:
    uv run manim -pql manim_scenes/synthesis_tour.py SynthesisTour
"""

from manim import *


class SynthesisTour(Scene):
    """A short tour: the dual-theory framework + the three findings."""

    def construct(self):
        # ---- Opening ----
        title = Tex(r"\textbf{The dual theory of relativity:}\\ \textbf{a tour}", font_size=44)
        subtitle = Tex(
            r"Verification campaign summary, branch \texttt{3-verify-equations-mathematica}",
            font_size=22,
        )
        subtitle.next_to(title, DOWN, buff=0.4)
        self.play(Write(title), Write(subtitle))
        self.wait(1.5)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ---- The central identity ----
        central_label = Tex(r"\textbf{Central identity:}", font_size=32).to_edge(UP)
        b_def = MathTex(
            r"b \;=\; \sqrt{c^2 + \mathbf{u}^2}, \qquad \mathbf{u} = \frac{d\mathbf{x}}{d\tau} = \gamma(\mathbf{w})\,\mathbf{w}",
            font_size=40,
            color=YELLOW,
        )
        duality = MathTex(
            r"\frac{\mathbf{w}}{c} \;=\; \frac{\mathbf{u}}{b}, \qquad \frac{1}{c}\frac{\partial}{\partial t} \;=\; \frac{1}{b}\frac{\partial}{\partial \tau}",
            font_size=36,
            color=YELLOW,
        )
        central = VGroup(b_def, duality).arrange(DOWN, buff=0.6).shift(DOWN * 0.3)
        self.play(Write(central_label))
        self.play(Write(b_def))
        self.play(Write(duality))
        self.wait(2.5)
        self.play(FadeOut(central_label), FadeOut(b_def), FadeOut(duality))

        # ---- The three pillars ----
        pillars_title = Tex(r"\textbf{Three pillars of the dual theory:}", font_size=30).to_edge(UP)

        pillar1 = VGroup(
            Tex(r"\textbf{Maxwell}", font_size=26, color=BLUE),
            Tex(r"$\nabla\!\times\!\mathbf{E} = -\tfrac{1}{b}\partial_\tau \mathbf{B}$", font_size=22),
            Tex(r"Eq. (3$'$): proper-time Maxwell", font_size=18),
        ).arrange(DOWN, buff=0.2)

        pillar2 = VGroup(
            Tex(r"\textbf{DRQM I}", font_size=26, color=GREEN),
            MathTex(r"K = \tfrac{H^2}{2mc^2} + \tfrac{mc^2}{2}", font_size=22),
            Tex(r"Eq. (I.6): canonical $\tau$-Hamiltonian", font_size=18),
        ).arrange(DOWN, buff=0.2)

        pillar3 = VGroup(
            Tex(r"\textbf{TCEP}", font_size=26, color=ORANGE),
            MathTex(r"\partial \tau' / \partial \tau = (\bar b / b)(r/s)", font_size=22),
            Tex(r"Eq. (3.36): retarded-time kinematics", font_size=18),
        ).arrange(DOWN, buff=0.2)

        pillars = VGroup(pillar1, pillar2, pillar3).arrange(RIGHT, buff=0.8).shift(DOWN * 0.3)
        self.play(Write(pillars_title))
        self.play(FadeIn(pillar1, shift=UP * 0.2))
        self.wait(0.8)
        self.play(FadeIn(pillar2, shift=UP * 0.2))
        self.wait(0.8)
        self.play(FadeIn(pillar3, shift=UP * 0.2))
        self.wait(2.0)
        self.play(FadeOut(pillars_title), FadeOut(pillars))

        # ---- The three findings ----
        findings_title = Tex(r"\textbf{Three findings flagged for author review:}", font_size=30).to_edge(UP)

        f1 = VGroup(
            Tex(r"\textbf{(1)} Maxwell Eq.\ (24)", font_size=24, color=BLUE),
            Tex(
                r"Two typos: missing factor of $c$ in $-e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}/(2m)$;\\"
                r"missing $+V^2/(2mc^2)$ term entirely.",
                font_size=20,
            ),
        ).arrange(DOWN, buff=0.15, aligned_edge=LEFT)

        f2 = VGroup(
            Tex(r"\textbf{(2)} DRQM I Eq.\ (III.22) -- the $g$-factor", font_size=24, color=GREEN),
            Tex(
                r"Stated $r_e = 0.499857150068631\,r_0$ gives $g = -2.0005714$,\\"
                r"\textbf{not} the experimental $-2.00231930436256$.  Required $r_e \approx 0.499420510\,r_0$.",
                font_size=20,
            ),
        ).arrange(DOWN, buff=0.15, aligned_edge=LEFT)

        f3 = VGroup(
            Tex(r"\textbf{(3)} TCEP Eq.\ (4.16) -- group velocity", font_size=24, color=ORANGE),
            Tex(
                r"Paper writes $v_g = v_g' - v$; algebra (and paper's own commentary) give $v_g = v_g' + v$.",
                font_size=20,
            ),
        ).arrange(DOWN, buff=0.15, aligned_edge=LEFT)

        findings = VGroup(f1, f2, f3).arrange(DOWN, buff=0.6, aligned_edge=LEFT).shift(DOWN * 0.3)
        self.play(Write(findings_title))
        self.play(FadeIn(f1, shift=LEFT * 0.2))
        self.wait(2.0)
        self.play(FadeIn(f2, shift=LEFT * 0.2))
        self.wait(2.0)
        self.play(FadeIn(f3, shift=LEFT * 0.2))
        self.wait(2.5)
        self.play(FadeOut(findings_title), FadeOut(findings))

        # ---- Methodology + closing ----
        method = VGroup(
            Tex(r"\textbf{Methodology:}", font_size=30),
            Tex(r"Every equation: grad-student-level expanded derivation +", font_size=22),
            Tex(r"single-line Mathematica check via the Wolfram MCP server.", font_size=22),
            Tex(r"\textit{Each finding above is reproducible from the recorded MCP input.}", font_size=20, color=YELLOW),
        ).arrange(DOWN, buff=0.25)
        self.play(FadeIn(method))
        self.wait(3.0)
        self.play(FadeOut(method))

        closing = VGroup(
            Tex(r"\textbf{11 physics papers verified.}", font_size=36, color=YELLOW),
            Tex(r"3 substantive findings, 1 prior open question resolved.", font_size=26),
            Tex(r"\texttt{See FINDINGS\_for\_author\_review.md}", font_size=22),
        ).arrange(DOWN, buff=0.4)
        self.play(Write(closing))
        self.wait(4.0)
