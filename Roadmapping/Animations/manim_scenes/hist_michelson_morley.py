"""
Michelson-Morley: predicted fringe shift under aether-wind vs observed null,
reframed by the dual framework as u=0 -> b=c -> no asymmetry.

Chapter reference:
  Roadmapping/History/02_classical_synthesis_1860_1900.md#23-1887----michelson-and-morleys-null-result

The scene:
  1. Sketches the L-shaped interferometer with beam splitter and two mirrors.
  2. Shows the parallel-arm round-trip time under aether-wind: t_par = 2L/(c(1-v^2/c^2)).
  3. Shows the perpendicular-arm round-trip: t_perp = 2L/(c sqrt(1-v^2/c^2)).
  4. Reads off the predicted fringe shift Delta n ~ L v^2 / (lambda c^2) ~ 0.4.
  5. Records the observed bound: ~0.01 fringes. Null.
  6. Switches to the dual-framework reframing: lab-rest source -> u=0 -> b=c, both arms
     symmetric, predicted shift identically zero -- no length contraction invoked.

Render (from Roadmapping/Animations):
    uv run manim -pql --media_dir rendered manim_scenes/hist_michelson_morley.py HistMichelsonMorley
    uv run manim -qh  --media_dir rendered manim_scenes/hist_michelson_morley.py HistMichelsonMorley
"""

from manim import (
    BLUE, Create, DOWN, FadeIn, FadeOut, GREEN, LEFT, Line, MathTex, ORIGIN, RED, RIGHT,
    Scene, Square, Tex, Transform, UP, VGroup, WHITE, Write, YELLOW, Dot,
)


class HistMichelsonMorley(Scene):
    """Predicted shift under aether wind vs observed null; dual reframing as u=0."""

    def construct(self):
        # ───── Title card ─────
        title = Tex(r"\textbf{Michelson--Morley: the null that broke the aether}", font_size=42)
        subtitle = Tex(r"Chapter 2, Classical Synthesis (1860--1900)", font_size=24)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ───── Interferometer schematic ─────
        # Beam splitter at center. Horizontal arm to the right; vertical arm up.
        beam_splitter = Square(side_length=0.3, color=WHITE).move_to(ORIGIN)
        bs_label = Tex(r"BS", font_size=18).next_to(beam_splitter, DOWN+LEFT, buff=0.1)

        horiz_arm = Line(start=[0.2, 0, 0], end=[3.5, 0, 0], color=BLUE, stroke_width=2)
        vert_arm = Line(start=[0, 0.2, 0], end=[0, 3.0, 0], color=BLUE, stroke_width=2)

        mirror_h = Line(start=[3.5, -0.3, 0], end=[3.5, 0.3, 0], color=WHITE, stroke_width=4)
        mirror_v = Line(start=[-0.3, 3.0, 0], end=[0.3, 3.0, 0], color=WHITE, stroke_width=4)

        arm_h_label = MathTex(r"L_\parallel", font_size=24).next_to(horiz_arm, UP, buff=0.15)
        arm_v_label = MathTex(r"L_\perp", font_size=24).next_to(vert_arm, RIGHT, buff=0.15)

        # Earth-velocity arrow indicating motion through the aether
        v_arrow_line = Line(start=[-4, -2.5, 0], end=[-2.5, -2.5, 0], color=RED, stroke_width=4)
        v_arrow_label = MathTex(r"\vec{v}_{\oplus} \approx 30\ \text{km/s}",
                                color=RED, font_size=22).next_to(v_arrow_line, DOWN, buff=0.15)

        schematic = VGroup(beam_splitter, bs_label, horiz_arm, vert_arm, mirror_h, mirror_v,
                           arm_h_label, arm_v_label, v_arrow_line, v_arrow_label)
        self.play(Create(schematic))
        self.wait(0.6)

        # ───── Predicted parallel-arm time ─────
        t_par = MathTex(
            r"t_\parallel \;=\; \frac{L}{c - v} + \frac{L}{c + v} \;=\; "
            r"\frac{2L}{c}\,\frac{1}{1 - v^2/c^2}",
            font_size=28,
        ).to_corner(DOWN + LEFT)
        self.play(Write(t_par))
        self.wait(1.5)

        # ───── Predicted perpendicular-arm time ─────
        t_perp = MathTex(
            r"t_\perp \;=\; \frac{2L}{c}\,\frac{1}{\sqrt{1 - v^2/c^2}}",
            font_size=28,
        ).to_corner(DOWN + RIGHT)
        self.play(Write(t_perp))
        self.wait(1.5)

        # ───── Predicted fringe shift ─────
        delta_n = MathTex(
            r"\Delta n \;\approx\; \frac{L \, v^2}{\lambda \, c^2} \;\approx\; 0.4 \text{ fringes}",
            color=YELLOW, font_size=32,
        ).to_edge(DOWN)
        self.play(Transform(VGroup(t_par, t_perp), delta_n))
        self.wait(2.0)

        # ───── Observed result ─────
        observed = MathTex(
            r"\text{Observed: } \Delta n_{\rm obs} \;<\; 0.01 \text{ fringes}",
            color=RED, font_size=32,
        ).next_to(delta_n, UP, buff=0.6)
        verdict = Tex(r"\textit{The aether's rest frame is undetectable.}",
                      color=RED, font_size=26).next_to(observed, UP, buff=0.4)
        self.play(Write(observed), Write(verdict))
        self.wait(2.5)

        # Clear schematic + predicted results before the dual reframing
        self.play(
            FadeOut(schematic), FadeOut(VGroup(t_par, t_perp)),
            FadeOut(observed), FadeOut(verdict),
        )

        # ───── Dual-framework reframing ─────
        dual_title = Tex(r"\textbf{Dual-framework reframing}", color=YELLOW, font_size=38).to_edge(UP)
        dual_step1 = Tex(
            r"In Gill's dual formulation, the relevant velocity for proper-time / observer-time "
            r"conversion is the \textit{source} velocity $u$.",
            font_size=24,
        ).next_to(dual_title, DOWN, buff=0.7)
        dual_step2 = MathTex(
            r"u \;=\; 0 \quad\Longrightarrow\quad b \;=\; \sqrt{c^2 + u^2} \;=\; c",
            color=YELLOW, font_size=36,
        ).next_to(dual_step1, DOWN, buff=0.7)
        dual_step3 = Tex(
            r"Both arms experience the same proper-time parameterisation. "
            r"No asymmetry. No need for length contraction.",
            font_size=24,
        ).next_to(dual_step2, DOWN, buff=0.7)
        dual_final = MathTex(
            r"\boxed{\;\Delta n_{\rm dual} \;=\; 0 \quad \text{(identically)}\;}",
            color=YELLOW, font_size=42,
        ).next_to(dual_step3, DOWN, buff=0.7)

        self.play(Write(dual_title))
        self.play(Write(dual_step1))
        self.wait(1.0)
        self.play(Write(dual_step2))
        self.wait(1.0)
        self.play(Write(dual_step3))
        self.wait(1.0)
        self.play(Write(dual_final))
        self.wait(3.0)

        # ───── Final note ─────
        note = Tex(
            r"Same null result. Different kinematics underneath. \\ "
            r"\textit{Standard and dual SR diverge only when the source moves at $u/c \gtrsim 10^{-5}$ --- "
            r"see Chapter 7 for GPS.}",
            font_size=22, color=WHITE,
        ).to_edge(DOWN)
        self.play(Write(note))
        self.wait(3.5)

        self.play(FadeOut(dual_title), FadeOut(dual_step1), FadeOut(dual_step2),
                  FadeOut(dual_step3), FadeOut(dual_final), FadeOut(note))
