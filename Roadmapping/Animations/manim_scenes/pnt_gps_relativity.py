"""
GPS satellite clock relativistic corrections: SR -7 us/day, GR +45 us/day,
net +38 us/day. Engineering solution: pre-detune the satellite clocks.
Dual-framework SR rederivation: same answer to within 10^-10.

Chapter reference:
  Roadmapping/History/Forward/07_PNT_GPS_SLR_QKD.md#b-gps-as-the-most-precise-everyday-test-of-relativity

Render:
    uv run manim -pql --media_dir rendered manim_scenes/pnt_gps_relativity.py PntGpsRelativity
    uv run manim -qh  --media_dir rendered manim_scenes/pnt_gps_relativity.py PntGpsRelativity
"""

from manim import (
    BLUE, DOWN, FadeOut, GREEN, LEFT, MathTex, RED, RIGHT, Scene, Tex, Transform,
    UP, VGroup, WHITE, Write, YELLOW,
)


class PntGpsRelativity(Scene):
    """GPS clock corrections: SR vs GR vs net vs dual framework."""

    def construct(self):
        title = Tex(r"\textbf{GPS clock corrections: SR vs GR vs dual}", font_size=38)
        subtitle = Tex(r"Chapter 7 \S B, PNT", font_size=24)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # SR contribution
        sr_label = Tex(r"\textbf{Special-relativistic time dilation}", color=BLUE, font_size=28).to_edge(UP)
        sr_setup = Tex(
            r"Satellite orbital velocity $v_{\rm sat} \approx 3.87$ km/s $\to$ "
            r"$v/c \approx 1.3\times 10^{-5}$",
            font_size=22,
        ).next_to(sr_label, DOWN, buff=0.4)
        sr_eq = MathTex(
            r"\Delta f / f \;=\; -\frac{v^2}{2 c^2} \;\approx\; -8.3 \times 10^{-11}",
            font_size=28,
        ).next_to(sr_setup, DOWN, buff=0.4)
        sr_result = MathTex(
            r"\Rightarrow\; -7\ \mu\text{s/day (satellite clock slower)}",
            color=BLUE, font_size=28,
        ).next_to(sr_eq, DOWN, buff=0.4)
        self.play(Write(sr_label))
        self.play(Write(sr_setup), Write(sr_eq))
        self.wait(0.5)
        self.play(Write(sr_result))
        self.wait(1.5)

        # GR contribution
        gr_label = Tex(r"\textbf{General-relativistic gravitational time dilation}",
                       color=GREEN, font_size=28).next_to(sr_result, DOWN, buff=0.7)
        gr_setup = Tex(
            r"Satellite altitude $h \approx 20{,}200$ km $\to$ $\Delta U/c^2 \approx +5.3\times 10^{-10}$",
            font_size=22,
        ).next_to(gr_label, DOWN, buff=0.4)
        gr_result = MathTex(
            r"\Rightarrow\; +45\ \mu\text{s/day (satellite clock faster)}",
            color=GREEN, font_size=28,
        ).next_to(gr_setup, DOWN, buff=0.4)
        self.play(Write(gr_label))
        self.play(Write(gr_setup))
        self.wait(0.5)
        self.play(Write(gr_result))
        self.wait(1.5)

        # Clear, then show net + engineering solution
        self.play(FadeOut(sr_label), FadeOut(sr_setup), FadeOut(sr_eq), FadeOut(sr_result),
                  FadeOut(gr_label), FadeOut(gr_setup), FadeOut(gr_result))

        net = MathTex(
            r"\boxed{\;\text{Net: } +38\ \mu\text{s/day}\;}",
            color=YELLOW, font_size=44,
        ).to_edge(UP)
        net_note = Tex(
            r"\textit{Without correction, GPS position error grows by} $\sim 10$ km/day.",
            font_size=24,
        ).next_to(net, DOWN, buff=0.5)
        engineering = Tex(
            r"\textbf{Engineering solution:} satellite atomic clocks are \\ "
            r"\textit{pre-detuned at manufacture} to tick at the GPS-Time rate \\ "
            r"after both relativistic corrections.",
            font_size=24,
        ).next_to(net_note, DOWN, buff=0.6)
        self.play(Write(net))
        self.play(Write(net_note))
        self.wait(1.0)
        self.play(Write(engineering))
        self.wait(2.5)
        self.play(FadeOut(net), FadeOut(net_note), FadeOut(engineering))

        # Dual-framework rederivation
        dual_title = Tex(r"\textbf{Dual-framework SR rederivation}",
                         color=YELLOW, font_size=32).to_edge(UP)
        dual_step1 = Tex(
            r"Use Gill's Maxwell paper Eq.\ (9): $t = (1/c)\int b(s)\,ds$",
            font_size=24,
        ).next_to(dual_title, DOWN, buff=0.5)
        dual_step2 = MathTex(
            r"b = \sqrt{c^2 + u^2}, \quad u/c \sim 1.3\times 10^{-5}",
            font_size=26,
        ).next_to(dual_step1, DOWN, buff=0.4)
        dual_step3 = MathTex(
            r"b/c - 1 \;\sim\; u^2/(2c^2) \;\sim\; 8.5 \times 10^{-11}",
            font_size=26,
        ).next_to(dual_step2, DOWN, buff=0.4)
        dual_result = MathTex(
            r"\Rightarrow\; \text{dual SR prediction} \;\approx\; -7\ \mu\text{s/day} \;\pm\; O(10^{-10})",
            color=YELLOW, font_size=28,
        ).next_to(dual_step3, DOWN, buff=0.4)
        dual_gr = Tex(
            r"\textbf{GR piece is \#gill-silent:} taken intact from standard GR. \\ "
            r"\textit{Net dual prediction: same $+38\ \mu$s/day as standard, to within $10^{-10}$.}",
            font_size=22, color=WHITE,
        ).next_to(dual_result, DOWN, buff=0.5)
        precision_note = Tex(
            r"\textit{Distinguishable in $\sim 2030$s by optical-lattice satellite clocks at $10^{-19}$.}",
            font_size=22, color=WHITE,
        ).next_to(dual_gr, DOWN, buff=0.4)

        self.play(Write(dual_title))
        self.play(Write(dual_step1))
        self.wait(0.4)
        self.play(Write(dual_step2))
        self.wait(0.4)
        self.play(Write(dual_step3))
        self.wait(0.4)
        self.play(Write(dual_result))
        self.wait(1.0)
        self.play(Write(dual_gr))
        self.wait(0.8)
        self.play(Write(precision_note))
        self.wait(3.0)
        self.play(FadeOut(dual_title), FadeOut(dual_step1), FadeOut(dual_step2),
                  FadeOut(dual_step3), FadeOut(dual_result), FadeOut(dual_gr),
                  FadeOut(precision_note))
