"""
Mercury's perihelion advance from GR: Schwarzschild metric -> geodesic equation
-> effective radial potential with extra 1/r^3 term -> orbit integral
-> 43''/century, matching Le Verrier's 1859 residual exactly.

Chapter reference:
  Roadmapping/History/Forward/07_PNT_GPS_SLR_QKD.md#a-mercurys-perihelion-the-pedagogical-bridge-to-gr-clocks

Render (from Roadmapping/Animations):
    uv run manim -pql --media_dir rendered manim_scenes/pnt_mercury_perihelion.py PntMercuryPerihelion
    uv run manim -qh  --media_dir rendered manim_scenes/pnt_mercury_perihelion.py PntMercuryPerihelion
"""

from manim import (
    BLUE, DOWN, FadeOut, GREEN, LEFT, MathTex, RIGHT, Scene, Tex, Transform, UP,
    VGroup, WHITE, Write, YELLOW,
)


class PntMercuryPerihelion(Scene):
    """Mercury perihelion advance from the Schwarzschild metric."""

    def construct(self):
        title = Tex(r"\textbf{Mercury's perihelion advance from GR}", font_size=40)
        subtitle = Tex(r"Chapter 7 \S A, PNT", font_size=24)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # Step 1: Le Verrier residual
        s1 = Tex(
            r"\textbf{Le Verrier (1859):} Mercury's perihelion advances at \\ "
            r"$\approx 43''/$century in excess of Newton + planetary perturbations.",
            font_size=26,
        ).to_edge(UP)
        self.play(Write(s1))
        self.wait(1.5)

        # Step 2: Schwarzschild metric
        s2 = MathTex(
            r"ds^2 = \left(1 - \frac{2GM}{c^2 r}\right) c^2 dt^2 "
            r"- \left(1 - \frac{2GM}{c^2 r}\right)^{-1} dr^2 - r^2 d\Omega^2",
            font_size=26,
        ).next_to(s1, DOWN, buff=0.6)
        self.play(Write(s2))
        self.wait(1.0)

        # Step 3: Effective radial potential
        s3 = MathTex(
            r"V_{\rm eff}(r) = -\frac{GM}{r} + \frac{L^2}{2 r^2} "
            r"- \frac{GM L^2}{c^2 r^3}",
            font_size=30,
        ).next_to(s2, DOWN, buff=0.6)
        s3_note = Tex(
            r"\textit{The $-GM L^2 / (c^2 r^3)$ term is the GR correction beyond Newton.}",
            color=YELLOW, font_size=22,
        ).next_to(s3, DOWN, buff=0.3)
        self.play(Write(s3))
        self.play(Write(s3_note))
        self.wait(2.0)

        self.play(FadeOut(s1), FadeOut(s2), FadeOut(s3), FadeOut(s3_note))

        # Step 4: Orbit integral
        s4 = Tex(
            r"\textbf{Integrate over one orbit:}",
            font_size=28,
        ).to_edge(UP)
        s4_eq = MathTex(
            r"\Delta\phi \;=\; \frac{6\pi G M}{c^2\,a\,(1 - e^2)}",
            color=YELLOW, font_size=40,
        ).next_to(s4, DOWN, buff=0.8)
        self.play(Write(s4))
        self.play(Write(s4_eq))
        self.wait(1.5)

        # Step 5: Numerical for Mercury
        s5 = Tex(
            r"\textbf{For Mercury:} $a = 5.79 \times 10^{10}$ m, $e = 0.2056$, $GM_\odot$ known.",
            font_size=24,
        ).next_to(s4_eq, DOWN, buff=0.6)
        s5_eq = MathTex(
            r"\Delta\phi \;\approx\; 5.0 \times 10^{-7}\ \text{rad/orbit} "
            r"\;\times\; 414\ \text{orbits/century} \;\approx\; 43''/\text{century}",
            font_size=24,
        ).next_to(s5, DOWN, buff=0.4)
        s5_match = Tex(
            r"\textit{Matches Le Verrier's residual exactly. First confirmation of GR.}",
            color=YELLOW, font_size=26,
        ).next_to(s5_eq, DOWN, buff=0.6)
        self.play(Write(s5))
        self.play(Write(s5_eq))
        self.wait(0.8)
        self.play(Write(s5_match))
        self.wait(2.5)

        # Final bridge to GPS
        bridge = Tex(
            r"\textit{Same physics --- clocks deeper in a gravitational potential tick slower ---} \\ "
            r"\textit{drives the GPS satellite clock correction. See \S C.}",
            font_size=24, color=WHITE,
        ).to_edge(DOWN)
        self.play(Write(bridge))
        self.wait(2.5)
        self.play(FadeOut(s4), FadeOut(s4_eq), FadeOut(s5), FadeOut(s5_eq),
                  FadeOut(s5_match), FadeOut(bridge))
