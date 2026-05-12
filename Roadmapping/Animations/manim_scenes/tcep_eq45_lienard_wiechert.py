"""
TCEP Eq. (3.45): the modified Lienard-Wiechert E-field with the three-term structure.

Derivation reference:
  Roadmapping/Equation_Verification/The_Classical_Electron_Problem.md#section-32--full-lienardwiechert-derivation

Geometry:
  - Source moves along a curve x'(tau').
  - Field point P (fixed) sees radiation from the retarded position x'(tau').
  - r = P - x'(tau'),  s = r - r.u/b,  r_u = r - (r/b) u

Three terms in E (paper Eq. 3.45):
  (1) velocity field:      q r_u (1 - u^2/b^2) / s^3
  (2) acceleration field:  q [r x (r_u x a)] / (b^2 s^3)
  (3) NEW dissipative:     q (u . a) [r x (u x r)] / (b^4 s^3)

The third term is the headline novelty: vanishes when a = 0, otherwise has a
longitudinal component (radiation reaction from the local clock).

Render:
    uv run manim -pql manim_scenes/tcep_eq45_lienard_wiechert.py TCEPLienardWiechert
"""

from manim import *
import numpy as np


class TCEPLienardWiechert(Scene):
    """Retarded-time geometry + the three-term E-field structure."""

    def construct(self):
        title = Tex(r"\textbf{Modified Li\'enard--Wiechert: three-term $\mathbf{E}$ field}", font_size=36)
        subtitle = Tex(
            r"Gill, Zachary, Lindesay (2001), \textit{The Classical Electron Problem}, Eq.\ (3.45)",
            font_size=22,
        )
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ---- Set up the spacetime geometry (2D x-vs-time-like) ----
        axes = Axes(
            x_range=[-5, 5, 1],
            y_range=[-1, 4, 1],
            x_length=10,
            y_length=5,
            tips=False,
            axis_config={"include_numbers": False},
        ).shift(DOWN * 0.5)
        x_label = MathTex(r"x", font_size=24).next_to(axes.x_axis.get_end(), RIGHT, buff=0.2)
        t_label = MathTex(r"t", font_size=24).next_to(axes.y_axis.get_end(), UP, buff=0.2)
        self.play(Create(axes), Write(x_label), Write(t_label))

        # ---- Source worldline (accelerated trajectory) ----
        def worldline(t):
            """Source position x'(t).  Slightly curved trajectory."""
            return -2.0 + 0.6 * t + 0.15 * t**2

        # Build worldline parametrically: time t runs along the y-axis, source position along x.
        worldline_pts = [axes.c2p(worldline(t), t) for t in np.linspace(0, 3.5, 50)]
        worldline_curve = VMobject(color=BLUE, stroke_width=4)
        worldline_curve.set_points_smoothly(worldline_pts)
        worldline_label = MathTex(r"\text{source worldline } \mathbf{x}'(\tau')", font_size=20, color=BLUE).move_to(
            axes.c2p(-2.5, 3.2)
        )
        self.play(Create(worldline_curve), Write(worldline_label))
        self.wait(0.5)

        # ---- Field point P (fixed in space, at large t) ----
        P_x, P_t = 3.0, 3.0
        P_pt = axes.c2p(P_x, P_t)
        P_dot = Dot(P_pt, color=YELLOW, radius=0.1)
        P_label = MathTex(r"P = (\mathbf{x}, \tau)", font_size=22, color=YELLOW).next_to(P_dot, RIGHT, buff=0.15)
        self.play(FadeIn(P_dot, scale=2.0), Write(P_label))
        self.wait(0.8)

        # ---- Retarded source position (light-cone intersection) ----
        # At retarded time tau', the source is at x'(tau') and the ray reaches (P_x, P_t) at speed c.
        # In our units, c = 1 along the +/-45 deg line.  Solve: P_t - tau' = |P_x - worldline(tau')|.
        # For simplicity pick a tau' that's roughly correct, then draw the ray.
        # We'll guess tau' ~ 0.7, then x' = -2 + 0.6*0.7 + 0.15*0.49 = -2 + 0.42 + 0.0735 = -1.5065
        # Check: P_x - x' = 3 - (-1.5) = 4.5; P_t - tau' = 3 - 0.7 = 2.3.  Not on light cone.
        # Need P_t - tau' = P_x - worldline(tau').  Solve numerically.
        def light_cone_eq(tp):
            return (P_t - tp) - (P_x - worldline(tp))
        # Crude bisection.  Bracket variable names deliberately avoid `a` (acceleration)
        # and `b` (collaborative speed of light), which are load-bearing physics symbols
        # elsewhere in this codebase.
        t_lo, t_hi = 0.0, 3.0
        for _ in range(40):
            mid = 0.5 * (t_lo + t_hi)
            if light_cone_eq(mid) > 0:
                t_lo = mid
            else:
                t_hi = mid
        tau_prime = 0.5 * (t_lo + t_hi)
        x_prime = worldline(tau_prime)

        ret_pt = axes.c2p(x_prime, tau_prime)
        ret_dot = Dot(ret_pt, color=GREEN, radius=0.09)
        ret_label = MathTex(r"\mathbf{x}'(\tau')", font_size=22, color=GREEN).next_to(ret_dot, DL, buff=0.15)
        self.play(FadeIn(ret_dot, scale=2.0), Write(ret_label))

        # Light ray from retarded position to P
        ray = Line(ret_pt, P_pt, color=WHITE, stroke_width=2)
        ray_label = MathTex(r"r = |\mathbf{x} - \mathbf{x}'|", font_size=20).move_to(
            (np.array(ret_pt) + np.array(P_pt)) / 2 + UP * 0.3
        )
        self.play(Create(ray), Write(ray_label))
        self.wait(0.8)

        # ---- Retarded-time integral relation ----
        light_cone_eq_box = MathTex(
            r"c(t - t') \;=\; \int_{\tau'}^{\tau} b(s)\,ds",
            font_size=28,
            color=GREEN,
        ).to_edge(DOWN, buff=0.5)
        self.play(Write(light_cone_eq_box))
        self.wait(2.0)

        # ---- Fade geometry, bring in formula ----
        self.play(
            FadeOut(axes), FadeOut(x_label), FadeOut(t_label),
            FadeOut(worldline_curve), FadeOut(worldline_label),
            FadeOut(P_dot), FadeOut(P_label),
            FadeOut(ret_dot), FadeOut(ret_label),
            FadeOut(ray), FadeOut(ray_label),
            FadeOut(light_cone_eq_box),
        )

        # ---- The three-term E-field formula ----
        intro = Tex(
            r"\textbf{After the (long) chain of kinematic identities (Eqs.\ 3.36, 3.40, 3.41), the field is:}",
            font_size=24,
        ).to_edge(UP)
        self.play(Write(intro))

        e_field = MathTex(
            r"\mathbf{E}(\mathbf{x}, \tau) \;=\;",
            r"\frac{q\,\mathbf{r_u}(1 - \mathbf{u}^2/b^2)}{s^3}",  # term 1
            r"\;+\;",
            r"\frac{q\,[\mathbf{r}\times(\mathbf{r_u}\times\mathbf{a})]}{b^2 s^3}",  # term 2
            r"\;+\;",
            r"\frac{q\,(\mathbf{u}\!\cdot\!\mathbf{a})\,[\mathbf{r}\times(\mathbf{u}\times\mathbf{r})]}{b^4 s^3}",  # term 3
            font_size=28,
        )
        self.play(Write(e_field))
        self.wait(1.5)

        # ---- Highlight each term ----
        term1_label = Tex(r"\textbf{Term 1:} velocity field", font_size=24, color=BLUE)
        term1_label.next_to(e_field[1], DOWN, buff=0.8)
        term1_box = SurroundingRectangle(e_field[1], color=BLUE, buff=0.1)
        self.play(Create(term1_box), Write(term1_label))
        self.wait(1.2)
        term1_note = Tex(
            r"\textit{Same form as Jackson's $E_v = e(n-\beta)(1-\beta^2)/(\kappa^3 R^2)$ with $\beta \to u/b$.}",
            font_size=20,
        ).next_to(term1_label, DOWN, buff=0.25)
        self.play(FadeIn(term1_note))
        self.wait(1.8)
        self.play(FadeOut(term1_box), FadeOut(term1_label), FadeOut(term1_note))

        term2_label = Tex(r"\textbf{Term 2:} radiation (acceleration) field", font_size=24, color=GREEN)
        term2_label.next_to(e_field[3], DOWN, buff=0.8)
        term2_box = SurroundingRectangle(e_field[3], color=GREEN, buff=0.1)
        self.play(Create(term2_box), Write(term2_label))
        self.wait(1.2)
        term2_note = Tex(
            r"\textit{Same form as Jackson's $E_a = (e/c)[n \times ((n-\beta)\times\dot\beta)]/(\kappa^3 R)$.}",
            font_size=20,
        ).next_to(term2_label, DOWN, buff=0.25)
        self.play(FadeIn(term2_note))
        self.wait(1.8)
        self.play(FadeOut(term2_box), FadeOut(term2_label), FadeOut(term2_note))

        term3_label = Tex(r"\textbf{Term 3:} \emph{new} dissipative term", font_size=24, color=RED)
        term3_label.next_to(e_field[5], DOWN, buff=0.8)
        term3_box = SurroundingRectangle(e_field[5], color=RED, buff=0.1)
        self.play(Create(term3_box), Write(term3_label))
        self.wait(1.2)
        term3_note1 = Tex(
            r"\textit{No counterpart in standard Li\'enard--Wiechert.}",
            font_size=20,
            color=RED,
        ).next_to(term3_label, DOWN, buff=0.25)
        term3_note2 = Tex(
            r"\textit{Vanishes when $\mathbf{a} = \mathbf{0}$; has a longitudinal component when $\mathbf{a} \neq \mathbf{0}$.}",
            font_size=20,
            color=RED,
        ).next_to(term3_note1, DOWN, buff=0.15)
        self.play(FadeIn(term3_note1))
        self.wait(1.5)
        self.play(FadeIn(term3_note2))
        self.wait(2.5)

        # ---- Closing: radiation reaction ----
        self.play(FadeOut(term3_box), FadeOut(term3_label), FadeOut(term3_note1), FadeOut(term3_note2), FadeOut(intro))
        closing = Tex(
            r"\textit{The dissipative term is the source of \textbf{radiation reaction}}\\"
            r"\textit{-- no Lorentz--Dirac equation, no self-energy divergence required.}",
            font_size=26,
            color=YELLOW,
        ).to_edge(UP)
        verified = Tex(
            r"\textit{Verified by Wolfram MCP: term~1 exactly matches Jackson (1D, $\mathbf{a}=\mathbf{0}$);}\\"
            r"\textit{term~3 vanishes identically when $\mathbf{u}\!\cdot\!\mathbf{a} = 0$.}",
            font_size=20,
        ).to_edge(DOWN)
        self.play(Write(closing))
        self.play(Write(verified))
        self.wait(3.0)
