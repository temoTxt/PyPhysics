"""
Maxwell paper Eqs. (10) and (11): proper-time boosts preserving the 4-velocity hyperboloid.

Derivation reference:
  Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-10--proper-time-boosts-of-position-velocity-acceleration

Eq. (10) velocity row:  u' = gamma (u - (v/c) b)
Eq. (11):                b' = gamma (b - u v/c)

Key invariant: b^2 - u^2 = c^2 is preserved by the boost.
Geometrically: (b, u) moves along the hyperboloid b^2 - u^2 = c^2 in (u, b) space.

Render:
    uv run manim -pql manim_scenes/maxwell_eq10_11_hyperboloid.py MaxwellEq10_11Hyperboloid
"""

from manim import *
import numpy as np


class MaxwellEq10_11Hyperboloid(Scene):
    """Animate (b, u) on the proper-velocity hyperboloid under a boost."""

    def construct(self):
        title = Tex(r"\textbf{The 4-velocity hyperboloid}", font_size=40)
        subtitle = Tex(
            r"Gill \& Zachary (2011), Eqs.\ (10)--(11): $(b, \mathbf{u})$ as a 4-vector",
            font_size=24,
        )
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ---- Step 1: the invariant ----
        step1 = Tex(r"\textbf{Step 1.}\ Recall $b = \sqrt{c^2 + \mathbf{u}^2}$, so:", font_size=28).to_edge(UP)
        eq1 = MathTex(
            r"b^2 - \mathbf{u}^2 \;=\; c^2",
            font_size=56,
            color=YELLOW,
        )
        self.play(Write(step1))
        self.play(Write(eq1))
        self.wait(1.5)
        sub_caption = Tex(
            r"\textit{The Minkowski length-squared of the 4-velocity $(b, \mathbf{u})$.}",
            font_size=24,
        ).next_to(eq1, DOWN, buff=0.6)
        self.play(Write(sub_caption))
        self.wait(1.5)
        self.play(FadeOut(eq1), FadeOut(sub_caption))

        # ---- Step 2: set up axes in (u, b) plane ----
        step2 = Tex(
            r"\textbf{Step 2.}\ Plot the hyperboloid in the $(u, b)$-plane (with $c=1$).",
            font_size=26,
        ).to_edge(UP)
        self.play(Transform(step1, step2))

        axes = Axes(
            x_range=[-3, 3, 1],
            y_range=[0, 4, 1],
            x_length=8,
            y_length=4,
            tips=False,
            axis_config={"include_numbers": False},
        ).shift(DOWN * 0.5)
        x_label = MathTex(r"u/c", font_size=28).next_to(axes.x_axis.get_end(), RIGHT, buff=0.2)
        y_label = MathTex(r"b/c", font_size=28).next_to(axes.y_axis.get_end(), UP, buff=0.2)

        # The hyperboloid b = sqrt(1 + u^2)  (in units c = 1)
        hyperboloid = axes.plot(
            lambda u: np.sqrt(1 + u**2),
            x_range=[-2.9, 2.9],
            color=BLUE,
        )
        hyperboloid_label = MathTex(r"b = \sqrt{c^2 + u^2}", font_size=28, color=BLUE).move_to(
            axes.c2p(2.0, 3.0)
        )

        self.play(Create(axes), Write(x_label), Write(y_label))
        self.play(Create(hyperboloid))
        self.play(Write(hyperboloid_label))
        self.wait(1.0)

        # ---- Step 3: place a point on the hyperboloid ----
        step3 = Tex(
            r"\textbf{Step 3.}\ A particle's state is a point $(u, b)$ on this curve.",
            font_size=24,
        ).to_edge(UP)
        self.play(Transform(step1, step3))

        u_initial = 0.5
        b_initial = np.sqrt(1 + u_initial**2)
        point_dot = Dot(axes.c2p(u_initial, b_initial), color=YELLOW, radius=0.1)
        point_label = MathTex(r"(u, b)", font_size=24, color=YELLOW).next_to(point_dot, UR, buff=0.1)
        self.play(FadeIn(point_dot), Write(point_label))
        self.wait(1.5)

        # ---- Step 4: apply the boost ----
        step4 = Tex(
            r"\textbf{Step 4.}\ Apply the boost (Eqs.\ 10, 11): $u' = \gamma(u - vb/c)$, $b' = \gamma(b - uv/c)$.",
            font_size=22,
        ).to_edge(UP)
        boost_label = MathTex(
            r"v/c = 0.3, \quad \gamma = 1/\sqrt{1 - v^2/c^2}",
            font_size=24,
        ).next_to(step4, DOWN, buff=0.3)
        self.play(Transform(step1, step4))
        self.play(Write(boost_label))

        # Animate the boost: parameterize from no boost to v/c = -0.7
        v_tracker = ValueTracker(0.0)

        def boosted_point():
            v_over_c = v_tracker.get_value()
            gam = 1.0 / np.sqrt(1.0 - v_over_c**2)
            u_b = gam * (u_initial - v_over_c * b_initial)
            b_b = gam * (b_initial - v_over_c * u_initial)
            return axes.c2p(u_b, b_b)

        # Replace point_dot with an updater
        new_dot = always_redraw(lambda: Dot(boosted_point(), color=YELLOW, radius=0.1))
        new_label = always_redraw(
            lambda: MathTex(r"(u', b')", font_size=24, color=YELLOW).next_to(new_dot, UR, buff=0.1)
        )
        self.play(FadeOut(point_dot), FadeOut(point_label))
        self.add(new_dot, new_label)
        self.wait(0.5)

        # Sweep the boost velocity (-0.7 .. +0.7)
        self.play(v_tracker.animate.set_value(-0.7), run_time=2.5)
        self.play(v_tracker.animate.set_value(0.7), run_time=4.0)
        self.play(v_tracker.animate.set_value(0.0), run_time=2.0)
        self.wait(0.5)

        # ---- Step 5: invariant preserved ----
        step5 = Tex(
            r"\textbf{Step 5.}\ At every boost: $(b')^2 - (u')^2 = c^2$. The point stays on the curve.",
            font_size=22,
        ).to_edge(UP)
        invariant_box = MathTex(
            r"\boxed{\;b^2 - \mathbf{u}^2 \;=\; (b')^2 - (\mathbf{u}')^2 \;=\; c^2\;}",
            font_size=36,
            color=YELLOW,
        ).next_to(boost_label, DOWN, buff=0.5)
        self.play(Transform(step1, step5))
        self.play(FadeOut(boost_label))
        self.play(Write(invariant_box))
        self.wait(2.0)

        verified = Tex(
            r"\textit{Verified by Wolfram MCP: $b'^2 - u'^2 - c^2$ FullSimplify $\to 0$.}",
            font_size=22,
        ).next_to(invariant_box, DOWN, buff=0.4)
        self.play(Write(verified))
        self.wait(3.0)
