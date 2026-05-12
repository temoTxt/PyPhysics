"""
Maxwell paper Eqs. (1) and (2): velocity and time-derivative duality.

Derivation reference:
  Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-1--velocity-duality

Eq. (1): w/c = u/b, with b = sqrt(c^2 + u^2).
Eq. (2): (1/c) d/dt = (1/b) d/dtau.

Render (from Roadmapping/Animations):
    uv run manim -pql manim_scenes/maxwell_eq01_02_duality.py MaxwellEq01_02Duality
"""

from manim import *


class MaxwellEq01_02Duality(Scene):
    """Show the algebraic seed of the dual theory: w/c = u/b and the time-derivative duality."""

    def construct(self):
        title = Tex(r"\textbf{Velocity and time-derivative duality}", font_size=44)
        subtitle = Tex(
            r"Gill \& Zachary (2011), \textit{Two Mathematically Equivalent Versions of Maxwell's Equations}, Eqs.\ (1)--(2)",
            font_size=24,
        )
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ---- Step 1: define proper time ----
        step1_label = Tex(r"\textbf{Step 1.}\ Proper time from the metric:", font_size=32).to_edge(UP)
        eq1 = MathTex(
            r"d\tau^2 \;=\; dt^2 - \frac{1}{c^2}\,d\mathbf{x}^2",
            font_size=44,
        )
        self.play(Write(step1_label))
        self.play(Write(eq1))
        self.wait(1.2)

        # ---- Step 2: introduce u = dx/dtau ----
        step2_label = Tex(
            r"\textbf{Step 2.}\ Define proper-time velocity $\mathbf{u} = d\mathbf{x}/d\tau$.", font_size=32
        ).to_edge(UP)
        eq2 = MathTex(
            r"dt^2 \;=\; d\tau^2 + \frac{1}{c^2} d\mathbf{x}^2 \;=\; d\tau^2 \left(1 + \frac{\mathbf{u}^2}{c^2}\right)",
            font_size=40,
        )
        self.play(Transform(step1_label, step2_label))
        self.play(Transform(eq1, eq2))
        self.wait(1.5)

        # ---- Step 3: introduce b ----
        step3_label = Tex(
            r"\textbf{Step 3.}\ Let $b = \sqrt{c^2 + \mathbf{u}^2}$ (collaborative speed of light).",
            font_size=30,
        ).to_edge(UP)
        eq3 = MathTex(
            r"\frac{dt}{d\tau} \;=\; \sqrt{1 + \frac{\mathbf{u}^2}{c^2}} \;=\; \frac{b}{c}",
            font_size=44,
        )
        self.play(Transform(step1_label, step3_label))
        self.play(Transform(eq1, eq3))
        self.wait(1.5)

        # ---- Step 4: derive w/c = u/b ----
        step4_label = Tex(
            r"\textbf{Step 4.}\ Then $\mathbf{u} = (dt/d\tau)\mathbf{w} = (b/c)\,\mathbf{w}$, so:",
            font_size=30,
        ).to_edge(UP)
        eq4 = MathTex(
            r"\boxed{\;\frac{\mathbf{w}}{c} \;=\; \frac{\mathbf{u}}{b}\;}",
            font_size=56,
            color=YELLOW,
        )
        self.play(Transform(step1_label, step4_label))
        self.play(Transform(eq1, eq4))
        self.wait(2.0)

        # ---- Step 5: time-derivative duality ----
        step5_label = Tex(
            r"\textbf{Step 5.}\ Chain rule: for any field $F(\mathbf{x},t)$,",
            font_size=30,
        ).to_edge(UP)
        eq5 = MathTex(
            r"\frac{\partial F}{\partial t} \;=\; \frac{d\tau}{dt}\,\frac{\partial F}{\partial \tau} \;=\; \frac{c}{b}\,\frac{\partial F}{\partial \tau}",
            font_size=40,
        )
        self.play(Transform(step1_label, step5_label))
        self.play(Transform(eq1, eq5))
        self.wait(1.5)

        # ---- Step 6: divide by c -> Eq (2) ----
        step6_label = Tex(
            r"\textbf{Step 6.}\ Divide both sides by $c$:",
            font_size=30,
        ).to_edge(UP)
        eq6 = MathTex(
            r"\boxed{\;\frac{1}{c}\frac{\partial}{\partial t} \;=\; \frac{1}{b}\frac{\partial}{\partial \tau}\;}",
            font_size=56,
            color=YELLOW,
        )
        self.play(Transform(step1_label, step6_label))
        self.play(Transform(eq1, eq6))
        self.wait(2.5)

        # ---- Closing: pair Eq (1) and Eq (2) ----
        closing_label = Tex(r"\textbf{The two duality identities, side by side:}", font_size=30).to_edge(UP)
        eq1_final = MathTex(r"\frac{\mathbf{w}}{c} \;=\; \frac{\mathbf{u}}{b}", font_size=44, color=YELLOW)
        eq2_final = MathTex(
            r"\frac{1}{c}\frac{\partial}{\partial t} \;=\; \frac{1}{b}\frac{\partial}{\partial \tau}",
            font_size=44,
            color=YELLOW,
        )
        eq1_final.shift(LEFT * 3)
        eq2_final.shift(RIGHT * 3)
        caption = Tex(
            r"\textit{Verified by Wolfram MCP: FullSimplify returns 0.}",
            font_size=24,
        ).shift(DOWN * 1.5)
        self.play(Transform(step1_label, closing_label))
        self.play(Transform(eq1, eq1_final), FadeIn(eq2_final))
        self.play(Write(caption))
        self.wait(3.0)
