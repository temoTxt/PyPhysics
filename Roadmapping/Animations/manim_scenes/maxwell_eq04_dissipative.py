"""
Maxwell paper Eq. (4): emergence of the dissipative term (u . a) / b^4 from curl-of-curl.

Derivation reference:
  Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-4--dual-wave-equations-with-dissipative-term

The key conceptual point: in the standard wave equation (1/c^2) d^2/dt^2 - laplacian = source,
nothing depends on time except the derivatives.  In the proper-time form, b = b(tau) is itself
time-dependent, which yields an extra first-tau-derivative term -- the source of radiation reaction.

Render:
    uv run manim -pql manim_scenes/maxwell_eq04_dissipative.py MaxwellEq04Dissipative
"""

from manim import *


class MaxwellEq04Dissipative(Scene):
    """The dissipative term -(u . a) / b^4 emerges from d/dtau (1/b)."""

    def construct(self):
        title = Tex(r"\textbf{Where does the dissipative term come from?}", font_size=40)
        subtitle = Tex(
            r"Gill \& Zachary (2011), Eq.\ (4) -- dual wave equations with dissipation",
            font_size=24,
        )
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.2)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ---- Step 1: take curl of Faraday ----
        step1 = Tex(r"\textbf{Step 1.}\ Take curl of $\nabla\times\mathbf{E} = -(1/b)\partial_\tau \mathbf{B}$:", font_size=28).to_edge(UP)
        eq1 = MathTex(
            r"\nabla\!\times\!\nabla\!\times\!\mathbf{E} \;=\; -\nabla\!\times\!\Big[\tfrac{1}{b}\,\partial_\tau \mathbf{B}\Big]",
            font_size=40,
        )
        self.play(Write(step1))
        self.play(Write(eq1))
        self.wait(1.5)

        # ---- Step 2: b commutes with grad (b depends on tau, not on x) ----
        step2 = Tex(
            r"\textbf{Step 2.}\ At the source worldline, $b = b(\tau)$ -- no spatial dependence.",
            font_size=26,
        ).to_edge(UP)
        eq2 = MathTex(
            r"\nabla\!\times\!\nabla\!\times\!\mathbf{E} \;=\; -\tfrac{1}{b}\,\partial_\tau (\nabla\!\times\!\mathbf{B})",
            font_size=40,
        )
        self.play(Transform(step1, step2))
        self.play(Transform(eq1, eq2))
        self.wait(1.5)

        # ---- Step 3: substitute Ampere ----
        step3 = Tex(
            r"\textbf{Step 3.}\ Substitute Amp\`ere: $\nabla\!\times\!\mathbf{B} = (1/b)\partial_\tau\mathbf{E} + (4\pi/b)\rho\mathbf{u}$.",
            font_size=24,
        ).to_edge(UP)
        eq3 = MathTex(
            r"\nabla\!\times\!\nabla\!\times\!\mathbf{E} \;=\; -\tfrac{1}{b}\,\partial_\tau\!\Big[\tfrac{1}{b}\partial_\tau\mathbf{E} + \tfrac{4\pi}{b}\rho\mathbf{u}\Big]",
            font_size=36,
        )
        self.play(Transform(step1, step3))
        self.play(Transform(eq1, eq3))
        self.wait(2.0)

        # ---- Step 4: the key derivative -- d/dtau(1/b) is NOT zero ----
        step4 = Tex(
            r"\textbf{Step 4.}\ The key derivative:",
            font_size=28,
        ).to_edge(UP)
        eq4a = MathTex(
            r"\partial_\tau\!\Big[\tfrac{1}{b}\partial_\tau\mathbf{E}\Big] \;=\; \tfrac{1}{b}\partial_\tau^2\mathbf{E} + (\partial_\tau \tfrac{1}{b})\,\partial_\tau\mathbf{E}",
            font_size=36,
        )
        eq4b = MathTex(
            r"\partial_\tau\!\!\Big[\tfrac{1}{b}\Big] \;=\; -\frac{\dot b}{b^2}, \qquad \dot b \;=\; \frac{\mathbf{u}\cdot\mathbf{a}}{b}",
            font_size=36,
            color=YELLOW,
        )
        eq4 = VGroup(eq4a, eq4b).arrange(DOWN, buff=0.5)
        self.play(Transform(step1, step4))
        self.play(Transform(eq1, eq4))
        self.wait(2.5)

        # ---- Step 5: combine, the dissipative term emerges ----
        step5 = Tex(
            r"\textbf{Step 5.}\ Therefore $\partial_\tau\!\Big[\tfrac{1}{b}\partial_\tau\mathbf{E}\Big] = \tfrac{1}{b}\partial_\tau^2\mathbf{E} - \tfrac{\mathbf{u}\cdot\mathbf{a}}{b^3}\partial_\tau\mathbf{E}$.",
            font_size=22,
        ).to_edge(UP)
        eq5 = MathTex(
            r"\boxed{\;\tfrac{1}{b^2}\partial_\tau^2\mathbf{E} \;-\; \frac{\mathbf{u}\cdot\mathbf{a}}{b^4}\,\partial_\tau\mathbf{E} \;-\; \nabla^2\mathbf{E} \;=\; \text{source}\;}",
            font_size=44,
            color=YELLOW,
        )
        self.play(Transform(step1, step5))
        self.play(Transform(eq1, eq5))
        self.wait(3.0)

        # ---- Closing: physical interpretation ----
        closing = Tex(
            r"\textit{The extra }$-(\mathbf{u}\cdot\mathbf{a})/b^4 \cdot \partial_\tau\mathbf{E}$\textit{ term has no counterpart in standard EM.}",
            font_size=24,
        ).to_edge(UP)
        interp1 = Tex(r"It vanishes when $\mathbf{a} = 0$ (no acceleration).", font_size=24)
        interp2 = Tex(r"When $\mathbf{a} \neq 0$, it acts like \textbf{damping} -- radiation reaction.", font_size=24)
        interp1.next_to(eq1, DOWN, buff=0.8)
        interp2.next_to(interp1, DOWN, buff=0.3)
        self.play(Transform(step1, closing))
        self.play(FadeIn(interp1))
        self.wait(1.5)
        self.play(FadeIn(interp2))
        self.wait(3.0)
