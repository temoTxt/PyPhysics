"""
Maxwell synthesis: Ampère + displacement current -> wave equation -> propagation at c.

Chapter reference:
  Roadmapping/History/02_classical_synthesis_1860_1900.md#4-key-derivations-worth-animating

The four-step assembly that gives electromagnetism = optics:
  1. Ampère's law (with displacement current).
  2. Take the curl; substitute Faraday's law on the right.
  3. Apply the vector triple-product identity on the left (with div B = 0).
  4. Read off the wave equation; identify v = 1/sqrt(mu0 eps0) = c.

Render (from Roadmapping/Animations):
    uv run manim -pql --media_dir rendered manim_scenes/hist_maxwell_synthesis.py HistMaxwellSynthesis
    uv run manim -qh  --media_dir rendered manim_scenes/hist_maxwell_synthesis.py HistMaxwellSynthesis
"""

from manim import (
    DOWN, FadeIn, FadeOut, LEFT, MathTex, RIGHT, Scene, Tex, Transform, UP, VGroup, WHITE,
    Write, YELLOW, Create, BLUE,
)


class HistMaxwellSynthesis(Scene):
    """Maxwell synthesis: displacement current + curl chain -> wave equation, v = c."""

    def construct(self):
        # ───── Title card ─────
        title = Tex(r"\textbf{Maxwell's synthesis: electromagnetism = optics}", font_size=42)
        subtitle = Tex(
            r"Chapter 2, Classical Synthesis (1860--1900)",
            font_size=24,
        )
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ───── Step 1: Ampère's law with displacement current ─────
        step1_label = Tex(r"Step 1. Amp\`ere's law (with displacement current):",
                          font_size=26).to_edge(UP)
        step1_eq = MathTex(
            r"\nabla \times \mathbf{B} \;=\; \mu_0 \mathbf{J} \;+\; \mu_0 \varepsilon_0 \, \partial_t \mathbf{E}",
            font_size=40,
        )
        displacement_box = Tex(r"$\mu_0 \varepsilon_0 \, \partial_t \mathbf{E}$ is Maxwell's addition", color=BLUE,
                               font_size=22).next_to(step1_eq, DOWN, buff=0.5)
        self.play(Write(step1_label), Write(step1_eq))
        self.play(FadeIn(displacement_box))
        self.wait(1.5)
        self.play(FadeOut(step1_label), FadeOut(displacement_box))

        # ───── Step 2: In vacuum (J = 0); take the curl of both sides ─────
        step2_label = Tex(r"Step 2. In vacuum ($\mathbf{J}=0$); take curl of both sides:",
                          font_size=26).to_edge(UP)
        step2_eq = MathTex(
            r"\nabla \times (\nabla \times \mathbf{B}) \;=\; \mu_0 \varepsilon_0 \, \partial_t (\nabla \times \mathbf{E})",
            font_size=40,
        )
        self.play(Transform(step1_eq, step2_eq), Write(step2_label))
        self.wait(1.5)
        self.play(FadeOut(step2_label))

        # ───── Step 3: Substitute Faraday's law on the right, identity on the left ─────
        step3_label = Tex(
            r"Step 3. Use Faraday: $\nabla\times\mathbf{E} = -\partial_t \mathbf{B}$; "
            r"identity $\nabla\times(\nabla\times\mathbf{B}) = \nabla(\nabla\cdot\mathbf{B}) - \nabla^2 \mathbf{B}$ "
            r"with $\nabla\cdot\mathbf{B}=0$:",
            font_size=22,
        ).to_edge(UP)
        step3_eq = MathTex(
            r"-\nabla^2 \mathbf{B} \;=\; -\mu_0 \varepsilon_0 \, \partial_t^2 \mathbf{B}",
            font_size=40,
        )
        self.play(Transform(step1_eq, step3_eq), Write(step3_label))
        self.wait(1.5)
        self.play(FadeOut(step3_label))

        # ───── Step 4: Wave equation, identify v = 1/sqrt(mu0 eps0) ─────
        step4_label = Tex(r"Step 4. Rearrange: a transverse wave equation.",
                          font_size=26).to_edge(UP)
        wave_eq = MathTex(
            r"\left(\nabla^2 - \mu_0 \varepsilon_0 \, \partial_t^2\right) \mathbf{B} \;=\; 0",
            font_size=40,
        )
        self.play(Transform(step1_eq, wave_eq), Write(step4_label))
        self.wait(1.5)
        self.play(FadeOut(step4_label), FadeOut(step1_eq))

        # ───── Final: speed identification + dual-framework note ─────
        final = MathTex(
            r"\boxed{\;v \;=\; \frac{1}{\sqrt{\mu_0 \varepsilon_0}} \;=\; c\;}",
            color=YELLOW, font_size=56,
        )
        self.play(Write(final))
        self.wait(1.0)

        identification = Tex(
            r"Weber \& Kohlrausch (1856) had measured $1/\sqrt{\mu_0 \varepsilon_0} \approx 3 \times 10^{10}$ cm/s. \\ "
            r"\textit{Light is an electromagnetic wave.}",
            font_size=26,
        ).next_to(final, DOWN, buff=0.5)
        self.play(Write(identification))
        self.wait(2.0)

        # ───── Dual-framework note ─────
        footnote = Tex(
            r"Dual framework: substitute $\partial_t = (b/c)\,\partial_\tau$ throughout "
            r"(Gill \& Zachary 2011, Eq.\ 2). \\ "
            r"For sources at rest in the lab frame, $u=0$, $b=c$, and the dual wave equation "
            r"is \textit{identical} to the standard one. \\ "
            r"Divergence appears only for accelerating or relativistic sources --- not in vacuum propagation.",
            font_size=22, color=WHITE,
        ).next_to(identification, DOWN, buff=0.6)
        self.play(Write(footnote))
        self.wait(3.0)

        self.play(FadeOut(final), FadeOut(identification), FadeOut(footnote))
