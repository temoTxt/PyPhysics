"""
Faraday's induction law visualised: a conducting loop translates through a non-uniform
B-field; the flux through the loop changes; the induced EMF drives a current whose own
magnetic moment opposes the change (Lenz's law).

Chapter reference:
  Roadmapping/History/01_early_electromagnetism_1800_1860.md#4-key-derivations-worth-animating

This is a *pedagogical* scene — no proper-time content. Chapter 1 is a null comparison
between standard and dual electromagnetism. The standard form

    EMF = -dPhi/dt

is the load-bearing equation. The animation:
  1. Shows a non-uniform B-field (denser arrows on the right).
  2. Translates a rectangular loop into the region of stronger field.
  3. Tracks Phi(t) as the loop crosses field-density gradient.
  4. Shows the induced current direction (CW from above, by right-hand-rule) and the
     induced field opposing the increase — Lenz's law.
  5. Closes on the boxed result EMF = -dPhi/dt with the dual-framework footnote:
     in the lab-frame, b ≈ c, so dPhi/dt = dPhi/dtau identically.

Render (from Roadmapping/Animations):
    uv run manim -pql --media_dir rendered manim_scenes/hist_faraday_induction.py HistFaradayInduction
    uv run manim -qh  --media_dir rendered manim_scenes/hist_faraday_induction.py HistFaradayInduction
"""

import numpy as np
from manim import (
    BLUE, DOWN, FadeIn, FadeOut, GREEN, GrowArrow, LEFT, MathTex, ORIGIN, RED, RIGHT,
    Rectangle, Scene, Tex, UP, VGroup, Vector, WHITE, Write, YELLOW, ArcBetweenPoints,
    Transform, np,
)


class HistFaradayInduction(Scene):
    """Faraday induction: loop translates through non-uniform B, Lenz sign visualised."""

    def construct(self):
        # ───── Title card ─────
        title = Tex(r"\textbf{Faraday's law of induction}", font_size=44)
        subtitle = Tex(
            r"Chapter 1, Early Electromagnetism (1800--1860) --- pedagogical",
            font_size=24,
        )
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ───── Set up a non-uniform magnetic field (denser → right) ─────
        # We render B as a grid of small "into-the-page" markers; density rises with x.
        field_markers = VGroup()
        for x in np.linspace(-5.5, 5.5, 12):
            density = max(1, int(2 + (x + 5.5) / 1.5))  # 1 → ~8 markers per column
            for y in np.linspace(-2.5, 2.5, density):
                m = MathTex(r"\otimes", color=BLUE, font_size=22).move_to([x, y, 0])
                field_markers.add(m)
        field_label = Tex(r"$\mathbf{B}$ into the page; density increases to the right",
                          font_size=22).to_corner(UP + LEFT)
        self.play(FadeIn(field_markers), Write(field_label))
        self.wait(0.5)

        # ───── Conducting loop starts on the left ─────
        loop = Rectangle(width=2.0, height=1.6, color=YELLOW, stroke_width=4)
        loop.move_to(LEFT * 4)
        loop_label = MathTex(r"C", color=YELLOW, font_size=28).next_to(loop, UP, buff=0.15)
        self.play(FadeIn(loop), Write(loop_label))
        self.wait(0.4)

        # ───── Flux equation, animated as the loop moves ─────
        flux_eq = MathTex(r"\Phi(t) \;=\; \int_S \mathbf{B}\cdot d\mathbf{A}",
                          font_size=32).to_corner(DOWN + LEFT)
        self.play(Write(flux_eq))
        self.wait(0.5)

        # Translate the loop to the right; flux through it grows.
        self.play(loop.animate.move_to(RIGHT * 2.5),
                  loop_label.animate.next_to(loop.copy().move_to(RIGHT * 2.5), UP, buff=0.15),
                  run_time=3)
        self.wait(0.5)

        # ───── Show induced EMF and current direction ─────
        emf_eq = MathTex(r"\mathcal{E} \;=\; \oint_C \mathbf{E}\cdot d\boldsymbol{\ell} \;=\; -\,\frac{d\Phi}{dt}",
                         font_size=32).to_corner(DOWN + RIGHT)
        self.play(Write(emf_eq))
        self.wait(0.5)

        # Induced current arrow — CCW viewed from front (B into page, increasing flux into-page
        # → induced current generates B out of page → CCW when viewed from outside the page).
        # Draw four short arrows along the loop perimeter in the CCW sense.
        loop_center = loop.get_center()
        arrow_specs = [
            (loop_center + [-0.9, -0.7, 0], loop_center + [0.9, -0.7, 0]),       # bottom: → right
            (loop_center + [0.9, -0.7, 0], loop_center + [0.9, 0.7, 0]),         # right: ↑
            (loop_center + [0.9, 0.7, 0], loop_center + [-0.9, 0.7, 0]),         # top: ← left
            (loop_center + [-0.9, 0.7, 0], loop_center + [-0.9, -0.7, 0]),       # left: ↓
        ]
        current_arrows = VGroup(*[
            Vector(end - start, color=RED).shift(start) for start, end in arrow_specs
        ])
        i_label = MathTex(r"I_{\rm ind}", color=RED, font_size=26).next_to(loop, RIGHT, buff=0.4)
        self.play(*[GrowArrow(a) for a in current_arrows], Write(i_label))
        self.wait(0.5)

        # Show the induced field as a dot in the centre of the loop (out of page).
        b_ind = MathTex(r"\odot \, \mathbf{B}_{\rm ind}", color=GREEN, font_size=28).move_to(loop_center)
        self.play(FadeIn(b_ind))
        self.wait(0.5)

        # Lenz callout
        lenz = Tex(r"Lenz: $\mathbf{B}_{\rm ind}$ opposes the increase in $\Phi$",
                   color=GREEN, font_size=28).next_to(loop, DOWN, buff=0.8)
        self.play(Write(lenz))
        self.wait(1.5)

        # ───── Final boxed result + dual-framework footnote ─────
        self.play(
            FadeOut(field_markers), FadeOut(loop), FadeOut(loop_label),
            FadeOut(current_arrows), FadeOut(i_label), FadeOut(b_ind),
            FadeOut(lenz), FadeOut(flux_eq), FadeOut(emf_eq), FadeOut(field_label),
        )
        final = MathTex(r"\boxed{\;\mathcal{E} \;=\; -\,\frac{d\Phi}{dt}\;}", color=YELLOW, font_size=56)
        self.play(Write(final))
        self.wait(1.2)

        footnote = Tex(
            r"Dual framework: $\partial_t = (b/c)\,\partial_\tau$ "
            r"with $b = \sqrt{c^2 + u^2}$. \\ "
            r"For sources at rest in the lab frame ($u=0$), $b = c$ exactly, and "
            r"$\frac{d\Phi}{dt} = \frac{d\Phi}{d\tau}$. \\ "
            r"\textit{Standard and dual induction laws are identical in this regime.}",
            font_size=24, color=WHITE,
        ).next_to(final, DOWN, buff=0.6)
        self.play(Write(footnote))
        self.wait(2.5)

        self.play(FadeOut(final), FadeOut(footnote))
