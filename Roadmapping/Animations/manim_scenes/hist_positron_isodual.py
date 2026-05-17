"""
Antimatter via two organisations: standard Dirac negative-energy sea + hole theory,
vs the dual Dirac + Santilli-isodual companion equation.

Chapter reference:
  Roadmapping/History/04_quantum_mechanics_1925_1948.md#4-key-derivations-worth-animating

The two organisations agree on every observable prediction for the positron --
cross-sections, decay widths, g-factor. They differ in mathematical *organisation*:
where standard Dirac requires a Pauli-exclusion-filled sea + hole reinterpretation,
the dual framework has an independent isodual companion equation whose solutions
*are* the positrons.

Render (from Roadmapping/Animations):
    uv run manim -pql --media_dir rendered manim_scenes/hist_positron_isodual.py HistPositronIsodual
    uv run manim -qh  --media_dir rendered manim_scenes/hist_positron_isodual.py HistPositronIsodual
"""

from manim import (
    BLUE, DOWN, FadeIn, FadeOut, GREEN, LEFT, Line, MathTex, ORIGIN, RED, RIGHT,
    Scene, Square, Tex, Transform, UP, VGroup, WHITE, Write, YELLOW, Dot,
)


class HistPositronIsodual(Scene):
    """Standard Dirac hole theory vs dual Dirac + isodual companion."""

    def construct(self):
        # ───── Title ─────
        title = Tex(r"\textbf{Antimatter: two mathematical organisations}", font_size=36)
        subtitle = Tex(r"Chapter 4, Quantum Mechanics (1925--1948)", font_size=24)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ───── Standard Dirac side ─────
        std_header = Tex(r"\textbf{Standard Dirac (1928--1931)}", color=BLUE, font_size=30).to_edge(UP).shift(LEFT * 3.2)
        std_eq = MathTex(r"(i\gamma^\mu \partial_\mu - mc/\hbar)\psi = 0", font_size=24).next_to(std_header, DOWN, buff=0.4)
        std_modes = Tex(
            r"Solutions in two continua: \\ "
            r"\quad $E > +m c^2$: electrons \\ "
            r"\quad $E < -m c^2$: \textit{negative-energy continuum}",
            font_size=22,
        ).next_to(std_eq, DOWN, buff=0.4)
        std_hole = Tex(
            r"Dirac (1930--31): \\ "
            r"\textit{negative-energy sea filled by Pauli exclusion}; \\ "
            r"a hole in the sea appears as a particle of \\ "
            r"electron mass, opposite charge. \\ "
            r"Anderson 1932: positron observed.",
            font_size=22, color=BLUE,
        ).next_to(std_modes, DOWN, buff=0.4)

        self.play(Write(std_header), Write(std_eq))
        self.play(Write(std_modes))
        self.wait(0.5)
        self.play(Write(std_hole))
        self.wait(2.0)

        # ───── Dual Dirac + isodual side ─────
        dual_header = Tex(r"\textbf{Dual Dirac + isodual}", color=GREEN, font_size=30).to_edge(UP).shift(RIGHT * 3.2)
        dual_eq = MathTex(
            r"K_D = \frac{(c\boldsymbol{\alpha}\cdot\mathbf{p} + \beta mc^2)^2}{2mc^2} + \frac{mc^2}{2}",
            font_size=22,
        ).next_to(dual_header, DOWN, buff=0.4)
        dual_pos = Tex(
            r"\textit{Positive-definite}; bounded below by $mc^2/2$. \\ "
            r"No negative-energy continuum at all.",
            font_size=22,
        ).next_to(dual_eq, DOWN, buff=0.4)
        dual_isodual = Tex(
            r"\textit{Santilli isodual} companion equation: \\ "
            r"reverse the metric sign, energy sign, time direction. \\ "
            r"$K_D^* = -\,K_D$ on the isodual sheet. \\ "
            r"Solutions of $K_D^*$ \textit{are} the positrons.",
            font_size=22, color=GREEN,
        ).next_to(dual_pos, DOWN, buff=0.4)

        self.play(Write(dual_header), Write(dual_eq))
        self.play(Write(dual_pos))
        self.wait(0.5)
        self.play(Write(dual_isodual))
        self.wait(2.5)

        # Clear both columns
        self.play(
            FadeOut(std_header), FadeOut(std_eq), FadeOut(std_modes), FadeOut(std_hole),
            FadeOut(dual_header), FadeOut(dual_eq), FadeOut(dual_pos), FadeOut(dual_isodual),
        )

        # ───── Same predictions ─────
        same = Tex(r"\textbf{Both organisations give the same observables:}",
                   color=YELLOW, font_size=32).to_edge(UP)
        same_list = VGroup(
            Tex(r"\quad\textbullet\ positron mass = electron mass", font_size=24),
            Tex(r"\quad\textbullet\ positron charge $= +|e|$", font_size=24),
            Tex(r"\quad\textbullet\ $e^+ e^- \to 2\gamma$ annihilation cross-section", font_size=24),
            Tex(r"\quad\textbullet\ positron g-factor $g_{e^+} = -g_{e^-}$", font_size=24),
            Tex(r"\quad\textbullet\ Bhabha scattering, pair production", font_size=24),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.25).next_to(same, DOWN, buff=0.6)

        self.play(Write(same))
        for item in same_list:
            self.play(Write(item))
            self.wait(0.2)
        self.wait(1.5)

        # ───── Closing ─────
        difference = Tex(
            r"\textit{The two organisations differ in mathematical structure, not in predictions.} \\ "
            r"\textit{The dual + isodual organisation avoids the Pauli-exclusion-filled-sea} \\ "
            r"\textit{interpretation in favour of a separate (sign-reversed) companion equation.}",
            color=YELLOW, font_size=24,
        ).next_to(same_list, DOWN, buff=0.6)
        self.play(Write(difference))
        self.wait(3.0)
        self.play(FadeOut(same), FadeOut(same_list), FadeOut(difference))
