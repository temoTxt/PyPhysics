"""
Bohr atom levels: standard derivation alongside the dual-framework redo
using Gill's positive-definite Hamiltonian K = H^2/(2mc^2) + mc^2/2.

Chapter reference:
  Roadmapping/History/03_old_quantum_theory_1900_1925.md#4-key-derivations-worth-animating

The narrative:
  1. Newton + Coulomb + L = n hbar -> r_n = n^2 a_0.
  2. Total energy E_n = -1/(2n^2) * (e^2 / 4 pi eps0 a_0).
  3. Side-by-side: the dual K classical Hamiltonian gives the same E_n at v/c << 1.
  4. Numerical agreement of order 10^-9 for hydrogen.

Render (from Roadmapping/Animations):
    uv run manim -pql --media_dir rendered manim_scenes/hist_bohr_proper_time.py HistBohrProperTime
    uv run manim -qh  --media_dir rendered manim_scenes/hist_bohr_proper_time.py HistBohrProperTime
"""

from manim import (
    BLUE, DOWN, FadeIn, FadeOut, GREEN, LEFT, MathTex, RIGHT, Scene, Tex, Transform,
    UP, VGroup, WHITE, Write, YELLOW,
)


class HistBohrProperTime(Scene):
    """Bohr atom: standard derivation plus dual-Hamiltonian redo, side by side."""

    def construct(self):
        # ───── Title ─────
        title = Tex(r"\textbf{The Bohr atom in standard and dual frameworks}", font_size=40)
        subtitle = Tex(r"Chapter 3, Old Quantum Theory (1900--1925)", font_size=24)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ───── Postulates ─────
        postulates = Tex(
            r"Bohr (1913) postulates: \\ "
            r"\quad 1. Circular orbits with $L = n\hbar$ (no radiation). \\ "
            r"\quad 2. Transitions: $h\nu = E_m - E_n$.",
            font_size=28,
        )
        self.play(Write(postulates))
        self.wait(1.5)
        self.play(FadeOut(postulates))

        # ───── Standard derivation block ─────
        std_header = Tex(r"\textbf{Standard derivation}", color=BLUE, font_size=32).to_edge(UP).shift(LEFT * 3)
        std_eqs = VGroup(
            MathTex(r"\frac{m_e v^2}{r} = \frac{e^2}{4\pi\varepsilon_0 r^2}", font_size=28),
            MathTex(r"m_e v r = n\hbar", font_size=28),
            MathTex(r"r_n = n^2 \, a_0,\ \ a_0 = \frac{4\pi\varepsilon_0 \hbar^2}{m_e e^2}", font_size=28),
            MathTex(r"\boxed{E_n^{\rm std} = -\frac{1}{2n^2}\,\frac{e^2}{4\pi\varepsilon_0 a_0}}", color=YELLOW,
                    font_size=30),
        ).arrange(DOWN, buff=0.3).next_to(std_header, DOWN, buff=0.4)

        self.play(Write(std_header))
        for eq in std_eqs:
            self.play(Write(eq))
            self.wait(0.4)
        self.wait(1.0)

        # ───── Dual-framework block (appears on the right) ─────
        dual_header = Tex(r"\textbf{Dual framework}", color=GREEN, font_size=32).to_edge(UP).shift(RIGHT * 3)
        dual_eqs = VGroup(
            MathTex(r"K = \frac{H^2}{2mc^2} + \frac{mc^2}{2}", font_size=28),
            MathTex(r"H = \frac{p^2}{2m_e} - \frac{e^2}{4\pi\varepsilon_0 r}", font_size=24),
            Tex(r"Stationary states of $K$ at $v/c \ll 1$:", font_size=22),
            MathTex(r"\boxed{E_n^{\rm dual} = E_n^{\rm std}\,\bigl[1 + O(v^4/c^4)\bigr]}", color=YELLOW,
                    font_size=28),
        ).arrange(DOWN, buff=0.3).next_to(dual_header, DOWN, buff=0.4)

        self.play(Write(dual_header))
        for eq in dual_eqs:
            self.play(Write(eq))
            self.wait(0.4)
        self.wait(1.5)

        # ───── Numerical comparison ─────
        # Clear the columns first
        self.play(
            FadeOut(std_header), FadeOut(std_eqs),
            FadeOut(dual_header), FadeOut(dual_eqs),
        )

        comparison_title = Tex(r"\textbf{Numerical comparison for hydrogen}", font_size=34).to_edge(UP)
        comparison = VGroup(
            Tex(r"Hydrogen Bohr velocity: $v/c \sim Z\alpha \approx 7 \times 10^{-3}$", font_size=26),
            Tex(r"Standard vs. dual relative difference: $\sim (v/c)^4 \sim 10^{-9}$", font_size=26),
            Tex(r"1913 spectroscopic precision: $\sim 10^{-4}$", font_size=26),
            Tex(r"Current hydrogen-spectroscopy precision: $\sim 10^{-14}$", font_size=26),
            Tex(r"At $10^{-14}$ precision, both framings are dominated by QED radiative corrections \\ "
                r"(Lamb shift, vacuum polarisation) computed equivalently in both.",
                font_size=22, color=WHITE),
        ).arrange(DOWN, buff=0.35).next_to(comparison_title, DOWN, buff=0.5)

        self.play(Write(comparison_title))
        for line in comparison:
            self.play(Write(line))
            self.wait(0.3)
        self.wait(2.5)

        self.play(FadeOut(comparison_title), FadeOut(comparison))

        # ───── Final boxed conclusion ─────
        conclusion = Tex(
            r"\textit{Bohr-atom predictions: identical in both framings at any practical precision. \\ "
            r"Dual framework's conceptual win: positive-definite $K$, no negative-energy modes,} \\ "
            r"\textit{no Klein--Gordon negative-probability problem at the relativistic level (see DRQM I, Ch 4).}",
            font_size=24, color=YELLOW,
        )
        self.play(Write(conclusion))
        self.wait(3.0)
        self.play(FadeOut(conclusion))
