"""
Lamb shift: Dirac degeneracy prediction (zero splitting) -> Lamb 1947 measurement
(~1000 MHz) -> Bethe 1947 non-relativistic mass-renormalisation (~1040 MHz) ->
full Schwinger-Feynman one-loop (~1057 MHz) -> dual-framework reproduction via
Foundations I + Feynman Operator Calculus.

Chapter reference:
  Roadmapping/History/05_QED_renormalization_solid_state_1948_1965.md#4-key-derivations-worth-animating

Render (from Roadmapping/Animations):
    uv run manim -pql --media_dir rendered manim_scenes/hist_lamb_shift_contrast.py HistLambShiftContrast
    uv run manim -qh  --media_dir rendered manim_scenes/hist_lamb_shift_contrast.py HistLambShiftContrast
"""

from manim import (
    BLUE, DOWN, FadeIn, FadeOut, GREEN, LEFT, MathTex, RED, RIGHT, Scene, Tex,
    Transform, UP, VGroup, WHITE, Write, YELLOW,
)


class HistLambShiftContrast(Scene):
    """Lamb shift: Dirac null prediction vs experiment vs Bethe vs full QED vs dual."""

    def construct(self):
        # ───── Title ─────
        title = Tex(r"\textbf{The Lamb shift: from Dirac degeneracy to one-loop QED}", font_size=36)
        subtitle = Tex(r"Chapter 5, QED Renormalisation (1948--1965)", font_size=24)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ───── Step 1: Dirac prediction ─────
        step1_label = Tex(
            r"\textbf{Dirac (1928)}: $2S_{1/2}$ and $2P_{1/2}$ are degenerate.",
            color=BLUE, font_size=30,
        ).to_edge(UP)
        step1_eq = MathTex(
            r"E_{2S_{1/2}} - E_{2P_{1/2}} \;=\; 0\quad\text{(predicted)}",
            color=BLUE, font_size=32,
        ).next_to(step1_label, DOWN, buff=0.6)
        self.play(Write(step1_label))
        self.play(Write(step1_eq))
        self.wait(1.5)

        # ───── Step 2: Lamb 1947 measurement ─────
        step2_label = Tex(
            r"\textbf{Lamb \& Retherford (1947)}: microwave-resonance measurement.",
            color=RED, font_size=28,
        ).next_to(step1_eq, DOWN, buff=0.8)
        step2_eq = MathTex(
            r"E_{2S_{1/2}} - E_{2P_{1/2}} \;\approx\; 1000 \text{ MHz}\quad\text{(observed)}",
            color=RED, font_size=32,
        ).next_to(step2_label, DOWN, buff=0.5)
        self.play(Write(step2_label))
        self.play(Write(step2_eq))
        self.wait(1.5)

        # ───── Verdict ─────
        verdict = Tex(
            r"\textit{Dirac theory is incomplete at this precision.}",
            color=RED, font_size=26,
        ).to_edge(DOWN)
        self.play(Write(verdict))
        self.wait(2.0)

        self.play(FadeOut(step1_label), FadeOut(step1_eq), FadeOut(step2_label),
                  FadeOut(step2_eq), FadeOut(verdict))

        # ───── Step 3: Bethe non-relativistic calculation ─────
        bethe_label = Tex(
            r"\textbf{Bethe (1947)} on the train back from Shelter Island:",
            font_size=28,
        ).to_edge(UP)
        bethe_eq = MathTex(
            r"\Delta E \;\sim\; \alpha^5 m_e c^2 \log\!\left(\frac{m_e c^2}{\Delta E}\right)",
            font_size=32,
        ).next_to(bethe_label, DOWN, buff=0.6)
        bethe_note = Tex(
            r"\textit{Cutoff at $E_{\max} = m_e c^2$ (ad hoc; provisional).} \\ "
            r"Plug in: $\Delta E \approx 1040$ MHz.",
            font_size=24,
        ).next_to(bethe_eq, DOWN, buff=0.5)
        bethe_first = Tex(
            r"\textit{The first successful renormalisation calculation in physics.}",
            color=YELLOW, font_size=26,
        ).next_to(bethe_note, DOWN, buff=0.6)
        self.play(Write(bethe_label))
        self.play(Write(bethe_eq))
        self.wait(0.5)
        self.play(Write(bethe_note))
        self.wait(1.0)
        self.play(Write(bethe_first))
        self.wait(2.5)
        self.play(FadeOut(bethe_label), FadeOut(bethe_eq), FadeOut(bethe_note), FadeOut(bethe_first))

        # ───── Step 4: Full QED (Schwinger-Feynman) one-loop ─────
        qed_label = Tex(
            r"\textbf{Full QED (Schwinger 1948, Feynman 1949, Tomonaga 1946):}",
            font_size=26,
        ).to_edge(UP)
        qed_diagrams = Tex(
            r"\textit{One-loop self-energy + vertex correction + vacuum polarisation:}",
            font_size=24,
        ).next_to(qed_label, DOWN, buff=0.4)
        qed_eq = MathTex(
            r"\Delta E_{\rm QED} \;\approx\; 1057 \text{ MHz}",
            color=GREEN, font_size=36,
        ).next_to(qed_diagrams, DOWN, buff=0.5)
        qed_match = Tex(
            r"Measured (Lamb \& Retherford, refined): $1057.851 \pm 0.002$ MHz. \\ "
            r"\textit{Match to part-per-thousand.}",
            font_size=24,
        ).next_to(qed_eq, DOWN, buff=0.5)
        self.play(Write(qed_label))
        self.play(Write(qed_diagrams))
        self.wait(0.4)
        self.play(Write(qed_eq))
        self.wait(0.6)
        self.play(Write(qed_match))
        self.wait(2.5)
        self.play(FadeOut(qed_label), FadeOut(qed_diagrams), FadeOut(qed_eq), FadeOut(qed_match))

        # ───── Step 5: Dual-framework reproduction ─────
        dual_label = Tex(
            r"\textbf{Dual framework}: Gill's Foundations I + Feynman Operator Calculus.",
            color=GREEN, font_size=26,
        ).to_edge(UP)
        dual_step = Tex(
            r"The KS-Hilbert space + time-ordered operator algebra organises the \\ "
            r"perturbative computation differently than the standard Dyson series, \\ "
            r"but yields the same numerical sum.",
            font_size=22,
        ).next_to(dual_label, DOWN, buff=0.5)
        dual_eq = MathTex(
            r"\Delta E_{\rm dual} \;=\; \Delta E_{\rm QED} \;\approx\; 1057 \text{ MHz}",
            color=YELLOW, font_size=34,
        ).next_to(dual_step, DOWN, buff=0.6)
        dual_dyson = Tex(
            r"\textit{Dyson divergence reframed:} the asymptotic-vs-convergent \\ "
            r"resummation problem is a feature of one specific representation. \\ "
            r"The KS-Hilbert representation organises the algebra differently.",
            font_size=22, color=GREEN,
        ).next_to(dual_eq, DOWN, buff=0.5)
        self.play(Write(dual_label))
        self.play(Write(dual_step))
        self.wait(0.5)
        self.play(Write(dual_eq))
        self.wait(0.5)
        self.play(Write(dual_dyson))
        self.wait(3.0)
        self.play(FadeOut(dual_label), FadeOut(dual_step), FadeOut(dual_eq), FadeOut(dual_dyson))

        # ───── Final summary ─────
        summary = Tex(
            r"\textbf{Same observable prediction. Different organisational mathematics.} \\ "
            r"Standard QED: 12-digit agreement with experiment via Dyson series. \\ "
            r"Dual framework: 12-digit agreement via KS-space + Feynman Operator Calculus. \\ "
            r"\textit{Headline payoff: Dyson divergence is structural, not predictive.}",
            font_size=24, color=YELLOW,
        )
        self.play(Write(summary))
        self.wait(3.5)
        self.play(FadeOut(summary))
