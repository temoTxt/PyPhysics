"""
Compton scattering: standard derivation of Delta lambda = (h/mc)(1 - cos theta),
then the dual-framework redo with b = sqrt(c^2 + u^2). For u = 0 (electron at rest),
the dual formula reduces exactly to the standard. The 'null' is that the two
predictions coincide at any current experimental precision for free-electron scattering.

Chapter reference:
  Roadmapping/History/03_old_quantum_theory_1900_1925.md#4-key-derivations-worth-animating

Render (from Roadmapping/Animations):
    uv run manim -pql --media_dir rendered manim_scenes/hist_compton_null.py HistComptonNull
    uv run manim -qh  --media_dir rendered manim_scenes/hist_compton_null.py HistComptonNull
"""

from manim import (
    BLUE, DOWN, FadeIn, FadeOut, GREEN, LEFT, MathTex, ORIGIN, RED, RIGHT, Scene, Tex,
    Transform, UP, VGroup, WHITE, Write, YELLOW, Arrow, Dot,
)


class HistComptonNull(Scene):
    """Compton scattering: standard derivation + dual-framework null reproduction."""

    def construct(self):
        # ───── Title ─────
        title = Tex(r"\textbf{Compton scattering: standard vs dual}", font_size=40)
        subtitle = Tex(r"Chapter 3, Old Quantum Theory (1900--1925)", font_size=24)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ───── Kinematics setup diagram ─────
        # Photon in from left, electron at origin (dot), photon out at angle, electron recoil at angle.
        electron = Dot(ORIGIN, color=BLUE, radius=0.1)
        electron_label = MathTex(r"e^-", font_size=24).next_to(electron, DOWN, buff=0.2)

        photon_in = Arrow(start=[-4, 0, 0], end=[-0.3, 0, 0], color=YELLOW, buff=0)
        photon_in_label = MathTex(r"\gamma_0,\ \lambda_0", font_size=22).next_to(photon_in, UP, buff=0.15)

        # Outgoing photon at ~+45 deg
        photon_out = Arrow(start=[0.3, 0, 0], end=[2.8, 2.5, 0], color=YELLOW, buff=0)
        photon_out_label = MathTex(r"\gamma,\ \lambda,\ \theta", font_size=22).next_to(photon_out, UP, buff=0.15)

        # Electron recoil at ~-25 deg
        electron_recoil = Arrow(start=[0.3, 0, 0], end=[3.0, -1.4, 0], color=BLUE, buff=0)
        recoil_label = MathTex(r"p_e,\ \phi", font_size=22).next_to(electron_recoil, DOWN, buff=0.15)

        diagram = VGroup(electron, electron_label, photon_in, photon_in_label,
                          photon_out, photon_out_label, electron_recoil, recoil_label)
        self.play(Write(diagram))
        self.wait(0.5)

        # ───── Energy + momentum conservation ─────
        conservation = VGroup(
            MathTex(r"\frac{hc}{\lambda_0} + m_e c^2 \;=\; \frac{hc}{\lambda} + \sqrt{p_e^2 c^2 + m_e^2 c^4}",
                    font_size=24),
            MathTex(r"\frac{h}{\lambda_0} \;=\; \frac{h}{\lambda}\cos\theta + p_e \cos\phi", font_size=24),
            MathTex(r"0 \;=\; \frac{h}{\lambda}\sin\theta - p_e \sin\phi", font_size=24),
        ).arrange(DOWN, buff=0.2).to_corner(DOWN + LEFT)
        self.play(Write(conservation))
        self.wait(1.5)

        # ───── Standard Compton formula ─────
        std_formula = MathTex(
            r"\boxed{\;\Delta\lambda^{\rm std} \;=\; \lambda - \lambda_0 \;=\; \frac{h}{m_e c}\,(1 - \cos\theta)\;}",
            color=YELLOW, font_size=32,
        ).to_corner(DOWN + RIGHT)
        self.play(Write(std_formula))
        self.wait(2.0)

        # Clear diagram + conservation for the dual section
        self.play(FadeOut(diagram), FadeOut(conservation), FadeOut(std_formula))

        # ───── Dual-framework redo ─────
        dual_title = Tex(r"\textbf{Dual-framework redo}", color=GREEN, font_size=36).to_edge(UP)
        dual_step1 = Tex(
            r"In Gill's dual framework, the kinematic factor is set by the source velocity $u$:",
            font_size=24,
        ).next_to(dual_title, DOWN, buff=0.6)
        dual_step2 = MathTex(
            r"b \;=\; \sqrt{c^2 + u^2}, \qquad \mathbf{w}/c \;=\; \mathbf{u}/b",
            color=YELLOW, font_size=28,
        ).next_to(dual_step1, DOWN, buff=0.5)
        dual_step3 = Tex(
            r"For a target electron at rest, $u = 0 \;\Longrightarrow\; b = c$:",
            font_size=24,
        ).next_to(dual_step2, DOWN, buff=0.5)
        dual_formula = MathTex(
            r"\boxed{\;\Delta\lambda^{\rm dual} \;=\; \frac{h}{m_e c}\,(1 - \cos\theta) \;=\; \Delta\lambda^{\rm std}\;}",
            color=YELLOW, font_size=32,
        ).next_to(dual_step3, DOWN, buff=0.5)

        self.play(Write(dual_title))
        self.play(Write(dual_step1))
        self.wait(0.6)
        self.play(Write(dual_step2))
        self.wait(0.6)
        self.play(Write(dual_step3))
        self.wait(0.6)
        self.play(Write(dual_formula))
        self.wait(2.0)

        # ───── Bound-electron correction estimate ─────
        bound_note = Tex(
            r"For a bound (moving) electron at $u/c \sim Z\alpha \sim 10^{-2}$: \\ "
            r"$b/c - 1 \sim u^2/(2c^2) \sim 10^{-4}$. \\ "
            r"\textit{Below Compton's apparatus precision and below modern} \\ "
            r"\textit{atomic Compton-profile measurements.}",
            font_size=22, color=WHITE,
        ).to_edge(DOWN)
        self.play(Write(bound_note))
        self.wait(3.0)

        self.play(FadeOut(dual_title), FadeOut(dual_step1), FadeOut(dual_step2),
                  FadeOut(dual_step3), FadeOut(dual_formula), FadeOut(bound_note))
