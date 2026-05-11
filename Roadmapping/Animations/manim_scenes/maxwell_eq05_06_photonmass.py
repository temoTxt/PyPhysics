"""
Maxwell paper Eqs. (5) and (6): scale transformation and effective photon mass.

Derivation reference:
  Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-5--kleingordonlike-form-after-rescaling

Eq. (5): substitution psi = (b/c)^{1/2} Psi_new converts the damped wave equation into
         a Klein-Gordon form with an effective potential coefficient.
Eq. (6): the effective "photon mass" mu, expressed in terms of u, u_dot, u_ddot.

Render:
    uv run manim -pql manim_scenes/maxwell_eq05_06_photonmass.py MaxwellEq05_06PhotonMass
"""

from manim import *


class MaxwellEq05_06PhotonMass(Scene):
    """Liouville substitution -> effective mass term -> photon mass formula."""

    def construct(self):
        title = Tex(r"\textbf{From dissipation to effective photon mass}", font_size=40)
        subtitle = Tex(r"Gill \& Zachary (2011), Eqs.\ (5) and (6)", font_size=24)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # ---- Step 1: damped wave equation (recap from Eq 4) ----
        step1 = Tex(r"\textbf{Step 1.}\ Recall the damped wave equation:", font_size=28).to_edge(UP)
        eq1 = MathTex(
            r"\tfrac{1}{b^2}\partial_\tau^2 \psi \;-\; \tfrac{\dot b}{b^3}\,\partial_\tau \psi \;-\; \nabla^2 \psi \;=\; \text{source}",
            font_size=40,
        )
        self.play(Write(step1))
        self.play(Write(eq1))
        self.wait(1.5)

        # ---- Step 2: introduce Liouville substitution ----
        step2 = Tex(
            r"\textbf{Step 2.}\ Substitute $\psi = (b/c)^{1/2}\,\Psi_{\rm new}$ to kill the first-order $\partial_\tau$.",
            font_size=24,
        ).to_edge(UP)
        eq2 = MathTex(
            r"\psi \;=\; \sqrt{\frac{b}{c}}\,\Psi_{\rm new}",
            font_size=44,
            color=YELLOW,
        )
        self.play(Transform(step1, step2))
        self.play(Transform(eq1, eq2))
        self.wait(2.0)

        # ---- Step 3: divide by (b/c)^{1/2}; effective mass emerges ----
        step3 = Tex(
            r"\textbf{Step 3.}\ Divide by $\sqrt{b/c}$. The first-derivative term cancels.",
            font_size=24,
        ).to_edge(UP)
        eq3 = MathTex(
            r"\tfrac{1}{b^2}\partial_\tau^2 \Psi_{\rm new} \;-\; \nabla^2 \Psi_{\rm new} \;+\; V_{\rm eff}\,\Psi_{\rm new} \;=\; \tilde{\text{source}}",
            font_size=38,
        )
        veff_box = MathTex(
            r"V_{\rm eff} \;=\; \frac{\ddot b}{2 b^3} \;-\; \frac{3\dot b^2}{4 b^4}",
            font_size=40,
            color=YELLOW,
        )
        veff_box.next_to(eq3, DOWN, buff=0.6)
        self.play(Transform(step1, step3))
        self.play(Transform(eq1, eq3))
        self.play(Write(veff_box))
        self.wait(2.5)

        # ---- Step 4: rewrite in terms of u, u_dot, u_ddot ----
        step4 = Tex(
            r"\textbf{Step 4.}\ Using $b^2 = c^2 + \mathbf{u}^2$, express $V_{\rm eff}$ via $\mathbf{u}, \dot{\mathbf{u}}, \ddot{\mathbf{u}}$:",
            font_size=22,
        ).to_edge(UP)
        v_eff_kinematic = MathTex(
            r"V_{\rm eff} \;=\; \frac{\mathbf{u}\!\cdot\!\ddot{\mathbf{u}} + \dot{\mathbf{u}}^2}{2 b^4} \;-\; \frac{5(\mathbf{u}\!\cdot\!\dot{\mathbf{u}})^2}{4 b^6}",
            font_size=40,
            color=YELLOW,
        )
        self.play(Transform(step1, step4))
        self.play(FadeOut(eq1))
        self.play(Transform(veff_box, v_eff_kinematic))
        self.wait(2.5)

        # ---- Step 5: identify as photon mass ----
        step5 = Tex(
            r"\textbf{Step 5.}\ Identify $V_{\rm eff}$ as (\textit{effective photon mass})$^2/\hbar^2$:",
            font_size=26,
        ).to_edge(UP)
        photon_mass = MathTex(
            r"\boxed{\;\mu \;=\; \frac{\hbar}{c}\sqrt{\frac{\mathbf{u}\!\cdot\!\ddot{\mathbf{u}} + \dot{\mathbf{u}}^2}{2 b^4} \;-\; \frac{5(\mathbf{u}\!\cdot\!\dot{\mathbf{u}})^2}{4 b^6}}\;}",
            font_size=36,
            color=YELLOW,
        )
        self.play(Transform(step1, step5))
        self.play(Transform(veff_box, photon_mass))
        self.wait(3.0)

        # ---- Closing: physical interpretation ----
        closing = Tex(
            r"\textit{Photon mass $\mu$ is dynamical and source-dependent}", font_size=26
        ).to_edge(UP)
        interp1 = Tex(r"$\mu = 0$ when the source is inertial ($\dot{\mathbf{u}} = 0$).", font_size=24)
        interp2 = Tex(r"$\mu > 0$ only during emission events with $\mathbf{a} \neq \mathbf{0}$.", font_size=24)
        interp1.next_to(photon_mass, DOWN, buff=0.6)
        interp2.next_to(interp1, DOWN, buff=0.3)
        self.play(Transform(step1, closing))
        self.play(FadeIn(interp1))
        self.wait(1.5)
        self.play(FadeIn(interp2))
        self.wait(3.0)
