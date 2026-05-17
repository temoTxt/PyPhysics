"""
BB84 quantum key distribution protocol: Alice's random-basis encoding,
channel transmission, Bob's random-basis measurement, public basis comparison,
key extraction + eavesdropper detection via QBER threshold.

Chapter reference:
  Roadmapping/History/Forward/07_PNT_GPS_SLR_QKD.md#d-quantum-key-distribution-bb84

Render:
    uv run manim -pql --media_dir rendered manim_scenes/pnt_qkd_bb84.py PntQkdBb84
    uv run manim -qh  --media_dir rendered manim_scenes/pnt_qkd_bb84.py PntQkdBb84
"""

from manim import (
    BLUE, DOWN, FadeOut, GREEN, LEFT, MathTex, RED, RIGHT, Scene, Tex, Transform,
    UP, VGroup, WHITE, Write, YELLOW,
)


class PntQkdBb84(Scene):
    """BB84 protocol walkthrough."""

    def construct(self):
        title = Tex(r"\textbf{BB84 quantum key distribution}", font_size=40)
        subtitle = Tex(r"Bennett \& Brassard 1984 --- Chapter 7 \S D, PNT",
                       font_size=22)
        subtitle.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), Write(subtitle))
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(subtitle))

        # Step 1: Alice's random encoding
        step1_title = Tex(r"\textbf{Step 1: Alice encodes random bits in random bases}",
                          color=BLUE, font_size=28).to_edge(UP)
        step1_table = VGroup(
            Tex(r"Random bit: 0 1 1 0 1 0 1 1 0", font_size=22),
            Tex(r"Random basis: + $\times$ $\times$ + + $\times$ + $\times$ +", font_size=22),
            Tex(r"$+$ basis: 0 $\to$ $\to$ (right); 1 $\to$ $\uparrow$ (up)", font_size=22),
            Tex(r"$\times$ basis: 0 $\to$ $\nearrow$ (diagonal); 1 $\to$ $\searrow$ (anti-diagonal)", font_size=22),
            Tex(r"Photon: $\to$ $\searrow$ $\searrow$ $\to$ $\uparrow$ $\nearrow$ $\uparrow$ $\searrow$ $\to$",
                font_size=22),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.3).next_to(step1_title, DOWN, buff=0.4)
        self.play(Write(step1_title))
        for line in step1_table:
            self.play(Write(line))
            self.wait(0.15)
        self.wait(1.0)
        self.play(FadeOut(step1_title), FadeOut(step1_table))

        # Step 2: Bob's random measurement
        step2_title = Tex(r"\textbf{Step 2: Bob measures in random bases}",
                          color=GREEN, font_size=28).to_edge(UP)
        step2_lines = VGroup(
            Tex(r"Bob's random basis: + + $\times$ $\times$ + $\times$ $\times$ $\times$ +", font_size=22),
            Tex(r"\textit{When basis matches Alice's: result deterministic.}", font_size=22),
            Tex(r"\textit{When basis mismatches: result is random (50/50).}", font_size=22),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.3).next_to(step2_title, DOWN, buff=0.4)
        self.play(Write(step2_title))
        for line in step2_lines:
            self.play(Write(line))
            self.wait(0.25)
        self.wait(1.0)
        self.play(FadeOut(step2_title), FadeOut(step2_lines))

        # Step 3: Basis sifting
        step3_title = Tex(r"\textbf{Step 3: Public basis comparison (sifting)}",
                          color=YELLOW, font_size=28).to_edge(UP)
        step3_lines = VGroup(
            Tex(r"Over a classical channel, Alice and Bob exchange the \textit{bases used},", font_size=22),
            Tex(r"\textit{not} the measurement results. Keep only events where bases matched.", font_size=22),
            Tex(r"Expected fraction kept: $\sim 1/2$ (in BB84). Result: the \textit{sifted key}.", font_size=22),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.3).next_to(step3_title, DOWN, buff=0.4)
        self.play(Write(step3_title))
        for line in step3_lines:
            self.play(Write(line))
            self.wait(0.25)
        self.wait(1.5)
        self.play(FadeOut(step3_title), FadeOut(step3_lines))

        # Step 4: Eavesdropper detection
        step4_title = Tex(r"\textbf{Step 4: Eavesdropper detection (QBER)}",
                          color=RED, font_size=28).to_edge(UP)
        step4_lines = VGroup(
            Tex(r"\textit{No-cloning theorem:} an eavesdropper Eve cannot copy unknown qubits.", font_size=22),
            Tex(r"Any in-channel measurement Eve performs collapses the state in some basis,", font_size=22),
            Tex(r"introducing errors when Bob's basis happens to match Alice's but not Eve's.", font_size=22),
            Tex(r"Alice and Bob \textit{publicly compare} a random subset of their sifted-key bits.", font_size=22),
            Tex(r"If $\text{QBER} > \sim 11\%$, abort. Otherwise: privacy amplification $\to$ final key.",
                font_size=22),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.3).next_to(step4_title, DOWN, buff=0.4)
        self.play(Write(step4_title))
        for line in step4_lines:
            self.play(Write(line))
            self.wait(0.25)
        self.wait(2.0)
        self.play(FadeOut(step4_title), FadeOut(step4_lines))

        # Dual-framework status
        dual = Tex(
            r"\textbf{Dual-framework status:} \#gill-silent on the protocol itself. \\ "
            r"The no-cloning theorem is a structural property of the Hilbert-space inner product, \\ "
            r"common to standard and dual QM. BB84 works identically in both framings. \\[0.5em] "
            r"\textit{Speculative extension (Ch 8): in Gill's KS-Hilbert space, the no-cloning} \\ "
            r"\textit{argument could in principle be re-examined --- whether the security guarantee} \\ "
            r"\textit{changes is an open question deferred to Chapter 8.}",
            font_size=22, color=YELLOW,
        )
        self.play(Write(dual))
        self.wait(4.0)
        self.play(FadeOut(dual))
