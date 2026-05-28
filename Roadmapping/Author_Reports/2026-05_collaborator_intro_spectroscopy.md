---
title: "Dual-theory reformulation of hydrogen-like precision spectroscopy"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-05-28"
subject: "Introductory note for a new experimental collaborator in optical/microwave atomic spectroscopy — unified method for hydrogen-like atoms + per-observable predictions and precision-measurement comparisons"
---

# Dual-theory reformulation of hydrogen-like precision spectroscopy — intro note for a new collaborator

**Date:** 2026-05-28.
**From:** Trey Morris (with Claude Opus 4.7).
**Re:** Draft introductory email to a new experimental collaborator in optical/microwave atomic spectroscopy. Tracks against [issue #88](https://github.com/temoTxt/PyPhysics/issues/88). Methods-focused framing; honest-scope paragraph is load-bearing and must be retained before sending. Source material: the Bethe–Salpeter precision campaign ([#50](https://github.com/temoTxt/PyPhysics/issues/50)), the Li$^{2+}$ Z-extension ([#78](https://github.com/temoTxt/PyPhysics/issues/78)), and the hydrogenic-ion g-factor Z-scan ([#82](https://github.com/temoTxt/PyPhysics/issues/82)).

---

## Draft email

*Replace `[Name]`, the salutation, and the sign-off before sending. Verify the numbers against the provenance table below.*

> **Subject:** Dual-theory reformulation of hydrogen-like precision spectroscopy — would value your read
>
> Dear [Name],
>
> I'm a co-author on a reformulation of relativistic quantum mechanics — Gill's "dual theory," built on a proper-time formulation — and we've been working through its predictions for precision atomic spectroscopy. Since your work is in optical/microwave spectroscopy of hydrogen-like systems, I'd value your read on the apparatus and your sense of which measurements are most diagnostic. This note summarizes the method and where it stands against measurement.
>
> **The method, in one paragraph.** The dual theory recasts relativistic dynamics through a proper-time canonical Hamiltonian, $K = H^2/(2mc^2) + mc^2/2$, with a dual-Dirac equation in place of the standard Dirac equation. For precision atomic observables this reduces to a clean structure: each $g_s$-dependent observable is the textbook (Sommerfeld–Dirac / Fermi-contact) result times $(g_s/{-2})^{n}$, with $n=1$ for spin–orbit and Fermi-contact contributions and $n=2$ for two-fermion spin–spin couplings. The electron spin g-factor $g_s$ itself enters through a single cutoff parameter in the framework's §III.D structure, and the whole apparatus scales across the hydrogen-like sequence with the expected $Z$-powers. The upshot is that one consistent prescription handles the Lamb shift, fine structure, and hyperfine splitting across H, He$^{+}$, Li$^{2+}$, and up the isoelectronic sequence.
>
> **Where it stands against measurement (optical/microwave regime).**
>
> | Observable | Measurement | Framework | Note |
> |---|---|---|---|
> | H Lamb shift $2S_{1/2}$–$2P_{1/2}$ | $1057.845(9)$ MHz | $\sim 1016$ MHz | Bethe-1947 estimate; $\sim 42$ MHz ($\sim 4\%$) residual at the estimate's precision floor |
> | H fine structure $2P_{3/2}$–$2P_{1/2}$ | $10\,969.13(10)$ MHz | $10\,962$ MHz | leading Dirac × anomalous-$g$ |
> | H 1S hyperfine (21-cm) | $1\,420.405\,751\,768(2)$ MHz | $1\,420.04$ MHz | Fermi contact at the framework cutoff; $\sim 0.4$ MHz residual |
> | Li$^{2+}$ Lamb shift $2S_{1/2}$–$2P_{1/2}$ | $62\,765(21)$ MHz (Schiffer 1995) | $\sim 61$–$63$ GHz | $Z^4$-scaled Bethe-estimate |
> | Li$^{2+}$ fine structure $2P_{3/2}$–$2P_{1/2}$ | — (no direct precision measurement located) | $\approx 888$ GHz | $Z^4 \times$ H |
> | Li$^{2+}$ 1s hyperfine | — (theory comparator: Pachucki 2023) | $\approx 29.8$ GHz | $Z^3 \times$ nuclear factors, $I=3/2$ |
>
> **Honest scope.** At the precision the apparatus currently delivers — the leading-$g_s$ / Bethe-estimate floor — the framework reproduces the standard QED results rather than predicting deviations from them; the agreement is by construction of the reformulation, not an independent corroboration. The frontier where a genuinely distinct prediction could appear is a full proper-time one-loop calculation, which we have not yet completed. I mention this up front because I'd rather be precise about what's reformulation and what would be new physics.
>
> **What would be most useful from you.** Two things: (1) an expert read on whether the proper-time / dual-Dirac reduction is faithful to how these observables are actually measured and extracted; and (2) your sense of which optical/microwave measurements — higher-$Z$ hydrogen-like fine structure or hyperfine, improved $2S$–$2P$ intervals, or something else — would most sharply test a reformulation of this kind, particularly any regime where the proper-time radiation-reaction structure might leave a measurable signature.
>
> Happy to send the full per-observable derivations (each reproduced symbolically in Mathematica) if useful. Thank you for considering it.
>
> Best regards,
> [Trey Morris]

---

## Backing provenance table (for Trey to verify before sending — not part of the email)

| Observable | Framework value | Source in repo | Measurement | Measurement source |
|---|---|---|---|---|
| H Lamb shift $2S_{1/2}$–$2P_{1/2}$ | $\sim 1016$ MHz | [`Bethe_Salpeter/10_CrossComparison.md`](../Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md) PR E row | $1057.845(9)$ MHz | CODATA-2018 |
| H fine structure $2P_{3/2}$–$2P_{1/2}$ | $10\,962$ MHz | [`03_FineStructure.md`](../Quantum_Mechanics/Bethe_Salpeter/03_FineStructure.md) BS-§14.2 | $10\,969.13(10)$ MHz | Hagley & Pipkin 1994 |
| H 1S hyperfine (21-cm) | $1\,420.04$ MHz | [`10_CrossComparison.md`](../Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md) §1 residual table | $1\,420.405\,751\,768(2)$ MHz | NIST / Karshenboim review |
| Li$^{2+}$ Lamb shift | $\sim 61$–$63$ GHz | `12_Li2plus_Spectroscopy.md` §12.2 (branch [#85](https://github.com/temoTxt/PyPhysics/pull/85)) | $62\,765(21)$ MHz | Schiffer et al., *PRL* **74**, 2188 (1995) |
| Li$^{2+}$ fine structure | $887{,}920$ MHz | `12_Li2plus_Spectroscopy.md` §12.3 (branch [#85](https://github.com/temoTxt/PyPhysics/pull/85)) | none located | — (theory: $Z^4$-Dirac) |
| Li$^{2+}$ 1s hyperfine | $\sim 29.8$ GHz | `13_Li2plus_Hyperfine.md` (branch [#86](https://github.com/temoTxt/PyPhysics/pull/86)) | none direct | Pachucki et al. 2023 (theory from Li$^{+}$) |

**Note on the Li$^{2+}$ rows:** the Li$^{2+}$ predictions live on the not-yet-merged PR branches [#85](https://github.com/temoTxt/PyPhysics/pull/85) (Lamb shift + fine structure) and [#86](https://github.com/temoTxt/PyPhysics/pull/86) (hyperfine). Verify against the merged versions once those PRs land, or against the branch docs in the interim. The Li$^{2+}$ fine-structure and hyperfine rows have **no direct precision measurement** — the email marks them "—" honestly; do not imply an experimental comparison exists for those two.

**Note on the honest-scope paragraph:** it is consistent with the campaign's own framing in [`10_CrossComparison.md`](../Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md): *"zero of the 28 recorded results constitute independent corroborations of the dual theory's content distinct from textbook QM."* The hydrogenic-ion g-factor Z-scan ([#82](https://github.com/temoTxt/PyPhysics/issues/82)) independently confirmed this on the Z-axis (the framework's cutoff is Z-trivial; it inherits QED's bound-state structure rather than predicting deviations). Keeping the honest-scope paragraph is non-negotiable for an external-scientist audience.

<!-- TODO: human reviews and fills in — Trey to (a) confirm the per-observable numbers against the merged campaign docs, (b) fill in [Name] / salutation / sign-off, (c) confirm the honest-scope paragraph is retained as worded, and (d) decide whether to attach the full per-observable Mathematica derivations or offer them on request. The methods-focused framing and the "—" (no measurement) marks on the two Li$^{2+}$ rows are substantive AI choices requiring sign-off per Crocco rule #1. -->
