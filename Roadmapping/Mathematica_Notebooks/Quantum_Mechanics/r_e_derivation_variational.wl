(* ::Package:: *)

(* r_e_derivation_variational.wl  --  Companion to issue #65 (candidate-2 variational route) *)
(*                                                                                            *)
(* Goal.  Determine r_e / r_0 first-principles by solving the renormalised dual-Dirac        *)
(* eigenvalue problem under a radial cutoff regulator and imposing a framework-internal      *)
(* mass-renormalisation closure condition.                                                    *)
(*                                                                                            *)
(* Source-of-record.                                                                          *)
(*   - DRQM I \[Section]III, Eqs. (III.1)-(III.23).  Cutoff r_e enters at (III.7) through    *)
(*     the small-(\[Lambda]-mc^2)/mc^2 approximation to (\[Lambda] - V_0 + mc^2)^(-1).       *)
(*   - FoundationsII-Classical \[Section]2.2 (Eq. 2.11).  r_0 is a critical point of the    *)
(*     dual Hamiltonian force F_K = -\[Del]V(1 + V/(mc^2)).                                  *)
(*   - .dev/research/STATE.md iterations 1-2: closure-condition inventory.                   *)
(*                                                                                            *)
(* Route X (this notebook).  Variational trial \[Psi]_1 = N exp(-r/aa) on radial domain      *)
(* [r_e, \[Infinity]).  Variational parameter aa (Bohr-radius-like).  Compute                *)
(* \[LeftAngleBracket]K_D\[RightAngleBracket] symbolically.  Closure condition #7            *)
(* (mass-renormalisation): \[LeftAngleBracket]K_D\[RightAngleBracket] = mc^2 + minimal      *)
(* binding-energy contribution.  Solve coupled system d/d(aa) -> 0 + closure for (r_e, aa). *)
(*                                                                                            *)
(* HONEST SCOPE.  This notebook is SUBSTANTIVE AI.  The choice of trial wavefunction         *)
(* (exp(-r/aa) on [r_e, \[Infinity])) is one of several plausible options; alternatives     *)
(* include hard-wall (\[Psi](r_e)=0) and soft regulator                                       *)
(* (1/r \[Rule] r/(r^2+r_e^2)).  The choice of "minimal binding-energy" prescription         *)
(* (\[CapitalDelta]E_bind = \[LeftAngleBracket]V_0\[RightAngleBracket] = -e^2                *)
(* \[LeftAngleBracket]1/r\[RightAngleBracket]) is a default reading of the framework's      *)
(* mass-renormalisation prescription and requires author endorsement to be terminal.         *)
(* Per CLAUDE.md Crocco-compliance: substantive-AI sections are flagged inline, and the     *)
(* human-acceptance section is left for the human reviewer.                                  *)
(*                                                                                            *)
(* WOLFRAM MCP gotchas (per CLAUDE.md).  Single-line code; potV (not V), ee (not e) for     *)
(* electron charge; Dot non-commutative -- handled by scalar surrogates here since the      *)
(* radial reduction is scalar.                                                                *)
(*                                                                                            *)
(* Author: Trey Morris with Claude Opus 4.7.  Date: 2026-05-26.  Branch: issue/#65.          *)

(* ============================================================ *)
(* Section 1.  Symbol setup and trial wavefunction.              *)
(* ============================================================ *)

ClearAll[m, c, hbar, ee, aa, re, r, r0, NN, psi1];

(* Length scales.  r0 = e^2/(mc^2) is the classical electron radius in Gaussian units.  *)
(* The variational length aa plays the role of a Bohr-radius-like trial parameter.       *)
r0Def = ee^2/(m c^2);

(* Trial upper spinor \[Psi]_1 (s-state, real, normalisable on [r_e, \[Infinity])).         *)
(* Choice: exponential decay exp(-r/aa).  No node at r_e (soft regulator: amplitude       *)
(* is just truncated below r_e, not forced to zero there).                                 *)
psi1[r_, aa_] := Exp[-r/aa]

Print["Section 1.  Trial wavefunction \[Psi]_1(r; aa) = exp(-r/aa) on [r_e, \[Infinity])."];

(* ============================================================ *)
(* Section 2.  Cutoff-restricted normalisation.                  *)
(* ============================================================ *)

(* Compute N^2 such that \[Integral]_{r_e}^{\[Infinity]} 4\[Pi] r^2 |\[Psi]_1|^2 dr = 1. *)
normInt = Integrate[r^2 psi1[r, aa]^2, {r, re, Infinity}, Assumptions -> {aa > 0, re > 0}];
Print["Cutoff-restricted radial norm integral: ", Simplify[normInt]];
Print["Limit re->0 (standard 1s normalisation, expected aa^3/4): ", Simplify[Limit[normInt, re -> 0, Assumptions -> {aa > 0}]]];
(* Wolfram MCP 2026-05-26: closed form = aa(aa^2 + 2 aa re + 2 re^2)/(4 Exp[2 re/aa]); limit re->0 = aa^3/4.  Confirmed. *)

Nsq = 1/(4 Pi normInt);
Print["N^2 = ", Simplify[Nsq]];
Print["Sanity (re->0): N^2 -> 1/(pi aa^3), expected match: ", Simplify[Limit[Nsq, re -> 0, Assumptions -> {aa > 0}] == 1/(Pi aa^3)]];

(* ============================================================ *)
(* Section 3.  PLACEHOLDER -- Compute \[LeftAngleBracket]K_D\[RightAngleBracket]_{r_e, aa}.  *)
(* ============================================================ *)
(*                                                                                            *)
(* Next iteration: build K_D term-by-term from DRQM I Eq. (III.4):                            *)
(*   K_D = pi^2/(2m) + \[Beta] V_0 + mc^2 - (e hbar \[CapitalSigma].B)/(2mc)                  *)
(*       + (V_0 \[Alpha].pi)/(mc) - (i hbar \[Alpha].\[Del]V_0)/(2mc) + V_0^2/(2mc^2).        *)
(* For the field-free s-state (no B, no spin-orbit), only the scalar terms survive:           *)
(*   <K_D>_{r_e,aa} = <pi^2/(2m)>_{r_e,aa} + mc^2 + <V_0>_{r_e,aa} + <V_0^2/(2mc^2)>_{r_e,aa} *)
(* (the spin-flip and chain-rule terms vanish on the spherically-symmetric trial).             *)
(*                                                                                            *)
(* Closure condition #7 (mass-renormalisation):                                                *)
(*   <K_D>_{r_e,aa} = mc^2 + <V_0>_{r_e,aa}   <-- the framework absorbs the bind-energy.      *)
(* Equivalently:                                                                               *)
(*   <pi^2/(2m)>_{r_e,aa} + <V_0^2/(2mc^2)>_{r_e,aa} = 0.                                      *)
(* This is one equation; the second comes from stationarity d/d(aa) of the LHS on the trial. *)
(*                                                                                            *)
(* TODO next iteration: implement <pi^2/(2m)>, <1/r>, <1/r^2> integrals symbolically and      *)
(* feed them into the coupled (aa, r_e) system.  Then cross-check r_e/r_0 against the         *)
(* triangulated 0.4994205099128317 and the Schwinger closed-form.                              *)

Print["Section 3: PLACEHOLDER. <K_D> computation queued for next iteration."];

(* ============================================================ *)
(* Human acceptance section (Crocco compliance).                  *)
(* ============================================================ *)
(*                                                                                            *)
(* Substantive-AI elements requiring human review before terminal acceptance:                  *)
(*   1. Choice of trial wavefunction (exp(-r/aa) soft cutoff vs hard-wall vs r/(r^2+r_e^2))   *)
(*   2. Reading of closure condition #7's <\[CapitalDelta]E_bind> as <V_0>                    *)
(*   3. Whether stationarity in aa is the appropriate companion condition (vs e.g.            *)
(*      stationarity in r_e directly, or virial-theorem constraint)                            *)
(*                                                                                            *)
(* <!-- TODO: human reviews and fills in --> *)
