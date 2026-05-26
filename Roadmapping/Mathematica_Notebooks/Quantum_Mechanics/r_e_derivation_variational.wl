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
(* Section 3.  Expectation values + closure equation assembly.   *)
(* ============================================================ *)
(*                                                                                            *)
(* K_D from DRQM I Eq. (III.4); for the field-free s-state (no B, no spin-orbit, no chain-   *)
(* rule term on a spherically symmetric trial), only the scalar terms survive:                *)
(*   <K_D>_{r_e,aa} = <T>_{r_e,aa} + mc^2 + <V_0>_{r_e,aa} + <V_0^2/(2mc^2)>_{r_e,aa}.        *)
(* <T> computed in the GRADIENT form <hbar^2 |grad psi|^2/(2m)> (positive-definite,           *)
(* variationally appropriate for a soft cutoff where \[Psi](r_e) != 0).                       *)

(* Wolfram MCP 2026-05-26 results (single-line cells for MCP-transport safety): *)
ClearAll[aa, re, r, mm, cc, hbar, ee, alpha];
mean1OverR = (aa + 2 re)/(aa^2 + 2 aa re + 2 re^2);
mean1OverRsq = 2/(aa^2 + 2 aa re + 2 re^2);
meanT = hbar^2/(2 mm aa^2);   (* gradient form -- independent of r_e on this trial *)
meanV0 = -ee^2 mean1OverR;
meanV0sqOver2mc2 = (ee^4/(2 mm cc^2)) mean1OverRsq;
meanKDminusMC2 = meanT + meanV0 + meanV0sqOver2mc2;
Print["<1/r>   = ", mean1OverR];
Print["<1/r^2> = ", mean1OverRsq];
Print["<T>     = ", meanT, "   (note: independent of r_e in gradient form)"];
Print["<K_D-mc^2> = ", meanKDminusMC2];

(* Dimensionless form: lengths in r_0 = e^2/(mc^2), energy in mc^2.                            *)
(* aHat = aa/r_0, rHat = r_e/r_0, alpha = e^2/(hbar c).                                        *)
(*   <T>/(mc^2)               =  1/(2 alpha^2 aHat^2)                                          *)
(*   <V_0>/(mc^2)             = -(aHat + 2 rHat)/(aHat^2 + 2 aHat rHat + 2 rHat^2)             *)
(*   <V_0^2/(2mc^2)>/(mc^2)   =  1/(aHat^2 + 2 aHat rHat + 2 rHat^2)                            *)
(*                                                                                              *)
(* CLOSURE (Route X, strong reading: <K_D> = mc^2, no separate binding subtraction):           *)
(*   Edim(aHat, rHat) = 1/(2 alpha^2 aHat^2) - (aHat + 2 rHat - 1)/(aHat^2 + 2 aHat rHat + 2 rHat^2) = 0. *)

ClearAll[Edim, aHat, rHat, alphaNum];
Edim[aH_, rH_, al_] := 1/(2 al^2 aH^2) - (aH + 2 rH - 1)/(aH^2 + 2 aH rH + 2 rH^2);
alphaNum = 1/137.035999;
Print["Edim at electron-radius scale (aHat=1, rHat=0.5): ", N[Edim[1, 0.5, alphaNum]]];
Print["Edim at Bohr scale (aHat=1/alpha^2, rHat=0.5):    ", N[Edim[1/alphaNum^2, 0.5, alphaNum]]];
Print["Edim at Bohr scale (aHat=1/alpha^2, rHat=0):      ", N[Edim[1/alphaNum^2, 0, alphaNum]]];
Print["Standard hydrogen binding -alpha^2/2 = ", N[-alphaNum^2/2]];

(* --------------------------------------------------------------------------------- *)
(* DIAGNOSTIC (substantive AI):                                                       *)
(*                                                                                    *)
(* (a) ELECTRON-RADIUS SCALE (aHat ~ 1, rHat ~ 0.5).  <T>/(mc^2) ~ 1/(2 alpha^2) ~    *)
(*     9400.  The non-relativistic kinetic-energy expression is invalid here          *)
(*     (localising the electron inside r_0 requires momentum ~ hbar/r_0 = mc/alpha,   *)
(*     super-relativistic).  DRQM I Eq. (III.4) is an EXPANSION valid only for        *)
(*     r >> hbar/(mc) (Compton wavelength) and V_0/(mc^2) << 1 -- both fail at r ~ r_e.*)
(*     So <K_D> as written is NOT the correct operator at the cutoff scale.            *)
(*                                                                                    *)
(* (b) BOHR SCALE (aHat ~ 1/alpha^2 ~ 18800, rHat ~ 0.5).  <K_D-mc^2>/(mc^2) ~        *)
(*     -alpha^2/2 (textbook hydrogen 1s).  But cutoff is INVISIBLE here -- Edim is     *)
(*     identical at rHat = 0 and rHat = 0.5 to 10 sig figs (cutoff is ~r_0; trial      *)
(*     concentrated at ~a_B = r_0/alpha^2; ratio ~alpha^2 ~ 5e-5).  So the closure     *)
(*     <K_D> = mc^2 at Bohr scale cannot pin rHat: dEdim/drHat ~ alpha^6 ~ 10^-13.    *)
(*                                                                                    *)
(* (c) NO INTERMEDIATE SCALE solves the closure: there is no aHat at which the         *)
(*     expansion (III.4) is valid AND the cutoff rHat couples to the closure.         *)
(*                                                                                    *)
(* CONCLUSION.  Closure condition #7 (mass-renormalisation) with the published         *)
(* expanded form of K_D and a non-relativistic exponential trial CANNOT determine     *)
(* rHat first-principles.  To make Route X tractable, one of:                          *)
(*   - The framework must supply an explicit \[CapitalDelta]E_SE^framework(r_e)        *)
(*     to act as the closure target (currently unspecified in DRQM I);                 *)
(*   - OR the calculation must be redone with the un-expanded full Dirac equation     *)
(*     H_D \[Psi] = lambda \[Psi] under a radial cutoff regulator (a 5-10 iteration   *)
(*     arc, not in this notebook's scope).                                            *)
(*                                                                                    *)
(* This is BLOCKED on author input: see Author_Reports/                                *)
(* 2026-05_re_derivation_candidates_for_gill.md Candidate 2 and issue #65.            *)

Print["Section 3: closure equation assembled.  See diagnostic comment block above:  *"];
Print["  closure has no solution at electron-radius scale (NR expansion invalid),    *"];
Print["  no rHat-coupling at Bohr scale (cutoff invisible to trial).                  *"];
Print["  BLOCKED on framework-internal \[CapitalDelta]E_SE(r_e) specification.        *"];

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
