(* ::Package:: *)

(* r_e_Zscan_fit.wl  --  Companion to Quantum_Mechanics/Bethe_Salpeter/14_HydrogenicIon_Zscan.md *)
(*                       and Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md       *)
(*                       Eqs. (III.21)-(III.23) + FINDINGS_for_author_review.md Finding 2.         *)
(*                                                                                                 *)
(* PURPOSE.  Issue #82.  Extend the Z=1 g-factor back-fit of r_e/r_0 (r_e_triangulation.wl) to a  *)
(* Z-scan across hydrogen-like ions.  For each ion of nuclear charge Z with a measured bound-      *)
(* electron g-factor g_meas(Z), test two readings of the DRQM I Eq. (III.22) cutoff prescription   *)
(*   g_r(x) = 2 (1 - 4/(2 x + 1)),  x = r/r_0 :                                                     *)
(*   (Z-i)  UNIVERSAL cutoff: x fixed at the Z=1 triangulated value 0.4994205099128317 -> the       *)
(*          framework predicts the SAME g = -2.00231930... at every Z.  chi^2 across all Z.        *)
(*   (Z-ii) Z-SCALED cutoff: invert g_r at each measured g to get x^(Z) = (2 - a)/(2 (2 + a))      *)
(*          with a = (|g_meas| - 2)/2; tabulate vs (Z alpha)^2 and fit a + b (Z alpha)^2.          *)
(*                                                                                                 *)
(* HONEST SCOPE (load-bearing).  The (Z-ii) per-Z back-fit is a one-to-one re-encoding of the      *)
(* measured g(Z): x^(Z) carries no information the measurement did not already supply.  The real   *)
(* falsifiable test is (Z-i) -- whether a SINGLE x reproduces all Z.  A clean (Z-ii) form-fit in   *)
(* (Z alpha)^2 only demonstrates that the back-fit inherits QED's bound-state g(Z alpha) = 2[1 -   *)
(* (1/3)(Z alpha)^2 - ...] structure; the framework supplies the g_r formula but NOT the binding   *)
(* coefficient -1/3, so a successful form-fit is Outcome C (QED inheritance), not B, unless the     *)
(* framework can derive that coefficient internally (it cannot, per PR #70 / DRQM I Sec III.D).    *)
(*                                                                                                 *)
(* CROCCO TAG.  SUBSTANTIVE AI on the Z-scaling functional-form choice (Section 4): the decision    *)
(* to model x^(Z) as a + b (Z alpha)^2 vs higher-order is an interpretive modelling choice.  See    *)
(* human-acceptance block at the foot of 14_HydrogenicIon_Zscan.md.                                 *)
(*                                                                                                 *)
(* MEASUREMENT PROVENANCE (full DOIs in 14_HydrogenicIon_Zscan.md):                                *)
(*   Z=2  3He+   g = 2.0021774158   Schneider 2022  Nature 606,878   10.1038/s41586-022-04761-7    *)
(*   Z=6  12C5+  g = 2.0010415902   Sturm 2014      Nature 506,467   10.1038/nature13026           *)
(*   Z=8  16O7+  g = 2.0000470254   Verdu 2004      PRL 92,093002    10.1103/PhysRevLett.92.093002 *)
(*   Z=14 28Si13+ g = 1.9953489587  Sturm 2011      PRL 107,023002   10.1103/PhysRevLett.107.023002*)
(*   Z=50 118Sn49+ g = 1.9105620590 Morgner 2023    Nature 622,53    10.1038/s41586-023-06453-2    *)
(*   (Z=3 7Li2+ imported from #78 Self-Energy branch -- NOT included here; brief placeholder        *)
(*    value 2.0000251707 is unphysical, see STATE.md iter 9.  Slot in when #78 delivers.)          *)
(*                                                                                                 *)
(* WOLFRAM MCP GOTCHAS (CLAUDE.md): single-line cells joined with ';'; no bare 'e' (Euler) or 'V'  *)
(* (Vanadium); use alpha / gMeas / etc.  No symbolic Dot here.                                      *)
(*                                                                                                 *)
(* RUN:  via Wolfram MCP, one cell at a time (each Section is one transport-safe line).             *)

(* ============================================================ *)
(* SECTION 0 -- constants + measured data                       *)
(* ============================================================ *)

ClearAll[alpha, gr, backfit, Zvals, gMeas, gSig, xUniversal, gFreeMag];

alpha = 7.2973525693*^-3; (* CODATA 2018 fine-structure constant *)

gFreeMag = 2.00231930436256; (* free-electron |g| (CODATA) *)

xUniversal = 0.4994205099128317; (* Z=1 triangulated r_e/r_0 (r_e_triangulation.wl Pass B) *)

Zvals = {2, 6, 8, 14, 50};

gMeas = {2.0021774158, 2.0010415902, 2.0000470254, 1.9953489587, 1.9105620590}; (* |g_meas| per ion *)

gSig = {4.5*^-10, 3.0*^-11, 4.6*^-9, 1.0*^-9, 1.0*^-9}; (* combined 1-sigma on g (absolute units) *)

(* ============================================================ *)
(* SECTION 1 -- g-formula (III.22) and its algebraic inverse    *)
(* ============================================================ *)

gr[x_] := 2 (1 - 4/(2 x + 1)); (* DRQM I Eq. (III.22); returns signed g_r (negative branch for electron via |.|) *)

backfit[a_] := (2 - a)/(2 (2 + a)); (* invert g_r(x) = -2(1+a) for x = r_e/r_0 given anomaly a *)

(* sanity: gr[1/2] == -2 (tree Dirac); gr[xUniversal] should give -2.00231930... *)
Print["S1 check: gr[1/2] = ", InputForm[gr[1/2]], " ; |gr[xUniversal]| = ", InputForm[Abs[N[gr[xUniversal], 18]]]];

(* ============================================================ *)
(* SECTION 2 -- per-Z anomalies and (Z-ii) back-fit cutoffs     *)
(* ============================================================ *)

aBound = (gMeas - 2)/2; (* a_e^bound(Z) = (|g_meas| - 2)/2, elementwise *)

xBackfit = backfit /@ aBound; (* per-Z r_e^(Z)/r_0 *)

Print["S2 (Z, a_bound, x_backfit):"]; Do[Print["  Z=", Zvals[[k]], "  a=", InputForm[N[aBound[[k]], 12]], "  x=", InputForm[N[xBackfit[[k]], 12]]], {k, 1, Length[Zvals]}];

(* ============================================================ *)
(* SECTION 3 -- (Z-i) UNIVERSAL-cutoff test: chi^2 with x fixed *)
(* ============================================================ *)

(* Prediction under (Z-i): |g| = |gr[xUniversal]| at EVERY Z (Z-independent). *)
gPredUniversal = Abs[N[gr[xUniversal], 18]];

residUniversal = gMeas - gPredUniversal; (* measured - predicted, elementwise *)

chi2Universal = Total[(residUniversal/gSig)^2];

Print["S3 (Z-i) universal-cutoff test:"]; Print["  g_pred (all Z) = ", InputForm[gPredUniversal]]; Do[Print["  Z=", Zvals[[k]], "  resid=", InputForm[N[residUniversal[[k]], 6]], "  (", InputForm[N[residUniversal[[k]]/gSig[[k]], 4]], " sigma)"], {k, 1, Length[Zvals]}]; Print["  chi^2 (Z-i, ", Length[Zvals], " ions, 0 free params) = ", InputForm[N[chi2Universal, 6]]];

(* ============================================================ *)
(* SECTION 4 -- (Z-ii) Z-scaling form-fit  x^(Z) = c0 + c2 (Z a)^2  [+ c4 (Z a)^4] *)
(* ============================================================ *)

za2 = (alpha Zvals)^2; (* (Z alpha)^2 per ion *)

(* Build {(Z alpha)^2, x_backfit} pairs and least-squares fit. *)
fitData2 = Transpose[{za2, xBackfit}];

ClearAll[u]; fitLin = Fit[fitData2, {1, u}, u]; (* x = c0 + c2 (Z a)^2 *)
fitQuad = Fit[fitData2, {1, u, u^2}, u]; (* + c4 (Z a)^4 *)

Print["S4 (Z-ii) Z-scaling fit of x_backfit vs (Z alpha)^2:"]; Print["  linear : ", InputForm[fitLin]]; Print["  quad   : ", InputForm[fitQuad]];

(* Compare intercept c0 to the Z=1 triangulated x and to the tree-Dirac 1/2. *)
Print["  intercept c0 (lin) = ", InputForm[N[fitLin /. u -> 0, 12]], "  vs xUniversal=", InputForm[xUniversal], "  vs 1/2"];

(* ============================================================ *)
(* SECTION 5 -- cross-check vs QED bound-state structure         *)
(* ============================================================ *)

(* QED leading bound-state: |g(Z a)| ~ 2(1 - (1/3)(Z a)^2 - (1/12)(Z a)^4) + alpha/Pi.  *)
(* The x_backfit curve should mirror this; the intercept-shift from 1/2 is the free anomaly. *)
gQEDlead[Z_] := 2 (1 - (1/3) (alpha Z)^2 - (1/12) (alpha Z)^4) + alpha/Pi;

Print["S5 QED-leading vs measured |g| (consistency, not a fit):"]; Do[Print["  Z=", Zvals[[k]], "  g_QEDlead=", InputForm[N[gQEDlead[Zvals[[k]]], 10]], "  g_meas=", InputForm[gMeas[[k]]]], {k, 1, Length[Zvals]}];

(* ============================================================ *)
(* MCP RESULTS -- executed 2026-05-27 (Wolfram MCP), STATE.md iter 11.                              *)
(*                                                                                                 *)
(* S1:  gr[1/2] = -2  (exact)                                                                      *)
(*      |gr[xUniversal]| = 2.00231930436256  (matches free |g|)                                    *)
(*                                                                                                 *)
(* S2  per-Z back-fit x^(Z) = r_e^(Z)/r_0  (and a_bound):                                           *)
(*      Z=2   a=+0.00108870790  x=0.49945594221                                                     *)
(*      Z=6   a=+0.00052079510  x=0.49973967024                                                     *)
(*      Z=8   a=+0.00002351270  x=0.49998824379                                                     *)
(*      Z=14  a=-0.00232552065  x=0.50116411391                                                     *)
(*      Z=50  a=-0.04471897050  x=0.52287086604                                                     *)
(*                                                                                                 *)
(* S3  (Z-i) UNIVERSAL-cutoff test (x fixed = xUniversal, 0 free params):                           *)
(*      g_pred (all Z) = 2.00231930436256                                                           *)
(*      residuals  = {-1.42e-4, -1.28e-3, -2.27e-3, -6.97e-3, -9.18e-2}                             *)
(*      sigma-cts  = {-3.15e5, -4.26e7, -4.94e5, -6.97e6, -9.18e7}                                  *)
(*      chi^2      = 1.0282e16   ==>  OUTCOME A DECISIVELY REJECTED.                                *)
(*                                                                                                 *)
(* S4  (Z-ii) Z-scaling form-fit of x^(Z) vs (Z alpha)^2:                                           *)
(*      linear : x = 0.4993835903 + 0.17639313 (Z a)^2                                              *)
(*      quad   : x = 0.4994206084 + 0.16627568 (Z a)^2 + 0.07415407 (Z a)^4                         *)
(*      quad-intercept c0 = 0.4994206084  vs  xUniversal = 0.4994205099  (match to 9.8e-8)          *)
(*      quad-slope    c2 = 0.16627568     vs  QED 1/6 = 0.16666667        (match to 2.3e-3)          *)
(*      linear residuals = {+3.5e-5, +1.8e-5, +3.5e-6, -6.1e-5, +4.3e-6} (structured -> needs quad) *)
(*                                                                                                 *)
(* INTERPRETATION (-> 14_HydrogenicIon_Zscan.md verdict):                                           *)
(*   * Outcome A (single universal cutoff) is excluded at chi^2 ~ 1e16 over 5 ions (Z=2..50).      *)
(*   * The (Z-ii) back-fit x^(Z) is well-described by x = c0 + c2(Za)^2 + c4(Za)^4 with             *)
(*       c0 = xUniversal (the Z=1 / free-electron triangulated cutoff, recovered as Z->0 to 1e-7), *)
(*       c2 ~ 1/6 (the QED bound-state leading coefficient = -dx/da * 1/3 = (1/2)(1/3)).           *)
(*   * The clean (Za)^2 form is a property of QED's bound-state g(Za) = 2[1-(1/3)(Za)^2-...],       *)
(*     INHERITED by the inversion -- NOT derived from the framework apparatus (which supplies only  *)
(*     gr[x] and leaves each state's cutoff free; the -1/3 coefficient is QED's, not the dual       *)
(*     framework's).  Therefore OUTCOME C, extending PR #70's lepton-axis verdict                   *)
(*     ("particle-specific through a_l") to the Z-axis ("Z-specific through a_e^bound(Za)").        *)
(* ============================================================ *)
