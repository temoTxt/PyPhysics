(* ::Package:: *)

(* r_e_schwinger_residual_test.wl  --  Candidate 3 empirical-test path *)
(*                                                                                                       *)
(* Companion to issue #66 (Candidate 3: Closed-form Schwinger identification of the triangulated r_e). *)
(*                                                                                                       *)
(* GOAL.  Test whether the residual between the triangulated r_e/r_0 = 0.4994205099128317 (PR #62)     *)
(* and the closed-form Schwinger one-loop value r_e/r_0 = (2 - alpha/(2 Pi))/(4 + alpha/Pi)             *)
(* systematically tracks the Karplus-Kroll two-loop QED correction to g_e (and the Laporta-Remiddi      *)
(* three-loop and Kinoshita-Fukuda-Sasaki-style four-loop terms at lower hierarchy).                    *)
(*                                                                                                       *)
(* STRUCTURAL CAVEAT (load-bearing -- read before interpreting the numerics).  The triangulated r_e    *)
(* from PR #62 is defined as the cutoff that makes g_r(r) = g_e^measured to 16 sig figs (Pass B fit). *)
(* Therefore the residual                                                                                *)
(*     Delta g_e^observed := g_e^measured - g_e^Schwinger-one-loop                                     *)
(* is EXACTLY EQUAL BY CONSTRUCTION to the all-orders QED contributions beyond Schwinger one-loop      *)
(* (Karplus-Kroll two-loop + Laporta-Remiddi three-loop + Kinoshita-Fukuda four-loop + ...).            *)
(* The "agreement" measured in Section 5 below is therefore not an independent test of whether the     *)
(* framework's cutoff prescription literally produces the closed-form expression -- it is a            *)
(* consistency check that the triangulated value lies on the all-orders-QED curve in r-space.  The     *)
(* test discriminates "framework cutoff = one-loop Schwinger + higher-order corrections from elsewhere *)
(* in the framework" from "framework cutoff = something else entirely that happens to coincide", but   *)
(* CANNOT distinguish "Schwinger one-loop encoded by design" from "all-orders QED encoded by design"   *)
(* without an independent observable not of the (g_s/-2)^n form.                                       *)
(*                                                                                                       *)
(* The substantive empirical-test outcome is therefore:                                                  *)
(*   - The triangulated r_e is consistent with the Schwinger-one-loop closed-form + standard QED       *)
(*     higher-order corrections to better than 10^-11 in g_e space.                                    *)
(*   - This is necessary for any intentional-encoding scenario but not sufficient to identify which    *)
(*     scenario holds.                                                                                  *)
(*   - The discriminating question is whether the framework's DRQM-I Section III.D derivation produces *)
(*     the closed-form (2 - alpha/(2 Pi))/(4 + alpha/Pi) as an identity (then higher-order QED must    *)
(*     come from loop corrections within the framework not yet derived) or produces the triangulated   *)
(*     value 0.499420509912831 7 (then the cutoff itself encodes all-orders QED, which requires        *)
(*     justification at the framework-derivation level).                                                *)
(*                                                                                                       *)
(* CROCCO COMPLIANCE.  Numerical evaluation, Wolfram MCP cell execution, and constant tabulation are   *)
(* PRAGMATIC AI.  The structural caveat above, the interpretation of the residual as not-independent, *)
(* and the formulation of the Tepper-discriminating question are SUBSTANTIVE AI.  The                  *)
(* human-acceptance section at the bottom is left for the human reviewer.                              *)
(*                                                                                                       *)
(* WOLFRAM MCP GOTCHAS.  No use of single-letter `e` or `V` (those collide with built-ins per          *)
(* CLAUDE.md).  Each Print/computation block is on a single line joined by `;`.                        *)
(*                                                                                                       *)
(* Author: Trey Morris with Claude Opus 4.7.  Date: 2026-05-26.  Iteration 2 of issue #66 loop.         *)


(* ============================================================ *)
(* Section 1.  Constants and the framework's g_r(r) formula.    *)
(* ============================================================ *)

ClearAll[gr, r, alpha, aPi, rClosedForm, rTriangulated, gSchwinger, gMeas];

(* CODATA 2018 fine-structure constant, to 30 sig figs working precision. *)
alpha = SetPrecision[7.2973525693*^-3, 30];

(* alpha/Pi -- the natural QED expansion parameter for a_e (electron anomaly). *)
aPi = alpha/Pi;

(* DRQM I Eq. (III.21) cutoff-to-g_r map.  Working in units r_0 = 1, so r here is r_e/r_0. *)
gr[r_] := 2 (1 - 4/(2 r + 1));

(* Schwinger one-loop QED anomalous moment:  g_e = -2 - alpha/Pi.  Bare one-loop, no higher-order. *)
gSchwinger = -2 - alpha/Pi;

(* CODATA 2018 measured g_e.  Uncertainty ~ 3 x 10^-13 in the last digit; this is all-orders QED. *)
gMeas = SetPrecision[-2.00231930436256, 20];

(* Inverting g_r(r) = g_Schwinger gives the closed-form r_e/r_0 from issue #66:               *)
(*   2 - 8/(2 r + 1) = -2 - alpha/Pi                                                            *)
(*   8/(2 r + 1) = 4 + alpha/Pi                                                                 *)
(*   2 r + 1 = 8/(4 + alpha/Pi)                                                                 *)
(*   r = (2 - alpha/(2 Pi))/(4 + alpha/Pi)                                                      *)
rClosedForm = (2 - alpha/(2 Pi))/(4 + alpha/Pi);

(* PR #62 triangulated value.  From r_e_triangulation.wl Pass B optimum.                       *)
rTriangulated = SetPrecision[0.4994205099128317, 20];

Print["===== Section 1: constants ====="];
Print["  alpha (CODATA 2018) = ", N[alpha, 20]];
Print["  alpha/Pi            = ", N[aPi, 20]];
Print["  r_closed-form       = ", N[rClosedForm, 20]];
Print["  r_triangulated      = ", N[rTriangulated, 20]];
Print["  g_Schwinger         = ", N[gSchwinger, 20]];
Print["  g_e (measured)      = ", N[gMeas, 20]];


(* ============================================================ *)
(* Section 2.  Verify the closed-form algebraic identity.       *)
(* ============================================================ *)

Module[{check}, check = N[gr[rClosedForm] - gSchwinger, 25]; Print["===== Section 2: closed-form verification ====="]; Print["  g_r(r_closed-form) - g_Schwinger = ", check, "  (should be 0 to working precision)"]];


(* ============================================================ *)
(* Section 3.  Residual r_triangulated - r_closed-form propagated to g_e. *)
(* ============================================================ *)

Module[{dr, dgrDr, dgeFromDr, dgeObserved}, dr = rTriangulated - rClosedForm; dgrDr = N[D[gr[r], r] /. r -> rClosedForm, 20]; dgeFromDr = dr*dgrDr; dgeObserved = gMeas - gSchwinger; Print["===== Section 3: r-residual and g_e-residual ====="]; Print["  dr = r_triangulated - r_closed-form = ", N[dr, 8]]; Print["  dg_r/dr at r_closed-form = 16/(2 r+1)^2 = ", N[dgrDr, 8]]; Print["  dg_e (from dr * dg/dr)           = ", N[dgeFromDr, 8]]; Print["  dg_e_observed = g_meas - g_Schwinger = ", N[dgeObserved, 8]]; Print["  (these two should agree if the framework's g_r(r) Taylor expansion in dr is dominated by the linear term)"]];


(* ============================================================ *)
(* Section 4.  Standard-QED higher-order coefficients for g_e.  *)
(* ============================================================ *)

(* a_e (the electron anomaly) = -(g_e + 2)/2.  Standard-QED expansion (mass-independent, dimensional regularisation): *)
(*   a_e = (alpha/2 Pi) - C2 (alpha/Pi)^2 + C3 (alpha/Pi)^3 - C4 (alpha/Pi)^4 + ...                        *)
(* with C2 = 0.328478965579193... (Karplus-Kroll 1950 + Sommerfield 1957 + Petermann 1957)                  *)
(*      C3 = 1.181241456587  ...  (Laporta-Remiddi 1996)                                                    *)
(*      C4 = 1.9106  ...           (Kinoshita-Fukuda-Sasaki numerical + Aoyama-Hayakawa-Kinoshita-Nio refinements) *)
(*                                                                                                            *)
(* Translating to g_e = -2 (1 + a_e):                                                                       *)
(*   delta g_e^two-loop   = -2 * (-C2 (alpha/Pi)^2) = +2 C2 (alpha/Pi)^2  (positive shift in g_e)           *)
(*   delta g_e^three-loop = -2 * (+C3 (alpha/Pi)^3) = -2 C3 (alpha/Pi)^3  (negative shift in g_e)           *)
(*   delta g_e^four-loop  = -2 * (-C4 (alpha/Pi)^4) = +2 C4 (alpha/Pi)^4  (positive shift in g_e)           *)

Module[{c2, c3, c4, kkTwo, lrThree, kfFour, dgePredAll, dgePredKK, dgeObserved, ratioKK, ratioAll, residual}, c2 = SetPrecision[0.328478965579193, 16]; c3 = SetPrecision[1.181241456587, 14]; c4 = SetPrecision[1.9106, 5]; kkTwo = 2*c2*aPi^2; lrThree = -2*c3*aPi^3; kfFour = 2*c4*aPi^4; dgePredAll = kkTwo + lrThree + kfFour; dgePredKK = kkTwo; dgeObserved = gMeas - gSchwinger; ratioKK = dgeObserved/dgePredKK; ratioAll = dgeObserved/dgePredAll; residual = dgeObserved - dgePredAll; Print["===== Section 4: higher-order QED predictions for g_e residual ====="]; Print["  C2 (Karplus-Kroll, two-loop) = ", N[c2, 16]]; Print["  C3 (Laporta-Remiddi, three-loop) = ", N[c3, 14]]; Print["  C4 (Kinoshita-Fukuda+, four-loop, approx) = ", N[c4, 5]]; Print["  dg_e^two-loop  = +2 C2 (alpha/Pi)^2 = ", N[kkTwo, 8]]; Print["  dg_e^three-loop = -2 C3 (alpha/Pi)^3 = ", N[lrThree, 8]]; Print["  dg_e^four-loop  = +2 C4 (alpha/Pi)^4 = ", N[kfFour, 8]]; Print["  dg_e^predicted (KK two-loop alone) = ", N[dgePredKK, 8]]; Print["  dg_e^predicted (KK + LR + KF) = ", N[dgePredAll, 8]]; Print["  dg_e^observed = g_meas - g_Schwinger = ", N[dgeObserved, 8]]; Print["  ratio dg_e^obs / dg_e^pred(KK alone) = ", N[ratioKK, 10]]; Print["  ratio dg_e^obs / dg_e^pred(KK+LR+KF) = ", N[ratioAll, 10]]; Print["  residual = dg_e^obs - dg_e^pred(KK+LR+KF) = ", N[residual, 6]]];


(* ============================================================ *)
(* Section 5.  Comparison + interpretation.                     *)
(* ============================================================ *)

Print["===== Section 5: interpretation ====="];
Print["The numerical observation: dg_e^obs = 3.515e-6 matches dg_e^pred(KK+LR+KF) to ~5e-12, i.e."];
Print["agreement at the 9-10 digit level, far below the framework's nominal 10^-6 precision floor."];
Print[""];
Print["This is CONSISTENT WITH (but does not uniquely identify) intentional Schwinger encoding."];
Print["Important caveat: the triangulated r_e is defined to reproduce g_e^measured.  Therefore the"];
Print["residual to g_e^Schwinger-one-loop is, by construction, exactly the all-orders QED correction"];
Print["beyond one-loop.  The match in Section 4 is a CONSISTENCY check, not an INDEPENDENT test of"];
Print["the framework's cutoff prescription."];
Print[""];
Print["For an independent test, observables NOT of the (g_s/-2)^n form are required (e.g. 1S-2S, muonic"];
Print["hydrogen Lamb shift, antiprotonic helium), and the framework's prediction formula for each must"];
Print["be derived.  This is queued for iteration 3 and beyond."];
Print[""];
Print["TENTATIVE OUTCOME-MATRIX BRANCH: B-conditional."];
Print["  B-confirmed iff Tepper indicates the DRQM I cutoff derivation produces the closed-form"];
Print["    (2 - alpha/(2 Pi))/(4 + alpha/Pi) as an identity (with higher-order QED coming from"];
Print["    framework loop corrections not yet derived in the campaign);"];
Print["  A-confirmed iff the framework's cutoff derivation produces directly the triangulated value"];
Print["    0.499420509912831 7 (requires the cutoff to encode all-orders QED at the derivation level,"];
Print["    which is suspicious without a derived loop structure);"];
Print["  D-confirmed iff the framework's cutoff derivation produces neither, in which case the"];
Print["    closed-form-Schwinger match is contingent."];


(* ============================================================ *)
(* Section 6.  Question for Tepper (queued via STATE.md).        *)
(* ============================================================ *)

Print["===== Section 6: question for Tepper queue ====="];
Print["Q.  In the DRQM I Section III.D derivation of the cutoff r_e/r_0 from the dual-Dirac"];
Print["    renormalisation prescription, does the derivation produce a closed-form expression"];
Print["    in alpha?  Specifically, is r_e/r_0 = (2 - alpha/(2 Pi))/(4 + alpha/Pi) the framework's"];
Print["    derived value (which would correspond to g_e = Schwinger one-loop, with higher-order QED"];
Print["    coming from elsewhere in the framework), or is it engineered to reproduce all-orders QED"];
Print["    directly (which would require the cutoff itself to encode the Karplus-Kroll, Laporta-"];
Print["    Remiddi, Kinoshita coefficients)?"];


(* ============================================================ *)
(* Section 8.  Muon cross-check -- breaking the back-fit         *)
(*             degeneracy via cross-particle test.               *)
(* ============================================================ *)

(* MOTIVATION.  Iter-2 found the electron's KK+LR+KF residual matches the closed-form Schwinger     *)
(* identification to ~10^-12, but iter 3 noted the agreement is back-fit-forced.  The cleanest test *)
(* that BREAKS the back-fit degeneracy is to apply the same framework formula -- DRQM I Eq. III.23: *)
(* g_mu^a = 2 (1 - 4 r_0^mu / (2 r_mu + r_0^mu)) -- to the muon and ask whether the per-particle    *)
(* back-fit cutoffs are equal (i.e. r_e/r_0^e = r_mu/r_0^mu), as a universal-closed-form-cutoff     *)
(* prediction would require.                                                                         *)
(*                                                                                                    *)
(* NB.  The published DRQM I does NOT make this universal-cutoff claim -- the paper writes the muon *)
(* and proton g-factor formulas in (III.23) but explicitly leaves r_mu, r_p as separate parameters  *)
(* without numerical values (verification doc lines 487-489).  This section tests a *stronger*      *)
(* prediction than the paper makes: IF the framework's cutoff prescription were universal-closed-   *)
(* form in alpha alone, then the same dimensionless cutoff r/r_0 must reproduce both a_e and a_mu.   *)
(* Since a_e and a_mu differ measurably (the muon has mass-dependent QED + hadronic + EW            *)
(* contributions), the universal-cutoff hypothesis is falsifiable.                                   *)

Module[{alpha, aPi, aSchwinger, aEexp, aMuExp, aMuSigma, rEbackFit, rMuBackFit, rClosedForm, deltaR, deltaA, falsifSigma}, alpha = SetPrecision[7.2973525693*^-3, 30]; aPi = alpha/Pi; aSchwinger = alpha/(2 Pi); aEexp = SetPrecision[1.15965218073*^-3, 18]; aMuExp = SetPrecision[1.16592059*^-3, 14]; aMuSigma = SetPrecision[2.2*^-10, 4]; rEbackFit = (2 - aEexp)/(4 + 2 aEexp); rMuBackFit = (2 - aMuExp)/(4 + 2 aMuExp); rClosedForm = (2 - aSchwinger)/(4 + 2 aSchwinger); deltaR = rEbackFit - rMuBackFit; deltaA = aMuExp - aEexp; falsifSigma = deltaR/(aMuSigma/4); Print["===== Section 8: muon cross-check ====="]; Print["  a_e^exp  (Fan 2023 / CODATA)        = ", N[aEexp, 14]]; Print["  a_mu^exp (FNAL Muon g-2 2023 final) = ", N[aMuExp, 12], "  +/- ", N[aMuSigma, 3]]; Print["  a_Schwinger one-loop = alpha/(2 Pi)  = ", N[aSchwinger, 14]]; Print["  delta_a := a_mu - a_e                = ", N[deltaA, 6], " = ", N[deltaA/aMuSigma, 5], " muon-sigmas"]; Print[""]; Print["  Framework's per-particle back-fit cutoffs (from g_r(r) = -2 - 2 a):"]; Print["    r_e/r_0^e   = (2 - a_e)/(4 + 2 a_e)    = ", N[rEbackFit, 18]]; Print["    r_mu/r_0^mu = (2 - a_mu)/(4 + 2 a_mu)  = ", N[rMuBackFit, 14]]; Print["    delta_r := r_e/r_0^e - r_mu/r_0^mu    = ", N[deltaR, 8]]; Print[""]; Print["  Schwinger one-loop universal closed-form: r/r_0 = ", N[rClosedForm, 18]]; Print[""]; Print["  Universal-cutoff falsification statistic:"]; Print["    delta_r / (sigma_{a_mu} / 4) = ", N[falsifSigma, 5], " sigma"]; Print["  -> UNIVERSAL-CLOSED-FORM-CUTOFF HYPOTHESIS RULED OUT AT >", N[falsifSigma/1000, 2], "k-sigma."]];

Print["===== Section 8: interpretation ====="];
Print["The per-particle back-fit cutoffs r_e/r_0^e and r_mu/r_0^mu differ by 3.13e-6, equivalent to"];
Print["~57,000 muon measurement-sigmas (using the FNAL 2023 a_mu uncertainty 2.2e-10).  Equivalently"];
Print["a_mu - a_e = 6.27e-6 differs from zero at ~28,500 muon-sigmas.  This rules out any prediction"];
Print["under which the dimensionless cutoff r/r_0 is universal across leptons."];
Print[""];
Print["The published DRQM I does NOT make a universal-cutoff prediction (the paper leaves r_mu free)."];
Print["So the section-8 test does not falsify the framework as written.  What it falsifies is the"];
Print["interpretation IF anyone were to claim (e.g. from issue #54's first-principles rederivation,"];
Print["if such a derivation produced a closed-form r_e/r_0 = (2 - alpha/(2 Pi))/(4 + alpha/Pi) without"];
Print["any particle-mass dependence) that the same closed form should hold for the muon.  Such a"];
Print["closed-form universal claim is now empirically forbidden."];
Print[""];
Print["The implication for outcome C-as-published is to STRENGTHEN it: the framework's cutoff is"];
Print["necessarily per-particle and empirically fit to each particle's measured anomaly.  The closed-"];
Print["form match seen for the electron is the back-fit cutoff that reproduces a_e^exp; since a_e^exp"];
Print["is dominated by Schwinger one-loop at the framework's nominal 10^-6 precision floor (the all-"];
Print["orders-QED-beyond-one-loop content of a_e is ~3.5e-6), the back-fit cutoff lies within ~10^-6"];
Print["of the closed-form Schwinger value automatically.  The same logic applied to the muon would"];
Print["give r_mu/r_0^mu = (2 - alpha/(2 Pi))/(4 + alpha/Pi) + delta_mu, with delta_mu propagated from"];
Print["a_mu^exp - alpha/(2 Pi) ~ 4.5e-6, of similar size to the electron's deviation.  Both are simply"];
Print["the all-orders QED + (for the muon, hadronic + EW) content of each particle's anomaly."];


(* ============================================================ *)
(* Section 9.  Human-acceptance section.                         *)
(* ============================================================ *)

(* <!-- TODO: human reviews and fills in -- confirms (a) the structural caveat in the header is the *)
(*      correct interpretation of the residual-tracks-KK observation (consistent with but not        *)
(*      independent of the closed-form identification), (b) the iter-3 finding that Section III.D    *)
(*      derives no closed-form r_e/r_0 in alpha (per Tepper's 2026-05-25 guidance that the published *)
(*      value was a uni-observable numerical search) settles the outcome-matrix at C-as-published,   *)
(*      (c) the iter-5 muon cross-check (Section 8) correctly frames the universal-closed-form-      *)
(*      cutoff hypothesis as a *stronger* prediction than the paper makes -- not a falsification of  *)
(*      the framework as written but a constraint on any future closed-form derivation (e.g. issue   *)
(*      #54) that would need to predict particle-mass-dependent cutoffs to be consistent with the    *)
(*      FNAL muon g-2 2023 measurement, (d) the (refined iter-3) Tepper question in STATE.md is the  *)
(*      right one to lift to a #66 comment, and (e) the verdict marker on Finding 2 stays at         *)
(*      [CHARACTERISED] with the structural argument now multi-step (back-fit-forcing + cross-       *)
(*      particle non-universality).  Note any choice the reviewer would have made differently. --> *)
