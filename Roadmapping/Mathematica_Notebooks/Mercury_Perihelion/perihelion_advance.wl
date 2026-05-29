(* ::Package:: *)

(* perihelion_advance.wl  --  Companion notebook to                              *)
(*   Roadmapping/Mercury_Perihelion/ (issue #81)                                 *)
(*                                                                               *)
(* Mercury perihelion advance under the Gill "Relativistic Newtonian Theory"     *)
(* (a.k.a. Dual Newtonian Theory) force law, compared with general relativity.   *)
(*                                                                               *)
(* HEADLINE RESULT                                                               *)
(*   The framework's force law  a = -(GM/r^2)(1 - GM/(c^2 r)) e_r  produces a     *)
(*   perihelion advance of EXACTLY -1/6 of the GR value:                         *)
(*       Delta phi_framework = -pi G M / (c^2 a (1-e^2))   per orbit             *)
(*       Delta phi_GR        = +6 pi G M / (c^2 a (1-e^2)) per orbit             *)
(*   i.e. OPPOSITE SIGN (retrograde, not prograde) and ONE-SIXTH the magnitude.  *)
(*   For Mercury: framework -7.17"/century  vs  GR +42.99"/century (= observed). *)
(*                                                                               *)
(* TWO METHODOLOGICAL POINTS THE NOTEBOOK MAKES EXPLICIT                          *)
(*   (1) The paper's section 2 computes the advance from a CIRCULAR orbit         *)
(*       (the omega_d . T_0 - 2pi construction). A circular orbit has no          *)
(*       perihelion. Cells 3-5 redo it the standard way: the apsidal angle of    *)
(*       an ECCENTRIC orbit. This is the change from "circular" to "eccentric".  *)
(*   (2) The framework's correction is a REPULSIVE 1/r^2 potential perturbation;  *)
(*       GR's (and the classic velocity-dependent-mass estimate's) correction is *)
(*       ATTRACTIVE. That sign difference is why the framework precesses the      *)
(*       WRONG WAY.                                                              *)
(*                                                                               *)
(* SYMBOL CONVENTIONS (avoid Wolfram built-in collisions)                         *)
(*   GG=G, MM=M (Sun), mm=m (Mercury), cc=c, aa=a, ecc=e, potV=potential.        *)
(*   ee/e avoided (Euler), V avoided (Vanadium) per repo convention.             *)
(*                                                                               *)
(* Each cell is a single logical line (re-runnable in Mathematica or via         *)
(* wolframscript, and against the Wolfram MCP). Confirmed 2026-05-29.            *)
(*                                                                               *)
(* Author: Trey Morris with Claude Opus 4.7.                                     *)


(* ========================================================================== *)
(* Cell 1 - Verify Eq. (h3): a = -(grad V / m)(1 + V/(m c^2)) from the kernel  *)
(*   K = pi^2/(2m) + m c^2 + V + V^2/(2 m c^2).                                *)
(* The radial Hamilton equation p-dot = -dK/dx gives the modified force.       *)
(* ========================================================================== *)

ClearAll[r, cc, mm, potV]; KK[r_] := potV[r] + potV[r]^2/(2 mm cc^2); Print["Cell 1 - Eq. (h3) residual (should be 0): ", FullSimplify[D[KK[r], r] - D[potV[r], r] (1 + potV[r]/(mm cc^2))]]

(* Expected: 0  (confirmed) *)


(* ========================================================================== *)
(* Cell 2 - Verify Eq. (h4): substitute V = -G M m / r (gravity) into (h3)     *)
(*   to recover  a = -(G M / r^2)(1 - G M/(c^2 r)) e_r.                        *)
(* ========================================================================== *)

ClearAll[r, cc, mm, GG, MM, potV]; potV[r_] := -GG MM mm/r; Print["Cell 2 - Eq. (h4) residual (should be 0): ", FullSimplify[(-D[potV[r], r]/mm (1 + potV[r]/(mm cc^2))) - (-(GG MM/r^2) (1 - GG MM/(cc^2 r)))]]

(* Expected: 0  (confirmed) *)


(* ========================================================================== *)
(* Cell 3 - The change from CIRCULAR to ECCENTRIC, and the EXACT -1/6.          *)
(*                                                                            *)
(* The force law's correction term is a +1/r^2 POTENTIAL perturbation:        *)
(*   Phi(r) = -G(M+m)/r + B/r^2,  B = G^2(M^2+m^2)/(2 c^2)  > 0 (repulsive).   *)
(* A 1/r^2 perturbation merges exactly with the centrifugal term, rescaling   *)
(* the angular momentum: ell^2 -> ell^2 + 2B. The apsidal angle is pi ell/    *)
(* Sqrt[ell^2+2B], so the precession per orbit is (to linear order in B):     *)
(*   Delta phi = -2 pi B / ell^2,  with ell^2 = G(M+m) a (1-e^2).             *)
(* Compare GR: Delta phi_GR = +6 pi G M /(c^2 a (1-e^2)).                      *)
(* ========================================================================== *)

ClearAll[BB, ellSq, GG, MM, mm, cc, aa, ecc]; dFramework = -2 Pi BB/ellSq /. {BB -> GG^2 (MM^2 + mm^2)/(2 cc^2), ellSq -> GG (MM + mm) aa (1 - ecc^2)}; dGR = 6 Pi GG MM/(cc^2 aa (1 - ecc^2)); Print["Cell 3 - framework precession per orbit: ", FullSimplify[dFramework]]; Print["Cell 3 - ratio framework/GR (general): ", FullSimplify[dFramework/dGR]]; Print["Cell 3 - ratio framework/GR (m << M): ", FullSimplify[dFramework/dGR /. mm -> 0]]

(* Expected: ratio = -(1/6)(M^2+m^2)/(M(M+m)) -> -1/6 as m->0.  Note the MINUS. *)


(* ========================================================================== *)
(* Cell 4 - Numerical values for Mercury.                                      *)
(*   Constants: NASA fact sheet + CODATA 2018 + IAU 2015.                     *)
(* ========================================================================== *)

ClearAll[Msun, mMerc, aMerc, eMerc, Torb, GG, cc, orbitsPerCentury, radToArcsec, dGRnum, dFWnum]; Msun = 1.98892*^30; mMerc = 3.3011*^23; aMerc = 5.7909*^10; eMerc = 0.20563; Torb = 87.969*86400.; GG = 6.6743*^-11; cc = 299792458.; orbitsPerCentury = 100*365.25*86400./Torb; radToArcsec = 180/Pi*3600.; dGRnum = 6 Pi GG Msun/(cc^2 aMerc (1 - eMerc^2)); dFWnum = -Pi GG (Msun^2 + mMerc^2)/(cc^2 aMerc (1 - eMerc^2) (Msun + mMerc)); Print["Cell 4 - GR        : ", N[dGRnum orbitsPerCentury radToArcsec, 6], " arcsec/century (matches observed ~43)"]; Print["Cell 4 - framework : ", N[dFWnum orbitsPerCentury radToArcsec, 6], " arcsec/century (negative = retrograde)"]; Print["Cell 4 - ratio     : ", N[dFWnum/dGRnum, 8]]

(* Expected: GR +42.99, framework -7.17, ratio -0.166667 *)


(* ========================================================================== *)
(* Cell 5 - INDEPENDENT CHECK: numerically integrate the orbit under the       *)
(*   framework force law and measure the perihelion shift DIRECTLY (no apsidal *)
(*   formula). Uses an exaggerated G M/(c^2 a) = 0.01 so the effect is large.  *)
(*   Confirms the sign is RETROGRADE and magnitude ~ -1/6 of GR.              *)
(* ========================================================================== *)

ClearAll[x, y, t, GMn, csqn, sol, xf, yf, rf, aao, eeo, rp, vp, Tn, p1, p2, angf]; GMn = 1.; csqn = 100.; aao = 1.; eeo = 0.3; rp = aao (1 - eeo); vp = Sqrt[GMn (1 + eeo)/(aao (1 - eeo))]; fxx[xx_, yy_] := With[{rr = Sqrt[xx^2 + yy^2]}, -(GMn/rr^2) (1 - GMn/(csqn rr)) xx/rr]; fyy[xx_, yy_] := With[{rr = Sqrt[xx^2 + yy^2]}, -(GMn/rr^2) (1 - GMn/(csqn rr)) yy/rr]; sol = First@NDSolve[{x''[t] == fxx[x[t], y[t]], y''[t] == fyy[x[t], y[t]], x[0] == rp, y[0] == 0, x'[0] == 0, y'[0] == vp}, {x, y}, {t, 0, 60}, MaxSteps -> 10^7]; xf[t_] := x[t] /. sol; yf[t_] := y[t] /. sol; rf[t_] := Sqrt[xf[t]^2 + yf[t]^2]; Tn = N[2 Pi]; p1 = tt /. FindRoot[rf'[tt] == 0, {tt, Tn}]; p2 = tt /. FindRoot[rf'[tt] == 0, {tt, 2 Tn}]; angf[tp_] := ArcTan[xf[tp], yf[tp]]; Print["Cell 5 - measured precession/orbit (rad): ", angf[p2] - angf[p1]]; Print["Cell 5 - apsidal-angle theory  -Pi/(csqn (1-e^2)): ", N[-Pi/(csqn (1 - eeo^2))]]; Print["Cell 5 - sign: ", If[angf[p2] - angf[p1] < 0, "RETROGRADE (wrong way vs GR/observation)", "prograde"]]

(* Expected: measured ~ -0.0342 rad/orbit, theory -0.0345; RETROGRADE.        *)
(* (The ~0.8% gap is the finite-eps higher-order term; for Mercury eps~1e-8.) *)


(* ========================================================================== *)
(* Cell 6 - Why an N-body planetary treatment cannot rescue the result.        *)
(*   The framework's relativistic factor G M_body/(c^2 r) scales with the      *)
(*   ATTRACTING mass. Jupiter (the largest perturber) contributes ~7e-5 of     *)
(*   the Sun-Mercury term; all planets summed ~1e-3 arcsec/century -- three    *)
(*   to four orders of magnitude too small to bridge -7 -> +43.               *)
(* ========================================================================== *)

ClearAll[GG, cc, Msun, aMerc, Mjup, aJup]; GG = 6.6743*^-11; cc = 299792458.; Msun = 1.98892*^30; aMerc = 5.7909*^10; Mjup = 1.898*^27; aJup = 7.785*^11; Print["Cell 6 - Sun-Mercury  G M/(c^2 r): ", ScientificForm[N[GG Msun/(cc^2 aMerc), 4]]]; Print["Cell 6 - Jupiter-Mercury G M/(c^2 r): ", ScientificForm[N[GG Mjup/(cc^2 aJup), 4]]]; Print["Cell 6 - Jupiter/Sun ratio: ", ScientificForm[N[(Mjup/aJup)/(Msun/aMerc), 4]]]

(* Expected: Sun 2.55e-8, Jupiter 1.81e-12, ratio 7.1e-5 *)


(* ========================================================================== *)
(* SUMMARY                                                                     *)
(*   - Eqs. (h3), (h4) verify exactly (Cells 1-2).                            *)
(*   - Proper (eccentric-orbit) perihelion advance = -1/6 of GR, EXACT        *)
(*     (Cell 3), confirmed by direct orbit integration (Cell 5).              *)
(*   - Numerically: framework -7.17"/cy vs GR/observed +42.99"/cy (Cell 4).   *)
(*   - The minus sign is RETROGRADE: the framework's V^2/(2 m c^2) term is a   *)
(*     repulsive 1/r^2 perturbation, opposite to GR's attractive correction.  *)
(*   - No N-body refinement closes the gap (Cell 6).                          *)
(* ========================================================================== *)
