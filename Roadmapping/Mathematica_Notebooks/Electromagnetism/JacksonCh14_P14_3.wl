(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh14_P14_3.wl  --  Non-relativistic Larmor + proper-time corr.  *)
(*  Companion to Ch14_Radiation_by_Moving_Charges.md Problem J3e-P14.3.    *)
(*                                                                         *)
(*  Verifies the proper-time correction to the Larmor formula:             *)
(*    P_PT = (2 e^2 a^2)/(3 b^3) = P_classical * [1 - (3/2)(u/c)^2 + ...]  *)
(*  with b = c sqrt(1 + u^2/c^2).                                          *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[ee, aa, cc, uu, bExpansion, larmorClassical, larmorPT, ratio,
   seriesExp];

bExpansion = cc Sqrt[1 + uu^2/cc^2];
larmorClassical = (2 ee^2 aa^2)/(3 cc^3);
larmorPT = (2 ee^2 aa^2)/(3 bExpansion^3);

Print["Classical Larmor: P = ", larmorClassical];
Print["Proper-time form: P_PT = ", FullSimplify[larmorPT]];

ratio = FullSimplify[larmorPT/larmorClassical];
Print["Ratio P_PT / P_classical = ", ratio];

seriesExp = Normal[Series[ratio, {uu, 0, 4}]];
Print["Series in u/c (through 4th order): ", seriesExp];

(* Quantitative estimates *)
Print[""];
Print["Fractional correction (1 - ratio) at various u/c:"];
Do[
  Print["  u/c = ", Indent /. N -> N,
     "  Delta P / P = ", N[1 - (ratio /. {uu -> cc r}) /. cc -> 1,
        Indent /. N -> N]],
  {r, {0.01, 0.06, 0.55}}];

(* For convenience, fixed correction formula: -(3/2)(u/c)^2 *)
Print[""];
Print["Leading-order correction: -(3/2)(u/c)^2"];
Print["  u/c = 0.01: ", N[-(3/2)(0.01)^2]];
Print["  u/c = 0.06: ", N[-(3/2)(0.06)^2]];
Print["  u/c = 0.55: ", N[-(3/2)(0.55)^2]];

(* Expected:                                                               *)
(*   Ratio = (1 + u^2/c^2)^(-3/2)                                          *)
(*   Series: 1 - 3 u^2/(2 c^2) + 15 u^4/(8 c^4) - ...                      *)
(*   Corrections: -1.5e-4 (chemistry), -5.4e-3 (1 keV), -0.45 (100 keV)   *)
