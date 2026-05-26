(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh11_P11_3.wl  --  Relativistic Doppler effect.                 *)
(*  Companion to Ch11_Special_Relativity.md Problem J3e-P11.3.             *)
(*                                                                         *)
(*  Verifies that the classical sqrt((1+beta)/(1-beta)) factor agrees     *)
(*  algebraically (after squaring) and numerically with the proper-time   *)
(*  forms (b+u)/c and c/(b-u), under b^2 = c^2 + u^2.                      *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[uu, cc, omegaPrime, bSubstituted, classicalDoppler,
   properTimeDoppler, altForm];

bSubstituted = Sqrt[cc^2 + uu^2];
classicalDoppler = omegaPrime Sqrt[(1 + uu/bSubstituted)/(1 - uu/bSubstituted)];
properTimeDoppler = omegaPrime cc/(bSubstituted - uu);
altForm = omegaPrime (bSubstituted + uu)/cc;

Print["Check 1 (squared): classical^2 - altForm^2 = ",
   FullSimplify[classicalDoppler^2 - altForm^2,
      Assumptions -> uu > 0 && cc > 0]];
Print["Check 2 (squared): properTime^2 - altForm^2 = ",
   FullSimplify[properTimeDoppler^2 - altForm^2,
      Assumptions -> uu > 0 && cc > 0]];

Print[""];
Print["Numerical check with uu = 1, cc = 1, omega' = 1:"];
Print["  classical  = ", N[classicalDoppler /. {uu -> 1, cc -> 1, omegaPrime -> 1}]];
Print["  properTime = ", N[properTimeDoppler /. {uu -> 1, cc -> 1, omegaPrime -> 1}]];
Print["  altForm    = ", N[altForm /. {uu -> 1, cc -> 1, omegaPrime -> 1}]];

(* Expected: Check 1 = 0, Check 2 = 0; all three numerical values ~ 2.414  *)
