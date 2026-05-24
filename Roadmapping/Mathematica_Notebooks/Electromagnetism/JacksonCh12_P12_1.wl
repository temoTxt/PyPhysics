(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh12_P12_1.wl  --  Free relativistic Lagrangian and action.    *)
(*  Companion to Ch12_Relativistic_Dynamics.md Problem J3e-P12.1.          *)
(*                                                                         *)
(*  Verifies that the classical Lagrangian L = -mc^2/gamma and the         *)
(*  proper-time Lagrangian density L_tau = -mc^2 give the same action     *)
(*  for a worldline of fixed endpoints, via the substitution dt = gamma   *)
(*  dtau.                                                                  *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[ww, cc, mm, gammaW, lClassical, dSdtau, bSubst, lInProperTimeVars,
   kKinetic];

(* Classical Lagrangian L = -mc^2 / gamma *)
gammaW = 1/Sqrt[1 - ww^2/cc^2];
lClassical = -mm cc^2/gammaW;

(* dt/dtau = gamma, so dS/dtau = L * gamma = -mc^2 *)
dSdtau = FullSimplify[lClassical gammaW];
Print["dS/dtau = L_classical * gamma = ", dSdtau];
Print["  Equals -mc^2?  ", dSdtau === -mm cc^2];

(* In proper-time variables: gamma = b/c, so 1/gamma = c/b, and
   L_classical = -mc^2 * (c/b) = -mc^3/b *)
bSubst = cc/Sqrt[1 - ww^2/cc^2];
lInProperTimeVars = -mm cc^2 cc/bSubst;
Print[""];
Print["L_classical in proper-time variables (= -mc^3/b) = ",
   FullSimplify[lInProperTimeVars]];
Print["  Matches L_classical?  ",
   FullSimplify[lInProperTimeVars - lClassical] === 0];

(* Gill's DRQM I proper-time kinetic energy K = (1/2) m b^2 *)
kKinetic = (1/2) mm bSubst^2;
Print[""];
Print["Gill's K = (1/2) m b^2 in observer-time variables = ",
   FullSimplify[kKinetic]];
Print["  At w -> 0: K -> ", FullSimplify[kKinetic /. ww -> 0]];
Print["  At w -> 0 expected: (1/2) m c^2"];

(* Expected output:                                                        *)
(*   dS/dtau = L_classical * gamma = -cc^2 mm                              *)
(*     Equals -mc^2?  True                                                 *)
(*   L_classical in proper-time variables = -mc^2 Sqrt[1 - ww^2/cc^2]      *)
(*     Matches L_classical?  True                                          *)
(*   Gill's K = m c^4 / (2 (cc^2 - ww^2))                                  *)
(*     At w -> 0: (1/2) m c^2                                              *)
