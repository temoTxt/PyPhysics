(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh14_P14_1.wl  --  Lienard-Wiechert potentials.                *)
(*  Companion to Ch14_Radiation_by_Moving_Charges.md Problem J3e-P14.1.    *)
(*                                                                         *)
(*  Confirms that the retarded-denominator s = R - R.v/c (classical)      *)
(*  equals s_PT = R - R.u/b (proper-time) under the velocity-duality      *)
(*  identity u/b = v/c.  The Lienard-Wiechert potentials are therefore   *)
(*  identical at the observable level in both formulations.               *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[rVec, vVec, uVec, capR, cVar, bVar];

classicalDenom = capR - (rVec . vVec)/cVar;
properTimeDenom = capR - (rVec . uVec)/bVar;

(* Apply velocity-duality identity: u/b = v/c, equivalently (rVec . uVec)/bVar = (rVec . vVec)/cVar *)
properTimeUnderDuality = properTimeDenom
   /. (rVec . uVec) -> (rVec . vVec) (bVar/cVar);

Print["s_classical = R - R . v/c = ", classicalDenom];
Print["s_proper-time (under duality) = ",
   FullSimplify[properTimeUnderDuality]];
Print["Match?  ",
   FullSimplify[properTimeUnderDuality - classicalDenom] === 0];

(* Expected:                                                               *)
(*   s_classical = capR - (rVec . vVec)/cVar                              *)
(*   s_proper-time = capR - (rVec . vVec)/cVar  (after duality)           *)
(*   Match?  True                                                          *)
