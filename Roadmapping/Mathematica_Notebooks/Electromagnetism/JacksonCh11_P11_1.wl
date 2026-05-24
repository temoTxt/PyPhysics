(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh11_P11_1.wl                                                   *)
(*                                                                         *)
(*  Companion notebook to Roadmapping/Electromagnetism/Jackson/            *)
(*  Ch11_Special_Relativity.md  ->  Problem J3e-P11.1.                     *)
(*                                                                         *)
(*  Relativistic velocity addition.  Verifies (i) the Minkowski-length    *)
(*  invariance of the 4-velocity (b^2 - u^2 = c^2 preserved under Lorentz *)
(*  boost), and (ii) the equivalence between the proper-velocity boost   *)
(*  (Eq. 11 of the Maxwell paper) and the classical observer-time         *)
(*  velocity-addition formula.                                             *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(*  Runnable independent of the Wolfram MCP.                               *)
(* ----------------------------------------------------------------------- *)

ClearAll[uMoving, vRel, cc, gammaV, bMovingExpr, uLab, bLab, labInvariant,
   wLab, wMovingFromProper, wLabClassical];

(* Set up: particle has proper-velocity uMoving in a frame S' that itself  *)
(* moves with lab-velocity vRel relative to lab frame S.  We transform    *)
(* (b, u) into the lab frame via the proper-time-group boost of Eq. (11). *)
gammaV = 1/Sqrt[1 - vRel^2/cc^2];
bMovingExpr = Sqrt[cc^2 + uMoving^2];

(* Eq. (11) and its complement: 4-velocity boost in (b, u) space.          *)
uLab = gammaV (uMoving + vRel bMovingExpr/cc);
bLab = gammaV (bMovingExpr + vRel uMoving/cc);

(* Check 1: Minkowski-length invariant b^2 - u^2 - c^2 should be zero.    *)
labInvariant = FullSimplify[bLab^2 - uLab^2 - cc^2,
   Assumptions -> 0 < vRel < cc && uMoving > 0];
Print["Lab-frame invariant b^2 - u^2 - c^2 = ", labInvariant];

(* Check 2: lab-frame observer-time velocity from proper-velocity         *)
(* composition vs from classical velocity-addition formula.                *)
wLab = cc uLab/bLab;
wMovingFromProper = cc uMoving/bMovingExpr;
wLabClassical = (wMovingFromProper + vRel)/(1 + wMovingFromProper vRel/cc^2);

Print["w_lab (from proper-velocity boost) = ",
   FullSimplify[wLab, Assumptions -> 0 < vRel < cc && uMoving > 0]];
Print["w_lab (from classical addition) = ",
   FullSimplify[wLabClassical,
      Assumptions -> 0 < vRel < cc && uMoving > 0]];
Print["Match?  ",
   FullSimplify[wLab - wLabClassical,
      Assumptions -> 0 < vRel < cc && uMoving > 0] === 0];

(* Expected output:                                                        *)
(*   Lab-frame invariant b^2 - u^2 - c^2 = 0                              *)
(*   w_lab (from proper-velocity boost) = cc (cc uMoving + bMovingExpr vRel) /  *)
(*     (cc bMovingExpr + uMoving vRel)                                    *)
(*   w_lab (from classical addition) = (same as above)                    *)
(*   Match?  True                                                         *)
