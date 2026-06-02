(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh11_P11_6.wl  --  Lorentz boost of E and B fields.            *)
(*  Companion to Ch11_Special_Relativity.md Problem J3e-P11.6.             *)
(*                                                                         *)
(*  Verifies the two Lorentz invariants of the EM field:                  *)
(*    I_1 = E^2 - B^2  (scalar)                                            *)
(*    I_2 = E . B      (pseudoscalar)                                      *)
(*  are preserved under the boost, in both classical (gamma, beta)        *)
(*  parametrisation and proper-time (b, u) parametrisation.                *)
(*                                                                         *)
(*  Structural observation: the proper-time form is LINEAR in (b, u)      *)
(*  whereas the classical form is non-linear in v (via gamma).            *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[capE1, capE2, capE3, capB1, capB2, capB3, vv, cc, gammaV, beta,
   eP, bP, uu, bb, ePProperTime, bPProperTime];

(* ----- Classical (gamma, beta) form ----- *)
gammaV = 1/Sqrt[1 - vv^2/cc^2];
beta = vv/cc;
eP = {capE1, gammaV (capE2 - beta capB3), gammaV (capE3 + beta capB2)};
bP = {capB1, gammaV (capB2 + beta capE3), gammaV (capB3 - beta capE2)};

diffI1 = FullSimplify[
   (eP . eP - bP . bP)
      - ((capE1^2 + capE2^2 + capE3^2)
         - (capB1^2 + capB2^2 + capB3^2)),
   Assumptions -> 0 < vv < cc];
diffI2 = FullSimplify[
   eP . bP - (capE1 capB1 + capE2 capB2 + capE3 capB3),
   Assumptions -> 0 < vv < cc];

Print["Classical form:"];
Print["  (E'^2 - B'^2) - (E^2 - B^2) = ", diffI1];
Print["  (E' . B') - (E . B) = ", diffI2];

(* ----- Proper-time (b, u) form ----- *)
ePProperTime = {capE1,
   (1/cc) (bb capE2 - uu capB3),
   (1/cc) (bb capE3 + uu capB2)};
bPProperTime = {capB1,
   (1/cc) (bb capB2 + uu capE3),
   (1/cc) (bb capB3 - uu capE2)};

diffI1PT = FullSimplify[
   (ePProperTime . ePProperTime - bPProperTime . bPProperTime)
      - ((capE1^2 + capE2^2 + capE3^2)
         - (capB1^2 + capB2^2 + capB3^2))
   /. bb -> Sqrt[cc^2 + uu^2],
   Assumptions -> uu > 0 && cc > 0];
diffI2PT = FullSimplify[
   ePProperTime . bPProperTime
      - (capE1 capB1 + capE2 capB2 + capE3 capB3)
   /. bb -> Sqrt[cc^2 + uu^2],
   Assumptions -> uu > 0 && cc > 0];

Print["Proper-time form (with b = Sqrt[c^2 + u^2]):"];
Print["  (E'^2 - B'^2) - (E^2 - B^2) = ", diffI1PT];
Print["  (E' . B') - (E . B) = ", diffI2PT];

(* Expected output:                                                        *)
(*   Classical form:                                                       *)
(*     (E'^2 - B'^2) - (E^2 - B^2) = 0                                    *)
(*     (E' . B') - (E . B) = 0                                            *)
(*   Proper-time form (with b = Sqrt[c^2 + u^2]):                         *)
(*     (E'^2 - B'^2) - (E^2 - B^2) = 0                                    *)
(*     (E' . B') - (E . B) = 0                                            *)
