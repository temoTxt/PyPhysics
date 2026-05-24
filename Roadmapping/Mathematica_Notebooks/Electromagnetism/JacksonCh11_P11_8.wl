(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh11_P11_8.wl  --  Lorentz invariants of the EM field.         *)
(*  Companion to Ch11_Special_Relativity.md Problem J3e-P11.8.             *)
(*                                                                         *)
(*  Identifies F^munu F_munu and F^munu dual(F)_munu as the two algebraic *)
(*  Lorentz invariants of the EM field tensor, in Gaussian units with     *)
(*  metric (+, -, -, -).                                                   *)
(*                                                                         *)
(*    F^munu F_munu        = 2 (B^2 - E^2) = -2 I_1                        *)
(*    F^munu dual(F)_munu  = -4 (E . B)    = -4 I_2                        *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[capE1, capE2, capE3, capB1, capB2, capB3, fUpper, fLower,
   fDualUpper, fDualLower, eta, i1Check, i2Check];

fUpper = {
   {0, capE1, capE2, capE3},
   {-capE1, 0, -capB3, capB2},
   {-capE2, capB3, 0, -capB1},
   {-capE3, -capB2, capB1, 0}};

eta = DiagonalMatrix[{1, -1, -1, -1}];
fLower = eta . fUpper . eta;

i1Check = Sum[
   fUpper[[m, n]] fLower[[m, n]],
   {m, 1, 4}, {n, 1, 4}];
Print["F^munu F_munu = ", FullSimplify[i1Check]];
Print["  Matches 2 (B^2 - E^2)?  ",
   FullSimplify[i1Check
      - 2 ((capB1^2 + capB2^2 + capB3^2)
         - (capE1^2 + capE2^2 + capE3^2))] === 0];

(* Dual tensor: dual(F)^{0i} = B^i, dual(F)^{ij} = epsilon^{ijk} E_k.     *)
fDualUpper = {
   {0, capB1, capB2, capB3},
   {-capB1, 0, capE3, -capE2},
   {-capB2, -capE3, 0, capE1},
   {-capB3, capE2, -capE1, 0}};
fDualLower = eta . fDualUpper . eta;

i2Check = Sum[
   fUpper[[m, n]] fDualLower[[m, n]],
   {m, 1, 4}, {n, 1, 4}];
Print["F^munu dual(F)_munu = ", FullSimplify[i2Check]];
Print["  Matches -4 (E . B)?  ",
   FullSimplify[i2Check
      - (-4 (capE1 capB1 + capE2 capB2 + capE3 capB3))] === 0];

(* Expected output:                                                        *)
(*   F^munu F_munu = 2 (B^2 - E^2)                                         *)
(*     Matches 2(B^2 - E^2)?  True                                         *)
(*   F^munu dual(F)_munu = -4 (E . B)                                     *)
(*     Matches -4(E . B)?  True                                            *)
