(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh14_P14_6.wl  --  Bremsstrahlung from linear deceleration.    *)
(*  Companion to Ch14_Radiation_by_Moving_Charges.md Problem J3e-P14.6.    *)
(*                                                                         *)
(*  Linear motion: u and a along same line, so u . a = +/- |u| |a|        *)
(*  is finite and the third term of Eq. (7) of the Maxwell paper engages. *)
(*                                                                         *)
(*  Computes the magnitude ratio of the third term to the classical       *)
(*  acceleration-field contribution:                                       *)
(*    |E_3rd| / |E_acc_classical| ~ (u/c)^2  at non-relativistic intensity *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[ee, aa, vv, cc, gammaW, lienardLinear];

gammaW = 1/Sqrt[1 - vv^2/cc^2];

(* Classical Lienard for LINEAR motion (w parallel to a):                 *)
(*   P = (2/3)(e^2/c^3) gamma^6 a^2  (gamma^6 scaling, vs gamma^4 for circular) *)
lienardLinear = (2 ee^2/(3 cc^3)) gammaW^6 aa^2;
Print["Classical Lienard for linear motion: P = ", lienardLinear];
Print["  Scales as gamma^6 (vs gamma^4 for circular motion in J3e-P14.5)"];

(* Order-of-magnitude estimate of third-term / classical ratio.            *)
(* Classical Larmor acceleration field at distance r:  |E_acc| ~ e a / (c^2 r) *)
(* Proper-time third-term field:  |E_3rd| ~ e |u| |u.a| / (b^4 r) ~ e u^2 a / (c^4 r) *)
(* Ratio: |E_3rd| / |E_acc| ~ u^2 / c^2 *)

Print[""];
Print["Ratio |E_3rd| / |E_acc_classical| at various u/c:"];
Do[
   Print["  u/c = ", N[r], "  ratio ~ (u/c)^2 = ", N[r^2]],
   {r, {0.01, 0.1, 0.5, 0.9}}];

(* At MeV electron energies (u/c ~ 0.94), the ratio is ~ 0.88 -- order unity. *)
(* At chemistry-scale (u/c ~ 0.01), ratio ~ 10^-4 -- below precision floor. *)
(* The MeV regime is the natural experimental target for the third-term test. *)

Print[""];
Print["Bremsstrahlung at MeV electron energies (u/c ~ 0.94):"];
Print["  Third-term ratio ~ 0.88  (order unity)"];
Print["  Precision bremsstrahlung-spectrum measurements achieve ~1-2% accuracy"];
Print["  -> Proper-time prediction is at the boundary of distinguishability"];

(* Expected output:                                                        *)
(*   Classical Lienard for linear motion: gamma^6 a^2 (2 e^2/3 c^3)         *)
(*   Ratio (u/c)^2:                                                       *)
(*     u/c = 0.01: 0.0001                                                  *)
(*     u/c = 0.1:  0.01                                                    *)
(*     u/c = 0.5:  0.25                                                    *)
(*     u/c = 0.9:  0.81                                                    *)
