(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh07_P7_3.wl  --  Fresnel reflection / refraction coefficients. *)
(*  Companion to Ch07_Plane_EM_Waves.md Problem J3e-P7.3.                  *)
(*                                                                         *)
(*  Encodes the Fresnel reflection (r_s, r_p) and transmission             *)
(*  (t_s, t_p) coefficients for a plane EM wave at a planar boundary     *)
(*  between two linear media (refractive indices n1, n2).                  *)
(*                                                                         *)
(*  Proper-time formulation: identical (boundary conditions depend on    *)
(*  field continuity, not on source motion).                              *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[n1, n2, thetaI, thetaT, rS, tS, rP, tP];

(* Snell's law relates thetaI and thetaT: n1 sin thetaI = n2 sin thetaT *)

rS = (n1 Cos[thetaI] - n2 Cos[thetaT])/(n1 Cos[thetaI] + n2 Cos[thetaT]);
tS = (2 n1 Cos[thetaI])/(n1 Cos[thetaI] + n2 Cos[thetaT]);
rP = (n2 Cos[thetaI] - n1 Cos[thetaT])/(n2 Cos[thetaI] + n1 Cos[thetaT]);
tP = (2 n1 Cos[thetaI])/(n2 Cos[thetaI] + n1 Cos[thetaT]);

Print["s-polarisation: r_s = ", rS, "   t_s = ", tS];
Print["p-polarisation: r_p = ", rP, "   t_p = ", tP];

(* Sanity check at normal incidence (thetaI = 0, thetaT = 0): *)
rSAtNormalIncidence = rS /. {thetaI -> 0, thetaT -> 0};
Print["At normal incidence: r_s = ", FullSimplify[rSAtNormalIncidence]];
Print["  Expected: (n1 - n2)/(n1 + n2)"];
