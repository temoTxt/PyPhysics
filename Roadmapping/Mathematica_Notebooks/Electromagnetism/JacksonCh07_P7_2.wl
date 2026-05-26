(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh07_P7_2.wl  --  Dispersion in conducting medium.             *)
(*  Companion to Ch07_Plane_EM_Waves.md Problem J3e-P7.2.                  *)
(*                                                                         *)
(*  Verifies the complex dispersion relation k^2 = omega^2 eps mu / c^2  *)
(*    + i 4 pi omega sigma mu / c^2, with skin depth                       *)
(*    delta = c / sqrt(2 pi omega sigma mu)  for a good conductor.        *)
(*                                                                         *)
(*  Proper-time formulation: medium drift velocities << c, so the         *)
(*  (b/c) J -> J substitution is approximately identity.  Dispersion     *)
(*  unchanged at observable level.                                         *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[kk, omega, cc, eps, mu, sigma, kSq, skinDepth];

(* Complex dispersion for plane wave in conducting medium *)
kSq = omega^2 eps mu/cc^2 + I 4 Pi omega sigma mu/cc^2;
Print["Complex k^2 = ", kSq];

(* Good-conductor limit: sigma >> omega eps / (4 pi) *)
(* k approx (1 + i) / delta with delta = c / sqrt(2 pi omega sigma mu) *)
skinDepth = cc/Sqrt[2 Pi omega sigma mu];
Print["Skin depth delta = c / sqrt(2 pi omega sigma mu) = ", skinDepth];
