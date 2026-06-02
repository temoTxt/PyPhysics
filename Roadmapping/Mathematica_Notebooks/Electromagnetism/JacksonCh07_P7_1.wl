(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh07_P7_1.wl  --  Plane EM wave in vacuum.                     *)
(*  Companion to Ch07_Plane_EM_Waves.md Problem J3e-P7.1.                  *)
(*                                                                         *)
(*  Verifies the wave-equation dispersion omega = k c for a plane wave    *)
(*  in vacuum.  Proper-time formulation reduces identically since there   *)
(*  is no source (u = 0 => b = c).                                         *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[capE0, kk, omega, cc, zz, tt, eField, residual];

eField = capE0 Cos[kk zz - omega tt];

residual = FullSimplify[
   (1/cc^2) D[eField, {tt, 2}] - D[eField, {zz, 2}]];

Print["Wave-equation residual: ", residual];
Print["  Solved when omega^2 = k^2 c^2  =>  omega = k c"];
