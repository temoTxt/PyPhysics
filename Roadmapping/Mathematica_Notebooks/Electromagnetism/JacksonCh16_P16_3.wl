(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh16_P16_3.wl  --  Pre-acceleration solution under step force. *)
(*  Companion to Ch16_Radiation_Damping.md Problem J3e-P16.3.              *)
(*                                                                         *)
(*  Solves the classical Abraham-Lorentz equation                          *)
(*    a - tau_0 da/dt = F_ext / m                                          *)
(*  for a step-function force F_ext(t) = F_0 theta(t), under the          *)
(*  asymptotic boundary condition a(infinity) = F_0/m.  The solution      *)
(*  exhibits PRE-ACCELERATION: a(t<0) = (F_0/m) exp(t/tau_0) > 0 BEFORE   *)
(*  the force is applied.                                                  *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[a, tt, tau0, ff0, mm, preAccSolNegT, preAccSolPosT];

(* For t < 0 (no external force): equation a - tau_0 da/dt = 0,           *)
(* with boundary condition a(0) = F_0/m (continuity with t > 0 solution). *)
preAccSolNegT = a[tt] /.
   DSolve[{a[tt] - tau0 D[a[tt], tt] == 0,
      a[0] == ff0/mm}, a[tt], tt][[1]];
Print["For t < 0:  a(t) = ", preAccSolNegT];
Print["  Expected: (F_0/m) exp(t/tau_0)"];

(* For t > 0 (constant force F_0): equation a - tau_0 da/dt = F_0/m,      *)
(* with boundary condition a(0) = F_0/m and a(infinity) = F_0/m (Newton). *)
preAccSolPosT = a[tt] /.
   DSolve[{a[tt] - tau0 D[a[tt], tt] == ff0/mm,
      a[0] == ff0/mm}, a[tt], tt][[1]];
Print["For t >= 0:  a(t) = ", preAccSolPosT];
Print["  Expected: F_0/m (Newton's law)"];

(* Pre-acceleration magnitude at various times *)
Print[""];
Print["Pre-acceleration magnitude (in units of F_0/m):"];
Do[
  Print["  t = ", n, " tau_0:  a(t) / (F_0/m) = ", Exp[n]],
  {n, {-3, -2, -1, -0.5, 0}}];

(* Conclusion: pre-acceleration is finite, with timescale tau_0 ~ 6e-24 s. *)
(* Observationally negligible at classical-physics timescales, but        *)
(* STRUCTURALLY UNPHYSICAL (violates causality).                           *)

Print[""];
Print["Proper-time treatment: first-order in partial_tau, standard causal"];
Print["IVP applies, no pre-acceleration.  Pathology is structurally absent."];
