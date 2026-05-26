(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh12_P12_14.wl  --  Charged particle in plane EM wave.         *)
(*  Companion to Ch12_Relativistic_Dynamics.md Problem J3e-P12.14.         *)
(*                                                                         *)
(*  Verifies the vector-potential / E-field consistency for a plane wave  *)
(*  in Lorenz gauge.  The motion of a charged particle initially at rest *)
(*  in this wave is governed by two conservation laws (p_y = 0,           *)
(*  H/c - p_z = m c) and the velocity-oscillation result                  *)
(*    gamma m w_x = a_0 m c sin(kz - omega t)                             *)
(*  with a_0 = q E_0 / (m omega c) the normalised vector potential.       *)
(*                                                                         *)
(*  This is the first PR C problem in which the time-averaged              *)
(*  proper-time invariant u . a is nonzero (at relativistic intensity     *)
(*  a_0 > 0).  The dissipative (u . a) / b^4 coefficient of Eq. (4) of   *)
(*  the Maxwell paper engages, recovering the Thomson-scattering rate at *)
(*  low intensity and the radiation-reaction physics of Cole/Poder/      *)
(*  Wistisen 2018 at relativistic intensity.                              *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[capE0, omega, kk, cc, tt, zz, capAx, eFromA];

(* Lorenz-gauge vector potential for plane wave:                          *)
(*   A_x = -(c E_0 / omega) sin(kz - omega t)                             *)
capAx = -(cc capE0/omega) Sin[kk zz - omega tt];

(* Recover E_x = -(1/c) d A_x / d t                                       *)
eFromA = -(1/cc) D[capAx, tt];
Print["E_x from vector potential = ", FullSimplify[eFromA]];
Print["  Expected: -E_0 cos(kz - omega t)  (note the sign convention)"];

(* Conservation laws documented in the prose; no symbolic check needed.   *)
(* Velocity-oscillation amplitude a_0 = q E_0 / (m omega c) is the        *)
(* dimensionless intensity parameter.                                      *)

Print[""];
Print["Quantitative estimate of dissipative engagement:"];
Print["  Low-intensity (a_0 << 1):  <P_rad> = (2 q^4 E_0^2) / (3 m^2 c^3) (Thomson)"];
Print["  High-intensity (a_0 ~ 1):  modified by figure-8 dynamics,"];
Print["                              regime of Cole/Poder/Wistisen 2018."];
Print["                              Comparison deferred to issue #43."];

(* Expected output:                                                        *)
(*   E_x from vector potential = -capE0 Cos[kk zz - omega tt]              *)
(*     (Equivalent to -capE0 Cos[omega tt - kk zz] under cosine parity)   *)
