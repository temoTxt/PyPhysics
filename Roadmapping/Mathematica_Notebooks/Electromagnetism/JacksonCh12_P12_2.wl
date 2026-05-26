(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh12_P12_2.wl  --  Cyclotron motion in uniform B.              *)
(*  Companion to Ch12_Relativistic_Dynamics.md Problem J3e-P12.2.          *)
(*                                                                         *)
(*  Verifies that the proper-time cyclotron angular frequency              *)
(*  omega_tau = q B / (m c)  is gamma-INDEPENDENT, and that the lab-frame *)
(*  frequency omega_lab = omega_tau / gamma matches the classical result *)
(*  omega_c = q B / (gamma m c) exactly.                                   *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[qq, mm, capB0, cc, vv, gammaW, omegaClassical, omegaTau,
   omegaLabFromProperTime];

gammaW = 1/Sqrt[1 - vv^2/cc^2];

(* Classical cyclotron frequency: omega_c = q B / (gamma m c)              *)
omegaClassical = qq capB0/(mm gammaW cc);

(* Proper-time cyclotron frequency: omega_tau = q B / (m c)                *)
(* This follows from m du/dtau = (q/c) u x B, in which no gamma appears.   *)
omegaTau = qq capB0/(mm cc);

(* Lab-frame angular frequency derived from omega_tau via time dilation:   *)
(*   dphi/dt = dphi/dtau * dtau/dt = omega_tau * (1/gamma) = omega_tau/gamma *)
omegaLabFromProperTime = omegaTau/gammaW;

Print["Classical cyclotron frequency omega_c = ", FullSimplify[omegaClassical]];
Print["Proper-time cyclotron frequency omega_tau = ", omegaTau];
Print["  (gamma-INDEPENDENT)"];
Print["Lab-frame frequency from omega_tau / gamma = ",
   FullSimplify[omegaLabFromProperTime]];
Print["Match classical?  ",
   FullSimplify[omegaLabFromProperTime - omegaClassical] === 0];

(* Numerical sanity check: take v = 0.6 c (gamma = 1.25), B = 1, q = 1, m = 1, c = 1 *)
sample = {qq -> 1, mm -> 1, capB0 -> 1, cc -> 1, vv -> 6/10};
Print[""];
Print["Numerical check at v = 0.6 c (gamma = 5/4):"];
Print["  omega_classical = ", N[omegaClassical /. sample]];
Print["  omega_tau       = ", N[omegaTau /. sample]];
Print["  omega_tau/gamma = ", N[omegaLabFromProperTime /. sample]];

(* Expected:                                                               *)
(*   omega_c = 4/5 (= 0.8)                                                 *)
(*   omega_tau = 1                                                         *)
(*   omega_tau / gamma = (1) / (5/4) = 4/5 = omega_c  ✅                   *)
