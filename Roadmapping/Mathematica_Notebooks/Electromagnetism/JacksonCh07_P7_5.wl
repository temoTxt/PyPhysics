(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh07_P7_5.wl  --  Polarisation states of plane wave.           *)
(*  Companion to Ch07_Plane_EM_Waves.md Problem J3e-P7.5.                  *)
(*                                                                         *)
(*  Stokes parameters (I, Q, U, V) for a general plane wave.              *)
(*  Polarisation states are field configurations, independent of source.  *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[capE0x, capE0y, delta, iStokes, qStokes, uStokes, vStokes];

(* Stokes parameters for plane wave with E_x = E_0x, E_y = E_0y * e^(i delta) *)
iStokes = capE0x^2 + capE0y^2;
qStokes = capE0x^2 - capE0y^2;
uStokes = 2 capE0x capE0y Cos[delta];
vStokes = 2 capE0x capE0y Sin[delta];

Print["Stokes parameters:"];
Print["  I = ", iStokes];
Print["  Q = ", qStokes];
Print["  U = ", uStokes];
Print["  V = ", vStokes];

(* Check fully-polarised constraint: I^2 = Q^2 + U^2 + V^2 *)
constraint = FullSimplify[iStokes^2 - (qStokes^2 + uStokes^2 + vStokes^2)];
Print[""];
Print["Constraint I^2 - (Q^2 + U^2 + V^2) = ", constraint];
Print["  Equals 0 for fully-polarised light?  ", constraint === 0];

(* Special cases: linear polarisation (delta = 0 or pi), circular (delta = pi/2, E_0x = E_0y) *)
Print[""];
Print["Linear polarisation (delta = 0):"];
Print["  V = ", vStokes /. delta -> 0];
Print["  (Zero, as expected for linear polarisation)"];
Print[""];
Print["Right circular polarisation (delta = pi/2, E_0x = E_0y = E_0):"];
Print["  Q = 0, U = 0, V = 2 E_0^2"];
