(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh01_P1_5.wl                                                    *)
(*                                                                         *)
(*  Companion notebook to Roadmapping/Electromagnetism/Jackson/            *)
(*  Ch01_Introduction_Electrostatics.md  ->  Problem J3e-P1.5.             *)
(*                                                                         *)
(*  Electrostatic self-energy of a uniformly charged solid sphere of       *)
(*  radius R and total charge Q, in Gaussian (CGS) units.                  *)
(*                                                                         *)
(*  Expected classical result:  W = (3/5) Q^2 / R.                         *)
(*                                                                         *)
(*  Per [§5 of the verified Maxwell paper]                                 *)
(*  (Roadmapping/Equation_Verification/...), the proper-time formulation   *)
(*  with u = 0 reduces b -> c identically, so the proper-time prediction   *)
(*  is identical to the classical result and no separate Mathematica       *)
(*  check is performed for branch (c).                                     *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(*  Runnable independent of the Wolfram MCP.                               *)
(* ----------------------------------------------------------------------- *)

(* Use capR and capQ to avoid any ambiguity with single-letter Wolfram     *)
(* protected symbols, in keeping with the discipline recorded in          *)
(* CLAUDE.md ("How equation verification works", gotchas list).            *)
ClearAll[r, capR, capQ];

(* Inside the sphere, the field is E = (capQ r)/capR^3 (Gaussian units).  *)
(* Energy density is E^2/(8 Pi); volume element is 4 Pi r^2 dr.            *)
wInside =
 Integrate[
  (capQ^2 r^2 / capR^6) / (8 Pi) * 4 Pi r^2,
  {r, 0, capR}, Assumptions -> capR > 0
 ];

(* Outside the sphere, the field is E = capQ/r^2 (Gaussian units).        *)
wOutside =
 Integrate[
  (capQ^2 / r^4) / (8 Pi) * 4 Pi r^2,
  {r, capR, Infinity}, Assumptions -> capR > 0
 ];

(* Total energy.  Expected: (3 capQ^2) / (5 capR).                        *)
wTotal = FullSimplify[wInside + wOutside];

Print["w_inside  = ", wInside];
Print["w_outside = ", wOutside];
Print["w_total   = ", wTotal];

(* Verification, returns True if the classical result holds.              *)
Print[
 "Matches (3/5) capQ^2 / capR?  ",
 FullSimplify[wTotal - 3 capQ^2 / (5 capR)] === 0
];

(* Expected output:                                                        *)
(*   w_inside  = capQ^2 / (10 capR)                                        *)
(*   w_outside = capQ^2 / (2 capR)                                         *)
(*   w_total   = (3 capQ^2) / (5 capR)                                     *)
(*   Matches (3/5) capQ^2 / capR?  True                                    *)
