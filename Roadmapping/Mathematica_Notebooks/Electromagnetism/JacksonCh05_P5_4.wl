(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh05_P5_4.wl                                                    *)
(*                                                                         *)
(*  Companion notebook to Roadmapping/Electromagnetism/Jackson/            *)
(*  Ch05_Magnetostatics_Faraday_Quasi_Static.md  ->  Problem J3e-P5.4.     *)
(*                                                                         *)
(*  Circular current loop of radius capR carrying steady current capI in   *)
(*  the xy-plane.  Magnetic field on the symmetry axis at distance capZ.   *)
(*                                                                         *)
(*  Expected classical result:                                             *)
(*    B_z(capZ) = 2 Pi capI capR^2 / (cc (capR^2 + capZ^2)^(3/2))  (CGS)    *)
(*                                                                         *)
(*  Proper-time formulation: steady current => partial_tau = 0;            *)
(*    Ampere source (4 Pi / b) rho u = (4 Pi / b) (b/c) capJ = (4 Pi/c) J  *)
(*  Exact cancellation, identical field equation, identical B field.       *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(*  Runnable independent of the Wolfram MCP.                               *)
(* ----------------------------------------------------------------------- *)

ClearAll[phi, capR, capZ, capI, cc, bb, capJ, dlCrossR, loopIntegral, bzField];

(* Check 1: Biot-Savart geometric integral over the loop.                  *)
(* dl x r = (capR capZ Cos[phi], capR capZ Sin[phi], capR^2) dphi.        *)
dlCrossR = {capR capZ Cos[phi], capR capZ Sin[phi], capR^2};
loopIntegral = Integrate[dlCrossR, {phi, 0, 2 Pi}];
Print["dl x r integrated over loop = ", loopIntegral];
Print["Transverse components zero?  ",
   loopIntegral[[1]] === 0 && loopIntegral[[2]] === 0];
Print["Axial component = 2 Pi capR^2?  ",
   FullSimplify[loopIntegral[[3]] - 2 Pi capR^2] === 0];

(* Check 2: Full Biot-Savart field on the axis in Gaussian units.          *)
(* B_z = (capI/cc) * loopIntegral[[3]] / (capR^2 + capZ^2)^(3/2)          *)
bzField = (capI/cc) loopIntegral[[3]]/(capR^2 + capZ^2)^(3/2);
Print["B_z(capZ) = ", FullSimplify[bzField]];
Print["Matches textbook 2 Pi capI capR^2 / (cc (capR^2 + capZ^2)^(3/2))?  ",
   FullSimplify[bzField - 2 Pi capI capR^2/(cc (capR^2 + capZ^2)^(3/2))] === 0];

(* Check 3: Proper-time Ampere-Maxwell source vs classical source.         *)
(* (4 Pi / bb) rho u = (4 Pi / bb) (bb/cc) capJ = (4 Pi / cc) capJ.       *)
ampMaxProperTime = (4 Pi/bb) (bb/cc) capJ;
ampMaxClassical = (4 Pi/cc) capJ;
Print["Proper-time source = ", FullSimplify[ampMaxProperTime]];
Print["Classical source = ", ampMaxClassical];
Print["Difference: ", FullSimplify[ampMaxProperTime - ampMaxClassical]];

(* Expected output:                                                        *)
(*   dl x r integrated over loop = {0, 0, 2 capR^2 Pi}                     *)
(*   Transverse components zero?  True                                     *)
(*   Axial component = 2 Pi capR^2?  True                                  *)
(*   B_z(capZ) = (2 Pi capI capR^2) / (cc (capR^2 + capZ^2)^(3/2))         *)
(*   Matches textbook ...?  True                                           *)
(*   Proper-time source = (4 Pi capJ) / cc                                 *)
(*   Classical source = (4 Pi capJ) / cc                                   *)
(*   Difference: 0                                                         *)
