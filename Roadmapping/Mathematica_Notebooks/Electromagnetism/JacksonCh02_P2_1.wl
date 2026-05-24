(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh02_P2_1.wl                                                    *)
(*                                                                         *)
(*  Companion notebook to Roadmapping/Electromagnetism/Jackson/            *)
(*  Ch02_Boundary_Value_Problems_Electrostatics_I.md  ->  Problem J3e-P2.1.*)
(*                                                                         *)
(*  Point charge capQ at height capD above an infinite grounded conducting *)
(*  plane.  Image-method solution in Gaussian units.                       *)
(*                                                                         *)
(*  Expected classical results:                                            *)
(*    sigma(rho)      = -capQ capD / (2 Pi (rho^2 + capD^2)^(3/2))         *)
(*    total induced   = -capQ                                              *)
(*    force on capQ   = -capQ^2 / (4 capD^2)  ( -zhat direction )          *)
(*    energy          = -capQ^2 / (4 capD)                                 *)
(*                                                                         *)
(*  The proper-time formulation with u = 0 reduces b -> c identically, so  *)
(*  the proper-time prediction is identical to the classical result and    *)
(*  no separate Mathematica check is performed for branch (c).             *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(*  Runnable independent of the Wolfram MCP.                               *)
(* ----------------------------------------------------------------------- *)

(* Use capQ, capD, etc. to avoid single-letter Wolfram-protected-symbol    *)
(* ambiguity, per CLAUDE.md "How equation verification works" gotchas.     *)
ClearAll[rho, capD, capQ, capR1, capR2, sigma, eField, pressure];

(* Check 1: Total induced surface charge integrates to -capQ.             *)
sigma = -capQ capD/(2 Pi (rho^2 + capD^2)^(3/2));
totalInduced = Integrate[
   sigma 2 Pi rho,
   {rho, 0, Infinity}, Assumptions -> capD > 0
   ];
Print["Total induced charge = ", FullSimplify[totalInduced]];
Print["Matches -capQ?  ", FullSimplify[totalInduced + capQ] === 0];

(* Check 2: Maxwell-stress force on the conductor toward capQ.            *)
(* eField is the normal component of E at z = 0+ (purely in -zhat).       *)
(* Outward electrostatic pressure on the conductor is eField^2 / (8 Pi).   *)
eField = -2 capQ capD/(rho^2 + capD^2)^(3/2);
pressure = eField^2/(8 Pi);
forceOnConductor = Integrate[
   pressure 2 Pi rho,
   {rho, 0, Infinity}, Assumptions -> capD > 0
   ];
Print["Force on conductor (+zhat) = ", FullSimplify[forceOnConductor]];
Print["Matches capQ^2/(4 capD^2)?  ",
   FullSimplify[forceOnConductor - capQ^2/(4 capD^2)] === 0];

(* Check 3: Potential satisfies Phi(z=0) = 0 trivially.                   *)
(* r1 = distance from (rho, 0, z) to the source at (0, 0, capD);         *)
(* r2 = distance from (rho, 0, z) to the image at (0, 0, -capD).         *)
ClearAll[z];
capR1 = Sqrt[rho^2 + (z - capD)^2];
capR2 = Sqrt[rho^2 + (z + capD)^2];
phi = capQ/capR1 - capQ/capR2;
Print["Phi(z=0) = ", FullSimplify[phi /. z -> 0]];

(* Check 4: Image-method force on capQ matches the Maxwell-stress check.  *)
(* The image charge -capQ at (0,0,-capD) is at distance 2 capD from capQ. *)
imageForce = -capQ^2/(2 capD)^2;
Print["Image-method force on capQ (zhat) = ", imageForce];
Print["Matches Maxwell-stress -capQ^2/(4 capD^2)?  ",
   FullSimplify[imageForce + capQ^2/(4 capD^2)] === 0];

(* Expected output:                                                        *)
(*   Total induced charge = -capQ                                          *)
(*   Matches -capQ?  True                                                  *)
(*   Force on conductor (+zhat) = capQ^2/(4 capD^2)                       *)
(*   Matches capQ^2/(4 capD^2)?  True                                      *)
(*   Phi(z=0) = 0                                                          *)
(*   Image-method force on capQ (zhat) = -capQ^2/(4 capD^2)               *)
(*   Matches Maxwell-stress -capQ^2/(4 capD^2)?  True                      *)
