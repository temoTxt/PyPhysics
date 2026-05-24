(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh06_P6_1.wl                                                    *)
(*                                                                         *)
(*  Companion notebook to Roadmapping/Electromagnetism/Jackson/            *)
(*  Ch06_Maxwell_Equations_Macroscopic_Media.md  ->  Problem J3e-P6.1.     *)
(*                                                                         *)
(*  Maxwell's equations with magnetic monopole sources, in Gaussian units, *)
(*  in SI, and in the proper-time reformulation.  Verification of the     *)
(*  electric-magnetic duality and of the source-term reductions to        *)
(*  classical for both electric and magnetic monopole species.            *)
(*                                                                         *)
(*  Structural observation surfaced in the per-problem document: the      *)
(*  Gill-Zachary framework as published specifies the proper time tau for *)
(*  a single source species; for problems with both electric and magnetic *)
(*  monopoles having distinct proper velocities u_e and u_m, the          *)
(*  displacement-current term in the Faraday equation acquires an         *)
(*  ambiguity (which species' proper-time governs partial_tau?).  This is *)
(*  not a flagged inconsistency, just a domain-of-validity observation.   *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(*  Runnable independent of the Wolfram MCP.                               *)
(* ----------------------------------------------------------------------- *)

ClearAll[bb, cc, capJe, capJm, eVec, bVec, rhoE, rhoM, jE, jM];

(* Check 1: Electric-source reduction in proper-time Ampere-Maxwell.       *)
electricSource = (4 Pi/bb) (bb/cc) capJe;
electricSourceClassical = (4 Pi/cc) capJe;
Print["Electric source = ", FullSimplify[electricSource]];
Print["  Reduces to classical?  ",
   FullSimplify[electricSource - electricSourceClassical] === 0];

(* Check 2: Magnetic-source reduction in proper-time Faraday.              *)
magneticSource = (4 Pi/bb) (bb/cc) capJm;
magneticSourceClassical = (4 Pi/cc) capJm;
Print["Magnetic source = ", FullSimplify[magneticSource]];
Print["  Reduces to classical?  ",
   FullSimplify[magneticSource - magneticSourceClassical] === 0];

(* Check 3: Duality transformation under {E -> B, B -> -E, rhoE -> rhoM,  *)
(*   rhoM -> -rhoE, jE -> jM, jM -> -jE} is order 4 (dual^4 = identity).  *)
dual = {eVec -> bVec, bVec -> -eVec,
   rhoE -> rhoM, rhoM -> -rhoE,
   jE -> jM, jM -> -jE};
dual2 = {eVec -> -eVec, bVec -> -bVec,
   rhoE -> -rhoE, rhoM -> -rhoM,
   jE -> -jE, jM -> -jM};

(* Original quantities *)
quantities = {eVec, bVec, rhoE, rhoM, jE, jM};

(* After applying duality twice: should equal -original (i.e., dual^2 = -1) *)
twiceTransformed = quantities /. dual /. dual;
Print["Dual^2 applied to {E, B, rhoE, rhoM, jE, jM} = ", twiceTransformed];
Print["  Equals -original (i.e., dual^2 = -1)?  ",
   twiceTransformed === (-quantities)];

(* Verify a Gauss-law transformation: 4 Pi rhoE -> 4 Pi rhoM under dual.   *)
gaussESource = 4 Pi rhoE;
gaussESourceDual = gaussESource /. dual;
Print["Gauss(E) source = ", gaussESource, "; under dual -> ", gaussESourceDual];

(* Expected output:                                                        *)
(*   Electric source = (4 Pi capJe) / cc                                   *)
(*     Reduces to classical?  True                                         *)
(*   Magnetic source = (4 Pi capJm) / cc                                   *)
(*     Reduces to classical?  True                                         *)
(*   Dual^2 applied to {E, B, rhoE, rhoM, jE, jM} = {-eVec, -bVec, -rhoE,  *)
(*     -rhoM, -jE, -jM}                                                    *)
(*     Equals -original (i.e., dual^2 = -1)?  True                         *)
(*   Gauss(E) source = 4 Pi rhoE; under dual -> 4 Pi rhoM                  *)
