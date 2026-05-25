(* ::Package:: *)

(* BetheSalpeter_S3.wl  --  Companion to Bethe_Salpeter/01_NonRelHydrogen.md *)
(*                                                                            *)
(* PR A scaffold result.  Verifies the load-bearing operator identity         *)
(*                                                                            *)
(*    K = H_0^2/(2 m c^2) + m c^2/2   with  H_0^2 = c^2 p^2 + m^2 c^4         *)
(*                                                                            *)
(* reduces *exactly* (not perturbatively) to                                  *)
(*                                                                            *)
(*    K = m c^2 + p^2/(2 m)                                                   *)
(*                                                                            *)
(* This is the campaign's load-bearing claim: the proper-time formulation     *)
(* recovers non-relativistic kinetic energy plus a constant rest-energy       *)
(* offset, with no (p/mc)^k correction at any order.  Coulomb potential V_0   *)
(* enters as a scalar addition unchanged from Bethe-Salpeter.                 *)
(*                                                                            *)
(* Author: Trey Morris with Claude Opus 4.7.  Date: 2026-05-25.               *)

ClearAll[m, c, p];
expr = ((c^2 p^2 + m^2 c^4)/(2 m c^2) + m c^2 / 2) - (m c^2 + p^2/(2 m));
Print["K - (m c^2 + p^2/(2m)) = ", FullSimplify[expr]];
Print["Expected: 0 (proper-time K reduces exactly to NR kinetic + rest energy)"];
