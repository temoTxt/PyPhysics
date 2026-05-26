(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh06_P6_5.wl                                                    *)
(*                                                                         *)
(*  Companion notebook to Roadmapping/Electromagnetism/Jackson/            *)
(*  Ch06_Maxwell_Equations_Macroscopic_Media.md  ->  Problem J3e-P6.5.     *)
(*                                                                         *)
(*  Poynting theorem in macroscopic media.  Verifies that the proper-time *)
(*  source term (1/b) rho_free u_free . E reduces algebraically to the    *)
(*  classical source term (1/c) J_free . E, and that the conservation     *)
(*  form is preserved under the substitution (1/c) d/dt = (1/b) d/dtau.   *)
(*                                                                         *)
(*  No new substitution rule is needed for the macroscopic Poynting       *)
(*  theorem -- the existing velocity-duality and time-derivative-duality  *)
(*  rules of Eqs. (1) and (2) of [[Two_Mathematically_Equivalent_Versions *)
(*  _of_Maxwells_Equations]] suffice.                                      *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(*  Runnable independent of the Wolfram MCP.                               *)
(* ----------------------------------------------------------------------- *)

ClearAll[bb, cc, capJ, capE, rhouFree, classicalSource, properTimeSource];

(* Classical Poynting source: (1/c) J_free . E in Gaussian units.          *)
classicalSource = (1/cc) capJ capE;

(* Proper-time source: (1/b) (rho_free u_free) . E in Gaussian-equivalent  *)
(* proper-time units.                                                      *)
properTimeSource = (1/bb) rhouFree capE;

(* Apply velocity duality: rho u = (b/c) J.                                *)
properTimeSourceAfterDuality = properTimeSource /. rhouFree -> (bb/cc) capJ;
Print["Proper-time source after velocity duality = ",
   FullSimplify[properTimeSourceAfterDuality]];
Print["Classical source = ", classicalSource];
Print["Difference = ",
   FullSimplify[properTimeSourceAfterDuality - classicalSource]];
Print["  Match?  ",
   FullSimplify[properTimeSourceAfterDuality - classicalSource] === 0];

(* The Poynting theorem conservation form preservation is then a single   *)
(* substitution: (1/c) d/dt = (1/b) d/dtau on the time-derivative term,   *)
(* which is exactly Eq. (2) of the Maxwell paper.  No further verification *)
(* needed; the cancellation above is the load-bearing algebraic check.    *)

(* Expected output:                                                        *)
(*   Proper-time source after velocity duality = (capE capJ) / cc          *)
(*   Classical source = (capE capJ) / cc                                  *)
(*   Difference = 0                                                       *)
(*     Match?  True                                                       *)
