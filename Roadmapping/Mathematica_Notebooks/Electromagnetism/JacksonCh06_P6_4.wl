(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh06_P6_4.wl                                                    *)
(*                                                                         *)
(*  Companion notebook to Roadmapping/Electromagnetism/Jackson/            *)
(*  Ch06_Maxwell_Equations_Macroscopic_Media.md  ->  Problem J3e-P6.4.     *)
(*                                                                         *)
(*  Electromagnetic momentum of a uniformly-moving point charge,           *)
(*  regularised as a uniformly charged sphere of radius capR with total    *)
(*  charge qq, moving with velocity vv along z-axis.                       *)
(*                                                                         *)
(*  Expected classical result:                                             *)
(*    P_field = (4/3) U_field v/c^2 = 4 qq^2 vv/(5 capR cc^2)              *)
(*  with U_field = 3 qq^2/(5 capR) (Gaussian, from PR 0 P1.5).             *)
(*                                                                         *)
(*  The famous "4/3 problem" of classical electron theory; the proper-time *)
(*  reformulation reproduces this result exactly via the velocity-duality  *)
(*  identity 1 - u^2/b^2 = 1 - w^2/c^2, without offering an alternative    *)
(*  resolution to the Abraham-Lorentz-Poincare puzzle.                     *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(*  Runnable independent of the Wolfram MCP.                               *)
(* ----------------------------------------------------------------------- *)

ClearAll[r, theta, phi, qq, capR, vv, cc, ww, uu, bb, bSolved, uSolved];

(* Check 1: integral of E^2 over space for uniformly charged sphere.       *)
eSquaredIn = Integrate[
   (qq r/capR^3)^2 4 Pi r^2,
   {r, 0, capR}, Assumptions -> capR > 0
   ];
eSquaredOut = Integrate[
   (qq/r^2)^2 4 Pi r^2,
   {r, capR, Infinity}, Assumptions -> capR > 0
   ];
totalESquared = FullSimplify[eSquaredIn + eSquaredOut];
Print["Integrate[E^2 dV] = ", totalESquared];
Print["  Matches 24 Pi qq^2/(5 capR)?  ",
   FullSimplify[totalESquared - 24 Pi qq^2/(5 capR)] === 0];

(* Check 2: field momentum P_z = (v/(4 Pi c^2)) (2/3) Integrate[E^2 dV].   *)
(* The (2/3) factor arises from the spherical symmetry of E:               *)
(* Integrate[E_x^2 + E_y^2 dV] = (2/3) Integrate[E^2 dV].                   *)
pField = (vv/(4 Pi cc^2)) (2/3) totalESquared;
Print["P_field = ", FullSimplify[pField]];

(* The result expressed in terms of U_field = 3 qq^2/(5 capR):             *)
uFieldEnergy = 3 qq^2/(5 capR);
expected = (4/3) uFieldEnergy vv/cc^2;
Print["(4/3) U_field v/c^2 = ", FullSimplify[expected]];
Print["  Match?  ", FullSimplify[pField - expected] === 0];

(* Check 3: velocity-duality identity 1 - u^2/b^2 = 1 - w^2/c^2.           *)
(* From w/c = u/b and b^2 = c^2 + u^2:                                      *)
(*   b = c/Sqrt[1 - w^2/c^2],  u = w/Sqrt[1 - w^2/c^2]                     *)
bSolved = cc/Sqrt[1 - ww^2/cc^2];
uSolved = ww/Sqrt[1 - ww^2/cc^2];
identityCheck = FullSimplify[
   1 - uSolved^2/bSolved^2 - (1 - ww^2/cc^2),
   Assumptions -> 0 < ww < cc
   ];
Print["(1 - u^2/b^2) - (1 - w^2/c^2) = ", identityCheck];
Print["  Identity holds?  ", identityCheck === 0];

(* This identity is the key insight: the proper-time Lienard-Wiechert      *)
(* formula uses (1 - u^2/b^2) where the classical formula uses             *)
(* (1 - w^2/c^2).  They are numerically identical for uniform motion;     *)
(* the proper-time framework does not introduce any new factor that       *)
(* alters the EM-momentum integral or the 4/3 puzzle.                      *)

(* Expected output:                                                        *)
(*   Integrate[E^2 dV] = 24 Pi qq^2/(5 capR)                              *)
(*     Matches 24 Pi qq^2/(5 capR)?  True                                 *)
(*   P_field = 4 qq^2 vv/(5 capR cc^2)                                    *)
(*   (4/3) U_field v/c^2 = 4 qq^2 vv/(5 capR cc^2)                        *)
(*     Match?  True                                                       *)
(*   (1 - u^2/b^2) - (1 - w^2/c^2) = 0                                    *)
(*     Identity holds?  True                                              *)
