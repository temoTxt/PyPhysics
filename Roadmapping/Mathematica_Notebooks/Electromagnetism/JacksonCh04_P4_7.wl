(* ::Package:: *)

(* JacksonCh04_P4_7.wl  --  Boundary conditions at dielectric interface.  *)
(*                                                                        *)
(* From Maxwell's equations + divergence/Stokes theorems:                *)
(*   D_1n = D_2n (continuity of normal D, no free surface charge)         *)
(*   E_1t = E_2t (continuity of tangential E)                              *)
(* With D_i = eps_i E_i:                                                  *)
(*   eps_1 E_1n = eps_2 E_2n  (refraction of field lines at interface)   *)
(* Author: Trey Morris with Claude Opus 4.7.  Date: 2026-05-24.           *)

ClearAll[capE1n, capE2n, capE1t, capE2t, eps1, eps2, capD1n, capD2n];

capD1n = eps1 capE1n;
capD2n = eps2 capE2n;

Print["Boundary conditions at dielectric interface (no free surface charge):"];
Print["  E_1t = E_2t  (tangential E continuous)"];
Print["  D_1n = D_2n  (normal D continuous)"];
Print["    Substituting D = eps E:"];
Print["    eps_1 E_1n = eps_2 E_2n"];
Print["    E-field bends across interface, analogous to Snell's law for refraction."];
