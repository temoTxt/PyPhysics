(* ::Package:: *)

(* JacksonCh04_P4_3.wl  --  Polarization in linear dielectric.            *)
(*                                                                        *)
(* P = chi_e E (Gaussian);  D = E + 4 pi P = eps E with eps = 1 + 4 pi chi_e *)
(* Bound charges: rho_b = -div P; sigma_b = P . nhat.                     *)
(* Stationary dielectric => u = 0, b = c; proper-time reduces identically.*)
(* Author: Trey Morris with Claude Opus 4.7.  Date: 2026-05-24.           *)

ClearAll[chi, capE, capP, capD, eps];

capP = chi capE;
capD = capE + 4 Pi capP;
eps = 1 + 4 Pi chi;

Print["Polarization: P = chi_e E = ", capP];
Print["Displacement: D = E + 4 pi P = ", capD];
Print["Permittivity: eps = 1 + 4 pi chi_e = ", eps];
Print["  D = eps E in Gaussian units"];
