(* ::Package:: *)

(* JacksonCh04_P4_5.wl  --  Energy in linear dielectric.                  *)
(*                                                                        *)
(* In Gaussian units: u = (1/(8 pi)) E . D = (eps/(8 pi)) E^2             *)
(* For free charge q at fixed configuration, presence of dielectric       *)
(* reduces stored energy by factor 1/eps (relative to vacuum).            *)
(* Author: Trey Morris with Claude Opus 4.7.  Date: 2026-05-24.           *)

ClearAll[capE, eps, uDens, uDensVac, ratio];

uDens = eps capE^2 /(8 Pi);
uDensVac = capE^2 / (8 Pi);
ratio = uDens / uDensVac;

Print["Energy density in dielectric: u = ", uDens];
Print["Energy density in vacuum:    u_vac = ", uDensVac];
Print["Ratio u / u_vac = ", FullSimplify[ratio]];
Print["  = eps (the relative permittivity)"];
Print["  For free-charge fixed configuration, total energy stored REDUCES by factor 1/eps"];
Print["  (the dielectric partially screens the field)"];
