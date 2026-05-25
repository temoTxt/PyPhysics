(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh16_P16_5.wl  --  Proper-time RR for Cole/Poder geometry.     *)
(*  Companion to Ch16_Radiation_Damping.md Problem J3e-P16.5.              *)
(*                                                                         *)
(*  Sketch-level parameter estimates for the proper-time RR prediction    *)
(*  in the Cole 2018 / Poder 2018 experimental geometry.                  *)
(*                                                                         *)
(*  Full quantitative comparison is the deliverable of issue #43:         *)
(*    Electromagnetism/Jackson/Experimental_Comparisons/                   *)
(*       radiation_reaction_2018.md                                        *)
(*                                                                         *)
(*  This notebook records only the input parameters and order-of-          *)
(*  magnitude estimates.                                                   *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[ebeam, mec2, gammaE, capIL, a0, chiQED];

(* Cole/Poder beam parameters *)
ebeam = 1.7 10^9;        (* 1.7 GeV electron beam, in eV *)
mec2 = 0.511 10^6;       (* electron rest mass-energy, in eV *)
gammaE = ebeam / mec2;
Print["Electron beam Lorentz factor: gamma_e = E/(m c^2) = ", N[gammaE]];

(* Laser parameters *)
capIL = 10^21;            (* 10^21 W/cm^2 *)
a0 = 50;                  (* dimensionless intensity at I ~ 10^21 W/cm^2 for optical laser *)
Print["Laser normalized vector potential: a_0 ~ ", a0, " (highly relativistic)"];

(* QED-relevant invariant chi *)
chiQED = 0.1;             (* order-of-magnitude estimate for chi ~ gamma * a_0 * (hbar omega)/(m c^2) *)
Print["QED invariant chi ~ ", chiQED, " (regime where RR + QED both matter)"];

(* Predicted radiation-reaction signatures *)
Print[""];
Print["RR predictions for Cole/Poder geometry:"];
Print["  Classical Landau-Lifshitz: Delta E/E ~ 10-20%"];
Print["  Quantum-corrected LL (Cole/Poder favoured): ~ 10% smaller than classical LL"];
Print["  Proper-time third-term correction: ~ chi * (LL prediction) ~ 10% correction to LL"];
Print["  Within Cole/Poder precision floor (10-20%): at boundary of distinguishability"];

Print[""];
Print["Full quantitative comparison is the deliverable of issue #43."];
Print["  See: Electromagnetism/Jackson/Experimental_Comparisons/radiation_reaction_2018.md"];
