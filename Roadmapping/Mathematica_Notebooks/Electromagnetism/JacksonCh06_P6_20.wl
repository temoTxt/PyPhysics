(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh06_P6_20.wl                                                   *)
(*                                                                         *)
(*  Companion notebook to Roadmapping/Electromagnetism/Jackson/            *)
(*  Ch06_Maxwell_Equations_Macroscopic_Media.md  ->  Problem J3e-P6.20.    *)
(*                                                                         *)
(*  Radiation pressure on a perfect conductor from a normally incident    *)
(*  plane wave.  Classical Gaussian computation.                           *)
(*                                                                         *)
(*  Expected result: P_rad = E_0^2 / (4 Pi) = 2 I / c, where               *)
(*    I = c E_0^2 / (8 Pi)  (incident intensity).                          *)
(*                                                                         *)
(*  Proper-time formulation: under the perfect-conductor idealisation     *)
(*  (infinite conductivity, instantaneous surface-charge response), the   *)
(*  dissipative (u . a) / b^4 term of Eq. (4) does NOT contribute.  The   *)
(*  proper-time prediction is identical to the classical result.  The     *)
(*  first non-null engagement of the dissipative term is deferred to     *)
(*  PR D (Ch. 14 Lienard-Wiechert).                                        *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(*  Runnable independent of the Wolfram MCP.                               *)
(* ----------------------------------------------------------------------- *)

ClearAll[t, omega, capE0, cc, kk, z];

(* Incident wave at z = 0:                                                 *)
(*   E_inc(0, t) = E_0 x cos(- omega t) = E_0 x cos(omega t)               *)
(*   B_inc(0, t) = E_0 y cos(- omega t) = E_0 y cos(omega t)               *)
(* Reflected wave at z = 0:                                                *)
(*   E_ref(0, t) = -E_0 x cos(omega t)                                     *)
(*   B_ref(0, t) =  E_0 y cos(omega t)                                     *)
(* Total at z = 0:                                                         *)
(*   E_total = 0  (boundary condition)                                     *)
(*   B_total = 2 E_0 y cos(omega t)                                        *)

eTotalAtZero = capE0 Cos[-omega t] + (-capE0) Cos[omega t];
bTotalAtZero = capE0 Cos[-omega t] + capE0 Cos[omega t];

Print["E_total at z = 0: ", FullSimplify[eTotalAtZero]];
Print["B_total at z = 0: ", FullSimplify[bTotalAtZero]];

(* Maxwell stress T_zz at z = 0:                                          *)
(*   T_zz = (1/(4 Pi))[E_z^2 + B_z^2 - (1/2)(E^2 + B^2)]                  *)
(*   Only E_x and B_y are nonzero, so T_zz = -B_total^2 / (8 Pi)          *)

tZZ = -bTotalAtZero^2/(8 Pi);
Print["T_zz at z = 0: ", FullSimplify[tZZ]];

(* Time-average over one period 2 Pi / omega.                              *)
timeAvgTZZ = Integrate[tZZ, {t, 0, 2 Pi/omega}] omega/(2 Pi);
timeAvgTZZ = FullSimplify[timeAvgTZZ];
Print["Time-averaged T_zz: ", timeAvgTZZ];

(* Radiation pressure = magnitude of the time-averaged stress.            *)
radiationPressure = -timeAvgTZZ;
Print["Radiation pressure: ", radiationPressure];

(* Incident intensity in Gaussian units: I = c <E^2>/(4 Pi).               *)
(* For monochromatic plane wave: <E^2> = E_0^2 / 2, so I = c E_0^2 / (8 Pi). *)
incidentI = cc capE0^2/(8 Pi);
Print["Incident intensity I: ", incidentI];
Print["2 I / c: ", FullSimplify[2 incidentI/cc]];

(* Verify P_rad = 2 I / c.                                                *)
Print["Match P_rad = 2 I / c?  ",
   FullSimplify[radiationPressure - 2 incidentI/cc] === 0];

(* Expected output:                                                        *)
(*   E_total at z = 0: 0                                                   *)
(*   B_total at z = 0: 2 capE0 Cos[omega t]                               *)
(*   T_zz at z = 0: -capE0^2 Cos[omega t]^2/(2 Pi)                         *)
(*   Time-averaged T_zz: -capE0^2/(4 Pi)                                   *)
(*   Radiation pressure: capE0^2/(4 Pi)                                    *)
(*   Incident intensity I: capE0^2 cc/(8 Pi)                              *)
(*   2 I / c: capE0^2/(4 Pi)                                              *)
(*   Match P_rad = 2 I / c?  True                                          *)
