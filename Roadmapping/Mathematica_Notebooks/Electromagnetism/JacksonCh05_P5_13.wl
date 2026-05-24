(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh05_P5_13.wl                                                   *)
(*                                                                         *)
(*  Companion notebook to Roadmapping/Electromagnetism/Jackson/            *)
(*  Ch05_Magnetostatics_Faraday_Quasi_Static.md  ->  Problem J3e-P5.13.    *)
(*                                                                         *)
(*  Magnetic dipole moment of a solid uniformly charged sphere of radius   *)
(*  capR, total charge capQ, rotating rigidly about z-axis with angular    *)
(*  velocity omega.                                                        *)
(*                                                                         *)
(*  Expected classical result:                                             *)
(*    m_z = capQ omega capR^2 / 5  (SI)                                    *)
(*    m_z = capQ omega capR^2 / (5 cc)  (Gaussian)                         *)
(*                                                                         *)
(*  Proper-time formulation: m_pt = (b/c) m_SI in the non-relativistic     *)
(*  limit (omega R / c -> 0) reduces to m_SI; first per-problem document   *)
(*  in the campaign where the proper-time formulation produces a non-zero  *)
(*  O((u/c)^2) deviation from the classical answer.                        *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(*  Runnable independent of the Wolfram MCP.                               *)
(* ----------------------------------------------------------------------- *)

ClearAll[r, theta, phi, capR, capRho, omega, capQ, mZ, mZSimplified];

(* Magnetic dipole moment integrand in spherical coordinates.              *)
(* m_z = (rho omega / 2) Integrate[r^2 sin^2(theta), dV] in SI;           *)
(* in Gaussian, divide by an additional cc.                                *)
(* dV = r^2 sin(theta) dr dtheta dphi.                                     *)

mZ = (capRho omega/2) Integrate[
   r^2 Sin[theta]^2 r^2 Sin[theta],
   {r, 0, capR}, {theta, 0, Pi}, {phi, 0, 2 Pi},
   Assumptions -> capR > 0
   ];
mZSimplified = FullSimplify[mZ /. capRho -> 3 capQ/(4 Pi capR^3)];
Print["m_z (SI) = ", mZSimplified];
Print["Matches Q omega R^2 / 5?  ",
   FullSimplify[mZSimplified - capQ omega capR^2/5] === 0];

(* In Gaussian units, divide by cc.                                        *)
mZGaussian = mZSimplified/cc;
Print["m_z (Gaussian) = ", mZGaussian];

(* Proper-time correction at leading order:                                *)
(* m_pt = (b/c) m_SI with b = c sqrt(1 + (omega R/c)^2 ... avg)           *)
(* The exact correction depends on the radial+latitudinal profile of u.   *)
(* At leading O((omega R / c)^2), the correction is positive (b > c).     *)
(* This notebook records only the leading-order classical computation;    *)
(* the proper-time correction is discussed in the prose of Ch05_*.md.     *)

(* Expected output:                                                        *)
(*   m_z (SI) = (capQ capR^2 omega) / 5                                    *)
(*   Matches Q omega R^2 / 5?  True                                        *)
(*   m_z (Gaussian) = (capQ capR^2 omega) / (5 cc)                         *)
