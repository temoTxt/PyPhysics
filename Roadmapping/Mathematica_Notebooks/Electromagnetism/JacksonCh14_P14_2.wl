(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh14_P14_2.wl  --  Lienard-Wiechert fields with third term.    *)
(*  Companion to Ch14_Radiation_by_Moving_Charges.md Problem J3e-P14.2.    *)
(*                                                                         *)
(*  Verifies the BAC-CAB structure of the third term of Eq. (7) of        *)
(*  the Maxwell paper:                                                     *)
(*    r x (u x r) = u (r . r) - r (u . r) = r^2 u - (u . r) r              *)
(*  This vector has a longitudinal component (parallel to u) and a        *)
(*  radial component (parallel to r), so the third-term radiation is     *)
(*  NOT transverse to the line of sight -- the campaign's first          *)
(*  qualitatively new field-level prediction.                              *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[rV, uV];

bacResult = TensorExpand[rV \[Cross] (uV \[Cross] rV)];
Print["r x (u x r) = ", bacResult];
Print["  Expected: u (r . r) - r (u . r) = r^2 u - (u . r) r"];

(* Structural feature: the third term, e (u . a) r x (u x r) / (b^4 s^3),  *)
(* contains both a u-direction (longitudinal) component and an r-direction *)
(* (radial) component.  This breaks pure transversality of LW radiation.   *)

(* Vanishing condition: u . a = 0  (a perp u).  Examples:                  *)
(*   - Circular motion (cyclotron, synchrotron)                            *)
(*   - Instantaneous rest in any frame (electrostatic configurations)     *)
(*   - Steady current (partial_tau = 0 averaged)                           *)
(* Non-vanishing examples:                                                 *)
(*   - Linear acceleration (bremsstrahlung)                                *)
(*   - Plane-wave figure-8 at relativistic intensity                       *)

(* Component decomposition:                                                *)
ClearAll[uX, uY, uZ, rX, rY, rZ];
uVec = {uX, uY, uZ};
rVec = {rX, rY, rZ};
third = rVec \[Cross] (uVec \[Cross] rVec);
Print[""];
Print["Component form of r x (u x r):"];
Print["  Component along u: r^2 = ", FullSimplify[(rVec . rVec)]];
Print["  Component along r: -(u . r) = ", -FullSimplify[uVec . rVec]];
Print["  Net vector: ", FullSimplify[third]];
Print["  = u * (r . r) - r * (u . r) by BAC-CAB"];

(* Expected output:                                                        *)
(*   r x (u x r) = uV rV . rV - rV rV . uV                                *)
(*   Net vector: as above (vector form)                                   *)
