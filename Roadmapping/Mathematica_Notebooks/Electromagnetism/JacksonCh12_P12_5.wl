(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh12_P12_5.wl  --  E x B drift in crossed uniform fields.      *)
(*  Companion to Ch12_Relativistic_Dynamics.md Problem J3e-P12.5.          *)
(*                                                                         *)
(*  Solves the drift-frame force-balance F = q [E + (1/c) w x B] = 0      *)
(*  for crossed fields E = E0 x_hat and B = B0 z_hat, recovering the      *)
(*  textbook drift w_d = -(c E0/B0) y_hat.                                 *)
(*                                                                         *)
(*  In the proper-time formulation, the same balance condition holds via *)
(*  the velocity-duality identity u/b = w/c; the dissipative term         *)
(*  (u.a)/b^4 vanishes throughout because the full motion (drift +        *)
(*  cyclotron) satisfies u . a = 0 at every instant.                       *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[capE0, capB0, cc, vY, gammaDrift, uDrift];

(* Classical force balance: E + (1/c) w_drift x B = 0                      *)
(* w_drift = vY y_hat; w_drift x B_0 z_hat = vY B_0 x_hat                  *)
balance = capE0 + (vY capB0)/cc;
Print["Force-balance condition: ", balance, " = 0"];

solution = Solve[balance == 0, vY];
Print["Drift velocity: ", solution];

vYDrift = vY /. solution[[1]];
Print["  w_drift = ", vYDrift, " y_hat  =  -(c E0 / B0) y_hat"];

(* Proper-velocity equivalent: u_drift = gamma_drift * w_drift             *)
gammaDrift = 1/Sqrt[1 - vYDrift^2/cc^2];
uDrift = gammaDrift vYDrift;
Print[""];
Print["Lorentz factor for drift: gamma_drift = ", FullSimplify[gammaDrift,
   Assumptions -> capE0 > 0 && capB0 > capE0 && cc > 0]];
Print["  (real and finite when |E0| < |B0|)"];
Print["Proper-velocity drift: u_drift = gamma_drift * w_drift = ",
   FullSimplify[uDrift,
      Assumptions -> capE0 > 0 && capB0 > capE0 && cc > 0]];

(* Expected output:                                                        *)
(*   Force-balance condition: capE0 + (capB0 vY)/cc = 0                   *)
(*   Drift velocity: {{vY -> -capE0 cc / capB0}}                          *)
(*     w_drift = -capE0 cc / capB0 y_hat                                  *)
(*   gamma_drift = 1/Sqrt[1 - capE0^2/capB0^2]                            *)
(*   u_drift = -capE0 cc / Sqrt[capB0^2 - capE0^2]                        *)
