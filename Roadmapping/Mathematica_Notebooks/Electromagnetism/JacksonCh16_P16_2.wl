(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh16_P16_2.wl  --  Radiation-reaction damped oscillator.        *)
(*  Companion to Ch16_Radiation_Damping.md Problem J3e-P16.2.              *)
(*                                                                         *)
(*  Damping rate of a harmonic oscillator from radiation reaction:         *)
(*    Gamma = tau_0 omega_0^2 = (2 e^2 omega_0^2) / (3 m c^3)               *)
(*  For visible-light atomic transitions (omega_0 ~ 3e15 rad/s):           *)
(*    Gamma ~ 5e7 s^-1  (the natural linewidth)                            *)
(*    Q = omega_0 / Gamma ~ 6e7  (high-Q oscillator)                       *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[mm, ee, cc, omega0, tau0, gamma1, qFactor];

tau0 = 2 ee^2/(3 mm cc^3);
gamma1 = tau0 omega0^2;
qFactor = omega0/gamma1;

Print["Radiation-reaction time: tau_0 = ", tau0];
Print["Damping rate: Gamma = tau_0 omega_0^2 = ", gamma1];
Print["Quality factor: Q = omega_0 / Gamma = ", FullSimplify[qFactor]];
Print["  = 3 m c^3 / (2 e^2 omega_0)"];

(* Numerical estimate for visible-light atomic transition.                *)
(* Use Gaussian-units numerical values:                                    *)
(*   e^2 ~ 2.3e-19 esu^2                                                  *)
(*   m ~ 9.1e-28 g                                                        *)
(*   c ~ 3e10 cm/s                                                        *)
(*   omega_0 ~ 3e15 rad/s (visible light, ~ 500 THz)                     *)
Print[""];
Print["Numerical estimate for visible-light atomic transition:"];
Print["  Gamma = (2/3) * 2.3e-19 * (3e15)^2 / ((3e10)^3 * 9.1e-28)"];
Print["       ~ (2/3) * 2.3e-19 * 9e30 / (2.7e+31 * 9.1e-28)"];
Print["       ~ (2/3) * 2.07e+12 / 2.46e+4"];
Print["       ~ 5.6e+7 s^-1"];
Print["  Q ~ 3e15 / 5.6e7 ~ 5.4e+7  (very high-Q)"];
Print[""];
Print["This is the classical analogue of spontaneous emission rate."];
Print["Einstein A-coefficient for atomic visible-light transitions matches"];
Print["this classical estimate within factors of order unity."];
