(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh16_P16_1.wl  --  Abraham-Lorentz runaway + dissolution claim. *)
(*  Companion to Ch16_Radiation_Damping.md Problem J3e-P16.1.              *)
(*                                                                         *)
(*  Confirms the classical runaway solution a(t) = a_0 exp(t/tau_0) of    *)
(*  the source-free Abraham-Lorentz equation, with                        *)
(*    tau_0 = 2 e^2 / (3 m c^3)  ~ 6.3 x 10^-24 s  for the electron.      *)
(*                                                                         *)
(*  The proper-time formulation's dissipative term is first-order in     *)
(*  partial_tau (not third-order like LAD), so the equation of motion    *)
(*  has no exp(t/tau_0) runaway branch and no pre-acceleration.  This is *)
(*  the structural content of the proper-time "dissolution claim".       *)
(*                                                                         *)
(*  This notebook verifies the classical pathology only.  The proper-    *)
(*  time dissolution is a STRUCTURAL claim, not a symbolic computation;  *)
(*  it follows from the first-order vs third-order distinction in the    *)
(*  equation of motion.                                                    *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[ee, mm, cc, tt, tau0, alpha, runawaySol, abrLorentz];

tau0 = 2 ee^2/(3 mm cc^3);
Print["Abraham-Lorentz radiation-reaction time:  tau_0 = 2 e^2/(3 m c^3) = ", tau0];

(* Source-free Abraham-Lorentz: m a = m tau_0 da/dt, with F_ext = 0       *)
(* equivalent to:               da/dt = a/tau_0                            *)
abrLorentz = D[a[tt], tt] == a[tt]/tau0;
Print["Source-free AL equation: ", abrLorentz];

runawaySol = DSolve[{abrLorentz, a[0] == alpha}, a[tt], tt];
Print["Solution with a(0) = alpha: ", runawaySol];

(* The exponential timescale is 1/tau_0 = 3 m c^3 / (2 e^2) ~ 1.6 x 10^23 s^-1 *)
(* Numerical estimate for electron in Gaussian units:                      *)
(*   e^2 ~ 2.3e-19 esu^2                                                   *)
(*   m c^3 ~ 9.1e-28 g * (3e10 cm/s)^3 ~ 2.5e+4 g cm^3 s^-3                *)
(*   tau_0 ~ (2/3) * 2.3e-19 / 2.5e+4 ~ 6.1e-24 s                          *)
(* For ANY initial acceleration alpha, a(t) grows exponentially with     *)
(* timescale 6e-24 s.  Unphysical.                                         *)

Print[""];
Print["Numerical estimate: tau_0 ~ 6.1 x 10^-24 s for electron"];
Print[""];

(* The structural dissolution claim of the proper-time formulation:        *)
Print["Proper-time formulation: dissipative term in Eq. (4) is first-order"];
Print["in partial_tau (NOT third-order in d/dt like LAD)."];
Print["  => Equation of motion has no exp(t/tau_0) runaway branch"];
Print["  => No pre-acceleration needed to enforce causality"];
Print["  => Classical Abraham-Lorentz pathologies are STRUCTURALLY ABSENT"];
Print[""];
Print["Whether the proper-time RR force is numerically correct at all"];
Print["energies is the open question of issue #43 (Cole/Poder/Wistisen 2018)."];

(* Expected output:                                                        *)
(*   Solution: a(t) = alpha Exp[t/tau_0]                                  *)
(*     = alpha Exp[3 c^3 m t / (2 e^2)]                                   *)
