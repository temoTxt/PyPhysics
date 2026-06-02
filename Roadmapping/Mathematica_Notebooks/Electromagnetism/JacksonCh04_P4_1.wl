(* ::Package:: *)

(* JacksonCh04_P4_1.wl  --  Multipole expansion of potential.            *)
(*                                                                        *)
(* phi(r) = q/r + (p . r_hat)/r^2 + Q_ij r_hat_i r_hat_j / (2 r^3) + ... *)
(* Static configuration => u = 0, b = c; proper-time reduces identically. *)
(* Author: Trey Morris with Claude Opus 4.7.  Date: 2026-05-24.           *)

ClearAll[rr, theta, phi, qq, pVec, capQTensor, rhatVec];

(* Field-point unit vector in spherical coords *)
rhatVec = {Sin[theta] Cos[phi], Sin[theta] Sin[phi], Cos[theta]};

(* Multipole terms *)
monopoleTerm = qq/rr;
dipoleTerm = (pVec . rhatVec)/rr^2;

Print["Multipole expansion:"];
Print["  Monopole:  ", monopoleTerm];
Print["  Dipole:    ", dipoleTerm];
Print["  Quadrupole: Q_ij rhat_i rhat_j / (2 r^3)  (tensor contraction)"];
