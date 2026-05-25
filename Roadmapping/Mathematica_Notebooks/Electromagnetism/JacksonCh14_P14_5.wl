(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh14_P14_5.wl  --  Synchrotron radiation from circular motion. *)
(*  Companion to Ch14_Radiation_by_Moving_Charges.md Problem J3e-P14.5.    *)
(*                                                                         *)
(*  Verifies the classical Lienard formula for circular motion reduces to *)
(*    P = (2/3)(e^2 a^2/c^3) gamma^4                                       *)
(*  via the identity w perp a => (w x a)^2 = w^2 a^2.                     *)
(*                                                                         *)
(*  Proper-time formulation: u . a = 0 for circular motion (a centripetal,*)
(*  u tangential), so the third term of Eq. (7) of the Maxwell paper      *)
(*  VANISHES identically.  The radiated power is identical to classical. *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[ee, aa, vv, cc, gammaW, lienardCircular, gamma4Form];

gammaW = 1/Sqrt[1 - vv^2/cc^2];

(* Lienard formula (Jackson 3e Eq. 14.26):                                *)
(*   P = (2/3)(e^2/c^3) gamma^6 [a^2 - (v x a)^2/c^2]                     *)
(* For circular motion w perp a: |v x a| = v a                            *)
lienardCircular = (2 ee^2/(3 cc^3)) gammaW^6 (aa^2 - vv^2 aa^2/cc^2);

gamma4Form = (2 ee^2 aa^2/(3 cc^3)) gammaW^4;

Print["Lienard for circular motion: ", FullSimplify[lienardCircular,
   Assumptions -> 0 < vv < cc]];
Print["Expected (2/3)(e^2 a^2/c^3) gamma^4: ", FullSimplify[gamma4Form,
   Assumptions -> 0 < vv < cc]];
Print["Match?  ", FullSimplify[lienardCircular - gamma4Form,
   Assumptions -> 0 < vv < cc] === 0];

(* Third-term vanishing: u . a = 0 by orthogonality.  Documented in prose. *)
Print[""];
Print["Proper-time third term of Eq. (7):"];
Print["  Coefficient: e (u . a) [r x (u x r)] / (b^4 s^3)"];
Print["  For circular motion: u . a = 0 identically (a centripetal, u tangential)"];
Print["  => Third term VANISHES."];
Print["  => Total radiated power equals classical Lienard formula."];

(* Numerical sample at v = 0.9 c (gamma ~ 2.29) *)
sample = {ee -> 1, aa -> 1, cc -> 1, vv -> 9/10};
Print[""];
Print["Numerical check at v = 0.9 c:"];
Print["  Lienard P = ", N[lienardCircular /. sample]];
Print["  (2/3) gamma^4 = ", N[gamma4Form /. sample]];
Print["  Expected gamma^4 = ", N[gammaW^4 /. sample]];
