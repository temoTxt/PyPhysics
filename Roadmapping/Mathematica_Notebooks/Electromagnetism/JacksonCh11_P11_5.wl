(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh11_P11_5.wl  --  Lorentz boost of position 4-vector.         *)
(*  Companion to Ch11_Special_Relativity.md Problem J3e-P11.5.             *)
(*                                                                         *)
(*  Verifies that the Minkowski interval c^2 t^2 - x^2 - y^2 - z^2 is     *)
(*  invariant under the standard Lorentz boost along x with velocity vv.   *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[tt, xx, yy, zz, vv, cc, gammaV, tPrime, xPrime, yPrime, zPrime,
   interval];

gammaV = 1/Sqrt[1 - vv^2/cc^2];

(* Lorentz boost along +x with velocity vv.                                *)
tPrime = gammaV (tt - vv xx/cc^2);
xPrime = gammaV (xx - vv tt);
yPrime = yy;
zPrime = zz;

(* Minkowski interval invariance check.                                    *)
interval = FullSimplify[
   cc^2 tPrime^2 - xPrime^2 - yPrime^2 - zPrime^2
      - (cc^2 tt^2 - xx^2 - yy^2 - zz^2),
   Assumptions -> 0 < vv < cc];

Print["Minkowski interval invariance: boost - original = ", interval];
Print["  Equals 0?  ", interval === 0];

(* Expected output:                                                        *)
(*   Minkowski interval invariance: boost - original = 0                  *)
(*     Equals 0?  True                                                    *)
