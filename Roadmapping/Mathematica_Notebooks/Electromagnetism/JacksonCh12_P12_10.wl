(* ::Package:: *)

(* ----------------------------------------------------------------------- *)
(*  JacksonCh12_P12_10.wl  --  Hamiltonian for charged particle.          *)
(*  Companion to Ch12_Relativistic_Dynamics.md Problem J3e-P12.10.         *)
(*                                                                         *)
(*  Free-particle Hamiltonian in observer-time variables (H = gamma m c^2) *)
(*  vs proper-time variables (H = m c b).  Both are equivalent under       *)
(*  b = gamma c.                                                           *)
(*                                                                         *)
(*  Eq. 24 of the Maxwell paper (Pauli-form Hamiltonian with flagged       *)
(*  spin-magnetic-moment and V^2/(2mc^2) terms) does NOT engage in this   *)
(*  problem because Jackson Ch. 12 P10 is classical (non-spin).  The     *)
(*  branched-treatment workflow per §5.1 of the campaign plan is deferred *)
(*  to spin/hyperfine problems in PR F+.                                  *)
(*                                                                         *)
(*  Author: Trey Morris with Claude Opus 4.7 (1M context).                 *)
(*  Date: 2026-05-24.                                                      *)
(* ----------------------------------------------------------------------- *)

ClearAll[uu, cc, mm, bb, ww, gammaW, hFree, hClassicalFree, diff];

(* Free Hamiltonian in proper-time variables *)
bb = Sqrt[cc^2 + uu^2];
hFree = mm cc bb;
Print["H_free (proper-time form) = ", hFree];

(* Free Hamiltonian in observer-time variables *)
gammaW = 1/Sqrt[1 - ww^2/cc^2];
hClassicalFree = gammaW mm cc^2;
Print["H_free (classical form) = ", FullSimplify[hClassicalFree]];

(* Substitute uu = ww/sqrt(1 - ww^2/cc^2) to convert proper-time -> classical *)
diff = FullSimplify[
   hFree /. uu -> ww/Sqrt[1 - ww^2/cc^2],
   Assumptions -> 0 < ww < cc && cc > 0 && mm > 0];
Print[""];
Print["H_free (proper-time) under u = gamma w substitution: ", diff];
Print["H_free (classical): ", FullSimplify[hClassicalFree]];
Print["Match?  ", FullSimplify[diff - hClassicalFree,
   Assumptions -> 0 < ww < cc && cc > 0 && mm > 0] === 0];

(* Hamilton's equation: dx/dt = dH/dp *)
(* For free particle, p = mu (proper momentum), so u = p/m, b = sqrt(c^2 + p^2/m^2) *)
ClearAll[bb];
hAsP = Sqrt[cc^2 + uu^2/1 (* p/m placeholder*) ];
ClearAll[pp];
hAsPRedone = Sqrt[cc^2 + (pp/mm)^2] mm cc;
vFromHam = D[hAsPRedone, pp];
Print[""];
Print["Hamilton's equation dx/dt = dH/dp (free particle):"];
Print["  dH/dp = ", FullSimplify[vFromHam]];
Print["  = p c^2 / H = p c / (m b)  (with p = mu, b = sqrt(c^2+u^2))"];
Print["  = u c / b = w  (observable velocity)  ✓"];

(* Eq. 24 engagement check: spelled out in companion document.  No new symbolic check needed; *)
(* Jackson Ch. 12 P10 is non-spin, so Eq. 24's spin-magnetic-moment and V^2/(2mc^2) terms     *)
(* do not enter.  Branched-treatment per §5.1 of plan deferred to PR F+.                       *)
Print[""];
Print["Eq. 24 of Maxwell paper:  Pauli-form Hamiltonian with flagged spin and V^2 terms."];
Print["This problem is classical (non-spin) -- Eq. 24 does NOT engage."];
Print["Branched-treatment workflow deferred to PR F+ (Pauli reduction problems)."];

(* Expected output:                                                        *)
(*   H_free (proper-time form) = cc mm Sqrt[cc^2 + uu^2]                  *)
(*   H_free (classical form) = (cc^2 mm) / Sqrt[1 - ww^2/cc^2]            *)
(*   ... after substitution, they match.                                  *)
