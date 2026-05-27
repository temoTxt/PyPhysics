(* ::Package:: *)

(* r_e_Li2plus_hyperfine.wl  --  Companion to Bethe_Salpeter/13_Li2plus_Hyperfine.md *)
(*                                                                            *)
(* Issue #78, observable #4: Li-7 1s hyperfine splitting in hydrogenic Li^2+. *)
(* Branch 78-li2plus-hyperfine.                                               *)
(*                                                                            *)
(* GOAL: extend the BS-S22 H 21-cm Fermi-contact apparatus (nuclear spin      *)
(* I=1/2) to Z=3 and the Li-7 nuclear spin I=3/2, then evaluate the           *)
(* dual-theory framework prediction = (g_s/-2)^1 * (textbook leading term)    *)
(* at the triangulated r_e/r_0 = 0.4994205099128317 (g_s = -2.00231930).      *)
(*                                                                            *)
(* SUBSTANTIVE-AI TAG (Crocco 2026 substantive use):                          *)
(* The I=1/2 -> I=3/2 extension of the Fermi-contact angular structure is a   *)
(* non-trivial substantive modeling choice, NOT a mechanical Z-rescaling.     *)
(* H (I=1/2) has a single two-level splitting F=0,1; Li-7 (I=3/2 + S=1/2)     *)
(* gives F=1,2, and the headline observable is the F=2<->1 interval.          *)
(* The angular factor is derived here from the Lande/interval rule and is     *)
(* flagged for human review in 13_Li2plus_Hyperfine.md.                       *)
(* <!-- TODO: human reviews and fills in the I=3/2 angular-structure choice -->*)
(*                                                                            *)
(* MEASUREMENT-PROVENANCE CAVEAT: no DIRECT experimental Li^2+ 1s HFS value   *)
(* exists. The authoritative comparator is the QED-theory value of Pachucki,  *)
(* Patkos, Yerokhin et al., "Hyperfine splitting in 6,7Li+", arXiv:2309.00436 *)
(* (2023), Table VII, which itself uses the experimental Li+ HFS as input.    *)
(*                                                                            *)
(* WOLFRAM-MCP GOTCHAS (per CLAUDE.md): each cell is a single line joined by  *)
(* ';'. No bare 'V' or 'e' symbols. All identifiers descriptive.              *)
(*                                                                            *)
(* Author: Trey Morris with Claude Opus 4.7.  Date: 2026-05-27.               *)


(* ---- Cell 1: H 1s sanity check -- reproduce the 21-cm 1420 MHz baseline ---- *)
(* Leading Fermi formula  dnu_H = (4/3) gP (me/Mp) alpha^4 (me c^2 / h)        *)
(* gP = mu_p/(I_p mu_N) = 2.7928473 / (1/2) = 5.5856946 (proton g-factor).     *)
ClearAll[alphaFS, meOverMp, gProton, mec2Hz, dnuHsanity]; alphaFS = 1/137.035999084; meOverMp = 1/1836.15267343; gProton = 5.5856946893; mec2Hz = 1.23558996028*10^20; dnuHsanity = (4/3)*gProton*meOverMp*alphaFS^4*mec2Hz/10^6; Print["Cell1: leading-Fermi H 1s HFS = ", dnuHsanity, " MHz  (expect ~1418-1421; measured 1420.4058 MHz)"];


(* ---- Cell 2: Li-7 / H scaling ratio of the textbook (g_s=-2) splitting ---- *)
(* Contact hyperfine constant a ~ (mu_I/I) |psi(0)|^2 ~ (mu_I/I) Z^3.          *)
(* Top-interval splitting dnu = a (I + 1/2)  =>  dnu ~ mu_I (2I+1)/(2I) Z^3.    *)
(* So  dnu_Li / dnu_H = [muLi (2 ILi+1)/(2 ILi)] / [muP (2 IP+1)/(2 IP)] * Z^3. *)
(* mu in nuclear magnetons; muP = 2.79284734, ILi = 3/2, IP = 1/2, Z = 3.      *)
ClearAll[muProton, muLi7, spinP, spinLi, zCharge, nucFactorP, nucFactorLi, scaleRatio, dnuHmeas, dnuLiTextbook]; muProton = 2.79284734; muLi7 = 3.256427; spinP = 1/2; spinLi = 3/2; zCharge = 3; nucFactorP = muProton*(2*spinP+1)/(2*spinP); nucFactorLi = muLi7*(2*spinLi+1)/(2*spinLi); scaleRatio = (nucFactorLi/nucFactorP)*zCharge^3; dnuHmeas = 1420.405751768; dnuLiTextbook = scaleRatio*dnuHmeas; Print["Cell2: scaleRatio (Li/H) = ", scaleRatio, " ;  textbook Li2+ 1s HFS (g_s=-2) = ", dnuLiTextbook, " MHz"];


(* ---- Cell 3: framework correction (g_s/-2)^1 at triangulated r_e ---------- *)
(* Hyperfine is LINEAR in g_s (n=1, per 10_CrossComparison.md S2).            *)
(* gsTri = -2.00231930436 at r_e/r_0 = 0.4994205099128317 (PR #62).           *)
ClearAll[gsTri, gsFactor, dnuLiFramework]; gsTri = -2.00231930436; gsFactor = gsTri/(-2); dnuLiFramework = dnuLiTextbook*gsFactor; Print["Cell3: (g_s/-2) = ", gsFactor, " ;  framework Li2+ 1s HFS prediction = ", dnuLiFramework, " MHz"];


(* ---- Cell 4: I=3/2 angular structure -- F levels and interval rule -------- *)
(* 1s_{1/2} (S=1/2) coupled to Li-7 I=3/2 gives F = 1, 2.                      *)
(* I.S eigenvalue = (1/2)[F(F+1) - I(I+1) - S(S+1)].                           *)
(* Splitting (F=2 <-> F=1) in units of the hyperfine constant a = a*(I+1/2).   *)
ClearAll[isEig, spinI, spinS, fUpper, fLower, splitInA]; spinI = 3/2; spinS = 1/2; isEig[ff_] := (1/2)*(ff*(ff+1) - spinI*(spinI+1) - spinS*(spinS+1)); fUpper = 2; fLower = 1; splitInA = isEig[fUpper] - isEig[fLower]; Print["Cell4: I.S(F=2)=", isEig[fUpper], " I.S(F=1)=", isEig[fLower], " ; F=2<->1 splitting = ", splitInA, " * a  (= (I+1/2) a = 2a, the headline interval)"];
