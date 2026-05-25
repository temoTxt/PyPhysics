(* ::Package:: *)

(* Companion Wolfram Language notebook for
   Roadmapping/Electromagnetism/Jackson/Experimental_Comparisons/
       bremsstrahlung_MeV_spectra.md

   Issue #48 — proper-time third-term predictions for MeV bremsstrahlung.
   Specialises Ch14 P14.6's general derivation to 1, 6, 18 MeV clinical-linac energies.

   Status: SCAFFOLD ONLY. The quantitative-prediction substantive blocks of §3 of the
   prediction document remain TODO-gated for Zorns-Lemmon's review; this notebook's
   computational cells correspondingly remain stubs until that gate clears.

   Cross-references:
     - Canonical derivation: Mathematica_Notebooks/Electromagnetism/JacksonCh14_P14_6.wl
     - Companion problem in chapter doc: Ch14_Radiation_by_Moving_Charges.md §J3e-P14.6
     - Predecessor (radiation-reaction comparison at GeV scale, issue #43):
       Roadmapping/Electromagnetism/Jackson/Experimental_Comparisons/radiation_reaction_2018.md
*)

(* ============================================================
   Section 0 — Universal constants and kinematic helpers
   ============================================================
   Pragmatic content. Reused unchanged from Ch14 P14.6 notebook conventions.
*)

ClearAll[gammaFromT, betaFromGamma, uOverC, properTimeB];

(* T in MeV, electron rest energy m_e c^2 = 0.5109989 MeV *)
gammaFromT[T_] := 1 + T / 0.5109989;
betaFromGamma[g_] := Sqrt[1 - 1/g^2];
uOverC[g_] := g * betaFromGamma[g];        (* u = gamma * v, so u/c = gamma * beta *)
properTimeB[g_] := g;                       (* b = gamma * c, so b/c = gamma *)

(* Verification table reproduced from Ch14 P14.6:
   echoes Section 2's table of (u/c) -> third-term-ratio for the three target energies. *)
verificationTable = Table[
  {T, gammaFromT[T], uOverC[gammaFromT[T]], (uOverC[gammaFromT[T]])^2},
  {T, {1, 6, 18}}
];
(* {{T, gamma, u/c, (u/c)^2 = approximate ratio E_3rd / E_acc,classical at non-rel.}, ...} *)

Print["=== Verification table — kinematics at 1, 6, 18 MeV ==="];
Print[Grid[
  Prepend[verificationTable, {"T (MeV)", "gamma", "u/c", "(u/c)^2"}],
  Frame -> All
]];

(* ============================================================
   Section 1 — TODO: longitudinal-polarisation prediction at 1 MeV
   ============================================================
   Substantive. To be filled in after Zorns-Lemmon engages on the
   prediction document's §3.1 quantitative-prediction block.

   Compute:
   - Angular distribution dW/dOmega for linear deceleration at T = 1 MeV
   - Decompose into transverse (classical) and longitudinal (third-term) components
   - Quantify the longitudinal-component magnitude as a function of emission angle theta
   - Compare against Seltzer-Berger 1986 / EGSnrc / PENELOPE Born-approximation
     at the same energy
*)

(* TODO: substantive computation — Zorns-Lemmon to author or to flag for routing *)

(* ============================================================
   Section 2 — TODO: longitudinal-polarisation prediction at 6 MeV
   ============================================================
   Same template as Section 1, evaluated at T = 6 MeV.
   Special attention: TG-51 dosimetry calibration operates at this energy with
   ~1% precision on integrated cross-section. The substantive question is whether
   the proper-time third term modifies the integrated cross-section at the percent
   level or only the angular-distribution / polarisation pattern.
*)

(* TODO: substantive computation *)

(* ============================================================
   Section 3 — TODO: longitudinal-polarisation prediction at 18 MeV
   ============================================================
   Same template, evaluated at T = 18 MeV.
   Note: 18 MeV is approaching the regime where the precision Born-approximation
   literature thins out. The relativistic beaming half-angle 1/gamma ~ 1.6° at this
   energy requires fine angular resolution for any polarisation-resolved test.
*)

(* TODO: substantive computation *)

(* ============================================================
   Section 4 — TODO: cross-energy summary
   ============================================================
   Compile the three per-energy results into a single comparison plot.
   Verdict per energy: ✅ / ⚠ / ❌ against Born-approximation comparators.
*)

(* TODO: substantive synthesis *)
