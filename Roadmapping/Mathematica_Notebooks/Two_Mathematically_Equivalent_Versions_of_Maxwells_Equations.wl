(* ::Package:: *)

(* ================================================================
   Equation verification for:
     Gill & Zachary (2011), "Two Mathematically Equivalent Versions
     of Maxwell's Equations", Foundations of Physics 41, 99-128.

   Companion to:
     ../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md

   Status: skeleton. Every cell marked (* PENDING *) is waiting on
   the Wolfram MCP / wolframscript license to be wired in. The cells
   are written to run standalone once a Mathematica kernel is
   available.

   Conventions:
     - Gaussian (c.g.s.) units throughout.
     - 3-vectors are symbolic; we use dot products explicitly so we
       can stay component-free where possible.
     - b[\[Tau]] := Sqrt[c^2 + u[\[Tau]].u[\[Tau]]] is the collaborative
       speed of light.
================================================================ *)


(* -- Eq. (1) ----------------------------------------------------- *)
(* Claim: w/c = u/b, with b = Sqrt[c^2 + u^2] and u = w/Sqrt[1-w^2/c^2]. *)

(* PENDING *)
ClearAll[w, u, c];
With[{b = Sqrt[c^2 + u^2]},
  Module[{wOfU},
    wOfU = w /. First @ Solve[u == w/Sqrt[1 - w^2/c^2] && w > 0 && c > 0 && u > 0, w];
    Print["Eq.(1): ", FullSimplify[wOfU/c - u/b, Assumptions -> {c > 0, u > 0}]]
  ]
];


(* -- Eq. (2) ----------------------------------------------------- *)
(* Claim: (1/c) \[PartialD]_t = (1/b) \[PartialD]_\[Tau]. *)

(* PENDING *)
With[{b = Sqrt[c^2 + u^2]},
  Print["Eq.(2): ", Simplify[(1/c)*(c/b) - 1/b]]
];


(* -- Eq. (3') ---------------------------------------------------- *)
(* Claim: substituting (1) and (2) into standard Maxwell's eqs gives
   the proper-time form. We confirm at the substitution-rule level. *)

(* PENDING: rule-rewriting check, no symbolic field manipulation needed. *)


(* -- Eq. (4) ----------------------------------------------------- *)
(* Claim: the wave equation for B (likewise E) in proper time has
   the form
     (1/b^2) \[PartialD]^2_\[Tau] B - (u.a / b^4) \[PartialD]_\[Tau] B - \[Del]^2 B = (1/b) [4 \[Pi] \[Del] x (\[Rho] u)].
   We verify the coefficient (u.a)/b^4 from d/d\[Tau](1/b). *)

(* PENDING *)
ClearAll[u, a, c, \[Tau]];
bF[\[Tau]_] := Sqrt[c^2 + u[\[Tau]] . u[\[Tau]]];
Print[
  "Eq.(4) coefficient check: ",
  FullSimplify[
    D[1/bF[\[Tau]], \[Tau]] /. {u'[\[Tau]] -> a[\[Tau]]}
      + (u[\[Tau]] . a[\[Tau]])/bF[\[Tau]]^3
  ]
];
(* Expected: 0. So d/d\[Tau](1/b) = -(u.a)/b^3, which contributes
   -(u.a)/b^4 after multiplication by 1/b in the curl-of-curl. *)


(* -- Eq. (5) ----------------------------------------------------- *)
(* Claim: \[Psi] = (b/c)^{1/2} \[CapitalPsi]_new kills the first-order \[PartialD]_\[Tau]
   in the damped wave equation and produces an effective mass term. *)

(* PENDING *)
ClearAll[b, \[Tau], \[Psi]new];
trans[\[Tau]_] := Sqrt[b[\[Tau]]/c] \[Psi]new[\[Tau]];
expr = (1/b[\[Tau]]^2) D[trans[\[Tau]], {\[Tau], 2}]
     - (b'[\[Tau]]/b[\[Tau]]^3) D[trans[\[Tau]], \[Tau]];
expanded = Expand[expr/Sqrt[b[\[Tau]]/c]];
massCoeffPredicted = b''[\[Tau]]/(2 b[\[Tau]]^3) - 3 b'[\[Tau]]^2/(4 b[\[Tau]]^4);
(* Print coefficient of \[Psi]new[\[Tau]] (no derivatives): *)
Print[
  "Eq.(5) effective-mass coefficient: ",
  Simplify[Coefficient[expanded, \[Psi]new[\[Tau]]] - massCoeffPredicted]
];
(* Expected: 0. *)


(* -- Eq. (6) ----------------------------------------------------- *)
(* Claim:
     b''/(2 b^3) - 3 b'^2/(4 b^4)
   equals
     (u . u'' + u' . u')/(2 b^4) - 5 (u . u')^2/(4 b^6).
   Use b = Sqrt[c^2 + u.u]. *)

(* PENDING *)
ClearAll[u, \[Tau], c];
bF2[\[Tau]_] := Sqrt[c^2 + u[\[Tau]] . u[\[Tau]]];
lhs = D[bF2[\[Tau]], {\[Tau], 2}]/(2 bF2[\[Tau]]^3) - 3 D[bF2[\[Tau]], \[Tau]]^2/(4 bF2[\[Tau]]^4);
rhs = (u[\[Tau]] . u''[\[Tau]] + u'[\[Tau]] . u'[\[Tau]])/(2 bF2[\[Tau]]^4)
       - 5 (u[\[Tau]] . u'[\[Tau]])^2/(4 bF2[\[Tau]]^6);
Print["Eq.(6): ", FullSimplify[lhs - rhs]];
(* Expected: 0. *)


(* -- Eq. (7) ----------------------------------------------------- *)
(* Modified Lienard-Wiechert E, B fields. Multi-page derivation in [18];
   here we will:
     (a) confirm that the third term vanishes when a = 0,
     (b) confirm Series expansion in u/c reproduces standard
         Lienard-Wiechert at leading order.
*)
(* PENDING: full derivation deferred. *)


(* -- Eqs. (8), (9) ----------------------------------------------- *)
(* (8) textbook. (9) integrates dt/d\[Tau] = b/c. *)

(* PENDING *)
ClearAll[bSym, s, \[Tau], c];
Print[
  "Eq.(9): ",
  FullSimplify[
    (1/c) Integrate[bSym[s], {s, 0, \[Tau]}]
      - ((1/\[Tau]) Integrate[bSym[s], {s, 0, \[Tau]}]) \[Tau]/c
  ]
];
(* Expected: 0. *)


(* -- Eq. (10) ---------------------------------------------------- *)
(* Position, velocity, acceleration transformations under proper-time
   group. Verify by chain-rule differentiation of Eq. (8) under Eq. (9). *)

(* PENDING: composite check; details in companion markdown. *)


(* -- Eq. (11) ---------------------------------------------------- *)
(* Differentiate bp \[Tau]/c = \[Gamma] (b \[Tau]/c - x.v/c^2) in \[Tau]
   to obtain bp = \[Gamma](b - u.v/c). *)

(* PENDING *)
ClearAll[bp, b, \[Tau], v, c, x, \[Gamma]];
(* Treat x as a function of \[Tau] with dx/d\[Tau] = u(\[Tau]) and b, bp as functions of \[Tau]: *)
ClearAll[bf, bpf, xf, uf];
relation = D[bpf[\[Tau]] \[Tau]/c, \[Tau]]
           - D[\[Gamma] (bf[\[Tau]] \[Tau]/c - xf[\[Tau]] . v/c^2), \[Tau]];
(* Substitute dx/d\[Tau] = uf[\[Tau]]: *)
relation2 = relation /. {xf'[\[Tau]] -> uf[\[Tau]]};
(* Solve for bpf[\[Tau]]: *)
sol = Solve[relation2 == 0, bpf[\[Tau]]];
Print["Eq.(11): bp[\[Tau]] = ", bpf[\[Tau]] /. First @ sol // Simplify];
(* Expected: \[Gamma] (b[\[Tau]] - u[\[Tau]] . v / c)  (after appropriate constants of integration). *)


(* -- Eqs. (12)-(15) ---------------------------------------------- *)
(* Charge / current density transformations.
   Algorithm:
     (a) Combine (11) with (13) to derive (14).
     (b) Substitute J/c = \[Rho] u/b into (14) to derive (15).
*)

(* PENDING *)
ClearAll[b, bp, \[Rho], \[Rho]p, J, u, v, c, \[Gamma]];
solRho = Solve[{bp == \[Gamma] (b - u . v/c),
                bp \[Rho]p == \[Gamma] (b \[Rho] - J . v/c)}, \[Rho]p][[1]];
rhoP14 = Simplify[\[Rho]p /. solRho];
target14 = (\[Rho] - (J . v)/(b c))/(1 - (u . v)/(b c));
Print["Eq.(14): ", FullSimplify[rhoP14 - target14]];

rhoP15 = Simplify[rhoP14 /. J -> \[Rho] u (c/b)];
target15 = \[Rho] (1 - (u . v)/b^2)/(1 - (u . v)/(b c));
Print["Eq.(15): ", FullSimplify[rhoP15 - target15]];
(* Expected: 0, 0. *)


(* -- Eqs. (16)-(18) ---------------------------------------------- *)
(* Canonical proper-time Hamiltonian K = H^2/(2 m c^2) + m c^2 / 2
   and the equations of motion that follow. Use Hamilton's equations
   directly. *)

(* PENDING: requires Poisson-bracket setup. Plan: define K as a
   symbolic function of pi, x, V and compute \[PartialD]K/\[PartialD]pi, \[PartialD]K/\[PartialD]x. *)


(* -- Eqs. (19)-(21) ---------------------------------------------- *)
(* Many-particle generalizations. PENDING. *)


(* -- Eqs. (22)-(24) ---------------------------------------------- *)
(* Quantum proper-time wave equations. PENDING; will use Dirac-matrix
   algebra (or FeynCalc) to verify Pauli/Dirac identities. *)


(* End of file *)
