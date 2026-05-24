# Jackson, *Classical Electrodynamics* — per-chapter status

Per-chapter status table for the Jackson canonical-problems campaign. Each row points to the chapter document and records the PR in which that chapter is worked, the unit-system regime per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), and the count of problems completed vs. planned.

The thread's README is at [`../README.md`](../README.md). The plan is at [`.dev/tasks/42-electromagnetism-jackson-proper-time.md`](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md).

## Status table

| Chapter | Title | Unit system | PR | Problems | Status |
|---|---|---|---|---|---|
| Ch. 1 | Introduction to Electrostatics | three-system (3e Chs. 1–10) | PR 0, then PR F+ | 1 / 2–6 | PR 0 in progress — see [`Ch01`](Ch01_Introduction_Electrostatics.md) |
| Ch. 2 | Boundary-Value Problems in Electrostatics I | three-system | PR 0, then PR F+ | 1 / 2–6 | PR 0 in progress — see [`Ch02`](Ch02_Boundary_Value_Problems_Electrostatics_I.md) |
| Ch. 3 | Boundary-Value Problems in Electrostatics II | three-system | PR F+ | 0 / 4–6 | planned |
| Ch. 4 | Multipoles, Electrostatics of Macroscopic Media | three-system | PR F+ | 0 / 4–6 | planned |
| Ch. 5 | Magnetostatics, Faraday's Law, Quasi-Static Fields | three-system | PR 0, then PR F+ | 2 / 2–6 | PR 0 drafted — see [`Ch05`](Ch05_Magnetostatics_Faraday_Quasi_Static.md) |
| Ch. 6 | Maxwell Equations + Macroscopic Media ⭐ | three-system | **PR A** | 5 / 5 | PR A drafted — see [`Ch06`](Ch06_Maxwell_Equations_Macroscopic_Media.md) |
| Ch. 7 | Plane Electromagnetic Waves, Wave Propagation | three-system | PR F+ | 0 / 4–6 | planned |
| Ch. 8 | Waveguides, Resonant Cavities, Optical Fibers | three-system | PR F+ | 0 / 4–6 | planned |
| Ch. 9 | Radiating Systems, Multipole Fields | three-system | PR F+ | 0 / 4–6 | planned |
| Ch. 10 | Scattering and Diffraction | three-system | PR F+ | 0 / 4–6 | planned |
| Ch. 11 | Special Theory of Relativity | CGS only (3e Chs. 11+) | **PR B** | 0 / 5–7 | planned |
| Ch. 12 | Dynamics of Relativistic Particles and EM Fields ⭐ | CGS only | **PR C** | 0 / 5–7 | planned |
| Ch. 13 | Collisions, Energy Loss, Scattering of Charged Particles | CGS only | PR F+ | 0 / 4–6 | planned |
| Ch. 14 | Radiation by Moving Charges ⭐ | CGS only | **PR D** | 0 / 5–7 | planned |
| Ch. 15 | Bremsstrahlung, Method of Virtual Quanta | CGS only | PR F+ | 0 / 4–6 | planned |
| Ch. 16 | Radiation Damping, Classical Models of Charged Particles ⭐ | CGS only | **PR E** | 0 / 3–5 | planned |

⭐ indicates a chapter where the proper-time vs classical contrast is most informative; see [§7 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list).

## Experimental comparisons

The [`Experimental_Comparisons/`](Experimental_Comparisons/) subdirectory holds cross-experiment analysis documents — predictions of the proper-time framework set against published experimental data. These are *not* per-problem documents in the sense of [`../_template_problem.md`](../_template_problem.md); they have their own structure (a §1 selection provenance, a §2 framework restatement, a §3 per-experiment block, and a §4 cross-experiment summary), but they follow the same voice and Crocco-substantive gating conventions as the chapter documents.

| Comparison document | Tracks | Status |
|---|---|---|
| [`Experimental_Comparisons/radiation_reaction_2018.md`](Experimental_Comparisons/radiation_reaction_2018.md) | issue [#43](https://github.com/temoTxt/PyPhysics/issues/43); Cole 2018, Poder 2018, Wistisen 2018 | draft — substantive blocks blank, pending Ch. 14 canonical derivation from PR D |

## Chapter contributions to PR 0

The PR 0 fluency-warm-up spans Chs. 1–2 and 5 with 3–4 problems where `u = 0` (electrostatics) or `∂_τ = 0` (steady currents), so that the proper-time formulation reduces to the classical answer up to at most a `J → (b/c) J` rescaling. Candidate problems are listed in [§7.3 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#73-pr-a-prequel-adopted); the final selection is made when PR 0 opens.

## Per-chapter document conventions

Each chapter document (`ChNN_*.md`) opens with a short header recording the chapter's PR, its unit-system regime, and the list of problems contained; the body is one entry per problem, structured by [`../_template_problem.md`](../_template_problem.md). Reviewers reading any single chapter document should not need to consult the plan to understand what they are looking at.
