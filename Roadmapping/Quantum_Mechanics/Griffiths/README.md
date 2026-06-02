# Griffiths, *Introduction to Quantum Mechanics* — per-chapter status

Per-chapter status table for the Griffiths canonical-problems campaign. Each row records the PR in which that chapter is worked and the count of problems completed vs. planned.

The thread's README is at [`../README.md`](../README.md). The plan is at [`.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md`](../../../.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md).

## Status table

| Chapter (3e) | Title | PR | Problems | Status |
|---|---|---|---|---|
| Ch. 1 | The Wave Function | **PR A** | 3 / 3 | ✅ drafted — see [`Ch01`](Ch01_Wave_Function.md) |
| Ch. 2 | Time-Independent Schrödinger Equation | **PR B** | 4 / 4 | ✅ drafted — see [`Ch02`](Ch02_TI_Schrodinger.md) |
| Ch. 3 | Formalism | **PR C** | 3 / 3 | ✅ drafted — see [`Ch03`](Ch03_Formalism.md) |
| Ch. 4 ⭐ | Quantum Mechanics in Three Dimensions (hydrogen) | **PR D** | 5 / 5 | ✅ drafted — see [`Ch04`](Ch04_QM_3D.md) |
| Ch. 5 | Identical Particles | **PR E** | 3 / 3 | ✅ drafted — see [`Ch05`](Ch05_Identical_Particles.md) |
| Ch. 6 | Symmetries and Conservation Laws | **PR F** | 3 / 3 | ✅ drafted — see [`Ch06`](Ch06_Symmetries_Conservation.md) |
| Ch. 7 ⭐ | Time-Independent Perturbation Theory (fine structure) | **PR G** | 5 / 5 | ✅ drafted — see [`Ch07`](Ch07_TI_Perturbation_Theory.md) |
| Ch. 8 | Variational Principle (helium) | **PR H** | 3 / 3 | ✅ drafted — see [`Ch08`](Ch08_Variational_Principle.md) |
| Ch. 9 | WKB Approximation | **PR I** | 3 / 3 | ✅ drafted — see [`Ch09`](Ch09_WKB_Approximation.md) |
| Ch. 10 | Scattering | **PR J** | 4 / 4 | ✅ drafted — see [`Ch10`](Ch10_Scattering.md) |
| Ch. 11 | Quantum Dynamics (3e renumber) | **PR K** | 4 / 4 | ✅ drafted — see [`Ch11`](Ch11_Quantum_Dynamics.md) |
| Ch. 12 | Afterword (3e renumber) | **PR L** | 2 / 2 | ✅ drafted — see [`Ch12`](Ch12_Afterword.md) |

⭐ indicates a chapter where the proper-time formulation has the most to say; see [§4 of the plan](../../../.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md#4-pr-sequencing-12-prs-al-per-issue-body).

## Notes on edition handling

Per [§2 of the plan](../../../.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md#2-repository-structure), the campaign uses 3e numbering as canonical with 2e cross-references in the source line. Chs. 11/12 are renumbered between editions and require explicit re-statement of which 2e chapter corresponds to which 3e chapter. The per-problem template is the canonical reference for this convention.

## Per-chapter document conventions

Each chapter document (`ChNN_*.md`) opens with a short header recording the chapter's PR, the relevant proper-time content level, and the list of problems contained. The body is one entry per problem, structured by [`../_template_problem.md`](../_template_problem.md).
