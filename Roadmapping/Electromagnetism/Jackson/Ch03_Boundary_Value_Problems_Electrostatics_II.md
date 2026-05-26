# Ch. 3 — Boundary-Value Problems in Electrostatics II

This chapter contains Jackson canonical problems on more advanced boundary-value electrostatics (separation of variables in spherical / cylindrical geometries, Green's function methods). Part of **PR K (final backfill)** per [§7 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list). All Ch. 3 problems are purely electrostatic; `u = 0`, `b = c`; proper-time reduces identically to classical. Compact format.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P3.1 — Conducting sphere in uniform external field](#problem-j3e-p31--conducting-sphere-in-uniform-external-field) | drafted (PR K) | fluency-builder |
| [Problem J3e-P3.4 — Dielectric sphere in uniform external field](#problem-j3e-p34--dielectric-sphere-in-uniform-external-field) | drafted (PR K) | fluency-builder |
| [Problem J3e-P3.6 — Spherical-harmonic expansion of a point charge](#problem-j3e-p36--spherical-harmonic-expansion-of-a-point-charge) | drafted (PR K) | fluency-builder |
| [Problem J3e-P3.8 — Green's function for grounded sphere](#problem-j3e-p38--greens-function-for-grounded-sphere) | drafted (PR K) | fluency-builder |

---

### Problem J3e-P3.1 — Conducting sphere in uniform external field

**Statement:** A conducting sphere of radius `R` is placed in a previously-uniform external field `\mathbf{E}_{0} = E_{0}\hat z`. Compute the induced surface charge distribution and the dipole moment of the polarised sphere.

**Classical solution:** Separation of variables with boundary condition `\phi(R) = 0`. Surface charge `\sigma(\theta) = (3E_{0}/(4\pi))\cos\theta`. Induced dipole moment `p = E_{0}R^{3}`. Field outside is `\mathbf{E}_{0}` plus pure-dipole field of `p = E_{0}R^{3}`.

**Proper-time:** Static, `u = 0`. Identical.

**Verdict:** ✅. **Companion:** [`JacksonCh03_P3_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh03_P3_1.wl).

---

### Problem J3e-P3.4 — Dielectric sphere in uniform external field

**Statement:** Replace the conductor of J3e-P3.1 with a linear dielectric of relative permittivity `\varepsilon`. Compute the field inside and the dipole moment.

**Classical solution:** The internal field is `\mathbf{E}_{\text{in}} = 3\mathbf{E}_{0}/(\varepsilon + 2)`. The induced dipole moment is `p = (\varepsilon - 1)/(\varepsilon + 2)\cdot R^{3}E_{0}` — the *Clausius–Mossotti* formula.

**Proper-time:** Static. Identical.

**Verdict:** ✅. **Companion:** [`JacksonCh03_P3_4.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh03_P3_4.wl).

---

### Problem J3e-P3.6 — Spherical-harmonic expansion of a point charge

**Statement:** Expand the potential `1/|\mathbf{r} - \mathbf{r}'|` in spherical harmonics centred at the origin.

**Classical solution:**

$$
\frac{1}{|\mathbf{r} - \mathbf{r}'|} = \sum_{\ell = 0}^{\infty}\frac{r_{<}^{\ell}}{r_{>}^{\ell+1}}P_{\ell}(\cos\gamma),
$$

with `r_{<}, r_{>}` the smaller / larger of `|r|, |r'|`, and `\gamma` the angle between them. Generalisable to `Y_{\ell m}` form via the addition theorem.

**Proper-time:** Static. Identical.

**Verdict:** ✅. **Companion:** [`JacksonCh03_P3_6.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh03_P3_6.wl).

---

### Problem J3e-P3.8 — Green's function for grounded sphere

**Statement:** Derive the Green's function `G(\mathbf{r}, \mathbf{r}')` for the Dirichlet problem `\phi = 0` on a grounded sphere of radius `R`.

**Classical solution:** Image-charge construction gives

$$
G(\mathbf{r}, \mathbf{r}') = \frac{1}{|\mathbf{r} - \mathbf{r}'|} - \frac{R/r'}{|\mathbf{r} - (R^{2}/r'^{2})\mathbf{r}'|},
$$

with the image charge at the inverse-radius point.

**Proper-time:** Static. Identical.

**Verdict:** ✅. **Companion:** [`JacksonCh03_P3_8.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh03_P3_8.wl).
