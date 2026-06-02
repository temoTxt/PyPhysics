# `flag.*` namespace — refresh log

## 2026-05-30 — initial seed (issue #98)

- **Cause:** initial implementation of the `flag.*` namespace per umbrella issue [#92](https://github.com/temoTxt/PyPhysics/issues/92) sub-issue [#98](https://github.com/temoTxt/PyPhysics/issues/98).
- **FLAG edition seeded:** FLAG 2024 (`flag_edition: "2024"`, *Eur. Phys. J. C* 84).
- **Quantities seeded:** `m_pi_charged`, `f_pi`, `f_K`, `f_K_over_f_pi`, `g_A`, `sigma_piN`.
- **Bib stub scaffolded:** `flag2024_averages` (issuing-body cite-key convention).
- **Smoke test:** `flag.get_quantity("f_pi")` returns 130.2(8) MeV with `flag_edition: "2024"` source-revision marker.

Subsequent entries:
- FLAG 2026/2027 edition bumps as they ship.
- Additional quantities (CKM elements, B-meson decay constants, quark masses by scheme) added by follow-up PRs.
