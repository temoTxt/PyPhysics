# `pdg.*` namespace — refresh log

Chronological record of refreshes, edition bumps, and citation diffs for the `pdg.*` namespace. Per the umbrella [`REFRESH_POLICY.md`](../../../REFRESH_POLICY.md) §7.

## 2026-05-30 — initial seed (issue #96)

- **Cause:** initial implementation of the `pdg.*` namespace per umbrella issue [#92](https://github.com/temoTxt/PyPhysics/issues/92) sub-issue [#96](https://github.com/temoTxt/PyPhysics/issues/96).
- **PDG edition seeded:** 2024 Review of Particle Physics (`pdg_edition: "2024"`).
- **Particles seeded:** electron (e⁻, MCID 11), muon (μ⁻, 13), tau (τ⁻, 15), proton (p⁺, 2212), neutron (n, 2112), charged pion (π⁺, 211), charged kaon (K⁺, 321), W boson (24), Z boson (23), Higgs boson (25). Anti-particles inferred by sign-flip from the corresponding particle entry; not yet seeded explicitly.
- **Anomalies seeded:** electron g-2 (single-value plus SM prediction) and muon g-2 (three values: Fermilab E989 experimental, BMW lattice-QCD, Theory Initiative 2020 data-driven SM).
- **Disagreement cases surfaced:** (a) neutron lifetime bottle-vs-beam puzzle at ~4σ; (b) muon g-2 lattice-vs-experiment-vs-data-driven tension at ~5σ as of 2024.
- **Bib stubs scaffolded:**
  - `pdg2024_review` (umbrella's non-paper-source cite-key convention)
  - `fan2023_electron_g_factor`
  - `fermilab_e989_2023`
  - `bmw2021_lattice_qcd`
  - `theory_initiative_2020`
- **Caller updates:** none in this PR; demonstrator refactor of [`Bethe_Salpeter/10_CrossComparison.md`](../../../../../Roadmapping/Quantum_Mechanics/Bethe_Salpeter/10_CrossComparison.md) deferred to issue [#97](https://github.com/temoTxt/PyPhysics/issues/97)'s demonstrator (which already touches that document for the `qed.*` migration; folding the `pdg.*` muon-mass reference into the same diff is cleaner).

Subsequent entries will track:
- PDG 2025 / 2026 edition bumps as they ship.
- Theory Initiative 2025+ update of the data-driven SM prediction for a_μ (expected to reduce the lattice-vs-data-driven tension).
- Additional particles added by follow-up PRs (anti-particles, neutral pion/kaon, J/ψ, Υ, top, etc.).
- Eventual switch from hand-curated seed to live PDG JSON API once acceptance §3 of [#96](https://github.com/temoTxt/PyPhysics/issues/96) is reached.
