# /loop prompt — Hydrogenic Z-scan g-factor (issue #82)

Copy the block below verbatim into a fresh Claude Code session in this folder, after starting `claude` from the VSCode extension. The session must remain open until the loop terminates.

---

```
/loop Advance the hydrogenic Z-scan g-factor branch (issue #82) one substantive step per the iteration protocol in .dev/research/brief.md. Read .dev/research/STATE.md to find the last iteration's pending-next; do exactly one step (read a source-of-record document not yet read, look up + record the precise g-factor measurement value + DOI/year for one ion in the Z-scan, scaffold or extend r_e_Zscan_fit.wl, execute/debug a Wolfram cell, draft a per-ion section in Bethe_Salpeter/14_HydrogenicIon_Zscan.md, or — once all ions are catalogued — execute the joint χ² + Z-scaling fit); append a dated entry to STATE.md recording what advanced, what's queued next, current ion focus (He⁺/Be³⁺/C⁵⁺/Si¹³⁺/Ca¹⁹⁺ or "joint fit"), outcome-matrix branch tentative (A/B/C/D), and any BLOCKED state with full measurement provenance; commit with a one-line message starting "zscan:"; push to origin/82-hydrogenic-z-scan-g-factor; then ScheduleWakeup for the next iteration with the aggressive cadence the user has requested — 60-180s (cache-warm; favor 60s, the runtime floor, unless the next step is explicitly idle-pending). DO NOT edit files owned by #78 branches: 11_Li2plus_HydrogenicIon.md, 12_Li2plus_Spectroscopy.md, 13_Li2plus_Hyperfine.md, r_e_Li2plus_*.wl all belong to parallel branches. Import the Li²⁺ (Z=3) g-factor value from the parallel Self-Energy branch's work; do not re-derive it. Stop the loop only when 5+ ions are catalogued + joint χ² + Z-scaling fit is reported + a Z-axis verdict is recorded, OR a hard BLOCKED state is recorded. Do not open PRs or post GitHub comments — orchestrator handles those after morning review.
```

---

## Notes

- Permission allowlist has `Bash(git *)` and `Bash(gh *)` from the broadened settings (2026-05-27). Loop should not be prompt-blocked on routine git/gh operations.
- Wolfram MCP available for symbolic + numerical work; obey CLAUDE.md gotchas.
- To peek at progress remotely: `git fetch origin && git log origin/82-hydrogenic-z-scan-g-factor --oneline -20`.
- To kill the loop cleanly: send `stop the loop` as a user message; the next iteration will not be scheduled.
