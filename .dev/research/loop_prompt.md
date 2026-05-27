# /loop prompt — Li²⁺ spectroscopy (Lamb shift + fine structure), issue #78

Copy the block below verbatim into a fresh Claude Code session in this folder, after starting `claude` from the VSCode extension. The session must remain open until the loop terminates.

---

```
/loop Advance the Li²⁺ spectroscopy branch (issue #78, Lamb shift + fine structure scope) one substantive step per the iteration protocol in .dev/research/brief.md. Prioritize observable #3 fine structure (real Z-axis discriminator) over #2 Lamb shift (weak discriminator per the parallel Self-Energy branch's iter-1 finding). Read .dev/research/STATE.md to find the last iteration's pending-next; do exactly one step (read a source-of-record document not yet read, scaffold or extend r_e_Li2plus_spectroscopy.wl, execute/debug a Wolfram cell, or draft a per-observable section in Bethe_Salpeter/12_Li2plus_Spectroscopy.md); append a dated entry to STATE.md recording what advanced, what's queued next, current observable focus (#2 or #3), outcome-matrix branch (A/B/C/D), (Z-i)/(Z-ii) reading differences if any, and any BLOCKED state with measurement provenance; commit with a one-line message starting "li2plus-spec:"; push to origin/78-li2plus-spectroscopy; then ScheduleWakeup for the next iteration with the aggressive cadence the user has requested — 60-180s (cache-warm; favor 60s, the runtime floor, unless the next step is explicitly idle-pending). DO NOT edit files owned by other branches: Bethe_Salpeter/11_Li2plus_HydrogenicIon.md, Bethe_Salpeter/13_Li2plus_Hyperfine.md, r_e_Li2plus_joint_fit.wl belong to Self-Energy and Hyperfine branches respectively. Stop the loop only when both #2 and #3 have framework predictions documented in 12_Li2plus_Spectroscopy.md with verdicts, OR a hard BLOCKED state is recorded. Do not open PRs or post GitHub comments — orchestrator handles those after morning review.
```

---

## Notes

- Permission allowlist has `Bash(git *)` and `Bash(gh *)` from the broadened settings (2026-05-27). Loop should not be prompt-blocked on routine git/gh operations.
- Wolfram MCP available for symbolic + numerical work; obey CLAUDE.md gotchas.
- To peek at progress remotely: `git fetch origin && git log origin/78-li2plus-spectroscopy --oneline -20`.
- To kill the loop cleanly: send `stop the loop` as a user message; the next iteration will not be scheduled.
