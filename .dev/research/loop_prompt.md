# /loop prompt — Li²⁺ hyperfine, issue #78

Copy the block below verbatim into a fresh Claude Code session in this folder, after starting `claude` from the VSCode extension. The session must remain open until the loop terminates.

---

```
/loop Advance the Li²⁺ hyperfine branch (issue #78, observable #4 only) one substantive step per the iteration protocol in .dev/research/brief.md. Read .dev/research/STATE.md to find the last iteration's pending-next; do exactly one step (read a source-of-record document not yet read, scaffold or extend r_e_Li2plus_hyperfine.wl, execute/debug a Wolfram cell, or draft the per-observable section in Bethe_Salpeter/13_Li2plus_Hyperfine.md); append a dated entry to STATE.md recording what advanced, what's queued next, observable focus (always #4), outcome-matrix branch (A/B/C/D), (Z-i)/(Z-ii) reading differences if any, and any BLOCKED state with measurement provenance; commit with a one-line message starting "li2plus-hfs:"; push to origin/78-li2plus-hyperfine; then ScheduleWakeup for the next iteration with the aggressive cadence the user has requested — 60-180s (cache-warm; favor 60s, the runtime floor, unless the next step is explicitly idle-pending). DO NOT edit files owned by other branches: Bethe_Salpeter/11_Li2plus_HydrogenicIon.md, Bethe_Salpeter/12_Li2plus_Spectroscopy.md, r_e_Li2plus_joint_fit.wl, r_e_Li2plus_spectroscopy.wl belong to Self-Energy and Spectroscopy branches respectively. Stop the loop only when #4 has a framework prediction documented in 13_Li2plus_Hyperfine.md with a verdict, OR a hard BLOCKED state is recorded. Pay extra attention to the Li-7 I=3/2 nuclear-spin extension of BS-§22's H I=1/2 apparatus — that's the substantive-AI move on this branch. Do not open PRs or post GitHub comments — orchestrator handles those after morning review.
```

---

## Notes

- Permission allowlist has `Bash(git *)` and `Bash(gh *)` from the broadened settings (2026-05-27). Loop should not be prompt-blocked on routine git/gh operations.
- Wolfram MCP available for symbolic + numerical work; obey CLAUDE.md gotchas.
- To peek at progress remotely: `git fetch origin && git log origin/78-li2plus-hyperfine --oneline -20`.
- To kill the loop cleanly: send `stop the loop` as a user message; the next iteration will not be scheduled.
