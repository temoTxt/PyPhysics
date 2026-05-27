# /loop prompt — Li²⁺ Z-extension (issue #78)

Copy the block below verbatim into a fresh Claude Code session in this folder, after starting `claude` from the VSCode extension. The session must remain open until the loop terminates.

---

```
/loop Advance the Li²⁺ Z-extension (issue #78) one substantive step per the iteration protocol in .dev/research/brief.md. Read .dev/research/STATE.md to find the last iteration's pending-next; do exactly one step (read a source-of-record document, scaffold or extend r_e_Li2plus_joint_fit.wl, execute/debug a Wolfram cell, draft a per-observable result section in Bethe_Salpeter/11_Li2plus_HydrogenicIon.md, or compute the joint χ² at Z=3 once all four predictions are in); append a dated entry to STATE.md recording what advanced, what's queued next, the current observable focus (#1 g-factor / #2 Lamb shift / #3 fine structure / #4 hyperfine), the current outcome-matrix branch (A/B/C/D), any (Z-i)/(Z-ii) reading differences, and any BLOCKED state with measurement provenance for any observable values added; commit with a one-line message starting "li2plus:"; push to origin/78-bethe-salpeter-z-extension-li2plus; then ScheduleWakeup for the next iteration with the aggressive cadence the user has requested — 60-180s (cache-warm; favor 60s, the runtime floor, unless the next step is explicitly idle-pending). Stop the loop only when all four observables are checked + a Z-axis verdict is recorded in STATE.md, OR a hard BLOCKED state requiring Tepper input is recorded. Do not open PRs or post GitHub comments — orchestrator handles those after morning review.
```

---

## Notes

- Permission allowlist now has `Bash(git *)` and `Bash(gh *)` (broadened 2026-05-27 to reduce friction). Loop should not be prompt-blocked on routine git/gh operations.
- Wolfram MCP available for symbolic + numerical work; obey CLAUDE.md gotchas (single-line cells; `ee` / `potV`; non-commutative `Dot`).
- To peek at progress remotely: `git fetch origin && git log origin/78-bethe-salpeter-z-extension-li2plus --oneline -20`.
- To kill the loop cleanly: send `stop the loop` as a user message; the next iteration will not be scheduled.
