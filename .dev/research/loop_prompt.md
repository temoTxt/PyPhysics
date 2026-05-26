# /loop prompt — Candidate 1 (issue #64)

Copy the block below verbatim into a fresh Claude Code session in this folder, after starting `claude` from the VSCode extension. The session must remain open overnight.

---

```
/loop Advance candidate-1 (issue #64) one substantive step per the iteration protocol in .dev/research/brief.md. Read .dev/research/STATE.md to find the last iteration's pending-next; do exactly one step; append a dated entry to STATE.md recording what advanced, what's queued next, the current outcome-matrix branch (A/B/C/D), and any BLOCKED state; commit with a one-line message starting "candidate-1:"; push to origin/64-theory-candidate-1-proper-time-self-energy-integral-derivation-of-r_e; then ScheduleWakeup for the next iteration (1200-1800s for substantive math, 60-270s for quick continuations). Stop the loop only when all acceptance criteria in issue #64 can be checked, OR STATE.md records a BLOCKED state requiring author input. Do not open PRs or post GitHub comments — orchestrator handles morning review.
```

---

## Notes

- If permission prompts interrupt the loop, run `/permissions` and add the offending command to the allowlist, then resume with a fresh `/loop` invocation.
- To peek at progress remotely: from any machine, `git fetch origin && git log origin/64-theory-candidate-1-proper-time-self-energy-integral-derivation-of-r_e --oneline -20`.
- To kill the loop cleanly: send `stop the loop` as a user message; the next iteration will not be scheduled.
