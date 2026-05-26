# /loop prompt — Candidate 3 (issue #66)

Copy the block below verbatim into a fresh Claude Code session in this folder, after starting `claude` from the VSCode extension. The session must remain open overnight.

---

```
/loop Advance candidate-3 (issue #66) one substantive step per the iteration protocol in .dev/research/brief.md. Focus is the empirical-test path (extending the PR #62 joint fit with higher-precision observables and testing whether the residual to the Schwinger closed-form tracks Karplus-Kroll two-loop QED corrections). Read .dev/research/STATE.md to find the last iteration's pending-next; do exactly one step; append a dated entry to STATE.md recording what advanced, what's queued next, the current outcome-matrix branch (A/B/C/D), every added observable's measurement provenance, and any BLOCKED state plus the running "Questions for Tepper" queue; commit with a one-line message starting "candidate-3:"; push to origin/66-theory-candidate-3-closed-form-schwinger-identification-of-the-triangulated-r_e; then ScheduleWakeup for the next iteration (1200-1800s for substantive math/data-curation work, 60-270s for quick continuations). Stop the loop only when the empirical-test path's acceptance criteria can be checked, OR STATE.md records a BLOCKED state. Do not open PRs or post GitHub comments — orchestrator handles morning review and lifts the Tepper-question queue to a #66 comment.
```

---

## Notes

- If permission prompts interrupt the loop, run `/permissions` and add the offending command to the allowlist, then resume with a fresh `/loop` invocation.
- To peek at progress remotely: from any machine, `git fetch origin && git log origin/66-theory-candidate-3-closed-form-schwinger-identification-of-the-triangulated-r_e --oneline -20`.
- To kill the loop cleanly: send `stop the loop` as a user message; the next iteration will not be scheduled.
