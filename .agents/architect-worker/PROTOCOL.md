# Architect-Worker Protocol

Protocol version: `architect-worker-v1`

## Purpose

Codex is the architect and manager. Claude Code is a bounded worker. The worker
may investigate, implement, and report, but it cannot approve its own work or
declare the user's task complete.

The protocol exists to preserve four things:

1. the user's objective without reinterpretation;
2. the repository's architectural laws and sealed results;
3. the complete evidence trail, including failures and withdrawn readings;
4. independent review between implementation and acceptance.

## Authority

```text
user       owns objective, taxonomy decisions, and requested scope
Codex      owns decomposition, delegation, review, verification, integration,
           escalation, and the claim of completion
Claude     owns only the execution of one bounded task packet and its report
```

Claude returns a candidate, never a verdict. A worker-reported `PASS` is a
claim requiring independent reproduction by Codex.

Claude must not:

- change the objective or acceptance laws;
- reinterpret a sealed tower result;
- invent taxonomy when the admitted vocabulary is insufficient;
- modify files outside its assigned worktree;
- commit, push, rebase, reset, merge, create branches, or manage worktrees;
- contact the user directly;
- mark a task `VERIFIED`, `INTEGRATED`, or `COMPLETE`.

## State machine

```text
DRAFT
  -> DELEGATED
  -> RUNNING
  -> CANDIDATE ---------> REVIEWING -> VERIFIED -> INTEGRATED -> COMPLETE
        |                     |             |          |
        |                     +-> REVISION_REQUESTED <-+
        |                               |
        +-------------------------------+-> RUNNING

RUNNING -> BLOCKED -> REVISION_REQUESTED or ESCALATED
RUNNING -> FAILED  -> REVISION_REQUESTED, REJECTED, or ESCALATED
any review state  -> REJECTED or ESCALATED
```

Only Claude may produce `CANDIDATE`, `BLOCKED`, or `FAILED` reports. Only Codex
may create every state to the right of `REVIEWING`.

For a read-only task, `VERIFIED -> COMPLETE` is permitted because there is no
candidate patch to integrate.

## Immutable task packet

Every delegation begins with a JSON task packet made from
`task.template.json`. It contains:

```text
task_id                 stable identity for the whole revision chain
objective               the bounded result Claude must produce
base_ref                commit from which the worker worktree is created
policy_files            architect-side policy snapshots Claude must read
sealed_facts            facts Claude may use but may not reopen
scope                    allowed and forbidden paths
                         and whether mathematical naming audit is mandatory
requirements             individually identified acceptance obligations
required_checks          exact worker-side checks and expected evidence
stop_conditions          conditions that require BLOCKED rather than guessing
deliverables             expected code, tests, documents, or analysis
execution                mode, model, effort, timeout, budget, and network law
```

The launcher copies the task and policy files into the run record and records
their hashes before Claude starts. They are never overwritten. A changed
objective requires a new task ID. A correction to execution uses a review and a
new turn under the same task ID.

The architect must translate a large user objective into requirements that are
independently checkable. A packet is not ready while it contains words such as
"correct", "clean", or "appropriate" without a law or observable consequence.

## Isolation

Each task runs in a detached Git worktree at the recorded `base_head`. The main
checkout is never Claude's working directory. The launcher refuses delegation
when an existing main-checkout change overlaps the packet's allowed paths.

The worktree protects the user's active checkout from ordinary editing errors.
It is not a security sandbox. Therefore the worker contract also forbids access
outside the worktree, and network tools are removed unless the packet explicitly
permits them.

Runtime records live under the repository's Git common directory:

```text
.git/architect-worker/runs/<task_id>/
```

Worker worktrees live in a repository-specific sibling directory. Neither is
part of the source diff.

## Communication record

Communication is file-based and append-only:

```text
worker-input/task.json            immutable original delegation
worker-input/policy/              immutable policy snapshots
worker-input/reviews/             copies of reviews exposed to Claude
events.jsonl                      state transitions
turn-001/request.md               exact prompt sent to Claude
turn-001/envelope.json            complete Claude CLI response
turn-001/report.json              schema-validated worker report
turn-001/candidate.json           actual changed-file hashes
turn-001/candidate.tar.gz         snapshot of changed files
turn-001/diff.patch               tracked-file patch
reviews/review-001.json           architect findings and decision
turn-002/...                      revision response; turn 1 is preserved
```

No RED result, failed attempt, mistaken interpretation, or review finding is
rewritten. Its resolution is appended in a later turn.

The launcher grants Claude explicit additional access only to `worker-input/`;
it does not supply the run record as worker context. Earlier envelopes, reports,
events, patches, and candidate archives remain architect-side evidence. Input
hashes are rechecked after every turn, and the canonical architect review
remains in `reviews/`, outside the worker-input directory. This separation is an
evidence boundary, not a hostile-process security sandbox.

Claude's report must satisfy `worker-report.schema.json`. The launcher also
checks facts Claude cannot certify about itself:

- the task ID and base commit agree;
- the worktree HEAD did not move;
- actual changed paths equal reported changed paths;
- every changed path is in scope and no forbidden path changed;
- every requirement and required check appears in the report;
- a candidate has every requirement satisfied and every required check passed;
- a candidate contains no unresolved decision question;
- read-only tasks made no changes;
- a code candidate contains an explicit naming audit.

A schema-valid report that fails these checks becomes `PROTOCOL_VIOLATION`, not
`CANDIDATE`.

## Architect review

Codex reviews in this order:

1. Restate the user's objective and compare it with the immutable packet.
2. Inspect the complete diff, including untracked and deleted files.
3. Check architecture, ownership, level direction, and sealed constraints.
4. Run the deterministic naming audit required by `AGENTS.md`.
5. Look for weakened tests, removed RED evidence, hidden fallbacks, and claims
   made true only by changing the experiment.
6. Reproduce required checks independently; do not reuse the worker's verdict.
7. Run broader regressions in proportion to the change's blast radius.
8. Record findings and evidence with `review.template.json`.

`VERIFIED` means the candidate satisfies the packet in its isolated worktree.
It does not mean the main checkout contains the change.

`INTEGRATED` means Codex transferred the reviewed change to the architect
checkout and inspected the resulting diff again.

`COMPLETE` means the integrated result passed the final relevant checks and no
packet requirement remains open. For a sealed tower, only Codex may recommend
sealing; the user retains any explicitly reserved sealing authority.

## Revision loop

A revision request identifies each defect separately:

```text
finding_id
severity
failed requirement or law
concrete evidence
required correction
checks that must be rerun
```

The original Claude session is resumed in the same worktree. Claude receives
the immutable task again plus the new review. It must address the findings or
return `BLOCKED`; it must not broaden scope to make the findings disappear.

Repeated failure is not erased. Codex decides whether another bounded revision
is useful, whether to decompose the task, or whether the obstacle requires user
judgment.

## Escalation to the user

Codex asks the user only when at least one of these is true:

- the objective admits materially different mathematical interpretations;
- a new canonical taxonomy word is required;
- a sealed architectural statement must be changed;
- scope or authority must expand;
- evidence contradicts an accepted requirement;
- external access or a destructive operation is required;
- the worker remains blocked after bounded, evidence-based revisions.

The escalation contains the exact decision, alternatives, evidence, and the
consequence of each choice. It does not forward an unfiltered worker transcript.

## Commands

From the repository root:

```bash
python3 .agents/architect-worker/pipeline.py doctor
python3 .agents/architect-worker/pipeline.py prepare path/to/task.json
python3 .agents/architect-worker/pipeline.py dispatch <task_id>
python3 .agents/architect-worker/pipeline.py status <task_id>
python3 .agents/architect-worker/pipeline.py review <task_id> path/to/review.json
python3 .agents/architect-worker/pipeline.py revise <task_id>
python3 .agents/architect-worker/pipeline.py retry <task_id>
python3 .agents/architect-worker/pipeline.py archive <task_id>
python3 .agents/architect-worker/pipeline.py discard <task_id>
```

`revise` resumes the existing Claude session. `retry` starts a new session after
an architect-recorded `REVISION_REQUESTED`; use it when transport failed before
Claude created a resumable session. The failed turn remains in the run record.

`archive` is allowed only after `COMPLETE` or `REJECTED`, and only for a clean
worker worktree. It removes the disposable worktree but preserves the complete
run record.

`discard` is allowed only after `REJECTED`. It verifies that the latest
candidate snapshot and digest exist, then removes the dirty disposable worktree.
The rejected patch remains recoverable from the append-only run record.

## Handoff to the user

Codex reports the work at architect altitude:

```text
current state
what Claude changed or discovered
architect findings
checks independently verified
remaining uncertainty or blocked decisions
next bounded step toward completion
```

The worker's prose is evidence, not the final answer. Codex is responsible for
every claim delivered to the user.
