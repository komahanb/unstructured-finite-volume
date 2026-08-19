# GTI roadmap: assumptions classified by core mathematical object

`doc/GTI.md` is the map of what exists. This document rules on what may
change: every assumption the current tower makes is attached to the
mathematical object it constrains and given one of two verdicts.

- **law** — a statement the tower is built on. Every suite may rely on
  it; drivers may assume it without checking twice; it changes only by
  deliberate amendment. Some laws are marked *amendable*: closed
  vocabularies designed to grow by amendment to one owning module,
  never by an open extension point.
- **scaffolding** — an assumption kept so that proofs stay small and
  suites stay exact. It may be retired without moving any public
  contract. Every scaffolding verdict names its retirement path and the
  invariant that survives retirement — retiring scaffolding that quietly
  bends a law is the failure mode this classification exists to prevent.

The sections follow the five core objects: the typed point a form acts
on, the derivative path algebra over it, the relation as a solvable
unit, the decision layer above the mechanisms, and the execution of the
mathematics beneath them.

# 1. Typed evaluation components

The core object is the evaluation point — the typed (U, xi, t) a
differentiable form acts on — together with its carriers: value
buffers, state/design bundles, and the form contract itself.

| Assumption | Where it lives | Verdict |
|---|---|---|
| One abstract form contract: name, input/output signature, max_degree, value, partial_action | `gti_form_interface` | law |
| Partial actions are multilinear contractions; no derivative tensor is ever stored | form contract, all drivers | law |
| Every partial_action self-validates through require_supported first | form contract | law |
| The argument vocabulary is closed: STATE, DESIGN, TIME, GEOM | `gti_form_interface` | law (amendable) |
| The state vocabulary is closed at differential orders 0..2: q, qdot, qddot | `gti_state_bundle` | law (amendable) |
| Output shape is [nentries, ncomp], declared by the form, re-checked on every call | `gti_form_evaluator` | law |
| Seats are not occupants: bundles split shape from occupancy, and absence is a lawful zero | bundles, chain assembly | law |
| Time is one scalar per evaluation point | `gti_evaluation_point` | law |
| Design is a flat list of field slots; design enters forms only through the bundle | `gti_design_bundle` | law |
| Value payloads are REAL64 vectors (GTI_VALUE_REAL is the only value kind) | `gti_value_buffer` | scaffolding |
| GEOM is a reserved seat with no occupant: no geometry bundle rides the point yet | `gti_evaluation_point` | scaffolding |
| The point carries identity seats window_id / step_id / stage_id that no law consumes yet | `gti_evaluation_point` | scaffolding |

Rulings on the scaffolding rows:

- **REAL64-only payloads.** Precision and payload kind (mixed precision,
  complex arithmetic) are representation, not contract. Retirement adds
  value kinds to the buffer module; the shape law and the
  seats-vs-occupants split survive unchanged, and no form signature
  moves.
- **The empty GEOM seat.** This is a reserved seat, not a missing law:
  the closed vocabulary already names geometry so that mesh/FV coupling
  fills an existing seat instead of amending the vocabulary under
  pressure. Retirement is a geometry bundle on the evaluation point plus
  forms that read it — no assembler or driver change.
- **Inert identity seats.** window_id / step_id / stage_id become
  load-bearing when windowing and checkpointing arrive (section 4).
  Until a law consumes them they are labels; nothing may *depend* on
  them before that law exists.

# 2. Degree-r path algebra

The core object is the path x(eps) through a form's arguments and the
total-derivative algebra over it: J_u q_u^(s) = -B^(s) with one J_u for
all s.

| Assumption | Where it lives | Verdict |
|---|---|---|
| The same J_u serves every degree; only the chain-rule right-hand side changes | degree-r driver, Newton/tangent/adjoint seats | law |
| B^(s) is defined by the suppression rule: unknown top seat zero, solved lower seats included, history included, design path as given | degree-r driver | law |
| One assembler owns all composition combinatorics; no driver re-derives a term | `gti_chain_rule_assembly` | law |
| A channel is one argument's path; seats above the requested degree are ignored; absence is zero; ordered tuples enumerate mixed terms | `gti_chain_rule_assembly` | law |
| User forms are smooth to their declared max_degree with commuting mixed partials — the pattern multiplicities count on it | form contract obligation | law |
| Directional derivatives live in caller-owned arrays indexed (degree, vertex); empty means zero | degree-r driver | law |
| The pattern table is enumerated by hand through degree 4, and r <= 4 | `gti_chain_rule_assembly`, degree-r options | scaffolding |
| The traversal's design path is affine: xi(eps) = xi + eps eta, higher seats absent | degree-r driver | scaffolding |
| Time is not a path variable: dt/deps = 0 in every traversal | all traversal drivers | scaffolding |

Rulings on the scaffolding rows:

- **The degree-4 table.** The bound is a table size, not an
  architecture. The degree-r driver already loops generically; the
  assembler's public verb takes the degree as data. Retirement replaces
  the hand table with a pattern generator over integer compositions —
  the channel laws, the public call shape, and every driver survive
  verbatim, and the hand table remains as the generator's degree <= 4
  oracle in the suite.
- **The affine design path.** The assembler already accepts occupied
  higher design seats; only the driver chooses to leave them absent.
  Retirement lets the caller supply a design path object with seats
  1..r. The invariant that survives: design derivative data enters the
  algebra only through channel seats, never through a side channel.
- **Frozen time.** GTI_ARG_TIME is a reserved lawful channel. The moment
  step sizes depend on design — h(xi), which any real adaptive
  controller induces — dt/deps is nonzero and the time channel must be
  occupied. The invariant that survives: time enters the algebra only as
  a channel, never as hidden state inside a driver.

# 3. Block implicit relations

The core object is the relation as a solvable unit. Today a relation is
an ordered sample tuple, a motif, one unknown sample position, and an
evaluation time; DIRK enters stage-by-stage, one relation per stage.

| Assumption | Where it lives | Verdict |
|---|---|---|
| A relation solves exactly one unknown sample | `gti_time_graph` | law |
| J_u is square: residual size equals unknown size | Newton and degree-r seats | law |
| Motif coefficient rows weight samples by scalars, one real per sample per row | `gti_time_local_scheme` | law (amendable) |
| Relations reference any arity of history samples | `gti_time_graph` | law |
| In stored order, every relation's non-unknowns are already solved (the history law) | forward driver | law |

This section's classification IS a ruling, and it is the decision this
document exists to record:

> **Block implicit schemes enter as product members, not as
> multi-unknown relations.** A coupled stage bundle — fully implicit
> Runge-Kutta, collocation — is ONE state instance whose q stacks the
> coupled stages. The residual form couples the stages internally; the
> relation still names one unknown sample; the member set absorbs the
> block structure.

Consequences, which is why the single-unknown relation is a law and not
scaffolding:

- Every driver works unchanged: J_u is already n x n for vector q, so a
  stacked stage vector is just a larger n. Newton, tangent, adjoint,
  reverse, and degree-r never learn the word "stage".
- Scalar coefficient rows suffice: a scalar weight on a stacked vector
  is exactly the Kronecker action w x I that block stencils need.
  Stage-internal coupling belongs to the form, not to the rows. (Should
  a genuinely matrix-weighted stencil ever be required, that is an
  amendment to the motif object — hence *amendable* — not a new relation
  species.)
- Simultaneous global-in-time solves (waveform iteration) are not block
  relations either: they are a policy-plus-backend composition over the
  same single-unknown relations, and they stay out of the relation
  object.

# 4. Policy/schedule layer

The core object is the decision: what to march next, what to accept,
what step or order to try. The standing law of the layer is the
mechanism/policy split — mechanisms below never decide.

| Assumption | Where it lives | Verdict |
|---|---|---|
| Drivers consume stored order; nothing sorts, schedules, or re-plans inside a driver | all traversal drivers | law |
| Acceptance is an external logical; the growth driver carries no estimator and no step law | adaptive growth driver | law |
| Non-convergence is a lawful reported answer; rejected or failed growth leaves no trace | forward and growth drivers | law |
| Functional terms anchor to vertices: vertex index plus evaluation time | functional seed driver | law |
| No automatic controller exists: no estimator, no step-size/order selection, no retry loop | absent | scaffolding (by absence) |
| BDF rows exist for uniform steps only | `gti_time_motif_builder` | scaffolding |
| No checkpointing/restart schedule exists | absent | scaffolding (by absence) |

Rulings on the scaffolding rows:

- **The missing controller.** It arrives as a client of
  `try_candidate`: estimator, step-size/order policy, retry loop,
  policy object — all above the drivers. The acceptance criterion for
  that future work is sharp: adding the controller must not add one
  line to any existing driver. The graph law and the driver verbs are
  already its complete interface.
- **Uniform BDF rows.** Variable-step BDF is a new builder minting
  different coefficient rows from the same types. The law that survives:
  named schemes exist only as row builders; no scheme name reaches a
  driver.
- **The missing schedule.** Windowing, checkpointing, and restart are
  schedule policies over the same graph; the inert identity seats of
  section 1 are their reserved hooks. Policies produce order and
  placement; mechanisms consume them.

# 5. Execution-plan/backend layer

The core object is the execution of the mathematics: storage, linear
algebra, and the scheduling of work. Everything here may change, and
the classification says so — this layer is where nearly all remaining
scaffolding lives. The one law is the boundary itself.

| Assumption | Where it lives | Verdict |
|---|---|---|
| Solves sit behind driver verbs; the graph never solves; the action interface is the contract | all drivers | law |
| Local solves are small dense private Gaussian eliminations, deliberately duplicated per driver | Newton/tangent/adjoint/degree-2/degree-r seats | scaffolding |
| J_u is materialized column-by-column from basis directions through the Jacobian action | Newton and directional seats | scaffolding |
| build_samples deep-copies relation-local state | `gti_time_graph` | scaffolding |
| Every Jacobian action call re-injects the trial q and rebuilds the evaluation point | Newton driver | scaffolding |
| The time graph stores plain arrays of vertices and relations holding value buffers | `gti_time_graph` | scaffolding as storage; law as control plane |
| Traversals are serial: one relation at a time, one degree at a time | all drivers | scaffolding |
| B^(s) is assembled fresh per degree; shared channel paths are recomputed | degree-r driver | scaffolding |

Rulings:

- **Dense private solves.** A production backend — sparse, iterative,
  matrix-free — replaces them behind the same verbs without changing
  the graph law. The duplicates are deliberate: no premature solver
  abstraction was allowed to harden while the laws were still being
  proven.
- **Column-built J_u.** The law is the action, not the matrix. A
  matrix-free Krylov backend consumes the identical `jacobian_action`
  and never forms J_u at all; the dense column harvest is the proof
  world's smallest correct implementation.
- **Deep copies and point rebuilds.** Copies enforce graph immutability
  in the proof world; an execution plan may replace them with a
  borrow/view discipline and hoisted point construction. The invariants
  that survive: traversals must not mutate the graph, and trial
  injection preserves field identity.
- **Array storage and serial order.** The time graph is a
  traversal/control plane, not hot-path storage; a backend owns any
  compiled or blocked storage tier. Batching independent relations and
  fusing degree loops are execution plans over the same stored order —
  independence is discovered by a policy, executed by a backend, and
  never inferred inside a driver.
- **Fresh assembly per degree.** The suppression rule fixes what the
  values are; the plan chooses when they are computed. Memoizing shared
  channel paths is plan-level caching with no semantic content.

# Retirement order

Laws never retire. The scaffolding retires in dependency order, each
step leaving every law intact:

1. **Pattern generator** for arbitrary r (section 2): local to the
   assembler, unlocks r > 4, keeps the hand table as its oracle.
2. **Caller-supplied design paths and time-channel occupancy**
   (section 2): curved paths and dt/deps, prerequisite for controller
   sensitivities under h(xi).
3. **Variable-step motif builders** (section 4): coefficient rows any
   real controller needs.
4. **The controller itself** (section 4): estimator, step/order policy,
   retry — a client of try_candidate, zero driver edits.
5. **Backend seats** (section 5): shared then replaceable local solves,
   matrix-free consumption of the action, borrow discipline, hoisted
   points, compiled storage, batched traversal.
6. **GEOM occupancy** (section 1): mesh/FV coupling fills the last
   reserved seat.

The order is chosen so that every step is testable against the step
before it, and so that the two objects users touch — the form contract
and the graph law — are never the thing that moves.
