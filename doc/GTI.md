# Graph Time Integrators

GTI is this repository's time-integration architecture: heterogeneous,
adaptive time marching expressed as traversals of a time graph, where
derivatives of every supported order arrive through the same machinery
that solves the primal problem. One differentiable-form contract sits at
the bottom, a time graph sits in the middle, and stateless traversal
drivers sit on top. The user supplies calculus; the graph layer owns
composition, ordering, and bookkeeping.

The core picture:

```
G_time = (S_time, R_time)

S_time : the member set of time/stage state instances
R_time : scheme relations between state instances

each relation owns:
    motif rows            one row per state component it builds
    coefficient weights   each coefficient row weights the relation's samples
    unknown sample index  which sample the relation solves for
    evaluation time       where the residual form is evaluated

fields/buffers carry:
    q                        primal state values on the vertices
    qdot / qddot             components built from motif rows, never stored
    adjoint seeds            in caller-owned seed arrays, outside the graph
    directional derivatives  q^(s) in caller-owned derivative arrays, outside the graph
```

Three facts are load-bearing and hold everywhere:

- Graph identity/topology is not the state values. What a vertex *is*
  never changes because of what it currently *holds*.
- Seeds and derivative arrays live outside the graph. Traversals read
  the graph and write caller-owned arrays.
- Graph q values are primal data. Only the forward primal traversal and
  transactional growth ever write them.

A relation is a scheme instance: an ordered tuple of vertices, a motif
(its coefficient rows), the unknown sample position, and an evaluation
time. A motif knows no scheme by name — BDF, DIRK, and ABM exist only as
builders that mint coefficient rows.

# User story

The public promise:

1. Implement the residual form R.
2. Implement the functional form F.
3. Build or choose a time graph.
4. Call the GTI driver.

```fortran
type(my_residual_form)   :: R
type(my_functional_form) :: F
type(gti_time_graph)     :: G
type(gti_time_functional):: J

call solve_gradient % solve(R, F, G, J, design, eta, options, result)

value      = result % value
dJ_dxi_eta = result % total_design_gradient_action
```

R and F share the same differentiable-form interface:

- `value` — evaluate the form at an evaluation point;
- `partial_action` — apply one exact mixed partial to given directions,
  as a multilinear contraction (no derivative tensor is ever built);
- `max_degree` — the highest partial order the form promises;
- `input_signature` / `output_signature` — which argument kinds the form
  reads, and the shape its value fills.

The two forms differ only in role: **R enforces equations** (each
relation demands R = 0 at its unknown), **F seeds and accumulates
functionals** (its terms produce a value, q-seeds, and an explicit
design action). Their derivative interface is the same, which is why
one chain-rule seat and one set of traversal drivers serve both.

# Module map

One file per seat under `src/`. A file `gti_X.f90` hosts the Fortran
module `gti_Xs` — the language denies a type its host module's name, so
modules speak in the plural while every public type keeps the singular.
The table names the seats by their files. Each seat has a law-and-refusal
suite under `test/gti-*`.

| Layer | Module | Public role | Does not own |
|---|---|---|---|
| form contract | `gti_form_interface` | differentiable form interface | no traversal, no graph |
| form evaluator | `gti_form_evaluator` | output-shape checked value/partial_action calls | no chain-rule combinatorics |
| chain-rule assembly | `gti_chain_rule_assembly` | degree 0..4 local higher-order chain-rule RHS assembly | no time graph, no Newton, no traversal |
| time-local scheme | `gti_time_local_scheme` | motif rows turn samples into q/qdot/qddot evaluation points | no named schemes, no graph traversal |
| named motif builders | `gti_time_motif_builder` | BDF/DIRK/ABM coefficient-row construction | no solve |
| unknown problem | `gti_time_local_unknown_problem` | inject trial unknown q into local samples | no Newton policy |
| local Newton | `gti_time_local_newton_driver` | local residual solve and J action | no time graph ownership |
| local tangent | `gti_time_local_tangent_driver` | local design tangent solve | no time traversal |
| local adjoint | `gti_time_local_adjoint_driver` | local adjoint/design action | no time traversal |
| time graph | `gti_time_graph` | G_time representation | no solving, no seeds, no derivatives |
| forward driver | `gti_time_forward_driver` | stored-order primal march | no reverse, no functional |
| reverse seed driver | `gti_time_reverse_driver` | reverse traversal from caller-owned vertex seeds | no functional seeding |
| functional seed driver | `gti_time_functional_seed_driver` | F_time -> value + q-seeds + explicit F_xi[eta] | no reverse call |
| solve-gradient driver | `gti_time_solve_gradient_driver` | forward + seed + reverse composition | no local calculus |
| adaptive growth driver | `gti_time_adaptive_growth_driver` | append/solve/accept-or-rollback graph growth | no estimator, no step-size policy |
| degree-2 directional driver | `gti_time_degree2_directional_driver` | compatibility/proof driver for degree 1/2 | superseded architecturally by degree-r |
| degree-r directional driver | `gti_time_degree_r_directional_driver` | q^(1)..q^(r), r = 1..4, over solved G_time | no reverse, no functional, no finite differences |

The import lists are the architecture: each seat may name only the seats
beneath it, and every boundary in the "does not own" column is enforced
by what the module is *allowed to import*, checked by the suites.

# Traversals

Four traversal modes exist. All of them run in the graph's stored
relation order (or its exact reverse); none of them sorts, schedules,
or re-plans.

**1. Forward primal traversal** — the forward traversal that solves the
graph:

```
for relation r in stored order:
    build samples
    solve R_r = 0 for unknown q
    write q back to unknown vertex
```

Non-unknown vertices of each relation must already carry solutions (the
history law). Non-convergence is a lawful answer, reported, never an
error stop; nothing is written back for a failed relation.

**2. Reverse seed traversal** — the adjoint sweep:

```
for relation r in reverse stored order:
    solve J_u^T lambda = seed_u
    accumulate -lambda^T R_xi[eta]
    propagate seed_h += -R_h^T lambda
```

Seeds enter from, and accumulate into, a caller-owned seed array — one
value buffer per vertex, an empty buffer meaning zero. The graph is
never written.

**3. Functional seeding** — how a functional becomes seeds:

```
for functional term k:
    value += F_k
    vertex_seed(v_k) += F_{k,q}^T
    explicit += F_{k,xi}[eta]
```

A time functional is a sum of terms, each naming one vertex and one
evaluation time, all evaluated through one functional form F. Duplicate
terms on the same vertex add.

**4. Degree-r directional traversal** — forward derivatives of every
order:

```
for relation r in stored order:
    build J_u once
    for s = 1,...,r:
        assemble B^(s)
        solve J_u q_u^(s) = -B^(s)
```

Derivatives land in a caller-owned derivative array indexed
(degree, vertex); later relations read earlier vertices' derivatives as
history.

The invariant all four share, and the architecture's central claim:
**the same J_u serves every degree — only B^(s) changes.** Newton
eliminates it once at degree 0, the tangent and adjoint seats reuse its
action at degree 1, and the directional traversal eliminates it r times.
Higher-degree terms are assembled locally by `gti_chain_rule_assemblies`;
no driver re-derives them.

The solve-gradient driver is not a fifth traversal: it is the pure
composition forward → seed → reverse → total = explicit + residual,
importing only the three time-level drivers it composes.

# Ownership laws

- `fractal_graph` identity is not modified by GTI.
- G_time stores primal time vertices and relations.
- q values live on time graph vertices.
- Reverse seeds live in caller-owned value-buffer arrays.
- Directional derivatives live in caller-owned value-buffer arrays.
- Motifs own coefficient rows, not solvers.
- Drivers are stateless verb sets.
- Forms own local value and partial_action logic.
- Chain-rule assembly owns combinatorics.
- Time drivers own traversal order.
- No scheme-specific sensitivity code is required.

And the negative laws:

- No parent pointer in graph.
- No hidden semantic map in graph.
- No adjoints stored in the graph.
- No derivative tensors materialized.
- No finite-difference sensitivity path in GTI production drivers.

Every law is loud: a violated precondition dies with
`error stop 'module: law statement'`, and each suite's refusal program
proves each message fires at its owning seat.

# Heterogeneous chaining

Proven by `test/gti-heterogeneous-chain` (test-only; no production code
was added to enable it): one G_time carrying a BDF1 relation, a DIRK
relation, and an ABM2 relation, marched and differentiated by the same
solve-gradient driver that serves homogeneous chains.

```
v1 --BDF1--> v2 --DIRK--> v3 --ABM2--> v4
```

- Heterogeneous chaining is graph composition. It is not scheme-pair
  glue.
- Each relation owns its motif rows.
- Vertices are shared across relation boundaries — sharing a vertex *is*
  the coupling.
- Drivers are scheme-blind: no driver names BDF, DIRK, or ABM.

The suite's crowning check: the reverse-propagated seed at v1 equals
the hand chain rule d(q4)/d(q1) = (2/3)(1/2)(2/3) = 2/9 straight across
both scheme boundaries.

# Adaptive graph growth

Current adaptive support is transactional graph growth:

```
candidate vertex + relation
  -> append to G_time
  -> solve appended relation
  -> accept keeps graph
  -> reject rolls back
```

A candidate names its own appended vertex by the sentinel index −1; an
accepted append persists, a rejected or non-converged append leaves no
trace. The accept decision is an external logical — the driver carries
no policy.

**This is not yet an automatic adaptive time-step controller.**

Current: transactional graph growth.

Future: estimator, step-size/order controller, retry loop, policy
object — all of them clients of `try_candidate`, none of them changes
to the graph law.

# Higher-degree directional traversal

`max_degree` r is currently 1..4. q^(s) is computed by stored forward
relation traversal: for each relation, J_u is built once from the primal
unknown Jacobian action, and q_u^(s) solves against a chain-rule
right-hand side —

```
J_u q_u^(s) = -B^(s),    s = 1,...,r

B^(s) = the total degree-s chain-rule residual
        excluding the unknown transport term R_u[q_u^(s)]
```

The suppression rule that defines B^(s):

```
unknown q_u^(s) is zero while assembling B^(s)
unknown q_u^(k), k < s, is included
history q_h^(k), k <= s, is included
xi^(1) = eta
xi^(k) = 0 for k >= 2
```

The last two lines are the affine design path xi + eps eta: the design
channel carries eta in its first seat and every higher seat is absent,
absence meaning zero.

B^(s) is produced by the higher-order chain-rule assembly seat
(internally, Faà di Bruno combinatorics over ordered channel tuples),
which owns every pattern and multiplicity through degree 4. A residual
whose own third and fourth partials vanish still has nonzero q^(3) and
q^(4): the right-hand side transports products of lower-degree
curvature, and the degree-r suite verifies this against closed forms at
every degree.

# Current status

Implemented and tested:

- local differentiable forms
- local chain-rule assembly through degree 4
- BDF/DIRK/ABM motifs
- local Newton/tangent/adjoint seats
- time graph representation
- forward primal traversal
- reverse seed traversal
- functional-to-seed construction
- solve-gradient composition
- heterogeneous chain proof
- adaptive graph growth
- degree-r directional traversal, r = 1..4

Not yet implemented:

- automatic adaptive time-step controller
- production public facade that hides all intermediate seats
- checkpointing/restart policy for long time graphs
- sparse/global linear algebra backend for large systems
- mesh/FV coupling
- documentation examples beyond toy forms

# API stability note

This is an internal architecture map. Some module names are still
implementation seats: they exist so each law has one owner, not because
a user should have to name them. The eventual public facade should be
smaller than this module list — on the order of

```fortran
call gti_solve(problem, graph, functional, design, options, result)
```

That facade does not exist yet and is deliberately not added here; the
solve-gradient driver is its closest current ancestor.

# Boundary audit

- **Kernel**: `src/fractal_graph.f90` must remain untouched by GTI work.
  Eighteen implementation phases have held this line.
- **Mesh/FV**: GTI currently has no mesh/FV dependency. No GTI module
  imports mesh, flux, or assembler seats.
- **Solver**: GTI local drivers currently use small dense private solves
  for tests/proofs. A production backend can replace these later without
  changing the graph law — the solve sits behind a driver verb, not in
  the graph.
- **Chain-rule assembly**: local form calculus only; no time graph
  dependency. The time drivers bring the graph to the assembler, never
  the reverse.
- **Time graph**: traversal/control plane; not hot-path storage for
  production large systems.
- **Degree-2 driver**: retained for compatibility/proof history. The
  degree-r driver is the architectural successor, and its suite proves
  the r = 2 case reproduces the degree-2 driver exactly.
