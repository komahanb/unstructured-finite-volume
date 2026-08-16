# Fractal graph kernel — feasibility spike report

Branch `tower-graph-as-sets-relations`. Prototype: `test/fractal-kernel/`
(kernel 227 code lines, views 291, suite 26 laws + 18 refusals, all
green, zero warnings, gfortran-15 `-std=f2023`). Production untouched.

## The verdict, first

**Feasible — as an ontology. Conditional — as an implementation.**

One primitive — `graph` with `branch(2)`, each branch NULL, UNKNOWN or
KNOWN — carries all four required inhabitations without a single
domain-specific change: the bare atom, the calculator `(2+3)×4 = 20`,
a finite-volume strip with boundary faces as NULL and a conserved
field, and the computational pair `G = (Q, R)` realized by solving.
The acceptance criterion held everywhere it was pushed: relational
graph and computational graph are one node; adjacency, incidence and
the ordinary neighbourhood are three readings of one pair list; the
calculator expression and the residual right-hand side are one
evaluator; `relational_graph` and `computational_graph` both die as
roots and survive as views.

The condition: the compression is **ontological, not operational**.
An adversarial performance review established that the semantic road
(handle-by-handle traversal) can never carry the hot path — at best
~24 ns/hop against production's 9 ns/query, degrading toward DRAM
latency at mesh scale, and the spike's naive walkers are O(E²). What
passes §66 is a storage tier of compiled contiguous blocks, offset
tables and borrowed slices — which is `csr_relation`'s arrays
readmitted behind the kernel. That is not a defeat; it is the
architecture: **one meaning above, the same fast arrays below.** The
kernel is the specification-and-construction layer; CSR becomes a
storage law it hides, exactly as the brief permitted ("separate
ontology from efficient representation if necessary"). But nobody
should migrate expecting the primitive itself to be fast.

---

## A. Proposed declarations

Semantically, as governed:

```fortran
type :: graph
   type(graph_branch) :: branch(2)
end type
```

Physically, as proven — because the conceptual
`class(graph), allocatable :: known` **fails three of the required
laws**: allocatables own, so sharing becomes deep copy (a twin, not
the same B); ownership is a tree, so cycles are unrepresentable; and
Fortran finalizes recursive allocatables by recursion, so deep
structures die of stack. The realization representation is therefore
a handle into an arena:

```fortran
type :: graph_arena                              ! the universe
   type(token), private :: identity              ! signs the roll once
   integer,     private :: n
   integer , allocatable, private :: status(:,:) ! (2, capacity)
   integer , allocatable, private :: target(:,:) ! (2, capacity)
   logical , allocatable, private :: carries(:)  ! value atoms
   real(dp), allocatable, private :: payload(:)  ! compressed numerals
end type

type :: graph                                    ! the handle
   class(graph_arena), pointer, private :: universe
   type(token),                 private :: stamp ! the signing that minted it
   integer,                     private :: seat_at
end type
```

`type(graph_branch)` as a public bundle was built, then **evicted by
the repository's own admission laws**: it is `status(side)` plus a
guarded `known(side)` — a composition (GENERATION), with zero
exercised citizens (INHABITATION), whose public components admitted
states the kernel's laws forbid. The ontology's `branch(2)` is the
pair of observations, not a type the kernel must sell.

Kernel public surface (the generators, and nothing else):

```text
constants     BRANCH_NULL  BRANCH_UNKNOWN  BRANCH_KNOWN
arena         graph_arena()  node(l, r)  value_atom(v)
              citizen(k)  population()
handle        status(side)  known(side)  realize(side, g)
              same_as(g)  carries_value()  value()  seat()
              universe_size()          [population re-exported across
                                        the capability wall - a reader
                                        never gains mint authority]
specs         branch_spec  null_branch()  unknown_branch()
              known_branch(g)          [ABSORPTION-clean: the 3-point
                                        status choice as a value; the
                                        overload road would need
                                        marker types]
```

## B. Identity, NULL, UNKNOWN, KNOWN

The three statuses are **values of one axis**, held in arrays — never
encoded as unallocated storage, never as subclasses. All nine 3×3
branch-state combinations are minted by the suite from one type.

The existing `graph_identity` token mechanism is **sufficient, used
at a different grain**: the universe signs the roll once; a node's
identity is `(universe token, seat)`, and every handle carries the
stamp of the signing that minted it. Two `(NULL, NULL)` atoms are two
atoms because two seats — the suite proves identity is the only thing
that tells them apart. No leaf ontology exists; an atom IS a graph.
The stamp is load-bearing, not ceremonial: a universe assigned over
is a different universe, and its orphaned handles are refused
("this universe is not the one that minted you"), not misread — a
hole the red team found and the stamp closed. An unhatched handle
equals nothing, itself included — the undeclared token's own law.

The two absences are distinct at the primitive: **UNKNOWN is ⊥**
(COMPUTATIONAL-GRAPH.md's epistemic bottom, not yet realized) and
**NULL is definite absence** (the boundary face's missing far side,
which needs no imaginary member). The naming pass legislated ⊥ ≠ ∅
one level up; the fractal kernel promotes it into the ground type.

## C. Sharing and cycles

Sharing: a realization is a seat; two branches holding one seat hold
**the same** realization — `A = (B, B)` answers `same_as` across its
branches. No duplication is representable, because branches hold
coordinates, not copies.

Cycles: the knot is tied by the kernel's **one transition** —
`realize`, UNKNOWN → KNOWN, exactly once. Mint A ignorant, mint B
knowing A, realize A's branch at B: two hops return home. The
representation is a graph, never a mere tree; particular views refuse
cycles by law where their questions demand it (`seq_len` walks under
a population budget and refuses "a spine that never ends"; the
evaluator carries the same budget after the red team hung it with a
cyclic expression).

The growth law: **knowledge only grows.** NULL never realizes
("absence is not ignorance"), KNOWN never re-realizes ("knowledge
grows once"), and arena growth reallocates arrays without
invalidating anyone, because a seat is an index, never an address —
proven by minting three hundred citizens past the initial capacity
and asking the first handle afterward.

From this the red team distilled the soundness theorem the views now
carry as the FINALITY law: *every status a lawful walker's answer
depends on is final at the moment of answering* — NULL and KNOWN
never change, and no walker answers through UNKNOWN; it refuses. So
no answer is ever invalidated by a later realize. (The first draft's
relation walkers silently read (KNOWN, UNKNOWN) pairs as boundary
occurrences — the only non-monotone observation in the system; they
now refuse: "an unlearned far side has no answer yet".)

## D. Recursive encoding of n-ary objects

One prescribed law, stated once in `fractal_views.f90`, agreed by
every reader:

```text
SEQUENCE  a cell's branch(1) is the element, KNOWN always;
          branch(2) is the continuation - KNOWN the next cell,
          NULL the end. The EMPTY sequence is a NULL branch at
          whoever holds it: there is no empty cell, so no second
          spelling of emptiness exists. An UNKNOWN tail is a
          sequence still being learned - no extent yet; walkers
          refuse it.

(a,b,c)   =  [a | *]->[b | *]->[c | NULL]     one spelling, one law
```

A tuple of arity k is a k-sequence; a relation is a sequence of
tuple nodes; `(a)` is a cell, never the bare atom `a`. Uniqueness is
**per reading**: within the sequence reading no second spelling
exists (the empty case is closed by law — there is no cell to spell
it with). Across readings, coincidence is the point, not a defect: a
singleton sequence cell and a boundary pair are the same shape
`(KNOWN, NULL)`, and only the reading decides — which is precisely
"two concepts that differ only by interpretation collapse onto the
same representation."

Two law-level findings from the adversarial pass, both closed:
the evaluator now reads the lawful empty spelling (a NULL argument
branch answers the identities: a sum over nobody is 0, a product
over nobody is 1); and the occurrence question is legislated — the
pair list is an **occurrence (edge-list) reading**: counts count
occurrences, `rel_has` quotients them; both lawful, neither the
other.

## E. Atom / value strategy

A value is a **compressed atom**: an atom may carry one number, which
stands for the numeral the structure could spell out and never
should. Not a second ontology — a storage annotation. Only atoms
carry values; meaning still arrives by interpretation: the
calculator's operations are bare atoms, and arithmetic attaches to
their **identity** through the reading's law table (`op same_as plus`)
— the constitution doctrine, fractally. A compressed numeral is never
a member of an environment (refused), so no binding can mean two
things to two readers.

Fields: value atoms minted consecutively form a block; a field is
the pair (domain sequence, value sequence), answered two ways — down
the semantic road cell by cell, and through the **contiguity
promise** (`block_value`: the k-th value in one step, seat
arithmetic). The suite proves both roads answer alike and that
interior face differences telescope to `q(4) − q(1)`. Production
generalizes the single `payload(:)` into typed pools (integer, real,
complex; SoA-segmented so structure nodes cost 8 bytes, not 28 — see
I below), with `seat order is mint order` stated as kernel law,
because `seat`/`citizen` make mint order observable and the storage
laws ride on it.

## F. View architecture

```text
              fractal kernel  (generators only)
                     |
        +------------+--------------------+----------------+
        |            |                    |                |
   encoding law   epistemic view      structural view   storage laws
   (sequences,    (G=(Q,R): the       (Gamma=(S,P):     (compiled
    pairs)        four state names     admission laws    blocks, CSR,
        |          read off two        over node         payload pools,
        |          branch statuses)    shapes)           borrowed slices)
        |
   +----+---------+-----------------+
   |              |                 |
 relation      directed          ordinary
 reading       reading           reading      ...  mesh schema, tags,
 (has/image)   (tail/head/deg)   (neighbours)      forms, constitution
```

Demonstrated in the spike, all over one storage:

- **Relation / operator unity**: `rel_has` reads a pair list as a
  binary relation; the residual `R` of example 4 is a sequence of
  (target, expression) pairs read by the same evaluator that computes
  the calculator — no separate inheritance roots for relation /
  residual / operator / differential operator. Differences are laws
  in readings.
- **Ordinary / directed as interpretations**: `tail_of/head_of/
  out_degree` and `neighbour_count` are two vocabularies over the one
  incidence list; neither duplicates storage; the boundary pair heads
  nowhere and neighbours nobody, by refusal rather than imaginary
  member.
- **The epistemic view**: the pair node's statuses read directly as
  COMPUTATIONAL-GRAPH.md's names — (UNKNOWN, KNOWN) is the operator
  graph; after the solve's one lawful `realize`, (KNOWN, KNOWN) is
  realized; a deliberately inconsistent pair still reads realized,
  because solved is another word. NULL rows add what the 2×2
  vocabulary could not say: *this problem has no such seat at all.*

The seam law (the Phase-2 verifier caught four rows violating it):
**the kernel answers status, realization, identity and enumeration —
nothing else.** `edge_tail`, `adjacent_vertices`, `num_vertices` and
every other schema word is view vocabulary; any design that lets the
kernel answer ordinary-graph questions has smuggled a schema back
into the primitive.

## G. Old-type → new-role deletion map

All nine competing roots die; exactly one ontology remains. What each
was for survives as a role:

| dying root | survives as |
|---|---|
| `graph_structure::relational_graph` | VIEW — the Γ=(S,P) admission laws (signature validity, each-domain-once, materialized-only ownership) over node shapes; the container type and its `held_*` seats die (branch chains carry heterogeneity natively) |
| `graph_state::computational_graph` | VIEW — the epistemic reading: four canonical names off two branch statuses; its two-seat SHAPE is what the kernel generalizes, and its single ⊥ is refined into NULL vs UNKNOWN |
| `graph_grammar::graph` (four-root grammar) | DELETE — admission creed survives as doctrine, not code |
| `interface_graph::graph`, `::digraph` | DELETE — algorithms consolidated into the relation-centered survivors (below); caches die with the module |
| `class_graph::stored_graph` | STORAGE donor — its four compressed lists are one of three rival CSRs; only ONE survives (below) |
| `class_stored_graph::{stored_graph, stored_digraph}` | DELETE |
| `class_graph_mesh::mesh` | VIEW (metric fields over Γ's carriers) + STORAGE (dense geometry) |
| `class_mesh::mesh` | mesh schema VIEW + geometry STORAGE; its graph plumbing dies |

Role totals from the 57-module classification (11 agents, verified by
an adversarial completeness check — census exact, every LOC figure
matches `wc -l`):

```text
KERNEL      1 module    101 LOC   graph_identity (the only production
                                  code the primitive absorbs; the
                                  fractal mechanics are greenfield)
VIEW       17 modules  4069 LOC
ALGORITHM  24 modules 11392 LOC
STORAGE     4 modules  2386 LOC
DELETE      6 modules  2127 LOC
UNTOUCHED   5 modules  1548 LOC   (string/file/verbosity/loaders)
```

The verifier's second-order findings — duties the migration inherits,
because killing rival roots is not yet killing rival implementations:

1. **Duplicate algorithm suites**: topological order, BFS, colouring,
   components, coarsen, refine and partition each survive TWICE (the
   relation-centered branch and the legacy `interface_graph`/
   `class_graph_walk` branch). One survivor per computation — keep
   the relation-centered one, retarget its suite, delete the twin.
   `interface_graph` (2,643 LOC) dies almost entirely under this rule.
2. **Three rival CSRs** (`csr_relation`, `class_graph`'s four lists,
   `interface_graph`'s xadj/adj): exactly ONE storage law survives —
   `csr_relation`'s build is the donor; the other two die with their
   hosts.
3. **Orphaned vocabulary needing homes**: tags/frames/boundary sets
   (a tag view: `P_tag ⊆ A×T`, AGENTS.md 35 — `class_robin_condition`
   consumes it, so it must exist before its supplier dies);
   `GRAPH_SIDE_*` landing constants (rehome to the operator view);
   the generic operation seam that newton/marcher/minimizers call
   through (becomes the action-view contract, named, not scattered).

Full per-module table in the appendix.

## H. Estimated production LOC after migration

Today: 21,623 total; 20,075 graph-related.

| tier | estimate | grounds |
|---|---|---|
| kernel (node/branch/arena/identity) | ~450 physical (~250 code) | the spike is 426/227 and complete |
| storage laws (SoA arena, ONE compiled-block/CSR law, payload pools, borrowed-slice API) | ~900 | donor: `csr_relation`'s ~600 |
| core views (encoding law, Γ admission, epistemic names, relation/binary vocabulary, ordinary+directed profile, tags, mesh schema) | ~2,600 | today's VIEW category is 4,069 with duplication |
| algorithms after pairwise consolidation | ~6,000 | 11,392 minus the legacy twins (interface_graph −2,300, walk/quotient/refine dupes −1,200, class_mesh plumbing −800), plus nothing new |
| physics/statement clients | ~1,000 | mostly untouched in role |
| **graph-related total** | **~11,000** | **≈ 45% reduction** |

The brief's compression target — "thousands of lines to approximately
O(1000) core lines" — is met **for the core ontology**: kernel +
storage laws + the structural/epistemic views ≈ 1,900 physical,
≈ 1,000 code lines. The spike's 518 code lines already carry all four
inhabitations. Algorithms and physics are clients of the core, not
core; they halve by deduplication, not by recursion.

## I. Migration order

Each step lands beside the tower with every certified suite as an
unchanged oracle; one deletion per commit; adapters live outside the
kernel and die last.

1. **Kernel to production** (`graph_kernel.f90` + arena), with the
   red team's hardening built in: token stamp, SoA segmentation
   (8-byte structure nodes, separate payload pools), the three
   lifetime laws in the header, cursor walkers instead of O(k)
   `seq_item`.
2. **The storage law**: compiled-sequence/CSR blocks + offset tables
   + borrowed-slice API, donor `csr_relation`; §66 benchmark gate
   before anything depends on it (targets: fibre query ≤ ~9 ns,
   field lookup ~1 ns, memory parity with CSR).
3. **Structural + epistemic views** over the kernel;
   `relational_graph` and `computational_graph` become adapters;
   calculator and learning towers green throughout.
4. **Carriers/relations/algebra** re-expressed (member sets = atom
   blocks; relations = tuple lists compiled to blocks); relation
   algebra retargets; graph-* suites green.
5. **Algorithm consolidation**, pairwise: retarget the survivor's
   suite, delete the legacy twin; `interface_graph`, `class_stored_
   graph`, `graph_grammar`, `graph_calculus` retire (Phase 11 debt).
6. **Fields and mesh**: payload pools, mesh schema view, tag view
   (before robin loses its supplier), loaders retarget.
7. **Inference through the action-view seam** (minimization, fitting,
   marching); the epistemic view's realize is the solve's public act.

## J. Unresolved laws and counterexamples

Open, recorded, deliberately not decided by the spike:

1. **Mutability doctrine.** `realize` is a lawful mutation of a
   structural object — AGENTS.md 54 says structural objects are
   immutable after construction. The reconciliation (the epistemic
   axis is single-assignment; structure never changes, knowledge
   grows once) needs an AGENTS.md amendment before production.
2. **Realize vs compiled blocks.** Growable knowledge and frozen
   contiguity conflict: a realize under a compiled storage block
   must either force recompilation or be refused post-compile. The
   generators audit proposes the cheapest oracle — an `epoch()`
   counter bumped on mint and realize — admitted only when the first
   real client (the FV adjacency freeze) exercises it; until then the
   discipline is an encoding law: *a transpose/compiled view is built
   after the last realize; a realize after it is the view's refusal.*
3. **Snapshot/restore.** Arena assignment over a populated arena is
   refused by the stamp for FRESH arenas, but `saved = a; a = saved`
   restores the SAME token at the same address — stale in-range
   seats answer silently. Adjoint checkpointing will want exactly
   this. Candidate: defined assignment that refuses a populated LHS
   ("a universe is never overwritten while it hosts"), making
   restore an explicit arena-level act with its own law.
4. **Cross-universe references.** Branches realize within their own
   universe only. Partition/assembly need structural mappings
   BETWEEN universes — a relation whose slots live in two arenas.
   Unrepresentable today; must be a first-class mapping object, not
   a weakened admit law.
5. **Reclamation.** The arena is append-only; nothing is ever freed.
   Regions/epochs (free a whole universe at once) fit the doctrine;
   per-node collection does not. Unresolved.
6. **Names.** The kernel is nameless by design (identity is not a
   name). Where the reader's names live — view-side name tables, or
   name atoms — is undecided; the naming pass's "the name is the
   reader's" points at views.
7. **Uninterpreted status pairs.** (NULL, KNOWN) and (UNKNOWN, NULL)
   are mintable and no canonical reading uses them; views refuse
   what they do not read. Legal-but-uninterpreted is the current
   law; a future reading may claim them.
8. **View-voiced refusals.** Some poisonous shapes still die with
   kernel diagnostics inside walkers (e.g. (NULL, KNOWN) fed to
   `rel_has`). Everything halts — nothing is misread — but the
   diagnostic attribution belongs to the view. Cosmetic, listed for
   the production pass. Likewise a `boundary_of(a, x)` constructor
   should name the shape the readers already accept.
9. **Deep-spine indexing.** `tail_of`/`head_of` answer occurrence k
   without establishing extent; on a cyclic spine occurrence k is
   indistinguishable from k mod cycle. Law needed: occurrence
   indices are meaningful only under a previously established
   `seq_len`, or the walkers pay the extent check.
10. **Partially learned environments.** A binding found before an
    UNKNOWN tail is final (shadowing + finality) — yet `env_value`
    currently refuses the whole environment because it measures
    extent first. Sound but over-strict; the refined walker answers
    on first match and refuses only past the match-free prefix.

## The red-team ledger

Fifteen findings, four reviewers, disposition of each:

```text
FIXED IN SPIKE   evaluate hangs on cyclic expression  -> depth budget
                 relation walkers answer through UNKNOWN -> refusal
                 empty application unreadable -> identities answer NULL
                 arena reassignment silently rebinds -> token stamp
                 5 nonconforming nested-mint statements -> hoisted
                 env accepts numeral members -> refused
                 block_value k<1 hole -> refused
                 solver rebind hazard -> bind-at-most-once guard
                 seq_item collapsed the two absences in a message
                 graph_branch/branch(): lawless, uninhabited -> deleted
                 branch_name in the kernel -> demoted to views
                 same_as totality -> law stated + tested
DEFERRED TO J    epoch/staleness oracle; snapshot/restore; view-voiced
                 refusals; boundary_of; deep-spine indexing; partial
                 environments; occurrence-vs-set counts (legislated,
                 kept)
```

Performance (i5-7500, gfortran-15, N=100k):

```text
                        -O0 suite flags     -O2
raw array walk               2.5 ns/el      1.8 ns/el
block road (storage law)    43.4 ns/el     15.0 ns/el   floor ~1 ns
                                                        with hoisted
                                                        guards
semantic road (spine)      113.7 ns/el     43.4 ns/el   ~2 ns/hop L1
                                                        floor; DRAM-
                                                        bound at mesh
                                                        scale
production references: csr borrowed view ~9 ns/query,
old stored_graph ~22 ns/query (doc/phase0-benchmark.md)
```

The spike's naive walkers are O(E²) by construction (`seq_item`
inside loops) — acceptable for a law suite, disqualifying for
production; cursor walkers and the compiled road are mandatory
(step 1–2 of I).

## Phase 3 — minimum generators

The audit (24 public symbols) closed at: **constants ×3, arena ×5
(ctor, node, value_atom, citizen, population), handle ×8 (status,
known, realize, same_as, carries_value, value, seat, universe_size),
specs ×4** — with `universe_size` honestly annotated as population's
answer re-exported across the capability wall (a reader must never
gain mint authority), and two evictions executed (`graph_branch` and
`branch()` deleted; `branch_name` demoted to views). Non-generability
proofs, one line each, are in the audit transcript; the load-bearing
ones: `status` is the only observation separating the two absences;
`known` is the only reader of the private target; `realize` is the
sole transition; `same_as` cannot be reproduced extensionally
(equal-extension atoms differ); `value_atom` is the only writer of
the private payload; `citizen` is the only re-entry to unreferenced
seats; `seat` is the only exit to the coordinate the storage laws
ride.

COMMUTATION holds where it must: `node` takes both branches
atomically (no first and second act); realizes on distinct
(node, side) pairs commute; the same pair twice is a lawful refusal;
mint/mint does NOT commute and is not supposed to — seat order is
mint order IS the contiguity promise, now stated at the kernel where
`seat`/`citizen` make it observable.

## Phase 4 — the four inhabitations

One kernel, no domain-specific changes, suite output:

```text
FRACTAL_CALCULATOR_RESULT = 20      (2+3)x4 off four atoms
conservation, fractally             interior differences telescope
the pair reads 'operator graph'     (UNKNOWN, KNOWN)
solve is the transition             bottom became knowledge, once
FRACTAL_REALIZED_RESULT = 20        q(c)=5, q(e)=20 - the old oracle
an inconsistent pair still reads 'realized' - not 'solved'
```

## Appendix — the full classification, module by module

| module | LOC | role | survives as |
|---|---|---|---|
| `graph_identity.f90` | 101 | **KERNEL** | token + mint_token -> the kernel's node identity itself: atom handles minted by the arena. matches/declared become handle equality and the undeclared-vs-declared (serial zero) distinction; the sign-once law becomes arena allocation discipline; the (image, serial) pair pre-fits a distributed arena. The type may dissolve into the arena h... |
| `class_graph.f90` | 1195 | **STORAGE** | class_graph::stored_graph DIES as a public root (its seat is taken by the fractal graph read through the ordinary-graph VIEW — resolving the two-stored_graphs name collision the audit flagged), but its body IS the surviving CSR STORAGE law: four compressed lists built once at construction (xinc/einc, xadj/vadj, xout/eout, xin/ein), hea... |
| `graph_binary_relation.f90` | 619 | **STORAGE** | csr_relation (xfwd/tgt/xbwd/src, counting-sort build, stamp-array collapse, O(degree) fibre slices, row-scan has) -> the canonical STORAGE law for binary relation readings, hidden behind the semantic interface exactly as the brief prescribes. Splits: binary_relation's arity-two vocabulary (source/target/image/preimage, view-vs-copy tie... |
| `class_graph_field.f90` | 395 | **STORAGE** | The dense realization behind the field view: the one-live-store rule (ivals/rvals/cvals/lvals/svals, exactly one allocated), the interleaving law position = (entry_position-1)*ncomp + component, and the size-check invariant stored = domain%size()*ncomp all survive as STORAGE laws of the arena's field payload. The public type 'field' as... |
| `class_graph_stencil.f90` | 177 | **STORAGE** | The compiled realization of R: weights-on-edges and constants-on-vertices become KNOWN branches of an operator node read through a CSR-style STORAGE law; the SpMV walk in stencil_apply is the traversal that law owns. Adjudication of roots it touches: class_graph::stored_graph (its `pattern` component) DIES — audit already marks it LEGA... |
| `graph_profile.f90` | 544 | **VIEW** | ordinary_graph_view -> VIEW: the T/H two-relation schema reading with its admission laws (E and V distinct, exactly one tail, at most one head) and the whole derived vocabulary (edge_tail/head, outgoing/incoming/incident_edges, adjacent/outgoing/incoming_vertices). Its 'boundary half-edge is an ABSENCE in H' becomes the kernel's NULL b... |
| `graph_carrier.f90` | 500 | **VIEW** | member_set contract (size/member/members/has/local_index, the two enumeration laws, is_subobject_of) -> the 'domain' reading of an atom list (right-nested chain of (NULL,NULL) nodes read as a set). Splits: counted_set's one-integer body -> STORAGE (an implicit 1..n range, no atoms materialized); subset_set's roll(:) -> STORAGE, its emb... |
| `graph_relation.f90` | 463 | **VIEW** | abstract relation contract (arity/domain/has/num_tuples/tuples/materialized) + slot signature -> the 'relation' reading of a node whose branch chain lists tuple nodes (tuples themselves right-nested chains); the construction refusals -> the schema-validity law of that reading; materialized() -> the safe-to-own law at the storage bounda... |
| `class_robin_condition.f90` | 351 | **VIEW** | robin_condition (tag, a, b, c) -> a parameter tuple node whose tag atom binds to the mesh VIEW's tagged-edge subset; robin/dirichlet/neumann constructors -> literal node shapes; wall_relation (the whole wall in two numbers w = a/denom, v = c/denom), operator_coefficients, boundary_values, and the four lhs/rhs coefficient families -> de... |
| `graph_state.f90` | 350 | **VIEW** | ADJUDICATION: computational_graph DIES as a root — uniquely among the competing roots, its SHAPE is what the kernel absorbs: two class(*) seats, each occupied-or-bottom, IS the two-branch node, so the kernel is this type generalized and made self-similar. Caveat the kernel improves on: graph_state's single bottom (unallocated seat) con... |
| `graph_structure.f90` | 344 | **VIEW** | ADJUDICATION: relational_graph DIES as a container root — the fractal node is the only container. held_set/held_relation seats -> DELETE outright (right-nested branch chains carry heterogeneous collections natively; the arena carries the storage). The create_graph refusal ladder (signature validity, each-domain-once, materialized-only ... |
| `class_graph_step.f90` | 280 | **VIEW** | step_operator -> VIEW: an operator interpretation of a node shape (branches: action = KNOWN subgraph pointing at the wrapped operation; a0/a1/a2/hs/reach and qold/qolder as atoms riding branches); step_dependencies' fan-in motif -> a schema reading answered off that shape — its stored_graph pattern build (competing-root call site) beco... |
| `class_graph_mesh.f90` | 231 | **VIEW** | ADJUDICATION: class_graph_mesh::mesh is the surviving mesh, but not as a type root — the `extends(stored_graph)` crossing onto class_graph::stored_graph dies with all the competing roots (interface_graph::graph/digraph, class_stored_graph::stored_graph/stored_digraph, class_mesh::mesh, and the class_graph/graph_structure/graph_state/gr... |
| `graph_field_calculus.f90` | 169 | **VIEW** | The field READING of a graph node: f : S -> V over one member-set domain becomes the canonical 'values live here' interpretation of a subset/carrier node shape. The abstract graph_field contract (name/units/domain/shape) survives as the field-view interface; the five GRAPH_FIELD_* kind tags survive as the view's payload vocabulary; the... |
| `class_conduction.f90` | 142 | **VIEW** | conduction's K tensor -> a parameter tuple node (nine atoms, or one tensor atom) in the one graph; isotropic/tensor constructors -> literal node shapes; normal_conductivity (keff = n^T K n) and edge_coefficients (keff * area, zero on headless faces) -> derived-vocabulary readings over the mesh VIEW's geometry STORAGE arrays. |
| `class_advection.f90` | 125 | **VIEW** | advection's velocity -> a parameter tuple node (three atoms); normal_speed (vn = v . n) and edge_coefficients (vn * area, zero on headless faces) -> derived-vocabulary readings over the mesh VIEW's normal/area STORAGE arrays; the upwind/central split stays out of it, as the banner insists — the scheme belongs to the calculus ALGORITHM ... |
| `graph_forms.f90` | 107 | **VIEW** | The form READING of a subset node: a basis family is a subset_set (a right-nested membership chain over a table of atoms) interpreted with three evaluation symbols — size_of, values(x,at), slopes(x,at,n). That interpretation survives verbatim as the forms/constitution view of the graph; restrict survives as an ordinary membership rebui... |
| `interface_mesh_loader.f90` | 104 | **VIEW** | The loader contract survives as the ingestion profile at the system boundary: 'a graph described by its edges' — member lists, edge/face/cell->vertex incidence, tag identity, coordinates as the only value payload. Its deliverable is retyped from ~25 out-arguments to one graph value (incidence lists as right-nested branch chains, coordi... |
| `class_harmonic_form.f90` | 101 | **VIEW** | The wave constitution of the form view: the law that a three-atom basis table reads as {1, sin(k.(x-at)), cos(k.(x-at))} with chain-rule slopes, parameterized by the wavenumber — which survives as an attribute (identity data) on the form node. As with the polynomial twin, the subset/counted_set construction dissolves into atoms plus a ... |
| `class_diffusion_statement.f90` | 99 | **VIEW** | diffusion_statement -> a constitution reading (VIEW): translates conduction + robin vocabulary into the assembly's neutral scales-per-face and wall-relation-per-tagged-face vocabulary, then delegates the act to fitted_balance_stencil (ALGORITHM); the polynomial_form default is part of the reading's schema. Root adjudication: its graph_... |
| `class_polynomial_form.f90` | 81 | **VIEW** | The polynomial constitution of the form view: the law that table entries 1..4 read as {1, x-at} with slopes {0, n} — the degree-one Taylor evaluation rule attached to a four-atom basis table. The subset/counted_set scaffolding in create() dissolves into plain arena atoms plus a membership branch; only the evaluation/slope formulas surv... |
| `module_solve_mode.f90` | 78 | **VIEW** | FORWARD/REVERSE (audit: INFERENCE) -> the direction-of-reading vocabulary for a relation node's two branches — tail-axis versus head-axis traversal of the same edges; WHOLE/DIAGONAL/LOWER_TRIANGLE/UPPER_TRIANGLE (audit: OPERATOR) -> subrelation selections, views filtering an operator-as-relation's edge nodes by row/column comparison; i... |
| `interface_graph.f90` | 2643 | **ALGORITHM** | The competing roots interface_graph::graph and interface_graph::digraph DIE, with vertex/edge: a vertex becomes an identity-bearing (NULL,NULL) node, an edge a tail/head branch pair, and directedness stops being a type — it becomes a property of the adjacency interpretation (a VIEW), not a second root. What dominates and survives is ~2... |
| `class_mesh.f90` | 1870 | **ALGORITHM** | ADJUDICATION: type class_mesh::mesh -> DELETE — it is the losing root of the mesh pair (extends interface_graph::graph, emits class_stored_graph::stored_graph from node_graph, straddling two dying graph stacks); class_graph_mesh's profile is the sole mesh reading. What survives is the module's bulk: invert_connectivities + create_mesh ... |
| `class_graph_differential_operator.f90` | 1199 | **ALGORITHM** | The S/G/D elementary steps, their reversed (adjoint) twins, and the parity chain runners (run_vertex_chain, apply_on_edges/vertices) -> ALGORITHM: the matrix-free interpreted evaluation of R over arena incidence. The declarative face — order/landing/coefficient/spacing/measure/boundary parameters and the named constructors gradient/int... |
| `class_graph_reduction.f90` | 773 | **ALGORITHM** | The four-step protocol (initialize/accumulate/combine/finalize) with its order-independence and finish-once guarantees, plus the real/complex/logical accumulation roads -> ALGORITHM: the partition-safe fold over a field's domain. The REDUCE_*/BROADCAST_* rule constants, the measure-seat-as-inner-product reading, and the adjoint-pair la... |
| `class_graph_partitioner.f90` | 605 | **ALGORITHM** | assign_owners_linear/assign_owners_breadth_first, gather_part, the edge-keep/tail-ownership rule, and carry_field transport -> ALGORITHM over the fractal graph; the part-to-whole record (vglobal/vowner/eglobal/eowner) it stamps on the part -> VIEW annotation (KNOWN branches from part nodes to their global identities); parent graph_calc... |
| `graph_minimization.f90` | 443 | **ALGORITHM** | type minimizer + solve_interface + delegation vocabulary (matvec/inner_product/norm/sweep_order/diagonal/constant) -> ALGORITHM base of the solver family; the operation face (solver_domain/solver_apply, 'a solver IS an operation') -> re-expressed as a VIEW reading of the solver node so it composes as an operator; seats on/coupling curr... |
| `class_graph_coarsener.f90` | 369 | **ALGORITHM** | Pairwise-matching walk (blocks_of), block-face redraw with dedup (coarsen_graph), and the add-vs-average lift (coarsen_data) -> ALGORITHM; the block_of aggregate map it computes -> a VIEW relation in the fractal graph (each fine node gets a KNOWN branch to its block node) consumed by multigrid; parent graph_calculus::graph_coarsener ->... |
| `class_graph_assembler.f90` | 361 | **ALGORITHM** | assemble_graph rename-home walk and gather_field owned-only collection -> ALGORITHM (the P-inverse of the partitioner, kept as its paired transform); the law assemble(partition(G))==G -> a contract test on that pair, not code; parent graph_calculus::graph_assembler -> DELETE; stored_graph construction re-targets the arena as with the p... |
| `class_graph_walk.f90` | 339 | **ALGORITHM** | colour (greedy proper colouring), breadth_first (visit order and depth), components (flood labelling) -> ALGORITHM, verbatim: these are the walks and orderings the ALGORITHM clause names explicitly, rewritten against arena adjacency for O(degree) neighbour reads. The WALK_* rule constants and the integers-per-vertex answer (colours, ra... |
| `class_graph_marcher.f90` | 306 | **ALGORITHM** | march / march_adjoint / state_seat / read_statement -> ALGORITHM: the forward walk, the explicit-rule adjoint reverse walk, and one governed solve per edge over the time-chain view; instants -> its chain construction retargets from class_graph::stored_graph (a competing root adjudicated with its own module, dies as a root) to building ... |
| `graph_fitting.f90` | 258 | **ALGORITHM** | fit's weight-finding minimization (B B^T lambda = r, w = B^T lambda, distance-priced, roster-masked) -> ALGORITHM over any point-constellation view of the one graph; form_optimizer abstract adapt contract -> ALGORITHM family head for model adaptation. Its parentage in graph_grammar::graph / graph_operation dies with that legacy root: d... |
| `graph_relation_algebra.f90` | 247 | **ALGORITHM** | restrict_slot / project_slots / compose_binary -> the relation-algebra algorithms: walks over tuple-node chains (filter by subobject membership, slot extraction, O(/R//S/) witness-search join) that build new relation-shaped nodes; their eager results land in the STORAGE laws (flat table for restrict/project, CSR for compose); the embed... |
| `class_graph_refiner.f90` | 244 | **ALGORITHM** | Split/sibling-join/face-carry-down computation (refine_graph) and parent-value injection (refine_data) -> ALGORITHM; the arithmetic child_of map ((v-1)*split+i) is a computed relation needing no stored map, so nothing of it goes to STORAGE; parent graph_calculus::graph_refiner -> DELETE; stored_graph output re-targets the arena. |
| `graph_algorithms.f90` | 242 | **ALGORITHM** | reachable (BFS walk) and topological_order (deterministic Kahn, refusal on cycles) -> ALGORITHM: walks and orderings over the directed-adjacency reading of kernel nodes, staying free procedures acting on a view from outside. sources/sinks straddle — empty-fibre scans are arguably derived vocabulary (the audit calls all four VIEW) — but... |
| `class_fitted_balance.f90` | 236 | **ALGORITHM** | fitted_balance_stencil + neighbourhood_of + grow -> ALGORITHM: the per-edge walk (structure), position gather (data), fit apply (algebra), incidence exchange, compiled into one operator. The compiled stencil_operator it emits is STORAGE (CSR law behind the operator view). Root adjudication: its throwaway class_graph::stored_graph const... |
| `class_graph_linearization.f90` | 234 | **ALGORITHM** | The differencing recipe — derivative_freeze, derivative_apply's two applies and quotient (S(q+eps v)-S(q))/eps, and the answered_on domain guard — -> ALGORITHM serving minimization/newton one rank up; the wrapper identity (of / at / base / step as branches of a tangent node) -> a thin VIEW shape so governors keep seeing an ordinary lin... |
| `class_graph_balance.f90` | 200 | **ALGORITHM** | The exactly-once incidence scatter in balance_apply (each edge gives to its head, takes from its tail; headless edges land on the tail alone) -> ALGORITHM: the generic edge-to-vertex fold that assembles R. The balance's shape — edge_terms(:) plus a source — -> VIEW: the constitution reading of a residual node whose branches are term-no... |
| `class_graph_multigrid.f90` | 193 | **ALGORITHM** | Two-grid cycle (smooth/restrict/solve-coarse/prolong/correct/smooth) and the Galerkin RAP fold over aggregates in setup() -> ALGORITHM in the minimizer family; the block stencil_operator it compiles -> a materialized STORAGE cache of the coarse operator; the aggregates(:) array it holds -> reads the coarsener's VIEW relation rather tha... |
| `class_graph_gmres.f90` | 160 | **ALGORITHM** | type gmres (restarted Arnoldi solve + restart parameter) -> ALGORITHM concretion; the Krylov basis/Hessenberg arrays are transient solve-local scratch, not STORAGE, and carry over untouched. |
| `class_graph_newton.f90` | 113 | **ALGORITHM** | type newton (solve + governed inner minimizer seat) -> ALGORITHM concretion: governance of inference. The tangent it freezes (difference_linearization) is not its own — it stays an OPERATOR citizen classified with class_graph_linearization, read as a VIEW-side R-candidate in the compressed ontology; newton keeps only the freeze/attach/... |
| `class_graph_conjugate_gradient.f90` | 98 | **ALGORITHM** | type conjugate_gradient (solve) -> ALGORITHM concretion; structure-free (never touches coupling), so it survives verbatim over the fractal graph via matvec/inner_product/norm alone. |
| `class_graph_gauss_seidel.f90` | 93 | **ALGORITHM** | type gauss_seidel (colour-swept solve + omega, SOR absorbed) -> ALGORITHM concretion; its sweep_order dependency is the colouring walk, itself ALGORITHM, applied to the coupling handle rather than any legacy root. |
| `class_graph_jacobi.f90` | 86 | **ALGORITHM** | type jacobi (solve + omega) -> ALGORITHM concretion of the minimizer family, unchanged in substance; only its inherited seats retype with the base. |
| `class_form_pruner.f90` | 80 | **ALGORITHM** | pruner % adapt -> ALGORITHM: computes per-member visibility norms over the constellation about its own centroid and calls shape % restrict; the threshold is one parameter atom on the pruner node. Nothing else in the module; it survives whole as a form-adaptation computation over the form VIEW. |
| `graph_grammar.f90` | 567 | **DELETE** | graph_grammar::graph DIES (one of the competing roots): its ordinary vertex/edge reading — edge_tail/edge_head/edge_has_head, named sets, the owned/borrowed/overlap frame — survives as the ordinary-graph VIEW over (NULL,NULL) identity nodes joined by tail/head branch pairs (the reading graph_profile::ordinary_graph_view already gives f... |
| `graph_calculus.f90` | 526 | **DELETE** | All ten abstract types die with their graph_grammar parents. What survives: the staged reduction discipline (initialize/accumulate/combine/finalize-once, order-independent combine) as an ALGORITHM contract on reductions; the round-trip laws (assemble(partition(G))=G, coarsen(refine(G))=G, average(copy(c))=c / sum(share(c))=c) stay in t... |
| `class_graph_functional.f90` | 490 | **DELETE** | Nothing survives as a type. Its own header concedes the compression: 'this is not a new kind of thing: it is the field at domain size one'. Under the one-graph ontology a functional is the field VIEW over a one-atom domain, and the dense field STORAGE already realizes it at n=1 — so the entire duplicated ten-adapter contract answered a... |
| `class_stored_graph.f90` | 269 | **DELETE** | class_stored_graph::stored_graph and ::stored_digraph DIE with their interface_graph roots — in the two-stored_graphs collision, this is the duplicate stack (its actual storage lives in the parent's vertices/edges/xadj arrays, and class_graph's CSR law is the one kept). Of its five constructors: create/create_directed collapse into the... |
| `class_array_mesh_loader.f90` | 159 | **DELETE** | Nothing of the type survives: under the one-graph ontology an in-memory mesh description IS a graph value, so a shim whose entire behavior is copying stored incidence arrays out to the contract's arguments is absorbed — derived meshes (class_mesh::refined's front door, its main consumer) become algorithms that return graphs directly. I... |
| `class_mesh_builder.f90` | 116 | **DELETE** | The module body is bridge code by its own banner ('the builder dies when the measurement machinery is ported'): it runs the dying class_mesh root once and rearranges its arrays into the new seat. Only the entry-point symbol mesh_from_gmsh (audit: DATA, 'the durable new-world entry point') survives, reseated as the load-pipeline ALGORIT... |
| `class_gmsh_loader.f90` | 674 | **UNTOUCHED** | Survives as-is: the MSH 4.1 parsing half of a mesh file loader (block-structured $Nodes/$Elements, entity->physical-tag lookup, sparse node-tag resolution). Only its deliverable retargets — get_mesh_data's output feeds the retyped ingestion VIEW, i.e. it hands over one graph value instead of 25 arrays — while the format-decoding body s... |
| `module_mesh_utils.f90` | 397 | **UNTOUCHED** | Survives whole as the leaf toolbox: distance/cross_product (tape-measure math used by the ported measurement ALGORITHMs), find/is_subset/isort (plain integer-array utilities used by the parser and wiring), elem_type_face_count/dimension/vertex_count and order_face_vertices (gmsh element-shape tables serving the parsing half). None of i... |
| `class_string.f90` | 218 | **UNTOUCHED** | string with tokenize / asinteger / asreal / equals / print stays a plain character-value utility feeding parsers; audit tags it DATA only as the nearest of six fits, tied to no carrier or graph domain. |
| `class_file.f90` | 218 | **UNTOUCHED** | file with open / close / get_unit / read_line / read_lines / get_num_lines stays the sequential line reader behind the mesh loaders' parsing halves; its banner's 'a file is a chain' is rhetoric — the module imports only class_string and no graph type. |
| `module_verbosity.f90` | 41 | **UNTOUCHED** | verbosity and set_verbosity stay exactly as they are — a protected global diagnostics knob and its single writer; audit files them VIEW only in the presentation sense, not as an interpretation over graph. |
