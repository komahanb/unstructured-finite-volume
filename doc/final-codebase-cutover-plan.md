# The final cutover: one ontology, many views

The target organization is

    ontology · views · maps · algorithms · representations · systems

This document is the evidence that must exist before any of it is moved.
It classifies every one of the 63 modules in `src/` into exactly one
bucket, and it does so from the **use-edge graph**, not from names. Where
the measurement contradicts the target picture, the measurement wins and
the contradiction is recorded here rather than quietly implemented.

Nothing is renamed, split or deleted by this document. It is the map that
tells the later phases which moves are legal.

---

## 1. How the edges were taken

The `use` statements of all 63 `src/*.f90` modules and all `test/**/*.f90`
programs were read, and the directed import graph built over them. Two
questions were then asked of every module:

    who imports it        src consumers, by MODULE name, not by path
    is it reachable       from the 49 modules the 32 test suites import
                          directly, following src->src edges to closure

A module with zero `src` consumers and many `test` consumers is **not**
dead. That is precisely the shape `graph_carrier` had before 133d29c, and
reading it as dead is the error this method exists to prevent.

**Result: the reachable closure is all 63 modules. Nothing in `src/` is
import-dead.** The `dead` bucket is empty and stays empty. Every deletion
the remaining phases want must therefore be *earned* by rewriting
consumers; none can be harvested by finding an orphan.

---

## 2. The bucket table

    bucket                 count   what the bucket means here
    ---------------------  -----   ------------------------------------
    core                       3   the kernel, and the laws every
                                   identified citizen draws on
    map                        3   identity -> answer, stored outside
    view                       6   an interpretation; stores nothing
    algorithm                  7   acts on structure, is not structure
    representation            11   how something is stored HERE
    system                    25   PDE, solve, discretization, statement
    legacy-compatibility       8   the ordinary-graph layer to retire
    dead                       0   (empty, and proven so)

### core — 3

| module | role | note |
|---|---|---|
| `fractal_graph` | \(G=(B_1,B_2)\), the single graph ontology | **byte-frozen**; 15 src and 149 test importers |
| `graph_identity` | `token`: mint once, copy whole, match by serial | infrastructure beneath the tower, not a level |
| `graph_relation` | a relation: identity + signature, the signature being a sequence of graph identities | `stored_relation` also holds `entry(:,:)` — a tuple table. **PR4 splits that half out**; the identity+signature half is core |

> **The bucket is not "ontology".** The ontology is one module and one
> pair of types — `fractal_graph`'s `graph` and `graph_branch` — and
> AGENTS.md admits nothing else to it. `graph_identity` says of itself
> "not a level"; `graph_relation` owns a type that is not built from
> branches. Both are here because `core` is what the tower stands on,
> not because they are ontology. If a later phase wants a bucket that
> means ontology exactly, it has one member.

### map — 3

| module | keyed on | answers |
|---|---|---|
| `graph_set_map` | set graph identity | which representation describes it (`size_of`, `member_of`, `index_in`, `extent_of`) |
| `graph_label_map` | set graph identity | what it is called |
| `graph_inclusion_map` | subobject identity | what it was declared into (`include_in`, `included`, `declared_into`) |

These three are already exactly what the target calls maps. No phase
below moves them.

### view — 6

| module | reads the graph as | publics |
|---|---|---|
| `graph_set_view` | a finite set | `set_size`, `set_member`, `set_has`, `set_local_index`, … |
| `graph_sequence_view` | a finite sequence | `sequence_size`, `sequence_element`, … |
| `graph_relational_view` | \((\mathcal S,\mathcal P)\) | `relational_binding`, `member_set_at`, `relation_at` |
| `graph_epistemic_view` | \((Q,R)\) — data and residual | `has_data`, `data_of`, `residual_of` |
| `graph_profile` | the ordinary directed graph as a **schema over two relations** \(T,H \subseteq E\times V\) | `ordinary_graph_view`, `directed_adjacency_view` |
| `graph_field_calculus` | a domain carrying values | `graph_field`, the five `GRAPH_FIELD_*` kinds, `set_graph` |

### algorithm — 7

`graph_algorithms` (`sources`, `sinks`, `reachable`, `topological_order`),
`graph_relation_algebra` (`restrict_slot`, `project_slots`,
`compose_binary`), `class_graph_partitioner`, `class_graph_assembler`,
`class_graph_coarsener`, `class_graph_refiner`, `class_graph_walk`.

The last five are the `graph_transform` and higher-order `graph_operation`
concretions: graph-to-graph verbs with no PDE semantics of their own.

### representation — 11

`graph_set_representation`, `graph_binary_relation`, `class_graph_field`,
`class_string`, `class_file`, `module_mesh_utils`, `module_verbosity`,
`module_solve_mode`, `interface_mesh_loader`, `class_gmsh_loader`,
`class_array_mesh_loader`.

> **Plan correction.** The eight-bucket list has no `io` bucket, but the
> final layout names `src/io`. The last five modules above are io and
> named-constant services, not representations of a mathematical object.
> They are bucketed `representation` because that is the bucket whose
> *rules* govern them correctly — no semantic identity, no map ownership,
> copyable into actions — and they are flagged here as PR6's `src/io`.

### system — 25

`class_graph_differential_operator`, `class_graph_balance`,
`class_graph_stencil`, `class_graph_step`, `class_graph_marcher`,
`class_graph_linearization`, `class_graph_reduction`,
`class_graph_functional`, `class_graph_multigrid`, `class_graph_newton`,
`class_graph_jacobi`, `class_graph_gauss_seidel`, `class_graph_gmres`,
`class_graph_conjugate_gradient`, `graph_minimization`, `graph_fitting`,
`graph_forms`, `class_polynomial_form`, `class_harmonic_form`,
`class_form_pruner`, `class_fitted_balance`, `class_diffusion_statement`,
`class_conduction`, `class_advection`, `class_robin_condition`.

### legacy-compatibility — 8

`graph_grammar`, `graph_calculus`, `class_graph`, `class_graph_mesh`,
`interface_graph`, `class_stored_graph`, `class_mesh`,
`class_mesh_builder`. Each has its own ledger entry in §4.

---

## 3. Six findings, four of which correct the plan

### 3.1 `graph_grammar` is 54% re-export

Every `use graph_grammar, only :` clause in `src/` and `test/` was tallied
by symbol. 183 symbol-imports across 68 files:

    OWNED by graph_grammar                    RE-EXPORTED from elsewhere
    ------------------------------            --------------------------
    graph               57  (21 src)          graph_field      55  (16 src)
    grammar_graph=>graph 2  ( 0 src)          set_graph        31  (21 src)
    graph_operation     24  ( 9 src)          GRAPH_FIELD_*    13  ( 8 src)
    graph_transform      1  ( 1 src)
    ------------------------------            --------------------------
                        84                                     99

`graph_field`, `set_graph` and the five kind constants are already public
in `graph_field_calculus` (which itself re-exports `set_graph` from
`fractal_graph`). **99 of the 183 imports can be redirected with no
semantic change of any kind** — the target module already exists and
already exports the symbol. This is the cheap, risk-free half of PR2 and
it should go first, alone, so the expensive half is measured against a
smaller number.

### 3.2 `graph_transform` has exactly one consumer

`graph_calculus`, and nobody else in `src/` or `test/`. It is extended
there four times (`graph_partitioner`, `graph_assembler`,
`graph_coarsener`, `graph_refiner`). A `graph_transform_view` module is a
one-file, one-consumer move — the smallest boundary in the whole cutover
and the right place to prove the pattern before touching `graph`.

### 3.3 The ordinary view carries a partition frame that is not a view question

`graph_grammar`'s abstract `graph` declares 36 deferred bindings.
Twenty-eight fall inside the five capabilities the plan allows a view.
Eight do not:

    num_parts  has_part_relation
    global_vertex_index  global_edge_index
    part_vertex_index    part_edge_index
    vertex_owner_part    edge_owner_part

In `class_graph` these are backed by per-graph *state* — `nparts`, `cut`,
`vglobal(:)`, `eglobal(:)`, `vowner(:)`, `eowner(:)`. Read for what they
are:

    vglobal(k)                the part's k-th vertex, named globally
                              = the EXTENSION of a subobject
    reverse_lookup(vglobal,v) where global v stands in the part
                              = the subobject's local index
    vowner(v)                 which part owns v
                              = an integer FIELD on the vertex set
    cut, nparts               whether this graph is a part, of how many

This is **not** a sixth capability, and PR2 does not have to stop. It is
inclusion/provenance and a field, stored inside the graph instead of in
the two homes that already answer them:

    part vertex set  --include_in-->  global vertex set   (inclusion_map)
    listed_set_representation of it   (set_map)
        member_of(k) == global_vertex_index(k)
        index_in(v)  == part_vertex_index(v, own part)
    owner            an integer field on the global vertex set

The frame therefore leaves the view in **PR3**, not PR2. PR2 must carry
all eight bindings into `graph_ordinary_view` unchanged, or it will
strand `class_graph_partitioner` and `class_graph_assembler` mid-flight.

### 3.4 Two unrelated `graph` hierarchies, and two types called `stored_graph`

The plan's PR3 reads `class_stored_graph` as the ordinary-graph
representation and `class_graph` as the possible abstract view. **The
code is the other way round, and they are not even in the same
hierarchy.**

    graph_grammar :: graph          (abstract, 33 deferred)
      `-- class_graph :: stored_graph          1303 lines, 9 src consumers
            `-- class_graph_mesh :: mesh        231 lines, 6 src consumers

    interface_graph :: graph, digraph  (abstract, 2650 lines, its own world)
      |-- class_stored_graph :: stored_graph, stored_digraph   269 lines
      `-- class_mesh :: mesh                                  1870 lines

Two modules export a type named `stored_graph`; the two types share no
ancestor, no contract and no consumer. This name collision is the single
largest ambiguity in `src/`, and removing it is the reduction PR3 must
deliver (global rule 4). `class_graph` is the concretion — it is not a
candidate for "keep only if it is the abstract ordinary view".

### 3.5 The second hierarchy is an island behind one function

Everything in the `interface_graph` world is reachable from exactly one
entry point:

    class_mesh_builder :: mesh_from_gmsh
      |-- class_graph_mesh   (the tower seat)
      `-- class_mesh  --> class_stored_graph --> interface_graph
             |-- class_array_mesh_loader, class_gmsh_loader,
             |   interface_mesh_loader, class_file, class_string,
             |   module_mesh_utils, module_verbosity
             `-- module_solve_mode

`class_mesh` has one src consumer (`class_mesh_builder`) and zero test
consumers. `interface_graph` has two src consumers, both inside the
island, and zero test consumers. The two test suites that reach any of it
— `graph-mesh` and `graph-minimization` — import `mesh_from_gmsh` and
nothing else from the island.

Those eleven modules are **6,678 lines** — 29% of `src/`'s 23,189 — and
not one of them is reachable except through `class_mesh_builder`'s own
116. A single function stands in front of all of it, and its header
already states its own deletion condition: *"The builder
dies when the measurement machinery is ported into the tower; its call
site will not change when that happens."* Global rule 5 is already
satisfied for this adapter. The blocker is not architecture — it is that
volumes, areas, normals, deltas and weights are computed in `class_mesh`
and nowhere else.

### 3.6 `graph_calculus` is a second legacy contract, and the plan does not name it

Its own header: *"THE LEGACY COMPATIBILITY CALCULUS."* Thirteen src
consumers, 9 test files across 7 suites. Its 29 symbol-imports are thin
and lopsided — `GRAPH_SIDE_VERTEX` alone is 14 of them, an absorbed
vertex/edge axis riding as a constant, exactly as the admission law
requires. It holds the nine legacy concretions of field, operation and
transform.

`graph_calculus` cannot be retired before `graph_grammar`, because it
extends `graph_transform` and `graph_operation`. It should be added to
PR2's scope as the module that *follows* the split rather than being
discovered during it.

---

## 4. The legacy-compatibility ledger

### `graph_grammar` — 661 lines

- **Consumers.** 24 src modules; 44 test files across 17 suites.
- **Role.** The ordinary binary-graph vocabulary, retyped onto the new
  ontology: abstract `graph`, `graph_operation`, `graph_transform`, plus
  re-exports of `graph_field`, `set_graph` and five kind constants.
- **Successors.** `graph_ordinary_view` ← `graph` (+ the frame, until
  PR3). `graph_operation_view` ← `graph_operation`.
  `graph_transform_view` ← `graph_transform`. `graph_field_view` **needs
  no new module** — it is `graph_field_calculus`, today. `set_graph` goes
  home to `fractal_graph`.
- **Blocker.** Nothing structural. 99 of 183 imports are re-exports and
  can move now; the remaining 84 are three type names across 68 files.
  The module dies when `grep -rl "use graph_grammar" src test` is empty.

### `graph_calculus` — 552 lines

- **Consumers.** 13 src modules; 9 test files across 7 suites.
- **Role.** The nine legacy concretions of field, operation and
  transform, plus the `GRAPH_SIDE_VERTEX`/`GRAPH_SIDE_EDGE` absorbed
  axis.
- **Successor.** Split with its parents: the transform concretions follow
  `graph_transform_view`, the operation concretions follow
  `graph_operation_view`. The two side constants belong wherever the
  ordinary view lands.
- **Blocker.** `graph_grammar` — it extends both abstracts.

### `class_graph` — 1303 lines

- **Consumers.** 9 src modules; 34 test files across 15 suites.
- **Role.** The concrete ordinary graph: edge list, derived adjacency,
  named/carved sets, neighbourhood queries, and the partition frame.
- **Successor.** `graph_ordinary_representation` (the name `stored_graph`
  cannot survive alongside `class_stored_graph`'s — §3.4).
- **Blocker.** The frame state (§3.3) must move to `inclusion_map` +
  `set_map` + a field first, or the representation carries maps, which
  global rule 6 and PR4 both forbid.

### `class_graph_mesh` — 231 lines

- **Consumers.** 6 src modules; 5 test files across 5 suites.
- **Role.** `mesh` extends `class_graph`'s `stored_graph` and seats the
  measurements as typed fields.
- **Successor.** `graph_mesh_view` (topology as a view) +
  `mesh_representation` (coordinates, volumes, areas, normals, deltas).
- **Blockers, both of them.** Its base class — it cannot be split before
  `class_graph` is. And its *supplier*: the measurement fields it seats
  are computed in `class_mesh` and carried across by
  `class_mesh_builder`, so `mesh_representation` cannot be given a home
  of its own until the island cut (§3.5) has one. It is the only legacy
  module that waits on both PR3 tracks, which is why it is not in the
  island's ledger despite being the thing the island exists to feed.

### `interface_graph` — 2650 lines

- **Consumers.** `class_mesh`, `class_stored_graph`. Zero test consumers.
- **Role.** The older abstract graph/digraph: retained CSR adjacency,
  traversal orders, orbits, power iteration, reverse-mode accumulation,
  the in-place partition stamp.
- **Successor.** None as a module. Its live algorithms are
  `graph_algorithms`' territory; its adjacency is a representation.
- **Blocker.** `class_mesh` (§3.5).

### `class_stored_graph` — 269 lines

- **Consumers.** `class_mesh`. Zero test consumers.
- **Role.** `stored_graph`/`stored_digraph` over `interface_graph`: the
  edge-list citizen, the quotient squint, the condensation.
- **Successor.** None. It dies with the island.
- **Blocker.** `class_mesh`.

### `class_mesh` — 1870 lines

- **Consumers.** `class_mesh_builder`. Zero test consumers.
- **Role.** Gmsh parsing plus **all** the geometric measurement:
  volumes from coordinates, areas by triangulation, normals, deltas,
  weights.
- **Successor.** `mesh_representation` for the arrays; the measurement
  code itself has no successor written yet.
- **Blocker.** The measurement machinery is not ported. This is the one
  blocker in the whole ledger that is real work rather than rewiring.

### `class_mesh_builder` — 116 lines

- **Consumers.** 0 src; `test/graph-mesh`, `test/graph-minimization`.
- **Role.** The declared bridge between the two worlds.
- **Successor.** Nothing — the call site `mesh_from_gmsh` is designed not
  to change.
- **Blocker.** Same as `class_mesh`, and stated in its own header.

---

## 5. What the evidence permits each remaining phase to do

    PR2  graph_grammar
         |-- FIRST, alone:  redirect the 99 re-export imports
         |                  (graph_field, set_graph, GRAPH_FIELD_*)
         |                  -> graph_field_calculus / fractal_graph
         |-- THEN:          graph_transform_view      1 consumer
         |-- THEN:          graph_operation_view     24 imports
         |-- LAST:          graph_ordinary_view      59 imports, and it
         |                  carries the frame UNCHANGED (finding 3.3)
         `-- graph_calculus follows, it does not lead (finding 3.6)

    PR3  the frame leaves the ordinary view: inclusion_map + set_map +
         an owner field. Only then can class_graph become
         graph_ordinary_representation, and only then does the
         stored_graph name collision resolve (finding 3.4).
         The island (finding 3.5) is a SEPARATE, later cut - its blocker
         is unported measurement code, not import structure.

    PR4  graph_relation           entry(:,:) -> table representation
         graph_binary_relation    csr_relation -> csr representation,
                                  binary_relation contract stays a view
         class_graph_field        -> graph_field_representation
         class_graph_mesh         -> view + mesh_representation

    PR5  systems consume views + representations. 25 modules, none of
         which define ontology today; the work is import rewriting that
         follows PR2-PR4, not new structure.

    PR6  directory layout, and the bucket table above is its input.
         Note that src/io has no bucket in the plan (see §2).

---

## 6. Verification record

Taken on the tree this document's parent commit produces — master plus
the benchmark repair — with one clean rebuild per suite:

    32 of 32 suites PASS, 0 FAIL

The reading survives either merge order, because this document adds no
code, changes no import, and is the only thing in its own commit.

`src/fractal_graph.f90` is byte-identical to master.
