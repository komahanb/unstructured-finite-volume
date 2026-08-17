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
    view                       7   an interpretation; stores nothing
    algorithm                  7   acts on structure, is not structure
    representation            11   how something is stored HERE
    system                    25   PDE, solve, discretization, statement
    legacy-compatibility       8   the ordinary-graph layer to retire
    dead                       0   (empty, and proven so)

The counts are as measured for PR1, over 63 modules. PR2 has since
carved two modules out of `graph_grammar` and deleted it:

    + graph_directed_view      view
    + graph_operation_view     view
    - graph_grammar            legacy-compatibility

so `view` is 8, `legacy-compatibility` is 7, and the total is 64. No
module has changed bucket.

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

### view — 8 (6 at PR1, plus the two carved from `graph_grammar`)

| module | reads the graph as | publics |
|---|---|---|
| `graph_set_view` | a finite set | `set_size`, `set_member`, `set_has`, `set_local_index`, … |
| `graph_sequence_view` | a finite sequence | `sequence_size`, `sequence_element`, … |
| `graph_relational_view` | \((\mathcal S,\mathcal P)\) | `relational_binding`, `member_set_at`, `relation_at` |
| `graph_epistemic_view` | \((Q,R)\) — data and residual | `has_data`, `data_of`, `residual_of` |
| `graph_profile` | the directed graph as a **schema over two relations** \(T,H \subseteq E\times V\) | `directed_incidence_view`, `directed_adjacency_view` |
| `graph_field_calculus` | a domain carrying values | `graph_field`, the five `GRAPH_FIELD_*` kinds, `set_graph` |
| `graph_directed_view` | the ordinary binary graph — vertices, edges, incidence, named sets, neighbourhoods | `directed_graph` (abstract). **Added by PR2**, carved from `graph_grammar`; renamed by PR3; carries the legacy partition frame, marked |
| `graph_operation_view` | the two verbs — within a graph, and between graphs | `graph_operation`, `graph_transform` (abstract). **Added by PR2**; the last of `graph_grammar` |

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

### legacy-compatibility — 7 (8 at PR1, less `graph_grammar`)

`graph_calculus`, `class_graph`, `class_graph_mesh`, `interface_graph`,
`class_stored_graph`, `class_mesh`, `class_mesh_builder`. Each has its
own ledger entry in §4.

~~`graph_grammar`~~ — **deleted in PR2**, once `use graph_grammar`
reached zero. It is the first module the cutover has retired, and it was
retired the way §1 says every deletion must be: by rewriting its
consumers until nothing imported it, never by finding it orphaned.

---

## 3. Six findings, four of which correct the plan

### 3.1 `graph_grammar` is 54% re-export

*(The reading that opened PR2. It has since been acted on — §6 measures
what is left.)*

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
all eight bindings into `graph_directed_view` unchanged, or it will
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
- **Successors.** `graph_directed_view` ← `graph` (+ the frame, until
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

    PR2  graph_grammar                                    [in progress]
         |-- DONE:  redirect the 99 re-export imports
         |          (graph_field, set_graph, GRAPH_FIELD_*)
         |          -> graph_field_calculus / fractal_graph
         |-- DONE:  measure what is left, and classify it   (6.1, 6.2)
         |-- DONE:  graph_directed_view, 59 imports, frame carried
         |          UNCHANGED and marked                    (6.4)
         `-- DONE:  graph_operation_view, holding BOTH verbs. The
                    measurement said one commit, and the ruling said
                    one MODULE: a transform is a graph-to-graph
                    operation contract, and two modules would have been
                    ceremony. graph_grammar reached zero and was
                    DELETED.                                    (6.5)

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

## 6. PR2 measurement, after the re-exports were redirected

Taken after the first draught (`graph_grammar` lends no name it did not
define). Every remaining symbol is classified into the five categories
the phase named.

### 6.1 What is left, and where each symbol goes

| symbol | files | src | test | successor |
|---|---|---|---|---|
| `graph` | 57 | 21 | 36 | `graph_directed_view`, carrying the frame marked |
| `graph_operation` | 25 | 9 | 16 | `graph_operation_view` |
| `graph_transform` | 1 | 1 | 0 | `graph_transform_view` |
| **true legacy forwarding** | **0** | **0** | **0** | — nothing left to forward |

**The fifth category is empty, and that is the result.** Before the
redirect, `graph_grammar` was 54% a forwarding shell. It is now 0%: all
three remaining symbols are abstract types the module itself defines. It
is no longer a shell to drain but three contracts awaiting three homes,
and each can leave in its own commit.

58 files still import it, down from 68. The ten that left named nothing
but re-exports.

*(§6.4 records what happened when the first of the three left.)*

### 6.2 The frame's real reach is eight files, not fifty-eight

Every call site of the 36 deferred bindings was measured across all 58
consumers. The eight frame relations — the bindings PR1 §3.3 found are
not a view question — are named in **8 files**, and only three of them
in `src/`:

    src/class_graph.f90              all eight        IMPLEMENTS them
    src/class_graph_assembler.f90    six             num_parts, has_part_relation,
    |                                                global_*_index, *_owner_part
    src/class_graph_partitioner.f90  two             global_vertex_index,
    |                                                global_edge_index
    test/graph-contract               five
    test/graph-partition              one
    partitioned-pde L4 / L5 / L6      six / one / one

    50 of 58 consumers touch no frame binding at all.

The six *frame sets* (`owned_`/`borrowed_`/`overlap_` vertices and edges)
are named in four files. They are carve-and-bind, which is an allowed
view capability, so they are **not** part of the PR3 move; only the eight
relations are.

This narrows PR3 sharply. Rehoming the frame is a change to one
implementer and **two** consuming algorithms, not to the fifty-eight
files that speak the ordinary graph.

### 6.3 The one clean seam

`test/visualization-tower/common/production_discretization_fixture.f90`
is the only consumer that imports `graph_operation` without also naming
`graph`. Every other operation consumer needs both, so
`graph_operation_view` cannot be split from `graph_directed_view` by
import surgery alone — the `apply` interface is written in `class(graph)`
and will import it from wherever the ordinary view lands.

`graph_transform`'s single consumer (`graph_calculus`) remains the
smallest boundary in the cutover.

### 6.4 After the ordinary view left

`graph_directed_view` now holds the abstract `graph`: 36 deferred
bindings and the 13 interfaces they name, moved verbatim. The eight
frame relations are carried **unchanged** and marked

    LEGACY PARTITION FRAME - MOVES IN PR3

with the successor for each written beside it. The six frame *sets* are
not marked, because §6.2 measured them as carve-and-bind, which a view
may do. `class_graph`, which implements the contract, changed its import
and nothing else — no rename, no frame move, both out of scope here.

What is left of `graph_grammar`:

    symbol            files  src  test   successor
    ----------------  -----  ---  ----   ---------------------------
    graph_operation      25    9    16   graph_operation_view
    graph_transform       1    1     0   graph_transform_view

    68 files -> 58 -> 25.  Two contracts, one module, no re-exports.

**The remaining split is a two-line decision, not two phases.** With the
ordinary view gone there is nothing left in `graph_grammar` to
disentangle the two verbs from: they share a file and nothing else. Both
are already written against `graph` imported from the view. Whoever takes
them can move both in one commit and delete the module in the same
breath, or leave them together under a truer name — that is a naming
question, and it no longer has a dependency in it.

The one thing PR2 did **not** resolve, and deliberately: the abstract
type is still called `graph`, which is why `set_graph => graph` is still
written at every door. The module name now says view; the type name does
not. Renaming it reaches every `class(graph)`, `type(graph)`,
`extends(graph)` and `import ::` site — redesign, and PR3's to make.

### 6.5 The grammar is deleted

`graph_operation_view` holds `graph_operation` and `graph_transform`
together, and that is a ruling, not an oversight. A transform *is* a
graph-to-graph operation contract; the two are named in the same files,
written against the same three imported names, and after the ordinary
view left there was nothing to disentangle them from. Two modules would
have been ceremony.

    graph_grammar, over PR2
    ---------------------------------------------------------------
    68 files  ->  58  ->  25  ->  0        and then deleted

    what it held                where it went
    --------------------------  -----------------------------------
    graph_field, GRAPH_FIELD_*  graph_field_calculus  (commit 1)
    set_graph                   fractal_graph         (commit 1)
    graph                       graph_directed_view   (commit 3)
    graph_operation             graph_operation_view  (commit 4)
    graph_transform             graph_operation_view  (commit 4)

Prose is not imported, so nothing forced its survival — but four blocks
of tower law lived only in that file's header and would have died with
it. They went where their subject went:

    WHAT A GRAPH IS MADE OF      -> graph_directed_view
    CAN A GRAPH CHANGE?          -> graph_directed_view
    THE ADMISSION LAW            -> graph_operation_view
    WHAT apply DOES TO A BUFFER  -> graph_operation_view
    THE FOUR ROLES               -> graph_operation_view, restated as
                                    a map of where each role now lives

One block was NOT carried: `GRAPH_FIELD. The carrier of values.`
`graph_field_calculus` already states the same shape invariant and the
same interleaving formula, so carrying it would have made two homes for
one law. That module's own header claimed the grammar "now re-exports
them" — no longer true after commit 1, and false in a stronger sense
now, so it says what is true instead.

**This is the cutover's first deleted module**, and it was earned the
way §1 requires: consumers rewritten until nothing imported it.

---

## 6.6 PR3: the last public non-ontology `graph`

> **Naming amended after PR3 merged.** `ordinary` names no mathematical
> role. The three names were renamed again, in a follow-up:
>
>     graph_ordinary_view    ->  graph_directed_view
>     ordinary_graph         ->  directed_graph
>     ordinary_stored_graph  ->  directed_stored_graph
>
> The hierarchy the vocabulary now states:
>
>     fractal_graph :: graph                 ontology, G=(B1,B2)
>     graph_directed_view :: directed_graph  the view D=(V,E,tail,head)
>     class_graph :: directed_stored_graph   one stored realization
>
> `directed` is what the structure IS — two finite domains and two maps
> tail, head : E -> V — where `ordinary` only said it was not exotic.
> The section below has been rewritten in the new names; the *evidence*
> it reports is unchanged, because a rename does not move an import.
>
> The last public name carrying the word went with them.
> `graph_profile :: ordinary_graph_view` is now
> `directed_incidence_view`, named for the relation it holds:
>
>     directed_incidence_view   T, H <= E x V   incidence, and every
>     |                                         neighbourhood question
>     |                                         as a composition of it
>     directed_adjacency_view   A <= V x V      one stored adjacency
>
> The pair names its own primitive now. The first answers
> D = (V, E, tail, head) read off **relations**; `class_graph` answers
> the same D read off **stored arrays** — two realizations of one
> contract, which is what the vocabulary should have said all along.
>
> `ordinary` survives in `src/` in exactly two comments that explain the
> ban and two that use the English word ("an ordinary linear question").
> The suite directory `test/graph-ordinary/` still carries it; renaming
> it would dangle two historical doc references, so it is flagged here
> rather than changed.

The abstract type `graph_directed_view :: graph` (then `ordinary_graph`) is now
`directed_graph`, and the concrete `class_graph :: stored_graph` is
`directed_stored_graph`.

### What made the rename safe

Three facts, established before a byte moved, each of which cut the
surface:

**The name resolves per program unit, not per file.** Every file was
split at `module`/`program` boundaries and its `graph` traced to the
module that supplies it:

    56 files   `graph` is graph_directed_view's      -> renamed
    62 files   `graph` is fractal_graph's            -> untouched
     2 files   `graph` is interface_graph's          -> untouched (island)
     0 files   two modules supply the name at once

Zero mixed files is the fact that made a file-wide rewrite legal.

**And the trap in that method, which caught this rename.** Keying on the
LOCAL name misses an import that renames at the door. Two files said

    use graph_directed_view, only : grammar_graph => graph

so their local name was `grammar_graph`, not `graph`, and they fell
outside the 56 — untouched by the rewrite, and broken by it, because the
module they import from no longer exported `graph`. `src` built clean
because both files are tests; two suites failed and named the line. The
alias existed only to dodge the collision the rename removes, so both
were collapsed to `only : directed_graph`.

The lesson is the method's, not the rename's: **classify by the symbol
imported, not by the name it is bound to locally.** The `stored_graph`
rename was checked for aliases before it ran, and had none.

**The ordinary type is abstract, so `type(graph)` can never name it.**
Of 370 `type(graph)` sites, **none** is in the 56 — they are all the
kernel, and a hit inside the 56 would have been a compile error today,
not work. The rename surface is `class(`, `extends(` and `import ::`
alone.

**Only one `extends(graph)` in `src/` is the ordinary one.** The four
sites split cleanly:

    class_graph          extends the ordinary abstract   -> renamed
    class_mesh           extends interface_graph's       -> untouched
    class_stored_graph   extends interface_graph's       -> untouched
    interface_graph      digraph extends its own         -> untouched

This is the one edit the compiler might not have caught, because both
names exist; it was enumerated rather than pattern-matched.

### The `stored_graph` collision is resolved

43 files import `stored_graph` from `class_graph`, one
(`class_mesh`) from `class_stored_graph`, and **no file imports both**.
The 43 were renamed — 125 declaration sites, 75 constructor calls, and
their imports. `class_stored_graph` still exports `stored_graph`, and
is now the only module that does.

The *module* rename was deferred: `class_graph` keeps its name. The
collision was between two public **types**, and that is what step 2
existed to remove.

### The goal, and the one thing that does not reach it

> after PR3, the only public type named `graph` should be the ontology

Two remain:

    fractal_graph  :: graph      the ontology - the intended survivor
    interface_graph:: graph      the mesh-loader island

The second is excluded by the ruling that the island is not to be
touched, and the exclusion is safe rather than merely permitted: no
file in the repository imports `graph` from two sources, so the
ambiguity is **dormant, not active**. It dies with the island, in the
mesh-measurement work — which is where the ruling put it.

### Not done, deliberately

The partition frame did not move. Its marker now reads

    LEGACY PARTITION FRAME - MOVES IN THE NEXT FRAME PR

as the ruling's correction instructs. Moving it stays scoped to
`class_graph`, `class_graph_partitioner`, `class_graph_assembler` and
the tests that assert frame laws — §6.2's measurement, unchanged.

---

## 6.7 PR4: the partition frame leaves the directed view

The eight bindings PR1 §3.3 measured as *not a view question* are gone
from `directed_graph`. They are `graph_partition_frame_representation`
now — a value a graph carries and an action may bind.

    D = (V, E, tail, head)      is what directed_graph says, and all
                                it says. 28 bindings, no frame.

### What the frame carries, and why it is more than coordinates

    part V, part E              identity  ]  so the frame can say
    whole V, whole E            identity  ]  WHICH graph it is about
    part / whole counts                   ]
    number, nparts, cut         provenance
    vglobal, eglobal            part-local -> whole
    vowner, eowner              who answers for each member

A bare numbering could not say which part it numbers, so a caller
holding two frames could hand the wrong one to the wrong graph and be
wrong in silence. `describes(g)` is the question that makes a captured
frame checkable, and it is what `defined_on_graph` now asks.

It holds **no** `set_map`, `label_map` or `inclusion_map`. What a set
means stays the caller's and arrives at the semantic boundary.

### How the two consumers get it

| | before | after |
|---|---|---|
| partitioner | stamped the frame into the piece | `partition_graph(g, part, frame)` hands it back; `partition_data` takes it |
| assembler | asked the piece | **binds** one at construction: `assembler(frame)` |

The assembler binding is admissible for the reason a CSR relation
copies its numbering: a representation carries no semantic identity and
cannot go stale under its holder. It still holds no map.

### One consequence the ruling did not anticipate

`transform_on_graph_interface` is **no longer `pure`**. Asking whether a
frame describes a graph means comparing set identities, and a set graph
carries a pointer component, so copying one out of an `intent(in)` dummy
is barred from a pure subprogram (F2018 C1594). The *signature* is
unchanged — the contract was not widened — and an impure interface still
permits a pure override, so coarsener, refiner and partitioner remain
pure.

### What the tests learned, which is the point

The migration's own failures taught the law better than a passing test
would have. A host-scope frame reused across two parts is *wrong*, and
five suites said so numerically:

    a check handed `part` must read the frame THAT part carries

which is exactly why the record is a value the graph carries rather
than a global. Every such site now reads `part % frame()` inside its
`select type`, and the shared partitioned fixture holds **one frame per
part** where it used to keep one.

### Not moved, deliberately

The six carved frame sets — `owned_`, `borrowed_`, `overlap_` vertices
and edges — stayed in the view. §6.2 measured them as carve-and-bind,
which a view may do. They now read ownership through the frame the graph
carries, but they are still the graph's own question.

---

## 7. Verification record

One clean rebuild per suite, at every commit recorded here:

    inventory (PR1)          32 of 32 suites PASS, 0 FAIL
    re-export redirect       32 of 32 suites PASS, 0 FAIL
    graph_directed_view      32 of 32 suites PASS, 0 FAIL   (as graph_ordinary_view)
    graph_grammar deleted    32 of 32 suites PASS, 0 FAIL
    directed_graph rename    32 of 32 suites PASS, 0 FAIL   (as ordinary_graph)
    ordinary -> directed     32 of 32 suites PASS, 0 FAIL
    partition frame moved    32 of 32 suites PASS, 0 FAIL

Seven tower import gates were re-asserted, not relaxed: a level that read
the field through the grammar is granted `graph_field_calculus` by name,
and adjoint level 9, which no longer imports `graph_grammar` at all, had
its grant withdrawn.

`src/fractal_graph.f90` is byte-identical to master.
