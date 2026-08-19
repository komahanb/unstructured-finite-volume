# The remaining migration to the fractal ontology

A re-measurement of `doc/final-codebase-cutover-plan.md` §4 against master
at `b411ce9`. The plan's ledger is authoritative on **roles and blockers**
and stale on **status**: it was written at PR2, and PRs #4–#10 have since
merged. This report carries only the delta — for each ledger entry, whether
the stated blocker is **discharged**, **standing**, or **replaced**.

## 1. Method

Every claim below is a grep or a line count against `src/` and `test/` at
`b411ce9`. Import censuses classify by the symbol **imported**, not by the
local name it is bound to — the rename method learned in PR #7, where
`only : alias => graph` made two files invisible to a local-name search.

## 2. The ledger, re-measured

    module                  LOC    blocker stated (PR2)      blocker today
    ---------------------------------------------------------------------
    graph_grammar           661    none structural           DELETED 72bbaec
    graph_calculus          612    graph_grammar             DISCHARGED
    class_graph            1250    the partition frame       DISCHARGED
    class_graph_mesh        231    base class + supplier     HALF-DISCHARGED
    interface_graph        2650    class_mesh                STANDING
    class_stored_graph      269    class_mesh                STANDING
    class_mesh             1870    measurement not ported    STANDING (real work)
    class_mesh_builder      116    same as class_mesh        STANDING

Two entries changed verdict since the plan was written. Both are now cuts
that can be taken, and neither was schedulable at PR2.

## 3. Findings

### 3.1 The frame blocker is discharged; `class_graph` is unblocked

The ledger required the partition frame to leave the type before
`class_graph` could become a representation. It has left. The type today:

    type, extends(directed_graph) :: directed_stored_graph
       integer              :: number, nv, ne
       integer, allocatable :: tail(:), head(:)
       integer, allocatable :: xinc(:), einc(:), xadj(:), vadj(:)
       integer, allocatable :: xout(:), eout(:), xin(:) , ein(:)
       character(len=:), allocatable :: vtag(:), etag(:)
       type(partition_relation) :: whole_rel
       type(set_graph)          :: vset, eset

`vglobal`, `vowner`, `eglobal`, `eowner` and `nparts` are **constructor
dummies only** — `grep "this % vglobal"` and its four siblings return
nothing. They are consumed at lines 306 and 320 to build `whole_rel`, and
never stored. PR #9/#10 discharged this: the frame became a value, and the
value is a relation.

What remains is an edge list, four derived adjacency tables, two tag
arrays, two set identities — and `whole_rel`, a `type(partition_relation)`
component handed out **by value** at line 1245.

The measurement stops there. That the five frame integers are gone is
proved; that a *relation-valued component* clears the representation bar of
global rule 6 and PR4 is **not** proved, and may be exactly what the
ledger's blocker meant to exclude. See §6.

### 3.2 `subset_set`, `member_set` and `graph_carrier` are historical

The set-re-rooting pass recorded a blocked deletion: `graph_forms.f90:52`
extending `subset_set`, 19 src constructions, 105 declarations across 29
modules. **None of it survives.** `grep -rn "subset_set\|member_set\|
graph_carrier" src/` returns 12 hits; every one is either prose in a
comment or the pair `num_member_sets` / `member_set_at`, which are the
relational view's accessors for the S of (S,P) — a different object that
keeps its name correctly.

Test-side hits are 41 files, all of them `README.md`, `NUCLEUS-
OBSERVATIONS.md`, or `check_imports.sh` lines that name the retired module
**in order to refuse it**. A gate naming a dead module is the gate working.

This entry must not appear on any future work list.

### 3.3 `graph_calculus` is the next unblocked atomic cut

Its documented blocker was `graph_grammar`, deleted at 72bbaec. Its nine
abstracts are refinements of three parents, and all three parents already
have view homes:

    graph_functional              <- graph_field         graph_field_calculus
    graph_reduction               <- graph_operation     graph_operation_view
    graph_broadcast               <- graph_operation     graph_operation_view
    discretization_operator       <- graph_operation     graph_operation_view
    linearization_operator        <- graph_operation     graph_operation_view
    graph_partitioner             <- graph_transform     graph_operation_view
    graph_assembler               <- graph_transform     graph_operation_view
    graph_coarsener               <- graph_transform     graph_operation_view
    graph_refiner                 <- graph_transform     graph_operation_view
    GRAPH_SIDE_VERTEX / _EDGE      (absorbed axis)       graph_directed_view

Nothing structural stands in the way. The cut is 13 src modules and 9 test
files, and it is the largest remaining piece of pure rewiring.

`GRAPH_SIDE_VERTEX`/`GRAPH_SIDE_EDGE` reach 5 src and 9 test files. The
ledger sends them "wherever the ordinary view lands"; that view landed as
`graph_directed_view`, which today publishes `directed_graph` alone.

### 3.4 The island is joined to the live tree at one file and one alias

The old-world graph stack is closed under three modules:

    interface_graph  (2650)  <- iso_fortran_env, module_solve_mode
      ^
    class_stored_graph (269)
      ^
    class_mesh (1870)        <- class_array_mesh_loader, class_string,
                                interface_mesh_loader, module_mesh_utils,
                                module_verbosity
      ^
    class_mesh_builder (116) <- ALSO class_graph_mesh  (the live tree)

`class_mesh_builder` is the only file that touches both worlds, and it does
so by importing the same public name from two sources:

    use class_graph_mesh, only : mesh
    use class_mesh      , only : legacy => mesh

This is the `mesh` counterpart of the `stored_graph` collision resolved in
PR #7/#8, with one difference: the rulings record that no file imports
`graph` from two sources, so that ambiguity is **dormant**. This one is
**live**, in one file, and that file is the bridge. `mesh_from_gmsh` is the
whole of the joint.

The loaders are **not** island. `interface_mesh_loader` depends on
`class_string` and `iso_fortran_env` only; `class_array_mesh_loader`,
`class_gmsh_loader`, `module_mesh_utils` and `class_file` never mention the
graph stack. They are I/O and geometry utilities and they survive the cut
unchanged.

**A suspected duplication is not real.** `class_gmsh_loader` carries 22
gmsh-format references; `class_mesh` carries 1. The parsing is not
duplicated — `class_mesh` delegates it. What `class_mesh` uniquely holds is
the measurement machinery: 10 procedures, ~445 lines.

    evaluate_vertex_weight            evaluate_cell_volumes
    evaluate_face_weight              evaluate_face_deltas
    evaluate_centroidal_vector        evaluate_cell_centers
    evaluate_face_centers_areas       evaluate_face_centers_areas_2d
    evaluate_face_tangents_normals    evaluate_face_tangents_normals_2d

This is the one blocker in the whole ledger that is work rather than
rewiring, and the binding ruling fixes its order: extract
`graph_mesh_measurement` / `mesh_geometry_representation` **first**; only
then may `mesh_from_gmsh` retire. The order does not invert.

Two notes on this entry. The ledger records `class_mesh_builder`'s
successor as "Nothing — the call site `mesh_from_gmsh` is designed not to
change"; the later binding ruling permits that call site to retire once the
extraction lands. **The ruling supersedes the ledger** — a reader following
the stale entry would preserve a bridge that is meant to go.

And the cut has a one-grep completion test, the same shape as the
emptiness gate that closed `graph_grammar`: the island is severed when no
file imports one public name from two sources. Today exactly one does.

### 3.5 Three rival CSR representations, as predicted

    graph_binary_relation :: csr_relation            (the designated donor)
    class_graph           :: directed_stored_graph   (xadj/vadj, xout/eout,
                                                      xin/ein, xinc/einc)
    interface_graph       :: graph                   (xadj/adj)

The third dies with the island. The first two must be reconciled, and that
reconciliation is priced, not free — see §4.

## 4. Two classes of work, and only one of them is free

The kernel spike measured semantic traversal at ~43 ns/hop against
production CSR at 9 ns/query, and concluded that fractal compression is
**ontological, not operational**. Every recommendation therefore falls into
one of two classes, and conflating them is what §66 exists to prevent.

**Class R — rewire and rename. No hot path is touched; cost is import
churn and one green suite run.**

    R1  graph_calculus       nine abstracts to their parents' view modules
                             GRAPH_SIDE_* to graph_directed_view
                             13 src + 9 test files
                             ORDER: the constants move first or together —
                             graph_forms and class_graph_reduction import
                             both from graph_calculus today, and a split
                             costs them two rewrites
    R2  class_graph          rename to a representation name
                             blocker discharged (3.1)
    R3  class_* prefix       37 modules; no public class_* name survives
                             mechanical, but must follow R1/R2
    R4  ordinary_graph_view  one inhabitant of a banned word (§6)

**Class P — re-root through the ontology. Priced, and each needs a storage
tier before it may proceed.**

    P1  CSR reconciliation   two survivors to one; csr_relation is the
                             donor; benchmark gate is mandatory
    P2  mesh measurement     10 procedures, ~445 lines, no successor
                             written; this is design, not migration
    P3  the island cut       4789 LOC retired, gated on P2

The island's 4789 lines are the largest number in the report and the
smallest amount of thinking: once P2 lands, P3 is deletion.

## 5. Cut order

    1.  R1  graph_calculus drained            unblocked NOW
    2.  R2  class_graph -> representation     unblocked NOW
    3.  R4  ordinary_graph_view ruling        a decision, not work
    4.  P2  mesh measurement extracted        real work; ruling fixes order
    5.  R3  class_* prefix retired            after 1 and 2
    6.  P3  island deleted                    after 4
    7.  P1  CSR reconciliation                after 6, with the benchmark

One deletion per commit. The seven towers stay unchanged oracles.

## 6. Decisions to surface, not schedule

- **`graph_profile :: ordinary_graph_view`** (581 lines). `ordinary` is
  banned from public names as of 2026-08-16 — it denotes no mathematical
  role. This is the last inhabitant, and it is a different object from the
  view that was renamed in PR #8: it holds the relational schema
  T, H ⊆ E×V. It awaits a ruling, not a commit.
- **Does a relation-valued component clear the representation bar?**
  `directed_stored_graph` holds `type(partition_relation) :: whole_rel`.
  The frame integers it replaced are gone (§3.1), but global rule 6 and
  PR4 forbid representations that carry maps. Whether a *relation* is a
  map for this purpose decides whether R2 is a rename or a fourth cut.
- **`class_graph_mesh` is half-discharged.** Its base-class blocker clears
  with R2; its supplier blocker clears with P2. It is the only entry
  waiting on both tracks.

## 7. Debt column

- **`similar_as`** (`fractal_graph.f90`, b411ce9) has **zero consumers** in
  `src/` and `test/`. It is kernel surface with no caller. Either a client
  earns it — adaptive refinement and coarsening termination are the two
  named candidates — or it comes back out.
- Legacy operation-host graphs: `class_fitted_balance`'s constellation, the
  minimizer host argument, the robustness and contract fixtures.
- Twenty `test/` directories hold nothing but stale build artifacts
  (`run`, `*.o`, `refusal`) from retired suites — `graph-structure`,
  `fractal-kernel`, `mesh` and seventeen others. All untracked, none with
  a Makefile or a source file. Free to remove.
- Phase 11 legacy consolidation, unstarted.
- Rectangular minimization family and natural join, both unearned.

## 8. Verification

**Baseline observed at `b411ce9`: 32 suites run, 32 green, 0 red.**
`test/*/run.sh` is the complete live suite — the 32 directories carrying a
`run.sh` are every suite there is.

    test/*/run.sh                    32/32 green at b411ce9
    test/graph-benchmark/run         against doc/phase0-benchmark.md

Two corrections to earlier records. `test/graph-benchmark` is **no longer
baseline-red** — 720537f recorded it red and 8d7eabc turned it green;
`doc/graph-benchmark-baseline-red.md` is stale. And a refusal binary is
**not a standalone test**: `test/fractal-graph/refusal` is a case
dispatcher that its `run.sh` drives as `./refusal "$case"`, checking each
error stop against an expected message. Invoked bare it answers
`refusal: no such case` and exits nonzero — correct behaviour, not a
failure. Any ad-hoc runner that calls these binaries directly will
misreport the suite.


Seven vertical towers stand as acceptance oracles, six of them carrying a
checked result marker:

    calculator-tower                 CALCULATOR_RESULT
    learning-tower                   LEARNING_RESULT
    derivative-action-tower          DERIVATIVE_RESULT
    adjoint-tower                    ADJOINT_RESULT
    partitioned-implicit-pde-tower   PARTITIONED_PDE_RESULT
    time-integration-tower           TIME_INTEGRATION_RESULT
    visualization-tower              (no result marker)

Each was built to require zero production change at every rung. Any cut
that moves a tower has changed the mathematics, not the packaging.
