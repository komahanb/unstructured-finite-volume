# Plan: one directed view, relation-seated

Status: PLANNED, not started. This document is deleted or absorbed
into ARCHITECTURE.md by the commit that completes phase 5.

## The duplication

The directed view D = (V, E, tail, head) is implemented twice:

    directed_stored_graph   [class_graph, ~1,230 lines]
        source of truth: flat arrays tail(:), head(:) plus derived
        CSR snapshots (xinc/einc, xadj/vadj, in/out lists), built
        through group_by_key. Carries the tag, ownership, and
        partition machinery. THE workhorse: every operation
        executes over it.

    directed_incidence_view, directed_adjacency_view
                            [graph_profile, ~580 lines]
        source of truth: two relations T, H <= E x V held as
        csr_relation; every question answered as a reading of
        their fibres. Consumers: graph_algorithms (sources, sinks,
        reachable, topological_order) and 17 test files across the
        five tower suites.

Both answer the same questions. Neither extends the other; the
profile views do not even extend the abstract `directed_graph`, so
they cannot be handed to an operation. AGENTS.md section 16 rules
the direction ("traditional graph objects should be expressed as
views over relational graphs") and section 66 forbids deleting old
storage before the hot paths are benchmarked.

## Facts the plan stands on (measured 2026-08-19)

1. All operation/solver consumers depend on the abstract
   `directed_graph` contract (28 deferred methods), not on either
   concrete type. The contract can stay frozen; only internals and
   the profile's consumers move.
2. `csr_relation % image_view / preimage_view` return POINTER
   slices into the compiled xfwd/tgt arrays - no allocation per
   query. The relation-backed reading is near array speed except
   for one `local_index` call per query.
3. The hot per-edge queries (edge_tail, edge_head, edge_has_head)
   run inside stencil_apply (every matvec of every solver
   iteration), the balance, and the differential operator's
   compile. Neighbourhood queries (incident/adjacent/in/out) run
   at setup time (walk, partitioner, coarsener, assembler).
4. The kernel ontology (AGENTS.md graph-ontology section) blesses
   compiled snapshots: "a CSR compiled from a graph stays valid
   when that graph changes". Keeping flat tail/head arrays as the
   compiled snapshot of T and H is lawful, so consolidation does
   not have to slow anything.
5. `test/graph-benchmark` exists and already times differential
   operators; it is the seed of the section-66 harness.

## The target

One concretion. `directed_stored_graph` keeps its name, its
constructor, and the frozen 28-method contract, but its source of
truth becomes the two relations:

    directed_stored_graph = (identity, V, E, T, H)
        T, H <= E x V         stored as csr_relation - the truth
        tail(:), head(:)      the compiled snapshot of T and H:
                              O(1) reads for the hot five
        xinc/einc, ...        derived snapshots, as today, built
                              by group_by_key from the tables
        tags, ownership,      unchanged
        partition seams

    graph_profile             DELETED. Its fibre-reading mechanics
                              move into the concretion; its two
                              view types retire.

    graph_algorithms          re-aimed at class(directed_graph):
                              sources/sinks/reachable/topological
                              need only in/out queries, which the
                              contract already declares.

    new constructor           directed_stored_graph from (T, H)
                              relations - the door the tower tests
                              that today build profile views walk
                              through instead.

This satisfies section 16 - the directed vocabulary is derived from
first-class relations - with one type instead of two, and satisfies
the snapshot doctrine so section 66 can pass.

## Phases

Each phase is one commit on master with the full battery green;
phases 3 and 5 additionally pass the benchmark gate.

    Phase 0   DONE (174ede7). The measuring stick.
              Extend test/graph-benchmark to time, on square-80 and
              box-3 scales: incident/adjacent/in/out queries,
              stencil matvec, differential operator apply, field
              assembly, partition construction. Commit the baseline
              numbers as a reference file the suite prints beside
              its results. Gate for later phases: no hot path slower
              than 1.10x baseline, else stop and write the
              optimization plan section 66 demands.

    Phase 1   DONE. The equivalence law.
              test/graph-ordinary's compare_topology engine asks
              the nine view-answerable questions of both
              constructions; phase 1 added a 400-vertex
              pseudo-random topology with walls, tuples shuffled.
              The set, tag, and ownership battery exists only on
              the stored graph; its characterization through the
              re-seat is the full 35-suite battery itself. Green
              BEFORE internals move and after every phase.

    Phase 2   DONE, by verification - no code change. Findings:
              counted_local_index is O(1) (range check plus
              identity); the benchmark's 23 ns fibre sweeps measure
              it in place; image_view/preimage_view cover the
              re-seat's needs. Requirement pinned for phase 3: the
              graph's T, H relations are built over COUNTED
              coordinate sets only - listed_local_index scans.
              The relation's borrow face.
              Verify csr_relation exposes everything the snapshot
              rebuild and the fibre reads need without copying:
              image_view/preimage_view (exist), whole-table reads,
              and O(1) local_index for counted coordinate sets.
              Measure local_index; if any coordinate path is O(n),
              fix it here, before anything depends on it.

    Phase 3   DONE - amended by the benchmark gate. The eager
              form (T, H built as components in every constructor)
              failed section 66: every pattern graph and part graph
              paid the build, costing stencil_of 1.66x, partition
              1.77x, assembly 1.65x. The landed form derives the
              relations ON REQUEST: tail_relation() and
              head_relation() build csr_relations from the stored
              table over counted coordinates when asked; a graph
              nobody reads relationally pays nothing, and every
              benchmark act is within noise of baseline. The
              equivalence law gained two accessor laws: T and H
              fibre every edge to its own ends, and their preimages
              are exactly the outgoing and incoming lists. The
              stored snapshots remain the compiled truth the hot
              paths read.

    Phase 4   One consumer contract.
              graph_algorithms takes class(directed_graph). The
              five tower level-4 suites and graph-ordinary /
              graph-algorithms tests construct the stored graph
              (directly or through the new from-relations
              constructor) instead of profile views. graph_profile
              gains no new consumers from here on.

    Phase 5   Retire graph_profile.
              Anything it alone still answers moves into
              class_graph or graph_algorithms; the module and its
              two types are deleted; a static (like the mesh one)
              keeps them deleted; ARCHITECTURE.md replaces the
              "known dual representation" section with "one
              directed view, relation-seated" in the same commit.
              Gates: equivalence law, battery, benchmark.

## Risks, named

    performance      the snapshot design should make phase 3 a
                     no-op on hot paths; the benchmark gate is the
                     proof, not the intention. If local_index or
                     fibre reads leak into a matvec, the gate
                     catches it at phase 3, not at phase 5.
    identity         vertex_set()/edge_set() must keep minting the
                     same identities across the re-seat; the
                     equivalence law checks same_as against maps
                     built before and after.
    purity           csr_relation faces are impure where the old
                     builders were pure; the earlier migration
                     ruled PURE loss acceptable, but each loss is
                     named in the phase-3 commit message.
    test surface     17 test files import graph_profile; phase 4
                     touches them mechanically but each tower suite
                     is its own oracle and must stay green
                     unmodified in what it asserts.
    memory shape     holding T, H, and snapshots stores the edge
                     list roughly twice more than today; acceptable
                     at current scales (largest suite mesh: 22k
                     faces), stated here so it is a decision and
                     not an accident.

## Order of execution

0 -> 1 -> 2 -> 3 -> 4 -> 5, no reordering: the measuring stick and
the equivalence law exist before anything moves (the mesh deletion
proved this discipline), the borrow face exists before the re-seat
depends on it, and the profile dies only after it has no consumer.
Phases 0-2 are small; phase 3 is the surgery; phases 4-5 are wide
but mechanical.
