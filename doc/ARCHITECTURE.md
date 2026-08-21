# The prime structure

This is the map of `src/`, generated from the code as it stands.
One line per idea; the constitution behind it is `AGENTS.md`. The
generating rule of the whole codebase is:

> A small set of prime types generates every higher type by
> composition. A type is prime when it cannot be expressed as a
> composition of the others. Everything that is not prime must say
> which primes it composes.

## The primes

Eight irreducible objects. Everything else in `src/` is a
composition, a view, or a concretion of these.

    graph            G = (B1, B2),  B in {NULL, UNKNOWN, KNOWN -> G}     graph_fractal
    token            object identity, minted once by declare             token_identity
    relation         P subset of A1 x ... x Ak                           relation_finitary
    field            f : A -> V on one set-graph domain                  field_calculus
    operation        (graph, fields) -> field                            operation_action
    transform        graph -> graph, data carried along                  transform_structure
    map              token -> {representation | name | inclusion | value}
                                                                         map_set, map_label,
                                                                         map_inclusion, map_value
    view             a reading of the kernel graph, never a new kind     view_directed,
                                                                         view_relational

The kernel graph is self-similar: each branch is NULL, UNKNOWN, or
another graph. Sets, tuples, relations, domains, and chains are all
this one object read differently. There is exactly one type named
`graph` in the living tower.

The namespace law (doc/coding-standards.md, 2026-08-19): a module
reads namespace order, `<prime>_<subject>` (`relation_binary`,
`field_stored`); a type reads english order (`binary_relation`,
`stored_field`); so no type can ever collide with a module name.
test/graph-ordinary enforces it statically.

## The composition tree

Indentation is `extends`; brackets name the file.

    relation                                   [relation_finitary]
    ├── stored_relation                        [relation_finitary]
    └── binary_relation                        [relation_binary]
        ├── csr_relation                       [relation_binary]
        └── transposed_relation  (the converse)  [relation_binary]

    group_by_key                               [relation_binary]
        the one grouping kernel (counting sort): the fibration of a
        stored relation over one slot. CSR builds, incidence lists,
        transpose_padded, and combine_triples are its callers; no
        other count/prefix/scatter exists.

    directed_graph  D = (V, E, tail, head)     [view_directed]
    └── stored_directed_graph                  [view_directed_stored]

    field                                      [field_calculus]
    ├── stored_field                           [field_stored]
    └── functional  (one entry; the scalar adapters, written once over the vector adapters)  [field_calculus]
        └── stored_functional                  [field_functional]

    operation                                  [operation_action]
    │     name, domain, apply; an argument space declared once by the
    │     constructor (declare_arguments: the readable input positions,
    │     not a required length); argument(k) names a position in it,
    │     and arguments of two operations never match. max_degree (0
    │     unless overridden) with partial_action over a list of
    │     variations - one (argument, direction) factor each - so
    │     differentiability is a capability of every operation;
    │     tangent_of picks the exact tangent when max_degree >= 1,
    │     the difference tangent otherwise
    ├── reduction   field -> functional  (sum, average, norm, ...)   [operation_reduction]
    ├── broadcast   functional -> field  (copy, share)               [operation_reduction]
    ├── discretization                             [operation_discretization]
    │   ├── stencil  space: matrix as weighted edges         [operation_stencil]
    │   │       constructed from triples, from a dense array, or
    │   │       compiled from an operation on the standard basis;
    │   │       transpose; its own tangent (a linear map, max_degree
    │   │       1); exports combine_triples (one entry per pair)
    │   └── scheme   time: a0 q_n + a1 q_(n-1) + a2 q_(n-2)        [operation_step]
    │           + h [theta S(q_n, xi) + (1-theta) S(q_(n-1), xi)]
    │           (theta 1 backward euler / BDF2, 0 forward euler,
    │           1/2 Crank-Nicolson); its argument space is the
    │           action's (state, auxiliaries) followed by history(k),
    │           k <= reach; every partial is derived from the action's
    │           by that formula, each factor restated on the action's
    │           argument first - the tangent of the step equation is
    │           exact when S's is. The history arrives as inputs
    │           m+1..m+reach, shape-checked against the state; the
    │           scheme stores no state of its own
    ├── linearization  D_a S(x)[v] at a frozen input tuple   [operation_linearization]
    │       the tangent of any statement in any of its arguments,
    │       derived by tangent_of(statement, wrt, at_inputs): exact
    │       through partial_action when max_degree >= 1, by differences
    │       otherwise - the primal is written once; dual_by_basis forms
    │       (D_a S)^T lambda under the Euclidean pairing on stored
    │       values, for rectangular blocks the square stencil cannot
    ├── minimizer                              [operation_minimization]
    │   │     attach holds the statement, the unknown domain and any
    │   │     held_inputs - inputs fixed during this solve, applied
    │   │     after the unknown, the affine part's included
    │   ├── jacobi, gauss_seidel, conjugate_gradient, gmres, multigrid
    │   ├── newton        Jacobian action = tangent_of(action)
    │   └── dense_direct  elimination on the attached operation's basis columns
    ├── differential_operator  d^n on graph sides, compiled     [operation_differential]
    │       the average, difference, and incidence steps as affine
    │       sparse maps, composed by parity into one matrix + constant;
    │       adjoint = transpose; stencil_of returns the square compiled
    │       form as a stencil
    ├── balance, fit, walk                     [clients and utilities]
    └── ...

    transform                            [transform_structure]
    ├── partitioner                            [transform_partitioner]
    ├── assembler                              [transform_assembler]
    ├── coarsener                              [transform_coarsener]
    └── refiner                                [transform_refiner]

    reversible_change  run_change: apply -> check -> keep | revert   [map_change_protocol]
    └── value_change                           [map_value_change]

    chain_rule + argument_path + path_derivative          [operation_chain_rule]
        total derivatives of S(x_1(s), ..., x_m(s)) by integer
        partitions of the degree, from partial actions alone

    marcher  time as a chain graph             [operation_marching]
        every rule is one scheme (forward euler theta 0, backward
        euler theta 1, BDF2 reach 2) and every verb one traversal:
        march and march_adaptive step by the residual at the zero
        state over -a0 when theta = 0, else by the held minimizer
        with the rest of the tuple held; march_adjoint traverses the
        converse chain, solves the transposed state block and
        subtracts the dual of each history block from the seed of the
        instant it reads; march_directional solves the state block
        against the chain rule over state, history and parameter
        paths. Nothing names a1 or a2. Every auxiliary argument of
        the action is supplied as a parameter, or the march is refused
    step_policy -> halving_policy              [operation_step_policy]

    form -> polynomial_form, harmonic_form     [field_forms]
    form_optimizer -> pruner;  fit (an operation)   [operation_fitting]

## Composition, stated as arithmetic

    field        = graph  x  values
    relation     = graph  x  tuples
    operation    = graph  x  fields  ->  field
    stencil      = operation with weights on edges          (a matrix)
    scheme       = operation x (a0, a1, a2, h, theta)       (a step rule)
    variation    = argument x direction                     (one derivative factor)
    tangent      = operation x argument, frozen at inputs   (a matrix action)
    derivative   = incidence o difference o ... o average   (one stencil)
    chain rule   = partitions x partial actions             (any derivative)
    minimizer    = operation + a stopping rule              (a solve)
    newton       = minimizer x tangent_of                   (nonlinear solve)
    marcher      = chain graph x scheme x minimizer         (integration)
    adjoint      = marcher read backward x transpose        (sensitivity)
    change       = mutation x memory x undo                 (transaction)
    mesh         = sets x relations x measurement fields    (see below)

## Level assignment (AGENTS.md section 2)

    0-3  kernel:      graph_fractal, token_identity, relation_finitary,
                      relation_binary, relation_partition,
                      map_set_representation, the four maps
    4    graph view:  view_directed, view_relational,
                      view_directed_stored, operation_walk, relation_algorithms,
                      transforms (partitioner/assembler/coarsener/refiner)
    5    field:       field_calculus, field_stored,
                      field_functional, operation_reduction
    6    calculus:    operation_discretization, operation_stencil,
                      operation_step, operation_step_policy,
                      operation_linearization, operation_chain_rule,
                      operation_differential, field_forms,
                      operation_fitting
    7    solve:       operation_minimization + members, marcher,
                      step_policy, change protocol + value map
    8-9  physics:     operation_balance, operation_fitted_balance,
                      operation_advection, operation_conduction,
                      operation_diffusion, operation_robin_condition,
                      view_mesh, view_mesh_geometry + the loaders

An import may only point downward in this list.

## The mesh, in the primes

A mesh file supplies member sets (vertices, cells, tagged boundary
faces), one relation (C2V, cell to vertex), one field (coordinates
on vertices), and tag names. Everything else is derived, in
`view_mesh_geometry`, by relation algebra and geometry:

    faces  = boundary faces + shared-vertex intersections of cell pairs
    V2C    = transpose(C2V)
    F2C    = {(f, c) : F2V(f) subset of C2V(c)}     (<= 2 cells per face)
    C2F    = transpose(F2C)
    measurements = pure functions of coordinates and the relations:
                   cell centres, face centres, one area vector S_f per
                   face (the scalar area is its norm, the unit normal
                   its direction, signed out of each cell), volumes by
                   the divergence theorem from the signed S_f,
                   centre-to-centre vectors, deltas |l_f.S_f|/|S_f|,
                   weights

`view_mesh_builder` seats the result as `view_mesh`: the
directed view whose vertices are cells and whose edges are the
two-cell faces (a one-cell face is an edge without a head), with
the measurements as fields and the tag names on the boundary edges.

The pre-tower graph stack (`interface_graph`, `class_stored_graph`,
`class_mesh`, `class_array_mesh_loader`, about 4,400 lines) was
deleted on 2026-08-19 after `view_mesh_geometry` reproduced its
derived arrays bitwise on all eleven sample meshes (2d quadrilateral
and triangle, 3d hexahedral and tetrahedral). A static in
test/graph-mesh keeps it deleted. AGENTS.md phase 11 is complete:
exactly one type named `graph` exists.

## The differential operator, in the primes

The differential operator's three elementary steps - average,
difference, incidence - are affine sparse maps (a matrix plus the
constant a boundary value leaves behind). The order-n operator is
their composition by parity, computed as one sparse triple product
with duplicates combined, so the operator IS a stencil: `stencil_of`
returns the square vertex-landing form as a `stencil`, and
the adjoint is the transpose of the composed matrix - no reversed
step kernels exist. The contract suite holds the agreement law:
applying the operator and applying its compiled stencil give the
same numbers, boundary constants and adjoint included.

## One directed view, relation-seated (consolidated 2026-08-19)

The directed view D = (V, E, tail, head) has one implementation:
`stored_directed_graph`. Its stored arrays are the compiled
snapshots the hot paths read; its `tail_relation()` and
`head_relation()` derive T <= E x V and H <= E x V as csr_relations
on request, over counted coordinates, so a graph nobody reads
relationally pays nothing (the section-66 benchmark rejected the
eager form at 1.7-2.2x on pattern and part construction; the lazy
form is within noise of baseline). `relation_algorithms` (sources,
sinks, reachable, topological_order) takes any `class(relation)`,
selects the binary reading, and refuses a non-binary or two-domain
adjacency. `graph_profile` and its two view types are deleted; a
static in test/graph-ordinary keeps them deleted, and that suite
holds the equivalence law: the stored graph against externally
built T and H relations, per edge and per vertex, lists compared
as sets, on topologies up to a 400-vertex pseudo-random one with
walls. AGENTS.md phases 9-11 are complete.

## Documentation rule

`doc/` holds current truth only: this file, `AGENTS.md` (the
constitution, at the repository root), `coding-standards.md`,
`GTI.md`, `GTI-roadmap.md`. Superseded plans and era reports live in
`doc/history/`. A document that stops being true moves there in the
same commit that falsifies it.
