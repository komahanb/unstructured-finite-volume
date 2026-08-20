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

    graph            G = (B1, B2),  B in {NULL, UNKNOWN, KNOWN -> G}     fractal_graph
    token            object identity, minted once by declare             graph_identity
    relation         P subset of A1 x ... x Ak                           graph_relation
    field            f : A -> V on one set-graph domain                  graph_field_calculus
    operation        (graph, fields) -> field                            graph_operation_view
    transform        graph -> graph, data carried along                  graph_operation_view
    map              token -> {representation | name | inclusion | value}
                                                                         graph_set_map, graph_label_map,
                                                                         graph_inclusion_map, graph_value_map
    view             a reading of the kernel graph, never a new kind     graph_directed_view,
                                                                         graph_relational_view, graph_profile

The kernel graph is self-similar: each branch is NULL, UNKNOWN, or
another graph. Sets, tuples, relations, domains, and chains are all
this one object read differently. There is exactly one type named
`graph` in the living tower.

## The composition tree

Indentation is `extends`; brackets name the file.

    relation                                   [graph_relation]
    ├── stored_relation                        [graph_relation]
    └── binary_relation                        [graph_binary_relation]
        ├── csr_relation                       [graph_binary_relation]
        └── transposed_view                    [graph_binary_relation]

    directed_graph  D = (V, E, tail, head)     [graph_directed_view]
    └── directed_stored_graph                  [class_graph]

    graph_field                                [graph_field_calculus]
    ├── field                                  [class_graph_field]
    └── graph_functional  (one entry)          [graph_field_calculus]
        └── functional                         [class_graph_functional]

    graph_operation                            [graph_operation_view]
    ├── graph_reduction   field -> functional  [class_graph_reduction]
    │   └── reduction   (sum, average, norm, ...)
    ├── graph_broadcast   functional -> field  [class_graph_reduction]
    │   └── broadcast   (copy, share)
    ├── discretization_operator                [graph_discretization]
    │   ├── stencil_operator  space: matrix as weighted edges   [class_graph_stencil]
    │   └── step_operator     time:  a0 q + a1 qold + a2 qolder + hs S(q)   [class_graph_step]
    ├── linearization_operator  J v at a frozen state           [graph_discretization]
    │   ├── difference_linearization           [class_graph_linearization]
    │   └── exact_linearization  + tangent_of chooser           [class_graph_exact_linearization]
    ├── differentiable_operation  exact partial actions to max_degree   [graph_discretization]
    ├── minimizer                              [graph_minimization]
    │   ├── jacobi, gauss_seidel, conjugate_gradient, gmres, multigrid
    │   ├── newton        Jacobian action = tangent_of(action)
    │   └── dense_direct  elimination; + dense_matrix_of, dense-array adapter
    ├── differential_operator  div/grad on graph sides          [class_graph_differential_operator]
    ├── balance, fit, walk                     [clients and utilities]
    └── ...

    graph_transform                            [graph_operation_view]
    ├── graph_partitioner -> partitioner       [class_graph_partitioner]
    ├── graph_assembler   -> assembler         [class_graph_assembler]
    ├── graph_coarsener   -> coarsener         [class_graph_coarsener]
    └── graph_refiner     -> refiner           [class_graph_refiner]

    reversible_change  apply -> check -> keep | revert   [graph_change_protocol]
    └── value_change                           [graph_value_change]

    chain_rule + argument_path + path_derivative          [class_graph_chain_rule]
        total derivatives of S(x_1(s), ..., x_m(s)) by integer
        partitions of the degree, from partial actions alone

    marcher  time as a chain graph             [class_graph_marcher]
        march, march_adjoint (backward substitution),
        march_directional (any order), march_adaptive
    step_policy -> halving_policy              [class_graph_step_policy]

    form -> polynomial_form, harmonic_form     [graph_forms and clients]
    form_optimizer -> pruner, fit              [graph_fitting]

## Composition, stated as arithmetic

    field        = graph  x  values
    relation     = graph  x  tuples
    operation    = graph  x  fields  ->  field
    stencil      = operation with weights on edges          (a matrix)
    step         = operation x (a0, a1, a2, hs)             (a scheme)
    tangent      = operation frozen at a state              (a matrix action)
    chain rule   = partitions x partial actions             (any derivative)
    minimizer    = operation + a stopping rule              (a solve)
    newton       = minimizer x tangent_of                   (nonlinear solve)
    marcher      = chain graph x step x minimizer           (integration)
    adjoint      = marcher read backward x transpose        (sensitivity)
    change       = mutation x memory x undo                 (transaction)
    mesh         = sets x relations x measurement fields    (see below)

## Level assignment (AGENTS.md section 2)

    0-3  kernel:      fractal_graph, graph_identity, graph_relation,
                      graph_binary_relation, graph_partition_relation,
                      graph_set_representation, the four maps
    4    graph view:  graph_directed_view, graph_relational_view,
                      graph_profile, class_graph, class_graph_walk,
                      transforms (partitioner/assembler/coarsener/refiner)
    5    field:       graph_field_calculus, class_graph_field,
                      class_graph_functional, class_graph_reduction
    6    calculus:    graph_discretization, stencil, step,
                      linearizations, chain rule, differential
                      operator, forms/fitting
    7    solve:       graph_minimization + members, marcher,
                      step_policy, change protocol + value map
    8-9  physics:     balance, fitted_balance, advection, conduction,
                      diffusion_statement, robin_condition,
                      mesh + loaders

An import may only point downward in this list.

## The mesh, in the primes

A mesh file supplies member sets (vertices, cells, tagged boundary
faces), one relation (C2V, cell to vertex), one field (coordinates
on vertices), and tag names. Everything else is derived, in
`graph_mesh_geometry`, by relation algebra and geometry:

    faces  = boundary faces + shared-vertex intersections of cell pairs
    V2C    = transpose(C2V)
    F2C    = {(f, c) : F2V(f) subset of C2V(c)}     (<= 2 cells per face)
    C2F    = transpose(F2C)
    measurements = pure functions of coordinates and the relations:
                   cell centers/volumes, face centers/areas/normals,
                   center-to-center vectors, deltas, weights

`class_mesh_builder` seats the result as `class_graph_mesh`: the
directed view whose vertices are cells and whose edges are the
two-cell faces (a one-cell face is an edge without a head), with
the measurements as fields and the tag names on the boundary edges.

The pre-tower graph stack (`interface_graph`, `class_stored_graph`,
`class_mesh`, `class_array_mesh_loader`, about 4,400 lines) was
deleted on 2026-08-19 after `graph_mesh_geometry` reproduced its
derived arrays bitwise on all eleven sample meshes (2d quadrilateral
and triangle, 3d hexahedral and tetrahedral). A static in
test/graph-mesh keeps it deleted. AGENTS.md phase 11 is complete:
exactly one type named `graph` exists.

## Review candidates (not duplication, but under watch)

- `class_graph_differential_operator` (1,208 lines, level 6): the
  original FV differential operator. Overlaps `stencil_operator`
  conceptually; consumed by `class_graph_balance` and ten suites.
  Candidate for re-expression over stencils.

## Documentation rule

`doc/` holds current truth only: this file, `AGENTS.md` (the
constitution, at the repository root), `coding-standards.md`,
`GTI.md`, `GTI-roadmap.md`. Superseded plans and era reports live in
`doc/history/`. A document that stops being true moves there in the
same commit that falsifies it.
