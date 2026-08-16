# Naming audit — every public symbol, one category

> **Superseded, 2026-08-15.** This document is premised on a split between a structural
> graph and a computational graph. That split is withdrawn: there is one graph ontology,
> `G=(B1,B2)` in `src/fractal_graph.f90`, and `(S,P)` and `(Q,R)` are views of it. See
> AGENTS.md, "The graph ontology". Retained as history; do not implement.

The §15 deliverable of the computational-graph naming pass
(`COMPUTATIONAL-GRAPH.md`): before any symbol is renamed, its category
must be explicit. Fifty-six source modules were read completely; every
symbol on every `public ::` list was classified into exactly one of:

```text
STRUCTURE   carriers, relations, algebra, the relational graph GAMMA,
            identity machinery, structural transforms of topology
DATA        values over domains - fields, functionals, measured
            numbers; candidate constituents of Q
OPERATOR    discrete residual/operator machinery - differential
            operators, stencils, balances, linearizations;
            candidate material for R
INFERENCE   processes between epistemic states - minimization,
            model discovery, marching
VIEW        interpretations reading a structural schema and
            answering derived vocabulary
LEGACY      compatibility contracts of the old four-root grammar,
            kept alive for migration and nothing else
```

No production symbol is renamed in this pass. The categories exist so
that the follow-up (releveling and any renames,
doc/computational-graph-releveling.md) starts from explicit ground.

## How it was verified

Ten independent auditors each read a disjoint slice of src/; an
eleventh swept tests and documents for reserved-letter and pair-notation
collisions. Two adversarial verifiers then re-derived the result:

- **Completeness**: a continuation-line-aware parse of every `public`
  statement across src/*.f90 found no omission. (The check was run
  against a snapshot taken before graph_state.f90 landed; that module's
  ten symbols are included in the table below by hand and are the
  vocabulary this pass introduces.) One incidental finding: 
  `module_mesh_utils.f90` carries no `private` statement, so its
  use-associated rename `dp => REAL64` is technically re-exported.
- **Categories**: all assignments were checked against the actual
  declarations; zero disputes survived. The two closest calls, both
  upheld on code evidence: `fit` is OPERATOR, not INFERENCE — its apply
  is deterministic stencil-weight generation (moment matrix, dual
  solve, weights), no epistemic state moves, and the module's own
  doctrine splits the coefficient sector (`fit`) from the form sector
  (`form_optimizer`, which IS model discovery and therefore
  INFERENCE); and `walk` is VIEW, not LEGACY — a live citizen that
  reads structure to answer derived vocabulary.

## The reserved letter, after the pass

The letter R now appears in exactly two capacities, and the audit
holds the line between them:

**Correct — residual semantics (stays R):**

```text
adjoint tower            R_q, R_q^T, R_qq, R_p     residual Jacobians
time-integration tower   R_BE(q) = q - q0 + h S(q) backward-Euler residual
learning tower           R(w) = 2w - 6             discovered residual
doc/coupled_diffusion_adjoint_fv_implementation_spec.md
                         R_u, R_v, R_i, R_h, R_D, R_alpha, R_beta
                                                   residual components
                                                   and derivatives
```

**Deferred — structural relations still wearing R (rename in the
releveling/rename commit, not here):**

```text
test/adjoint-tower       R_dep <= T x V            the structural
                                                   dependency relation;
                                                   ~24 comment/doc
                                                   sites, no code
                                                   identifiers.
                                                   Proposed: T_dep,
                                                   in the adjoint
                                                   tower's own commit
                                                   so its Gate B work
                                                   is not disturbed.
```

Everything else was renamed in this pass: the structural pair is
GAMMA = (S, P) in every contract document, tower README, banner and
module header; the calculator flow relation is T_flow (identifiers
t_flow/t_out3/t_in3 in the calculator tower tests); the compose law
reads compose_binary(P_AB, P_BC) = P_BC o P_AB; and the field
value-space is written `-> reals`, never `-> R`.

## Confusables the auditors flagged (no action, recorded)

- **materialized vs realized**: `materialized()` (graph_relation and
  the algebra) means self-contained storage safe to own; `realized`
  (graph_state) means an epistemic seat is occupied. Different words
  for different facts — keep both, never blur them.
- **state**: `graph_state`'s state is epistemic (which seats are
  occupied); a solver's state vector and a physical state field are Q
  material. Comments should say "epistemic state" when near solver
  vocabulary.
- **operator**: `graph_operation` is the legacy generic verb;
  `operator graph` is the epistemic state (bot, R); a discrete
  operator is OPERATOR-category machinery. Bare "operator" in new
  text must say which.
- The full per-line list of auditor-flagged ambiguities follows the
  symbol table.

## graph_state.f90 — the vocabulary this pass introduces

| module | header claim | symbol | kind | category |
|---|---|---|---|---|
| `graph_state.f90` | NOT YET A LEVEL . THE COMPUTATIONAL GRAPH | `computational_graph` | type | **(the pair itself)** |
| `graph_state.f90` |  | `GRAPH_STATE_VOID` | constant | **(epistemic state)** |
| `graph_state.f90` |  | `GRAPH_STATE_DATA` | constant | **(epistemic state)** |
| `graph_state.f90` |  | `GRAPH_STATE_OPERATOR` | constant | **(epistemic state)** |
| `graph_state.f90` |  | `GRAPH_STATE_REALIZED` | constant | **(epistemic state)** |
| `graph_state.f90` |  | `state_name` | procedure | **(diagnostics)** |
| `graph_state.f90` |  | `void_graph` | procedure | **(constructor)** |
| `graph_state.f90` |  | `data_graph` | procedure | **(constructor)** |
| `graph_state.f90` |  | `operator_graph` | procedure | **(constructor)** |
| `graph_state.f90` |  | `realized_graph` | procedure | **(constructor)** |

These ten deliberately refuse the six categories: the computational
graph G = (Q, R) is the object the categories describe FROM — its data
seat will hold DATA citizens, its residual seat OPERATOR citizens, and
INFERENCE moves it between states. Classifying the classifier into its
own boxes would be a category error, and the completeness verifier said
as much.

## Symbol table

| module | header claim | symbol | kind | category |
|---|---|---|---|---|
| `class_advection.f90` | LEVEL 3 OF THE STRATIFICATION. An advection law holds one velocity | `advection` | type | **DATA** |
| `class_array_mesh_loader.f90` | none | `array_mesh_loader` | type | **STRUCTURE** |
| `class_conduction.f90` | LEVEL 3 OF THE STRATIFICATION. A conduction law holds one tensor, | `conduction` | type | **DATA** |
| `class_diffusion_statement.f90` | LEVEL 3 OF THE STRATIFICATION. This is where the physics words | `diffusion_statement` | procedure | **OPERATOR** |
| `class_file.f90` | none | `file` | type | **DATA** |
| `class_fitted_balance.f90` | LEVEL 3 - the statements' shared machinery. Level 2 holds | `fitted_balance_stencil` | procedure | **OPERATOR** |
| `class_form_pruner.f90` | none | `pruner` | type | **INFERENCE** |
| `class_gmsh_loader.f90` | none | `gmsh_loader` | type | **STRUCTURE** |
| `class_gmsh_loader.f90` | none | `gmsh_loader (generic constructor interface -> create)` | interface | **STRUCTURE** |
| `class_graph.f90` | none | `stored_graph` | type | **LEGACY** |
| `class_graph_assembler.f90` | none | `assembler` | type | **STRUCTURE** |
| `class_graph_balance.f90` | none | `balance` | type | **OPERATOR** |
| `class_graph_coarsener.f90` | none | `coarsener` | type | **STRUCTURE** |
| `class_graph_coarsener.f90` | none | `COARSEN_PAIRWISE` | constant | **STRUCTURE** |
| `class_graph_coarsener.f90` | none | `COARSEN_ADOPTED` | constant | **STRUCTURE** |
| `class_graph_conjugate_gradient.f90` | none | `conjugate_gradient` | type | **INFERENCE** |
| `class_graph_differential_operator.f90` | none | `differential_operator` | type | **OPERATOR** |
| `class_graph_differential_operator.f90` | none | `edge_differential_operator` | procedure | **OPERATOR** |
| `class_graph_differential_operator.f90` | none | `vertex_differential_operator` | procedure | **OPERATOR** |
| `class_graph_differential_operator.f90` | none | `gradient` | procedure | **OPERATOR** |
| `class_graph_differential_operator.f90` | none | `interpolation` | procedure | **OPERATOR** |
| `class_graph_differential_operator.f90` | none | `divergence` | procedure | **OPERATOR** |
| `class_graph_differential_operator.f90` | none | `laplacian` | procedure | **OPERATOR** |
| `class_graph_field.f90` | none | `field` | type | **DATA** |
| `class_graph_functional.f90` | none | `functional` | type | **DATA** |
| `class_graph_gauss_seidel.f90` | none | `gauss_seidel` | type | **INFERENCE** |
| `class_graph_gmres.f90` | none | `gmres` | type | **INFERENCE** |
| `class_graph_jacobi.f90` | none | `jacobi` | type | **INFERENCE** |
| `class_graph_linearization.f90` | LEVEL 1 OF THE STRATIFICATION. The first concretion of the | `difference_linearization` | type | **OPERATOR** |
| `class_graph_marcher.f90` | LEVEL 2 OF THE STRATIFICATION. Time is not a special dimension | `marcher` | type | **INFERENCE** |
| `class_graph_marcher.f90` | LEVEL 2 OF THE STRATIFICATION. Time is not a special dimension | `MARCH_FORWARD` | constant | **INFERENCE** |
| `class_graph_marcher.f90` | LEVEL 2 OF THE STRATIFICATION. Time is not a special dimension | `MARCH_BACKWARD` | constant | **INFERENCE** |
| `class_graph_marcher.f90` | LEVEL 2 OF THE STRATIFICATION. Time is not a special dimension | `MARCH_BDF2` | constant | **INFERENCE** |
| `class_graph_mesh.f90` | LEVEL 1 OF THE STRATIFICATION - for measurements are calculus | `mesh` | type | **DATA** |
| `class_graph_mesh.f90` | LEVEL 1 OF THE STRATIFICATION - for measurements are calculus | `mesh (generic constructor interface -> create)` | interface | **DATA** |
| `class_graph_multigrid.f90` | none | `multigrid` | type | **INFERENCE** |
| `class_graph_newton.f90` | none | `newton` | type | **INFERENCE** |
| `class_graph_partitioner.f90` | none | `partitioner` | type | **STRUCTURE** |
| `class_graph_partitioner.f90` | none | `PARTITION_LINEAR` | constant | **STRUCTURE** |
| `class_graph_partitioner.f90` | none | `PARTITION_BREADTH_FIRST` | constant | **STRUCTURE** |
| `class_graph_partitioner.f90` | none | `PARTITION_ADOPTED` | constant | **STRUCTURE** |
| `class_graph_reduction.f90` | none | `reduction` | type | **OPERATOR** |
| `class_graph_reduction.f90` | none | `REDUCE_SUM` | constant | **OPERATOR** |
| `class_graph_reduction.f90` | none | `REDUCE_AVERAGE` | constant | **OPERATOR** |
| `class_graph_reduction.f90` | none | `REDUCE_MINIMUM` | constant | **OPERATOR** |
| `class_graph_reduction.f90` | none | `REDUCE_MAXIMUM` | constant | **OPERATOR** |
| `class_graph_reduction.f90` | none | `REDUCE_NORM` | constant | **OPERATOR** |
| `class_graph_reduction.f90` | none | `REDUCE_COUNT` | constant | **OPERATOR** |
| `class_graph_reduction.f90` | none | `REDUCE_ALL` | constant | **OPERATOR** |
| `class_graph_reduction.f90` | none | `REDUCE_ANY` | constant | **OPERATOR** |
| `class_graph_reduction.f90` | none | `broadcast` | type | **OPERATOR** |
| `class_graph_reduction.f90` | none | `BROADCAST_COPY` | constant | **OPERATOR** |
| `class_graph_reduction.f90` | none | `BROADCAST_SHARE` | constant | **OPERATOR** |
| `class_graph_refiner.f90` | none | `refiner` | type | **STRUCTURE** |
| `class_graph_stencil.f90` | none | `stencil_operator` | type | **OPERATOR** |
| `class_graph_step.f90` | LEVEL 1 OF THE STRATIFICATION. One step of a time recurrence, as | `step_operator` | type | **OPERATOR** |
| `class_graph_step.f90` | LEVEL 1 OF THE STRATIFICATION. One step of a time recurrence, as | `backward_euler` | procedure | **OPERATOR** |
| `class_graph_step.f90` | LEVEL 1 OF THE STRATIFICATION. One step of a time recurrence, as | `bdf` | procedure | **OPERATOR** |
| `class_graph_walk.f90` | none | `walk` | type | **VIEW** |
| `class_graph_walk.f90` | none | `WALK_COLOURING` | constant | **VIEW** |
| `class_graph_walk.f90` | none | `WALK_VISIT_ORDER` | constant | **VIEW** |
| `class_graph_walk.f90` | none | `WALK_COMPONENT` | constant | **VIEW** |
| `class_graph_walk.f90` | none | `WALK_DEPTH` | constant | **VIEW** |
| `class_harmonic_form.f90` | none | `harmonic_form` | type | **OPERATOR** |
| `class_mesh.f90` | none | `mesh` | type | **LEGACY** |
| `class_mesh.f90` | none | `mesh (generic constructor interface -> create_mesh)` | interface | **LEGACY** |
| `class_mesh_builder.f90` | none | `mesh_from_gmsh` | procedure | **DATA** |
| `class_polynomial_form.f90` | none | `polynomial_form` | type | **OPERATOR** |
| `class_robin_condition.f90` | LEVEL 3 OF THE STRATIFICATION. The constitution says what the | `robin_condition` | type | **DATA** |
| `class_robin_condition.f90` | LEVEL 3 OF THE STRATIFICATION. The constitution says what the | `robin` | procedure | **DATA** |
| `class_robin_condition.f90` | LEVEL 3 OF THE STRATIFICATION. The constitution says what the | `dirichlet` | procedure | **DATA** |
| `class_robin_condition.f90` | LEVEL 3 OF THE STRATIFICATION. The constitution says what the | `neumann` | procedure | **DATA** |
| `class_stored_graph.f90` | none | `stored_graph` | type | **LEGACY** |
| `class_stored_graph.f90` | none | `stored_digraph` | type | **LEGACY** |
| `class_string.f90` | none | `string` | type | **DATA** |
| `graph_algorithms.f90` | LEVEL 4 OF THE NEW TOWER . THE GRAPH ALGORITHMS | `sources` | procedure | **VIEW** |
| `graph_algorithms.f90` | LEVEL 4 OF THE NEW TOWER . THE GRAPH ALGORITHMS | `sinks` | procedure | **VIEW** |
| `graph_algorithms.f90` | LEVEL 4 OF THE NEW TOWER . THE GRAPH ALGORITHMS | `reachable` | procedure | **VIEW** |
| `graph_algorithms.f90` | LEVEL 4 OF THE NEW TOWER . THE GRAPH ALGORITHMS | `topological_order` | procedure | **VIEW** |
| `graph_binary_relation.f90` | LEVEL 1 OF THE NEW TOWER . THE BINARY SPECIALIZATION | `binary_relation` | type | **STRUCTURE** |
| `graph_binary_relation.f90` | LEVEL 1 OF THE NEW TOWER . THE BINARY SPECIALIZATION | `csr_relation` | type | **STRUCTURE** |
| `graph_binary_relation.f90` | LEVEL 1 OF THE NEW TOWER . THE BINARY SPECIALIZATION | `transposed_view` | type | **STRUCTURE** |
| `graph_binary_relation.f90` | LEVEL 1 OF THE NEW TOWER . THE BINARY SPECIALIZATION | `transpose_of` | procedure | **STRUCTURE** |
| `graph_binary_relation.f90` | LEVEL 1 OF THE NEW TOWER . THE BINARY SPECIALIZATION | `inclusion_of` | procedure | **STRUCTURE** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `graph_functional` | type | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `graph_reduction` | type | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `graph_broadcast` | type | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `discretization_operator` | type | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `linearization_operator` | type | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `graph_partitioner` | type | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `graph_assembler` | type | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `graph_coarsener` | type | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `graph_refiner` | type | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `GRAPH_SIDE_VERTEX` | constant | **LEGACY** |
| `graph_calculus.f90` | THE LEGACY COMPATIBILITY CALCULUS | `GRAPH_SIDE_EDGE` | constant | **LEGACY** |
| `graph_carrier.f90` | LEVEL 0 OF THE NEW TOWER . THE CARRIERS | `member_set` | type | **STRUCTURE** |
| `graph_carrier.f90` | LEVEL 0 OF THE NEW TOWER . THE CARRIERS | `counted_set` | type | **STRUCTURE** |
| `graph_carrier.f90` | LEVEL 0 OF THE NEW TOWER . THE CARRIERS | `subset_set` | type | **STRUCTURE** |
| `graph_field_calculus.f90` | LEVEL 5 OF THE NEW TOWER . THE FIELD CALCULUS | `graph_field` | type | **DATA** |
| `graph_field_calculus.f90` | LEVEL 5 OF THE NEW TOWER . THE FIELD CALCULUS | `GRAPH_FIELD_INTEGER` | constant | **DATA** |
| `graph_field_calculus.f90` | LEVEL 5 OF THE NEW TOWER . THE FIELD CALCULUS | `GRAPH_FIELD_REAL` | constant | **DATA** |
| `graph_field_calculus.f90` | LEVEL 5 OF THE NEW TOWER . THE FIELD CALCULUS | `GRAPH_FIELD_COMPLEX` | constant | **DATA** |
| `graph_field_calculus.f90` | LEVEL 5 OF THE NEW TOWER . THE FIELD CALCULUS | `GRAPH_FIELD_LOGICAL` | constant | **DATA** |
| `graph_field_calculus.f90` | LEVEL 5 OF THE NEW TOWER . THE FIELD CALCULUS | `GRAPH_FIELD_CHARACTER` | constant | **DATA** |
| `graph_fitting.f90` | LEVEL 2 OF THE STRATIFICATION . THE FITTING FAMILY | `fit` | type | **OPERATOR** |
| `graph_fitting.f90` | LEVEL 2 OF THE STRATIFICATION . THE FITTING FAMILY | `form_optimizer` | type | **INFERENCE** |
| `graph_forms.f90` | LEVEL 1 . THE FORMS | `form` | type | **OPERATOR** |
| `graph_grammar.f90` | THE LEGACY ORDINARY-GRAPH COMPATIBILITY CONTRACT | `graph` | type | **LEGACY** |
| `graph_grammar.f90` | THE LEGACY ORDINARY-GRAPH COMPATIBILITY CONTRACT | `graph_field` | reexport | **LEGACY** |
| `graph_grammar.f90` | THE LEGACY ORDINARY-GRAPH COMPATIBILITY CONTRACT | `graph_operation` | type | **LEGACY** |
| `graph_grammar.f90` | THE LEGACY ORDINARY-GRAPH COMPATIBILITY CONTRACT | `graph_transform` | type | **LEGACY** |
| `graph_grammar.f90` | THE LEGACY ORDINARY-GRAPH COMPATIBILITY CONTRACT | `GRAPH_FIELD_INTEGER` | reexport | **LEGACY** |
| `graph_grammar.f90` | THE LEGACY ORDINARY-GRAPH COMPATIBILITY CONTRACT | `GRAPH_FIELD_REAL` | reexport | **LEGACY** |
| `graph_grammar.f90` | THE LEGACY ORDINARY-GRAPH COMPATIBILITY CONTRACT | `GRAPH_FIELD_COMPLEX` | reexport | **LEGACY** |
| `graph_grammar.f90` | THE LEGACY ORDINARY-GRAPH COMPATIBILITY CONTRACT | `GRAPH_FIELD_LOGICAL` | reexport | **LEGACY** |
| `graph_grammar.f90` | THE LEGACY ORDINARY-GRAPH COMPATIBILITY CONTRACT | `GRAPH_FIELD_CHARACTER` | reexport | **LEGACY** |
| `graph_identity.f90` | THE IDENTITY . INFRASTRUCTURE BENEATH THE TOWER | `token` | type | **STRUCTURE** |
| `graph_identity.f90` | THE IDENTITY . INFRASTRUCTURE BENEATH THE TOWER | `mint_token` | procedure | **STRUCTURE** |
| `graph_minimization.f90` | LEVEL 7 OF THE NEW TOWER . THE MINIMIZATION | `minimizer` | type | **INFERENCE** |
| `graph_profile.f90` | LEVEL 4 OF THE NEW TOWER . THE ORDINARY GRAPH PROFILE | `ordinary_graph_view` | type | **VIEW** |
| `graph_profile.f90` | LEVEL 4 OF THE NEW TOWER . THE ORDINARY GRAPH PROFILE | `directed_adjacency_view` | type | **VIEW** |
| `graph_relation.f90` | LEVEL 1 OF THE NEW TOWER . THE RELATIONS | `relation` | type | **STRUCTURE** |
| `graph_relation.f90` | LEVEL 1 OF THE NEW TOWER . THE RELATIONS | `stored_relation` | type | **STRUCTURE** |
| `graph_relation.f90` | LEVEL 1 OF THE NEW TOWER . THE RELATIONS | `slot` | type | **STRUCTURE** |
| `graph_relation_algebra.f90` | LEVEL 2 OF THE NEW TOWER . THE RELATION ALGEBRA | `restrict_slot` | procedure | **STRUCTURE** |
| `graph_relation_algebra.f90` | LEVEL 2 OF THE NEW TOWER . THE RELATION ALGEBRA | `project_slots` | procedure | **STRUCTURE** |
| `graph_relation_algebra.f90` | LEVEL 2 OF THE NEW TOWER . THE RELATION ALGEBRA | `compose_binary` | procedure | **STRUCTURE** |
| `graph_structure.f90` | LEVEL 3 OF THE NEW TOWER . THE RELATIONAL GRAPH | `relational_graph` | type | **STRUCTURE** |
| `graph_structure.f90` | LEVEL 3 OF THE NEW TOWER . THE RELATIONAL GRAPH | `held_set` | type | **STRUCTURE** |
| `graph_structure.f90` | LEVEL 3 OF THE NEW TOWER . THE RELATIONAL GRAPH | `held_relation` | type | **STRUCTURE** |
| `interface_graph.f90` | none | `graph` | type | **LEGACY** |
| `interface_graph.f90` | none | `digraph` | type | **LEGACY** |
| `interface_graph.f90` | none | `vertex` | type | **LEGACY** |
| `interface_graph.f90` | none | `edge` | type | **LEGACY** |
| `interface_graph.f90` | none | `counting_sort` | procedure | **STRUCTURE** |
| `interface_graph.f90` | none | `transpose_adjacency` | procedure | **STRUCTURE** |
| `interface_graph.f90` | none | `power_iteration` | procedure | **INFERENCE** |
| `interface_mesh_loader.f90` | none | `mesh_loader` | type | **STRUCTURE** |
| `module_mesh_utils.f90` | none | `distance` | interface | **DATA** |
| `module_mesh_utils.f90` | none | `distanceAB` | procedure | **DATA** |
| `module_mesh_utils.f90` | none | `cross_product` | procedure | **DATA** |
| `module_mesh_utils.f90` | none | `find` | procedure | **STRUCTURE** |
| `module_mesh_utils.f90` | none | `is_subset` | procedure | **STRUCTURE** |
| `module_mesh_utils.f90` | none | `isort` | procedure | **STRUCTURE** |
| `module_mesh_utils.f90` | none | `elem_type_face_count` | procedure | **STRUCTURE** |
| `module_mesh_utils.f90` | none | `elem_type_dimension` | procedure | **STRUCTURE** |
| `module_mesh_utils.f90` | none | `elem_type_vertex_count` | procedure | **STRUCTURE** |
| `module_mesh_utils.f90` | none | `order_face_vertices` | procedure | **STRUCTURE** |
| `module_solve_mode.f90` | none | `FORWARD` | constant | **INFERENCE** |
| `module_solve_mode.f90` | none | `REVERSE` | constant | **INFERENCE** |
| `module_solve_mode.f90` | none | `WHOLE` | constant | **OPERATOR** |
| `module_solve_mode.f90` | none | `DIAGONAL` | constant | **OPERATOR** |
| `module_solve_mode.f90` | none | `LOWER_TRIANGLE` | constant | **OPERATOR** |
| `module_solve_mode.f90` | none | `UPPER_TRIANGLE` | constant | **OPERATOR** |
| `module_solve_mode.f90` | none | `is_valid_mode` | procedure | **INFERENCE** |
| `module_solve_mode.f90` | none | `is_valid_part` | procedure | **OPERATOR** |
| `module_verbosity.f90` | none | `verbosity` | constant | **VIEW** |
| `module_verbosity.f90` | none | `set_verbosity` | procedure | **VIEW** |

## Counts

```text
STRUCTURE   41
DATA        22
OPERATOR    36
INFERENCE   17
VIEW        13
LEGACY      29
total      158
```

## Rationales (per symbol)


### class_advection.f90

- `advection` (type) — **DATA**: A flow parameter payload — one velocity vector — that supplies per-face coefficient arrays (v.n, vn*area); a parameter constituent of Q that owns no operator and no balance; the same-named generic interface is its constructor.

### class_array_mesh_loader.f90

- `array_mesh_loader` (type) — **STRUCTURE**: An in-memory source honoring the loader contract: it stores and copies out the structural incidence lists, tag identity, and coordinate payload from which a mesh's structure is wired — a carrier of structural description, not of computed values over a realized domain.

### class_conduction.f90

- `conduction` (type) — **DATA**: A material parameter payload — one conductivity tensor K — that supplies per-face coefficient arrays (n^T K n, keff*area); a parameter constituent of Q that owns no operator and no balance; the same-named generic interface is its isotropic/tensor constructor.

### class_diffusion_statement.f90

- `diffusion_statement` (procedure) — **OPERATOR**: Translates a conduction law and robin conditions into scales and wall relations and delegates to the fitted-balance assembly, returning the compiled discrete diffusion operator — construction of R material.

### class_file.f90

- `file` (type) — **DATA**: A sequential reader delivering stored lines/values from disk into memory — an ingestion conduit for measured/stored numbers (Q-candidate inputs like mesh files), and not itself structure, operator, or inference; its same-name generic interface is the constructor.

### class_fitted_balance.f90

- `fitted_balance_stencil` (procedure) — **OPERATOR**: Assembles per-edge exactly-fitted flux weights, exchanged through incidence, into one compiled stencil_operator — a discrete balance/discretization action producing candidate R material.

### class_form_pruner.f90

- `pruner` (type) — **INFERENCE**: Concrete form_optimizer that strikes basis members the point constellation cannot see by mutating the form's roster (shape % restrict) — model adaptation/discovery, not a VIEW-style pruner since it changes the model rather than re-presenting it.

### class_gmsh_loader.f90

- `gmsh_loader` (type) — **STRUCTURE**: A parser of MSH 4.1 files that satisfies the loader contract by producing the structural incidence lists, element-type codes, and the entity-to-physical-tag identity table — it delivers structural description, measuring nothing.
- `gmsh_loader (generic constructor interface -> create)` (interface) — **STRUCTURE**: Constructor that points the structural-description parser at a .msh file and sizes its line buffer; purely part of the structural source's identity.

### class_graph.f90

- `stored_graph` (type) — **LEGACY**: Concrete citizen of graph_grammar's legacy graph root, storing incidence/adjacency lists, tags, part frame, and the two carrier sets; semantically STRUCTURE.

### class_graph_assembler.f90

- `assembler` (type) — **STRUCTURE**: Concrete inverse structural transform (part back to whole-graph order via recorded index/ownership maps), acting purely on topology and its carried indexing, with no residual or solver content.

### class_graph_balance.f90

- `balance` (type) — **OPERATOR**: The discrete conservation balance — face terms reduced onto cells through incidence plus a source — which is exactly the residual assembly a solver drives; core candidate for R (the same-named generic interface is its constructor).

### class_graph_coarsener.f90

- `coarsener` (type) — **STRUCTURE**: Concrete structural transform of topology (fine cells glued into blocks, faces redrawn between blocks); the aggregate map it owns is a structural assignment, and its data carry rides that structure.
- `COARSEN_PAIRWISE` (constant) — **STRUCTURE**: Rule selector for how the structural gluing of cells into blocks is chosen.
- `COARSEN_ADOPTED` (constant) — **STRUCTURE**: Rule selector saying the cell-to-block map (a structural assignment) is taken from outside.

### class_graph_conjugate_gradient.f90

- `conjugate_gradient` (type) — **INFERENCE**: A concrete Krylov solution process (CG) minimizing the residual over conjugate directions — solver machinery, not structure, data, or the operator itself.

### class_graph_differential_operator.f90

- `differential_operator` (type) — **OPERATOR**: A matrix-free discrete differential operator of arbitrary order built from the S/G/D elementary steps — prime candidate material for the residual R.
- `edge_differential_operator` (procedure) — **OPERATOR**: Constructor producing the edge-landing differential operator; it configures operator machinery, nothing else.
- `vertex_differential_operator` (procedure) — **OPERATOR**: Constructor producing the vertex-landing differential operator, including its adjoint (transposed-walk) form.
- `gradient` (procedure) — **OPERATOR**: Named constructor for the order-1 edge operator G q — a discrete differential operator by definition.
- `interpolation` (procedure) — **OPERATOR**: Named constructor for the order-0 edge-average operator S q, an elementary discretization action.
- `divergence` (procedure) — **OPERATOR**: Named constructor for the order-1 vertex operator (incidence step on an edge field) — the discrete divergence.
- `laplacian` (procedure) — **OPERATOR**: Named constructor for the order-2 vertex operator D G q — the discrete second derivative.

### class_graph_field.f90

- `field` (type) — **DATA**: The one concrete store of values over a member-set domain (name/units/domain/width plus one live typed store); a direct constituent of Q, and the same public name also serves as its constructor generic interface.

### class_graph_functional.f90

- `functional` (type) — **DATA**: A single reduced value (total, objective, norm, yes-or-no answer) — the field at domain size one, so a scalar constituent of Q, with reduction workspace (tally, weight) that is still stored numbers.

### class_graph_gauss_seidel.f90

- `gauss_seidel` (type) — **INFERENCE**: A concrete colour-swept relaxation process (Gauss-Seidel/SOR) driving the residual to zero — a solution process, with omega as absorbed governance.

### class_graph_gmres.f90

- `gmres` (type) — **INFERENCE**: A concrete restarted Krylov solution process (GMRES/Arnoldi) minimizing the residual norm — inference machinery independent of any structural coupling.

### class_graph_jacobi.f90

- `jacobi` (type) — **INFERENCE**: A concrete solution process (damped Jacobi iteration) that drives the residual toward zero — pure epistemic-state movement, no structure or data of its own.

### class_graph_linearization.f90

- `difference_linearization` (type) — **OPERATOR**: The finite-difference tangent J v ~ (S(q+eps v) - S(q))/eps of a statement, wrapped as an operation: a linearization of residual machinery, hence material for R (its same-name generic interface is the constructor).

### class_graph_marcher.f90

- `marcher` (type) — **INFERENCE**: Time marching as forward inference (and its explicit adjoint as reverse inference): it moves the data q along the instants chain, governing an inner minimizer for implicit rules.
- `MARCH_FORWARD` (constant) — **INFERENCE**: Tag selecting the explicit rule of the marching inference process; it configures how the forward walk realizes states, not the operator itself.
- `MARCH_BACKWARD` (constant) — **INFERENCE**: Tag selecting the implicit backward-euler rule, i.e. one governed solve per edge of the marching process.
- `MARCH_BDF2` (constant) — **INFERENCE**: Tag selecting the second-order implicit rule of the marching process (implicit solve per edge, started by one backward step).

### class_graph_mesh.f90

- `mesh` (type) — **DATA**: Its own content is seven typed fields (volumes, centers, areas, deltas, normals, weights) — values over the cell/face carriers, candidate Q-material — riding an inherited stored_graph rather than being structure itself.
- `mesh (generic constructor interface -> create)` (interface) — **DATA**: Gates raw measurement arrays against the structure and seats them as fields on the graph's own carriers, producing the data-bearing mesh; the structural half is delegated to the parent's constructor.

### class_graph_multigrid.f90

- `multigrid` (type) — **INFERENCE**: A two-grid solution process (a minimizer governing a smoother and a coarse solve) that moves a problem from posed to answered; it is epistemic-state-moving machinery, not the operator itself.

### class_graph_newton.f90

- `newton` (type) — **INFERENCE**: A concrete nonlinear solution process that freezes a linearization and delegates each linear question to a governed inner minimizer — governance of inference, owning no operator mathematics itself.

### class_graph_partitioner.f90

- `partitioner` (type) — **STRUCTURE**: Concrete structural transform of topology (whole graph to part graph with ownership/index maps), acting on carriers and relations, not on values or residuals.
- `PARTITION_LINEAR` (constant) — **STRUCTURE**: Rule selector for how the structural cut of the member set is made; it parameterizes a topology transform.
- `PARTITION_BREADTH_FIRST` (constant) — **STRUCTURE**: Rule selector for a connectivity-following structural cut; parameterizes the same topology transform.
- `PARTITION_ADOPTED` (constant) — **STRUCTURE**: Rule selector saying the ownership map (a structural assignment of members to parts) is taken from outside.

### class_graph_reduction.f90

- `reduction` (type) — **OPERATOR**: Discrete integration/aggregation operator mapping a field (and optional measure) to one functional — a discretization action (sum, integral, norm, inner product), candidate observation-operator material for R.
- `REDUCE_SUM` (constant) — **OPERATOR**: Rule selector defining which aggregation operator the reduction applies.
- `REDUCE_AVERAGE` (constant) — **OPERATOR**: Rule selector for the measure-weighted average operator, finalized by one division.
- `REDUCE_MINIMUM` (constant) — **OPERATOR**: Rule selector for the ordering (minimum) aggregation operator over real fields.
- `REDUCE_MAXIMUM` (constant) — **OPERATOR**: Rule selector for the ordering (maximum) aggregation operator over real fields.
- `REDUCE_NORM` (constant) — **OPERATOR**: Rule selector for the p-norm operator, whose root is taken once at finalize.
- `REDUCE_COUNT` (constant) — **OPERATOR**: Rule selector for the entry-counting aggregation operator.
- `REDUCE_ALL` (constant) — **OPERATOR**: Rule selector for the logical conjunction reduction over a logical field.
- `REDUCE_ANY` (constant) — **OPERATOR**: Rule selector for the logical disjunction reduction over a logical field.
- `broadcast` (type) — **OPERATOR**: The transpose of a reduction — one functional value filled into a field (copy transposes a sum, share an average) — an adjoint operator action, not data itself.
- `BROADCAST_COPY` (constant) — **OPERATOR**: Rule selector making broadcast the transpose of a sum.
- `BROADCAST_SHARE` (constant) — **OPERATOR**: Rule selector making broadcast the transpose of an average, pinned by reduce(broadcast(J)) = J.

### class_graph_refiner.f90

- `refiner` (type) — **STRUCTURE**: Concrete structural transform of topology (each cell opened into a joined family of children, faces carried down); injection of values is transport along that structure, not an operator or inference.

### class_graph_stencil.f90

- `stencil_operator` (type) — **OPERATOR**: The compiled sparse linear operator (weights on a dependency pattern plus affine constants) — directly R material; the same-named generic interface is its constructor from (row, column, weight) triples.

### class_graph_step.f90

- `step_operator` (type) — **OPERATOR**: The temporal discretization stencil a0 q + a1 qold + a2 qolder + h S(q): it forms a discrete residual per instant, i.e. candidate material for R.
- `backward_euler` (procedure) — **OPERATOR**: Constructor stamping the reach-1 [1,-1] table onto a step residual operator; it manufactures operator material, not a solve.
- `bdf` (procedure) — **OPERATOR**: Constructor for the k-step BDF residual operator (tables for k=1,2); order rides as an argument on the same operator family.

### class_graph_walk.f90

- `walk` (type) — **VIEW**: Reads only the structure of the host graph and re-presents it as derived vocabulary (colours, BFS order, components, depths) — an interpretation over the structural schema, not an operator on values.
- `WALK_COLOURING` (constant) — **VIEW**: Rule selector for the greedy proper-colouring interpretation of the structure.
- `WALK_VISIT_ORDER` (constant) — **VIEW**: Rule selector for the breadth-first visit-order labelling, a derived reading of structure.
- `WALK_COMPONENT` (constant) — **VIEW**: Rule selector for connected-component labelling, purely a structural interpretation.
- `WALK_DEPTH` (constant) — **VIEW**: Rule selector for BFS distance-from-seed labelling, another derived structural vocabulary.

### class_harmonic_form.f90

- `harmonic_form` (type) — **OPERATOR**: The wave basis {1, sin(k.(x-at)), cos(k.(x-at))} at a fixed wavenumber, whose evaluation and chain-rule slopes let a fit differentiate waves exactly — concrete operator shape material for R; the public name also serves as its constructor generic interface.

### class_mesh.f90

- `mesh` (type) — **LEGACY**: The old-world mesh extending the abstract graph root of interface_graph, retained as the parsing-and-measuring machine behind the bridge (class_mesh_builder aliases it 'legacy'); semantically it would be DATA (geometry measurements over the structural cell graph) fused with the STRUCTURE it derives (face->cell, cell->face relations).
- `mesh (generic constructor interface -> create_mesh)` (interface) — **LEGACY**: Old-front-door constructor that wires and measures a mesh from a mesh_loader; kept alive for the bridge, and would semantically be DATA (producing measured values on a derived structure).

### class_mesh_builder.f90

- `mesh_from_gmsh` (procedure) — **DATA**: Constructs the tower's level-1 mesh — measurements seated as typed fields on a structural carrier (Q-material realized from a file) — by running the legacy machinery once at load; the symbol itself is the durable new-world entry point, so it is not classified LEGACY even though its current body is.

### class_polynomial_form.f90

- `polynomial_form` (type) — **OPERATOR**: The degree-one Taylor basis {1, x-at} with exact values and slopes — concrete shape material from which linear reconstruction/differentiation operators (candidate material for R) are assembled; the public name also serves as its constructor generic interface.

### class_robin_condition.f90

- `robin_condition` (type) — **DATA**: A boundary parameter payload — a tag and three numbers a, b, c of a*phi + b*dphi/dn = c — that computes per-face coefficient and wall-relation arrays; parameters and derived numbers (Q constituents) that only downstream assembly turns into R material.
- `robin` (procedure) — **DATA**: Constructor of the general robin_condition parameter object from a tag and the three coefficients.
- `dirichlet` (procedure) — **DATA**: Constructor specialization producing the a=1, b=0 parameter payload (held value).
- `neumann` (procedure) — **DATA**: Constructor specialization producing the a=0, b=1 parameter payload (held gradient/flux).

### class_stored_graph.f90

- `stored_graph` (type) — **LEGACY**: Concrete stored citizen of interface_graph's legacy graph root (edge-list, quotient, and refinement constructors); semantically STRUCTURE.
- `stored_digraph` (type) — **LEGACY**: Concrete stored citizen of the legacy digraph (directed edge list plus the condensation constructor); semantically STRUCTURE.

### class_string.f90

- `string` (type) — **DATA**: A self-contained character value with parsing (tokenize, asinteger, asreal) whose products become parameters and stored numbers; a plain utility datum not tied to any structural domain, and DATA is the nearest of the six categories (its same-name generic interface is the constructor).

### graph_algorithms.f90

- `sources` (procedure) — **VIEW**: Reads the interpretation's predecessor fibres and re-presents the empty-preimage members as a subset_set subobject; a derived-vocabulary query that creates no new topology and moves no epistemic state.
- `sinks` (procedure) — **VIEW**: The dual reading over successor fibres, likewise only re-presenting members of the existing domain as a subobject.
- `reachable` (procedure) — **VIEW**: Answers directed-path existence by breadth-first reading of the adjacency interpretation's fibres; a pure combinatorial question about structure, not an inference between epistemic states.
- `topological_order` (procedure) — **VIEW**: Derives a linear order consistent with the interpreted arrows (deterministic Kahn) or refuses on a cycle; it re-presents the interpretation's order structure and neither transforms topology nor realizes any Q or R.

### graph_binary_relation.f90

- `binary_relation` (type) — **STRUCTURE**: The abstract arity-two structural relation adding source/target/image/preimage — canonical binary-relation structure over carriers.
- `csr_relation` (type) — **STRUCTURE**: The indexed (CSR-backed) concretion of a binary structural relation — the same mathematical set of pairs with O(degree) fibres.
- `transposed_view` (type) — **STRUCTURE**: Mathematically the transpose relation R^T — a canonical slot permutation of a structural relation; the 'view' is only a borrowing mechanism, not an ordinary-graph interpretation profile, so STRUCTURE rather than VIEW.
- `transpose_of` (procedure) — **STRUCTURE**: Constructs the canonical slot permutation of a binary structural relation — a structural transform delivered lazily.
- `inclusion_of` (procedure) — **STRUCTURE**: Builds the inclusion relation I_S <= S x A of a subobject — the relational face of embedding, pure structural machinery.

### graph_calculus.f90

- `graph_functional` (type) — **LEGACY**: Legacy-compat abstract of the one-entry field, a single stored value J; semantically DATA (a functional's value, candidate constituent of Q).
- `graph_reduction` (type) — **LEGACY**: Legacy-compat abstract of the many-to-one map field->functional (sum, integral, norm, inner product with a measure seat); semantically OPERATOR.
- `graph_broadcast` (type) — **LEGACY**: Legacy-compat abstract of the transpose map functional->field (copy/share fills); semantically OPERATOR.
- `discretization_operator` (type) — **LEGACY**: Legacy-compat abstract that binds a continuous statement to a graph's arithmetic and owes its dependency pattern; semantically OPERATOR, direct material for the residual R.
- `linearization_operator` (type) — **LEGACY**: Legacy-compat abstract of the tangent J v at a frozen state, the linear question newton governs; semantically OPERATOR (a linearization of R).
- `graph_partitioner` (type) — **LEGACY**: Legacy-compat abstract of P, cutting the whole into parts with matching data maps; semantically STRUCTURE (a structural transform of topology).
- `graph_assembler` (type) — **LEGACY**: Legacy-compat abstract of P^-1, restoring whole-graph order collecting only owned values; semantically STRUCTURE.
- `graph_coarsener` (type) — **LEGACY**: Legacy-compat abstract of the coarsening transform (fewer, larger cells for a multigrid level); semantically STRUCTURE.
- `graph_refiner` (type) — **LEGACY**: Legacy-compat abstract of the refinement transform (one cell becomes several, one-sided law); semantically STRUCTURE.
- `GRAPH_SIDE_VERTEX` (constant) — **LEGACY**: Explicitly documented as the legacy differential operator's output-landing choice and nothing else; semantically OPERATOR metadata (where a discrete operator's answer lands).
- `GRAPH_SIDE_EDGE` (constant) — **LEGACY**: The edge-side twin of the legacy output-landing choice; semantically OPERATOR metadata.

### graph_carrier.f90

- `member_set` (type) — **STRUCTURE**: The abstract member set A = {a_1..a_n}, the ground carrier of the structural graph's S — a member set by definition.
- `counted_set` (type) — **STRUCTURE**: The concrete carrier enumerating members 1..n — a member set, hence a carrier of the structural tower.
- `subset_set` (type) — **STRUCTURE**: A subobject S embedded in an ambient carrier, itself a declared member set — carrier machinery with the embedding order.

### graph_field_calculus.f90

- `graph_field` (type) — **DATA**: The abstract contract for a function f : S -> V over one member-set domain — the archetypal constituent of Q; not LEGACY because this module is the rehomed owner and the old grammar now re-exports from here.
- `GRAPH_FIELD_INTEGER` (constant) — **DATA**: Value-kind tag naming which payload type a field's stored values carry — an attribute of Q's values, not of structure or operators.
- `GRAPH_FIELD_REAL` (constant) — **DATA**: Value-kind tag for real-valued field payloads, an attribute of the values in Q.
- `GRAPH_FIELD_COMPLEX` (constant) — **DATA**: Value-kind tag for complex-valued field payloads (the seat of complex-step perturbations), still an attribute of stored values.
- `GRAPH_FIELD_LOGICAL` (constant) — **DATA**: Value-kind tag for logical field payloads, an attribute of stored values.
- `GRAPH_FIELD_CHARACTER` (constant) — **DATA**: Value-kind tag for character field payloads (e.g. boundary names), an attribute of stored values.

### graph_fitting.f90

- `fit` (type) — **OPERATOR**: Despite the INFERENCE-flavoured name, a fit computes exact-reproduction weights for a directional-derivative functional over a point constellation (generalized finite-difference stencil generation, fed straight into stencil_operator entries), a discretization action whose internal CG minimization is purely instrumental; the same-named generic interface is its constructor.
- `form_optimizer` (type) — **INFERENCE**: Abstract governor that adjusts which basis members of a form stand active based on data (positions) and residual demand — model discovery/adaptation, machinery that changes the model rather than the topology.

### graph_forms.f90

- `form` (type) — **OPERATOR**: An abstract family of basis functions with evaluation (values) and directional-derivative (slopes) actions — the shape material discrete reconstruction and differentiation operators are built from, hence candidate material for R; its subset_set parentage and restrict only administer which basis members stand, so the STRUCTURE reading is secondary to what the symbol mathematically does.

### graph_grammar.f90

- `graph` (type) — **LEGACY**: Abstract root of the old four-root grammar answering ordinary vertex/edge structure; kept as the compatibility contract, semantically STRUCTURE.
- `graph_field` (reexport) — **LEGACY**: Re-export of the field abstraction owned by graph_field_calculus, kept here only for remaining legacy consumers; semantically DATA (values over a domain, candidate constituent of Q).
- `graph_operation` (type) — **LEGACY**: Abstract 'verb within' root of the old grammar (name/domain/apply); a generic data-to-data verb explicitly not automatically the residual, semantically OPERATOR-shaped.
- `graph_transform` (type) — **LEGACY**: Abstract 'verb between' root of the old grammar holding the two admissibility questions for partition/assemble/coarsen/refine; semantically STRUCTURE (transforms of topology).
- `GRAPH_FIELD_INTEGER` (reexport) — **LEGACY**: Re-exported value-kind tag of the legacy field's absorbed axis; semantically DATA (labels the payload kind of a field).
- `GRAPH_FIELD_REAL` (reexport) — **LEGACY**: Re-exported value-kind tag for real-valued field payloads; semantically DATA.
- `GRAPH_FIELD_COMPLEX` (reexport) — **LEGACY**: Re-exported value-kind tag for complex-step-derivative payloads; semantically DATA.
- `GRAPH_FIELD_LOGICAL` (reexport) — **LEGACY**: Re-exported value-kind tag for mask payloads; semantically DATA.
- `GRAPH_FIELD_CHARACTER` (reexport) — **LEGACY**: Re-exported value-kind tag for boundary/material name payloads; semantically DATA.

### graph_identity.f90

- `token` (type) — **STRUCTURE**: The opaque identity stamp that declares a structural domain or relation to BE itself — identity/token machinery, explicitly STRUCTURE.
- `mint_token` (procedure) — **STRUCTURE**: The sole minting act of the identity law — fresh, unrepeatable stamps for structural citizens — identity machinery.

### graph_minimization.f90

- `minimizer` (type) — **INFERENCE**: Abstract base of the residual-driving solver family — machinery that moves a statement from (Q?,R) toward a realized state; it merely wears the legacy graph_operation face for composition.

### graph_profile.f90

- `ordinary_graph_view` (type) — **VIEW**: An interpretation that reads the two-relation T/H schema off a relational graph and answers the ordinary directed-graph vocabulary (tails, heads, incidence, adjacency) with no storage of its own; constructor rides the same-named interface.
- `directed_adjacency_view` (type) — **VIEW**: An interpretation that reads one graph-owned same-domain binary relation as a directed adjacency (successors/predecessors), borrowing storage and computing nothing; constructor rides the same-named interface.

### graph_relation.f90

- `relation` (type) — **STRUCTURE**: The abstract finite-arity structural relation over a signature of carriers — the P of Gamma=(S,P) itself.
- `stored_relation` (type) — **STRUCTURE**: The materialized tuple-table concretion of the structural relation — a stored set of tuples over carrier domains.
- `slot` (type) — **STRUCTURE**: One signature seat holding a carrier concretion — schema machinery that types a relation's cartesian product.

### graph_relation_algebra.f90

- `restrict_slot` (procedure) — **STRUCTURE**: Restriction R|_S at one slot — a primitive of the relation algebra acting on structural relations and subobjects.
- `project_slots` (procedure) — **STRUCTURE**: Projection pi(R) onto chosen slots — relation-algebra primitive producing a new structural relation with the selected signature.
- `compose_binary` (procedure) — **STRUCTURE**: Binary composition S o R over a shared middle domain — the relation algebra's generation of structural relations from structural relations.

### graph_structure.f90

- `relational_graph` (type) — **STRUCTURE**: The structural container Gamma = (S, P) itself: owned member sets plus typed relations under the signature validity law, with identity; a same-named generic interface provides its constructor.
- `held_set` (type) — **STRUCTURE**: A polymorphic seat wrapping one carrier (member set) so heterogeneous carriers can sit in one array; pure carrier-storage machinery, with a same-named constructor interface.
- `held_relation` (type) — **STRUCTURE**: A polymorphic seat wrapping one structural relation for the graph's heterogeneous relation array; relation-storage machinery, with a same-named constructor interface.

### interface_graph.f90

- `graph` (type) — **LEGACY**: Abstract carrier of the old vertex/edge world (retained adjacency, traversal, coloring, partition bookkeeping, gather/scatter/dot), a compatibility contract of the pre-relational ontology; semantically STRUCTURE.
- `digraph` (type) — **LEGACY**: Abstract directed extension of the old graph carrying dependency order, strong components, condensation and adjoint accumulation; semantically STRUCTURE with INFERENCE-flavoured riders (accumulate_adjoint).
- `vertex` (type) — **LEGACY**: Plain member token of the old ontology (label plus owning part); semantically STRUCTURE (identity/token machinery).
- `edge` (type) — **LEGACY**: Plain tail/head pair token of the old ontology; semantically STRUCTURE (one instance of a structural relation).
- `counting_sort` (procedure) — **STRUCTURE**: A count/prefix-sum/scatter kernel that materializes CSR relation storage from (key,value) pairs; relation-storage machinery in its own right, not old vocabulary.
- `transpose_adjacency` (procedure) — **STRUCTURE**: Computes the converse of a stored bipartite adjacency (for each value, the keys touching it) - the transpose of a structural relation.
- `power_iteration` (procedure) — **INFERENCE**: Iterated matvec-and-dot estimation of a linear action's dominant stretch (spectral radius); a solution process that infers a number, not an operator or a structure.

### interface_mesh_loader.f90

- `mesh_loader` (type) — **STRUCTURE**: An abstract contract whose one deferred procedure (get_mesh_data) delivers the structural description of the mesh — member lists and incidence relations (edge/face/cell -> vertex) plus tag identity — with coordinates as the only value payload; its deliverable is a graph described by its edges, i.e. structure before any measurement.

### module_mesh_utils.f90

- `distance` (interface) — **DATA**: Generic exposing the two-point Euclidean distance — measurement arithmetic that produces geometric values (Q-feed material), not structure or an operator on fields. Note: this module has no private statement, so every symbol below is public by default rather than via a 'public ::' list.
- `distanceAB` (procedure) — **DATA**: The specific implementation behind the distance generic; computes a measured number from two coordinate points.
- `cross_product` (procedure) — **DATA**: Vector cross product used as the tape measure for areas and normals; produces geometric values rather than acting as a discrete residual operator or a structural transform.
- `find` (procedure) — **STRUCTURE**: Identity/token machinery: maps a global member number to its index in a member list (used to resolve sparse gmsh node tags).
- `is_subset` (procedure) — **STRUCTURE**: A membership (subset) test between vertex token sets — relation algebra on member sets.
- `isort` (procedure) — **STRUCTURE**: Orders integer member-token lists in support of the membership test; utility on member-set tokens, not on values.
- `elem_type_face_count` (procedure) — **STRUCTURE**: A per-element-type table of how many faces bound each cell shape — local incidence schema used to size the structural face set algebraically.
- `elem_type_dimension` (procedure) — **STRUCTURE**: Table giving the topological dimension of each gmsh element type; it classifies members into cells/faces/edges, a schema-level distinction.
- `elem_type_vertex_count` (procedure) — **STRUCTURE**: A per-element-type table of corner counts — the local incidence schema of each cell shape.
- `order_face_vertices` (procedure) — **STRUCTURE**: Rewinds a face's vertex tokens into the cell type's outward-wound wiring-table order — an orientation of an incidence relation (its geometric consequence, outward normals, is computed elsewhere).

### module_solve_mode.f90

- `FORWARD` (constant) — **INFERENCE**: Direction tag for the primal solve A x = b; it parametrizes an inference process, not the operator.
- `REVERSE` (constant) — **INFERENCE**: Direction tag for the adjoint/transpose solve A^T x = b, i.e. the reverse direction of the same inference.
- `WHOLE` (constant) — **OPERATOR**: Names the full discrete operator as the part a jacobian-vector product acts on; the referent is a portion of R's linearization.
- `DIAGONAL` (constant) — **OPERATOR**: Names the diagonal part of the discrete operator used by stationary splittings; an operator-decomposition tag.
- `LOWER_TRIANGLE` (constant) — **OPERATOR**: Names the strictly lower triangle of the discrete operator matrix; an operator-decomposition tag.
- `UPPER_TRIANGLE` (constant) — **OPERATOR**: Names the strictly upper triangle of the discrete operator matrix; an operator-decomposition tag.
- `is_valid_mode` (procedure) — **INFERENCE**: Entry validator for solve-direction tags; guards the inference machinery's configuration.
- `is_valid_part` (procedure) — **OPERATOR**: Entry validator for operator-part tags; guards which portion of the operator a product may act on.

### module_verbosity.f90

- `verbosity` (constant) — **VIEW**: A protected module variable (read-only knob outside its setter) governing how much diagnostic printing happens; pure re-presentation machinery, no mathematics.
- `set_verbosity` (procedure) — **VIEW**: The single writer of the diagnostics knob; configures presentation, touches no structure, data, operator, or inference.

## Ambiguities recorded by the auditors

- `class_advection.f90:2` — `! The advection law: the flow's answer to a state.`
  - (c) 'state' here means the physical transported state, not one of the four epistemic states (void/data/operator/realized) of the computational graph; the banner sentence uses the contested word unqualified.
- `class_diffusion_statement.f90:2` — `! The diffusion statement: the constitution speaks, an operator`
  - (e) 'operator' means the compiled stencil (residual-operator material), but in the legacy grammar the operation role (graph_operation) uses near-identical vocabulary; the two senses coexist in this file via the fit it indirectly employs.
- `class_diffusion_statement.f90:28` — `use graph_grammar        , only : graph`
  - (a) imports the LEGACY abstract graph root (unused in the body), inviting the reading that a statement acts on 'the graph' without saying which of the three graph meanings — while in the new ontology a statement contributes to the computational graph's R.
- `class_file.f90:3` — `! vertices, "next line" is the only edge, and the reader can do just`
  - Confusion (a): the header borrows the legacy ordinary-graph vertex/edge vocabulary as a metaphor for a text file; amid three live meanings of 'graph', a reader could take this chain for a structural relational_graph or a computational graph when it is neither.
- `class_fitted_balance.f90:5` — `! optimizers only: things that drive a residual or govern one that`
  - (b) 'residual' here is the level-2 fit misfit that form optimizers watch, not the reserved R of G=(Q,R); the sentence defines a level boundary using the contested word.
- `class_fitted_balance.f90:141` — `! Algebra: one apply, aimed along the normal at the face.`
  - (e) 'apply' is the legacy graph_operation verb, yet the thing being applied (the fit) is producing entries of a genuine residual operator; the legacy role vocabulary and OPERATOR sense cross in one line.
- `class_form_pruner.f90:7` — `! runs. What used to hide inside a solve as a pivot trick is now a`
  - (d) 'a solve' means an internal linear-algebra pivot situation, not an epistemic move to a realized state; sitting inside inference-family machinery the word invites the stronger reading.
- `class_graph.f90:63` — `public :: stored_graph`
  - The name stored_graph is also the public concrete of class_stored_graph on the other legacy branch; one name denotes two unrelated types extending two different abstract graphs, so a reader cannot tell which structural world a 'stored_graph' inhabits.
- `class_graph.f90:48` — `! the graph accumulates state, and its answers come to depend on the`
  - 'state' here is accumulated mutable object/solver state, a sense distinct from the four epistemic states; the sentence argues immutability but the word overlaps the new vocabulary.
- `class_graph.f90:35` — `! Colouring, traversal order, partitioning and the rest are operations`
  - 'operations' is used in the legacy graph_operation sense for things (colouring, partitioning) that in the new taxonomy are VIEW/STRUCTURE work, inviting conflation with the residual-operator sense of 'operator'.
- `class_graph_assembler.f90:4` — `! P inverse, and only that. It puts a piece back into whole-graph`
  - P here means the partition transform, colliding with P as the structural-relations component of Gamma = (S, P).
- `class_graph_assembler.f90:16` — `!         assemble( partition( G ) )     ==  G`
  - G here denotes the structural graph being cut and reassembled, but G is now the reserved letter for the computational graph G = (Q, R); the round-trip law reads as if it were about the computational object.
- `class_graph_balance.f90:20` — `! up as a crash - it shows up as a solution that is quietly not the`
  - Confusable (d): 'solution' names the computed balance value; under the new ontology this is merely a realized quantity, and 'solution' smuggles in the solver-stage sense the header itself elsewhere disclaims.
- `class_graph_balance.f90:39` — `! A solver calls the result a residual. That word names a stage in a`
  - Confusable (b): the balance IS the natural residual, and the letter R is now reserved for exactly that; this comment banishes 'residual' from the library at the moment the ontology promotes it to a first-class object, a direct terminology collision.
- `class_graph_balance.f90:67` — `type, extends(graph_operation) :: balance`
  - Confusable (e): the residual-candidate balance is typed through the LEGACY generic graph_operation root, conflating the old grammar role with the residual-operator sense of 'operator'.
- `class_graph_coarsener.f90:74` — `! Add the fine values or average them. A residual adds, because`
  - "residual" here is a numeric field being restricted, easily read as the operator R of G = (Q, R); the word names data in a file about structure.
- `class_graph_coarsener.f90:75` — `! it is a total. A state averages, because it is a level.`
  - "state" means a physical/solver state variable, not one of the four epistemic states (void/data/operator/realized) the repo now also calls states.
- `class_graph_coarsener.f90:241` — `! this to build its coarse operator; the map is the coarsener's to`
  - "operator" here means the residual/discrete operator (Galerkin coarse operator), but in this repo "operator" also echoes the legacy graph_operation root; the sentence sits in a structural module and could be read either way.
- `class_graph_coarsener.f90:306` — `! Getting that choice wrong is quiet: a residual that gets averaged`
  - Again "residual" as a numeric vector being coarsened, confusable with R the residual operator of the computational graph.
- `class_graph_conjugate_gradient.f90:4` — `! For a symmetric positive operator, the best correction in every`
  - (e) 'operator' here means the attached residual operator (the matvec action), but the same word names the legacy graph_operation contract the action is typed as.
- `class_graph_conjugate_gradient.f90:22` — `module class_graph_conjugate_gradient`
  - (a) The 'graph' prefix reads as structural/legacy graph vocabulary while the module is inference on the computational graph.
- `class_graph_differential_operator.f90:2` — `! Differential operators on a graph.`
  - Confusable (a): 'graph' here is the legacy abstract structural host, but in a module that manufactures candidate R material the bare word invites reading it as the computational graph G=(Q,R).
- `class_graph_differential_operator.f90:219` — `type, extends(graph_operation) :: differential_operator`
  - Confusable (e): a genuine residual-operator candidate is spelled through the LEGACY generic root graph_operation, blurring 'operation' (old grammar role) with 'operator' (R material).
- `class_graph_differential_operator.f90:745` — `!      three parts of a velocity, the five of a conserved state. The`
  - Confusable (c): 'state' here means the physical conserved-variable vector, which a reader of the new ontology could misread as one of the four epistemic states of G.
- `class_graph_differential_operator.f90:1147` — `select type (state => input_data(1))`
  - Confusable (c): the associate-name 'state' labels an incoming data field (candidate Q constituent), not an epistemic state; the same choice recurs at line 1182.
- `class_graph_field.f90:60` — `module class_graph_field`
  - Every consumer writes 'use class_graph_field' although the type inside is the graph-free 'field' whose header insists it needs only a domain; the module name keeps 'graph' alive and is readable as either structural-graph attachment or computational-graph (Q) membership.
- `class_graph_functional.f90:26` — `! Logical is here so a question such as "is this graph acyclic" comes`
  - 'this graph' means the structural/ordinary graph being interrogated, but a functional is itself a Q constituent of the computational graph G=(Q,R), so a new-ontology reader can misattach the question to G.
- `class_graph_functional.f90:39` — `use graph_grammar      , only : graph`
  - Imports the LEGACY abstract graph root of the old four-root grammar (apparently unused in the body) into a data module; the bare name 'graph' is confusable with both the structural relational_graph and the computational graph G.
- `class_graph_gauss_seidel.f90:4` — `! Jacobi corrects every cell from the old state; gauss-seidel lets`
  - (c) 'old state' means the previous iterate, not an epistemic state of the computational graph.
- `class_graph_gauss_seidel.f90:5` — `! each correction see the ones already made. On a graph the safe`
  - (a) 'On a graph' is ambiguous: the colouring is of the structural coupling graph, not of the execution host nor the computational graph G=(Q,R).
- `class_graph_gauss_seidel.f90:19` — `module class_graph_gauss_seidel`
  - (a) The 'graph' prefix reads as structural/legacy graph vocabulary while the module is inference on the computational graph.
- `class_graph_gmres.f90:4` — `! For any operator, symmetric or not: build an orthonormal basis of`
  - (e) 'operator' means the residual operator behind matvec, but that seat is typed as a legacy graph_operation — the two senses share the word.
- `class_graph_gmres.f90:22` — `module class_graph_gmres`
  - (a) The 'graph' prefix reads as structural/legacy graph vocabulary while the module is inference on the computational graph.
- `class_graph_jacobi.f90:16` — `module class_graph_jacobi`
  - (a) The 'graph' prefix suggests the structural relational_graph or legacy graph vocabulary, but the content is inference over the computational statement G=(Q,R).
- `class_graph_linearization.f90:7` — `! standing state q, computed as a difference,`
  - Confusion (c): 'standing state' (and 'frozen state' at line 165, freeze's 'standing state' at line 82) is the physical linearization point q, not an epistemic state of the computational graph.
- `class_graph_linearization.f90:126` — `! THE DIFFERENCE IS TAKEN WHERE THE OPERATION LIVES.`
  - Confusion (e): 'the operation' is the legacy graph_operation role of the differenced statement, yet mathematically the object is the residual operator whose tangent is being formed; the word does double duty across ontologies.
- `class_graph_marcher.f90:5` — `! here - it is one more graph: instants are vertices, steps are`
  - Confusion (a): 'graph' here is the structural chain of instants spoken in legacy vertex/edge vocabulary; a reader adopting G=(Q,R) could read 'time is one more graph' as the computational graph, which this module never builds.
- `class_graph_marcher.f90:12` — `! and the state moves against it. Three rules, absorbed:`
  - Confusion (c): 'state' throughout this module (state, state_domain, state_seat, 'THE EVOLVING STATE' at line 229) means the marched solution vector q, i.e. solver/physical state, not one of the four epistemic states of G.
- `class_graph_marcher.f90:24` — `! it one step operator per edge - the level-1 citizen from the time`
  - Confusion (e): the 'step operator' is mathematically residual material for R, but it is held and passed through the legacy graph_operation role (the 'action' arguments), so 'operator/operation' straddles the old role and the new residual sense.
- `class_graph_mesh.f90:6` — `! This is the tower's one inheritance crossing: the mesh IS a graph`
  - Bare 'graph' in a new-tower header: it means the structural stored graph, but with G=(Q,R) being introduced a reader could take a data-carrying mesh that 'IS a graph' for the computational graph (confusion a).
- `class_graph_mesh.f90:41` — `module class_graph_mesh`
  - The public name 'graph_mesh' reads equally as 'mesh over the computational graph' and 'mesh extending the stored (structural) graph'; only the latter is meant (confusion a).
- `class_graph_mesh.f90:30` — `! one feeds; an operator receives those numbers at construction and`
  - 'operator' here is the residual/differential operator (R material), but the legacy graph_operation root still owns that word elsewhere in the repo (confusion e).
- `class_graph_multigrid.f90:8` — `! down, so the residual travels to the blocks, is answered there,`
  - "residual" is the numeric defect vector r = rhs - Ax, confusable with R the residual operator of G = (Q, R).
- `class_graph_multigrid.f90:22` — `! already accumulates. That is R A P with summing restriction and`
  - In R A P, R means restriction and P prolongation - both letters collide with the ontology's reserved meanings (R as residual/structural relation, P as structural relations of Gamma).
- `class_graph_multigrid.f90:25` — `!      solve_coarse( R(A(P e)) ) = e        the coarsened statement`
  - Same R and P collision inside the commutation law; the equation also uses "solve" for what is a realization of the coarse statement, blurring solved vs realized.
- `class_graph_multigrid.f90:17` — `! THE GALERKIN ROAD. The coarse operator is not re-derived: it is`
  - "operator" here is the residual-operator sense (a stencil_operator), but the word also names the legacy graph_operation root; the migration reserves OPERATOR for R-material and the prose does not mark which is meant.
- `class_graph_newton.f90:16` — `! citizen - the linearization operator - and newton merely governs:`
  - (e) 'the linearization operator' is genuinely OPERATOR-category material (candidate for R), yet it is seated and passed around as a legacy graph_operation — the word covers both.
- `class_graph_newton.f90:18` — `! freeze the linearization at the standing state, hand the linear`
  - (c) 'standing state' means the current iterate q, easily confused with the epistemic states of the computational graph.
- `class_graph_newton.f90:25` — `module class_graph_newton`
  - (a) The 'graph' prefix reads as structural/legacy graph vocabulary while the module is inference on the computational graph.
- `class_graph_partitioner.f90:4` — `! P cuts a graph into parts. One concrete type holds the rule, so a`
  - The letter P names the partition transform, but in the two-object ontology P is the structural-relations component of Gamma = (S, P); a reader can mistake the transform for the relation set.
- `class_graph_partitioner.f90:192` — `! P. Work out who owns what, gather this part's cells, and rebuild`
  - Same single-letter P collision: here P is the partition operator, not the P of Gamma = (S, P).
- `class_graph_reduction.f90:212` — `class(graph_functional), allocatable, intent(inout)     :: state`
  - Confusable (c): 'state' names the running accumulator payload (a tally-and-weight data value) throughout initialize/accumulate/combine/finalize, not an epistemic state of G=(Q,R).
- `class_graph_reduction.f90:139` — `! The operation face: the field reduced over its own domain,`
  - Confusable (e): 'operation face' invokes the LEGACY graph_operation vocabulary for what is mathematically a functional-evaluation operator, blurring the old grammar role with the R-sense of operator.
- `class_graph_refiner.f90:4` — `! R is the other way from C. One cell becomes several:`
  - R names the refine transform, directly colliding with the RESERVED letter R (historically structural relations, now the residual of G = (Q, R)).
- `class_graph_refiner.f90:128` — `! R. Split every cell, join each family together, and carry every`
  - Same collision: a bare "R." heading reads as the residual or the relation set, but means the refiner.
- `class_graph_stencil.f90:3` — `! A sparse matrix is a graph with a number on every edge. This type`
  - Confusable (a): 'graph' means the stored structural dependency pattern inside a compiled operator; a reader of the new ontology could take 'a matrix is a graph' as a claim about the computational graph G=(Q,R).
- `class_graph_stencil.f90:165` — `! The contract's answer: the pattern IS a graph, handed out whole.`
  - Confusable (a): the pattern handed out is the LEGACY abstract class(graph) — an ordinary-graph/structural object — yet the emphatic 'IS a graph' says nothing about which of the three graph senses is meant.
- `class_graph_step.f90:128` — `!      R_n = a0 q_n + a1 q_(n-1) + a2 q_(n-2) + h S(q_n)`
  - Confusion (b): R_n here is the step residual — consistent with the computational-graph R, but the letter R is being reserved after its historical use for structural relations, so an unlabelled R in a comment is exactly the collision the reservation exists to prevent (and S(q) simultaneously reuses the letter of Gamma's member sets S).
- `class_graph_step.f90:48` — `class(graph_operation), allocatable :: action`
  - Confusion (e): the held 'action' is mathematically the residual operator S being discretized, but its type is the legacy graph_operation root, so the operator-vs-operation vocabulary crosses old and new ontologies inside an OPERATOR type.
- `class_graph_step.f90:219` — `! the answer lands on. The state's width is read from the state -`
  - Confusion (c): 'state' (also in the line-248 refusal message) means the physical unknown vector q, not an epistemic state of G=(Q,R).
- `class_graph_walk.f90:6` — `! every cell. That makes every walk a vertex field operation: the`
  - Confusable (e): a walk is mathematically a structural view, yet the comment classifies it as an 'operation' — the LEGACY graph_operation role — inviting confusion with residual-operator material.
- `class_graph_walk.f90:68` — `type, extends(graph_operation) :: walk`
  - Confusable (e): the type inherits the LEGACY generic operation root, so a VIEW wears OPERATOR-flavoured vocabulary; its semantic category remains VIEW.
- `class_mesh.f90:2` — `! The mesh IS a graph that happens to live in space. Cells are the`
  - Bare 'graph' at the top of the file: it means the ordinary/structural cell graph, but the word now has three meanings in the repo and this header does not say which (confusion a).
- `class_mesh.f90:83` — `type, extends(graph) :: mesh   ! The mesh IS the cell graph: cells are the graph's vertices and interior faces are its edges.`
  - The public type extends a root literally named 'graph' — the LEGACY abstract graph of interface_graph — while the comment says 'cell graph'; three graph senses collide at one declaration (confusion a).
- `class_mesh.f90:479` — `! interior face. This adjacency is also the nonzero pattern of the`
  - The sentence continues '...Jacobian': it identifies a structural relation (adjacency, historically R, now P) with the sparsity of the residual's linearization — the exact seam where structural-relation R and residual R could be conflated (confusions a and b).
- `class_mesh.f90:1079` — `! The length of an edge as the diffusion operator feels it: project`
  - 'operator' means the discrete residual/differential operator (candidate R), not the legacy graph_operation role that still owns the word as a class name (confusion e).
- `class_mesh.f90:1265` — `! ends - each cell watches flux leave through its own side of the`
  - The next line calls the result 'the assembled operator' — again the residual operator sense of 'operator', in a file whose inheritance vocabulary comes from the old four-root grammar where graph_operation means something else (confusion e).
- `class_robin_condition.f90:33` — `! THE OPERATOR ROAD, for a > 0: the diffusive part rides the`
  - (e) 'OPERATOR' names the stencil/residual-operator route through the calculus, while the legacy grammar's operation role (graph_operation) shares the vocabulary; a reader could take the road as a graph_operation pathway.
- `class_robin_condition.f90:78` — `procedure :: operator_coefficients`
  - (e) public naming choice: 'operator' in this method name means coefficients destined for the residual operator, not anything to do with the legacy graph_operation role — the word carries both senses in this repo.
- `class_stored_graph.f90:35` — `public :: stored_graph`
  - Identical public type name to class_graph's stored_graph on the other legacy branch; the same word names two unrelated concrete graphs, and neither is the relational_graph or the computational G.
- `class_stored_graph.f90:30` — `use interface_graph, only : graph, digraph`
  - The bare name 'graph' is imported from the legacy interface while a different abstract 'graph' lives in graph_grammar; inside this file only the use statement discloses which structural branch's graph is in scope.
- `graph_algorithms.f90:10` — `! caller - the calculator's dependency walk - has earned:`
  - A calculator's dependency walk is precisely the evaluation-scheduling concern of the coming computational graph G=(Q,R), so a reader can take these algorithms as computational-graph machinery when they are views over the structural interpretation.
- `graph_algorithms.f90:227` — `error stop 'graph_algorithms: a topological order needs an acyclic graph'`
  - Bare 'graph' in a user-facing refusal, in the one routine a computational-graph evaluator would call first; with three meanings of 'graph' in the repo the message does not say it refuses a cyclic INTERPRETATION (ordinary-graph view), not a cyclic G=(Q,R).
- `graph_binary_relation.f90:6` — `!      R  <=  A x B`
  - The letter R names a structural binary relation while R is reserved for the computational residual; the clash repeats at line 88 with '(R^T)^T = R'.
- `graph_binary_relation.f90:74` — `!      the graph OWNS stable relations;`
  - The unqualified 'the graph' in this forward-declared ownership policy is the coming structural relational_graph, but a reader could read it as the computational graph G=(Q,R) owning its R.
- `graph_binary_relation.f90:153` — `! The CSR citizen: both directions materialized once, every query`
  - 'Materialized' here means precomputed index storage, yet it neighbors the computational graph's 'realized' epistemic state — a reader could hear an epistemic-state transition where only storage is meant.
- `graph_calculus.f90:332` — `! GRAPH_REFINER. R: the other way, one cell becomes several. The`
  - The refiner is named with the bare letter R, which is being RESERVED for structural relations (historically) and now denotes the residual in G=(Q,R); this is a third competing reading of R.
- `graph_calculus.f90:59` — `!    exact        assemble( partition( G ) ) = G       both ways`
  - G names the ordinary structural graph in the round-trip laws while G is reserved for the computational graph (Q,R); the laws could be misread as statements about computational realization.
- `graph_calculus.f90:357` — `! The reduction's four steps and its one-call form. The state is`
  - 'state' here is the reduction's running accumulator (an unfinished functional) - a third sense of the word beside epistemic state and solver/physical state.
- `graph_calculus.f90:21` — `!                by freezing a state]-----> linearization_operator`
  - 'freezing a state' means the physical/solver state at the linearization point, not an epistemic state; near the new vocabulary a reader could hear a transition between epistemic states.
- `graph_calculus.f90:38` — `! THE DERIVED OPERATORS. Two families of operations are built FROM`
  - 'operators' and 'operations' interleave in one sentence: the derived types are genuine residual-operator material (R) yet are typed as legacy graph_operations, blurring confusion class (e) in both directions.
- `graph_field_calculus.f90:12` — `! field needs a domain; it does not need a graph container.`
  - Here 'graph' means the structural relational_graph, but fields are the canonical constituents of Q in the computational graph G=(Q,R); a reader of the new ontology could take 'does not need a graph' as denying the field's membership in G rather than its detachment from Gamma.
- `graph_field_calculus.f90:43` — `public :: graph_field`
  - The 'graph_' prefix on a domain-only value carrier suggests attachment to a graph — structural or computational, the name does not say which — and it is the same word the legacy graph_grammar root vocabulary uses, though this symbol is the rehomed non-legacy owner.
- `graph_fitting.f90:12` — `!                         on the form's span. A fit is an OPERATION`
  - (e) OPERATION here means the legacy graph_operation role, yet the fit's product is residual-operator (stencil) material; the legacy role word and OPERATOR status are easy to conflate.
- `graph_fitting.f90:13` — `!                         on the point set, and its solve is a`
  - (d) 'its solve' is an instrumental least-squares minimization, not a transition of the computational graph toward a realized state; 'solve'/'solver' vocabulary makes the fit read like a solution process.
- `graph_fitting.f90:20` — `!                         admitting what the residual demands. A`
  - (b) 'residual' means the fit's misfit, not the reserved R (the residual operator of G=(Q,R)); a reader could think the form optimizer consults the computational graph's R.
- `graph_fitting.f90:30` — `!      B(m,j) = basis_m(x_j)        r(m) = scale * d(basis_m)/dn |at`
  - (b) the symbol r names the moment right-hand-side vector while the letter R is being reserved (historically structural relations, now the residual); the doctrinal formula uses the contested letter for a third thing.
- `graph_fitting.f90:149` — `class(graph), intent(in)                       :: input_graph`
  - (a) 'graph' is the LEGACY abstract root with ordinary-graph vocabulary, but paired with input_data(:) it reads like the computational graph G with its Q, and could equally be mistaken for the structural relational_graph.
- `graph_forms.f90:9` — `!      graph -> graph_support -> support -> form -> polynomial | wave`
  - The root 'graph' in this lineage is the LEGACY abstract graph of graph_grammar, but the unqualified word is readable as the structural relational_graph or the computational graph G=(Q,R).
- `graph_forms.f90:14` — `! member set - a graph act, owned by the form sector one level up.`
  - 'a graph act' means a structural transform (building a smaller member set), but a reader could take 'graph' as the computational graph and misread pruning as an epistemic move on G.
- `graph_forms.f90:37` — `! operator holds coefficients - as data about shape.`
  - 'operator' here means the residual/differential operator (R material) holding its coefficient data, but the legacy role graph_operation also claims the word, so the sentence is confusable between the two operator senses.
- `graph_grammar.f90:361` — `!        G'  =  P^-1 ( A ( P ( G ) ) )`
  - The central law names the ordinary structural graph G, but G is now reserved for the computational graph G=(Q,R); a reader could take this law as acting on (Q,R) rather than on topology.
- `graph_grammar.f90:72` — `!    graph_operation .. verb within   how data becomes other data`
  - The generic legacy verb role reads as 'operator' vocabulary; a reader could conflate graph_operation with the residual/operator R of the computational graph, which the migration context explicitly warns against.
- `graph_grammar.f90:170` — `! the ordinary state, complex for a complex-step derivative,`
  - 'state' here means the physical/solver state a real field carries, not one of the four epistemic states (void/data/operator/realized) of the new ontology.
- `graph_minimization.f90:5` — `! residual map R : U -> Y, vary the values on the UNKNOWN domain U`
  - (b) The letter R names the residual map here while R historically meant structural relations; a reader must know which reservation is in force on this line.
- `graph_minimization.f90:8` — `! vertices; the graph argument survives only as the legacy`
  - (a)/(e) 'the graph argument' is the legacy abstract graph host, not the structural relational_graph nor the computational G=(Q,R); three graph meanings collide in one clause.
- `graph_minimization.f90:29` — `!                       at zero state; the assembled right hand side`
  - (c) 'zero state' means the zero-valued iterate, but in a repo defining epistemic states (void = bot) it can be misread as the void epistemic state.
- `graph_minimization.f90:126` — `! The operation face: a solver IS an operation - the one that`
  - (e) 'operation' here is the legacy graph_operation contract, yet the attached action a solver drives is the residual operator R — two senses of operator in one banner.
- `graph_minimization.f90:207` — `! earn R^n -> R^m later; it has not yet.`
  - (b) A third use of the letter R (real coordinate spaces) alongside the reserved structural R and the residual R.
- `graph_minimization.f90:399` — `! The solver's answer is a solution on U.`
  - (d) Calls the output a 'solution', but solver_apply returns x = 0 unsolved when no input_data is present — the field is merely realized, not solved.
- `graph_minimization.f90:418` — `! x IS a state on the unknown domain; say so.`
  - (c) 'state' means the solver iterate/field of values, easily confused with the epistemic states of the computational graph.
- `graph_minimization.f90:435` — `out = field('solution', this % unknown_domain, ncomp=this % ncomp)`
  - (d) The output field is named 'solution' unconditionally, including the branch where solve() was never called and x is all zeros — realized data labelled as solved.
- `graph_profile.f90:153` — `class(relational_graph), target, intent(in) :: g`
  - The dummy argument (keyword-addressable, hence a public naming choice, repeated at line 468) names the STRUCTURAL graph with the single letter g, while G is being reserved for the computational graph G=(Q,R); call sites reading create_view(g, ...) can mistake which of the two objects the profile interprets.
- `graph_relation.f90:7` — `!      R  <=  A_1 x A_2 x ... x A_k ,        k >= 1`
  - Uses the letter R for a structural relation while R is now reserved for the residual/operator of the computational graph G=(Q,R); the same clash recurs at lines 24-25 and 31 (R_CC, R_CF, R_end).
- `graph_relation.f90:12` — `! arrives on a higher level, will CONTAIN relations; a relation`
  - Unqualified 'a graph ... will CONTAIN relations' means the structural relational_graph, but a reader carrying G=(Q,R) — a graph that also contains an R — could take it as the computational graph.
- `graph_relation.f90:122` — `procedure :: materialized`
  - The public naming choice 'materialized' (self-contained, safe-to-own storage) is one synonym away from the computational graph's 'realized' epistemic state, inviting the realized-vs-merely-stored confusion.
- `graph_relation_algebra.f90:10` — `!      restrict_slot     R|_S       keep the tuples whose i-th part`
  - R denotes the structural relation being restricted, colliding with the reserved residual letter of G=(Q,R); the whole primitive table (lines 10-19) and line 196's 'compose_binary(R_AB, R_BC) = R_BC o R_AB' share the clash.
- `graph_relation_algebra.f90:24` — `!      R_flow restricted to the output port, projected to O x X,`
  - 'R_flow' names a structural relation with the reserved letter, and next to 'flow' a skimming reader could take it for a flow residual/operator rather than a relation of the calculator schema.
- `graph_relation_algebra.f90:34` — `! Every result here is MATERIALIZED, and says so: restriction and`
  - 'MATERIALIZED' (eagerly stored result) sits one synonym from the computational graph's 'realized' epistemic state, risking a stored-vs-realized conflation.
- `graph_structure.f90:6` — `!      G = ( S, R )`
  - The header banner writes the STRUCTURAL graph with both reserved symbols at once: G is now the computational graph G=(Q,R) and R is being reserved for the residual, so a reader can take this line as the computational object; under the new ontology this line should read Gamma = (S, P).
- `graph_structure.f90:152` — `!     the same relation twice         R_i /= R_j, by identity`
  - R_i, R_j denote structural relations, but the letter R is reserved for the residual/operator of G=(Q,R); a reader scanning for residual machinery could misread this uniqueness law.
- `interface_graph.f90:71` — `public :: graph, digraph, vertex, edge`
  - A second public abstract type named 'graph', unrelated to graph_grammar's graph, the relational_graph, or the computational G=(Q,R); the bare name cannot tell a reader which of the repo's three graph meanings is in scope.
- `interface_graph.f90:66` — `use module_solve_mode, only : FORWARD, REVERSE`
  - Visit-order direction constants arrive from a module named 'solve_mode', dressing a purely structural traversal axis in solver vocabulary; a reader could take FORWARD/REVERSE as solution passes rather than visit directions.
- `interface_graph.f90:34` — `! graph-shaped data.`
  - 'A partition of a graph is still graph-shaped data' describes structural bookkeeping with the word 'data', which in the new ontology names Q constituents; partition records are STRUCTURE, not candidate Q.
- `interface_graph.f90:2235` — `! Which components are knots - the ones a walk must SOLVE at, not`
  - 'SOLVE' is solver-process (inference) vocabulary stated as a property of structure; near the new epistemic vocabulary it could be misread against 'solved'/'realized' as an epistemic transition.
- `interface_graph.f90:2333` — `! the caller with access to the required system state.`
  - 'system state' means the physical/solver state edge_apply reads to evaluate derivatives, not one of the four epistemic states of the computational graph.
- `interface_mesh_loader.f90:5` — `! physical tags. That is a graph described by its edges, delivered`
  - Bare 'a graph' in the contract header: it means the structural incidence description, but in a repo where 'graph' also names the computational G=(Q,R) and a legacy abstract root, the sense is unqualified (confusion a).
- `module_solve_mode.f90:3` — `! integrator and the graph.`
  - Confusion (a): 'the graph' is unqualified among the repo's three graph meanings; here it means the legacy abstract graph the solvers act through, but a reader could take it as the structural relational_graph or the computational G=(Q,R).
