# Releveling report — the tower under G=(Q,R)

Status: **proposal only**. The naming pass (this branch, see
`COMPUTATIONAL-GRAPH.md`) reserved the vocabulary and renumbered nothing.
This report is the separate deliverable it owes: for every level-bearing
module, the current claim, the proposed seat under the computational-graph
architecture, the reason, and the dependencies a move would touch.
Naming and releveling are separate commits; this document belongs to the
naming commit, the moves belong to the next one.

## The problem being solved

Two numbering systems currently coexist in module headers:

```text
"LEVEL n OF THE NEW TOWER"        the relation-centered tower (AGENTS.md 2)
"LEVEL n OF THE STRATIFICATION"   the older physics stratification
```

A reader meeting `graph_fitting` ("LEVEL 2 OF THE STRATIFICATION") beside
`graph_relation_algebra` ("LEVEL 2 OF THE NEW TOWER") sees one number
meaning two unrelated altitudes. The releveling commit must retire the
old stratification labels entirely and seat every citizen in one tower.

## The proposed tower

The computational graph earns a seat by the placement question of
AGENTS.md 3: *what becomes true here that was not true below?* Below it,
data exists (fields, level 5) and operators exist (discretizations,
level 6), but no object states **which knowledge is realized and which is
bottom**. The epistemic pair is that statement, and inference is defined
as movement between its states — so the pair must sit above both of its
constituents' levels and below everything that solves, fits or marches.

```text
level 0   carriers                              (unchanged)
level 1   relations                             (unchanged)
level 2   relation algebra                      (unchanged)
level 3   relational graph      Gamma = (S, P)  (unchanged)
level 4   graph calculus        views, walks    (unchanged)
level 5   field calculus        Q's material    (unchanged)
level 6   discretization        R's material    (unchanged)
level 7   computational graph   G = (Q, R)      <- inserted
level 8   inference             solve | discover | march
level 9   constitution                          (was 8)
level 10  statement                             (was 9)
```

Level 8 unifies what the old stratification scattered: minimization
(solution, \((\bot,R)\to(Q,R)\)), fitting (operator discovery,
\((Q,\bot)\to(Q,R)\)) and time marching (forward inference) are one
mathematical stratum — processes between epistemic states — with three
canonical directions, not three levels.

## The table

| Module | Current claim | Proposed | Reason | Dependencies affected |
|---|---|---|---|---|
| `graph_identity` | not a level (service) | unchanged | identity is a law, not an altitude | none |
| `graph_carrier` | 0, new tower | 0 | unchanged | none |
| `graph_relation` | 1, new tower | 1 | unchanged | none |
| `graph_binary_relation` | 1, new tower | 1 | unchanged | none |
| `graph_relation_algebra` | 2, new tower | 2 | unchanged | none |
| `graph_structure` | 3, new tower | 3 | unchanged; now written \(\Gamma=(\mathcal S,\mathcal P)\) | none |
| `graph_profile` | 4, new tower | 4 | interpretation of structural schemas; VIEW, not structure | `graph_algorithms` |
| `graph_algorithms` | 4, new tower | 4 | unchanged | none |
| `graph_field_calculus` | 5, new tower | 5 | unchanged; reframed as Q's raw material — a field is a candidate *constituent* of Q, never Q itself | `class_graph_field`, `graph_grammar` |
| `graph_state` | not yet a level (naming pass) | **7** | the epistemic pair; today it imports only level 3, but its seats will hold level-5 data and level-6 operators, and level-8 inference consumes it — the seat is semantic, taken only when Q/R storage decisions land | none yet (no consumers) |
| `graph_minimization` | 7, new tower | **8** | the solution direction of inference, \((\bot,R)\to(Q,R)\); shifts up one to sit above the pair it acts on | `class_graph_jacobi`, `class_graph_conjugate_gradient`, `class_graph_gauss_seidel`, `class_graph_gmres`, `class_graph_newton`, `class_graph_multigrid`, `class_graph_marcher` |
| `graph_fitting` | 2, stratification | **6 / 8, split** | the audit split the module along its own doctrine: `fit` is OPERATOR — deterministic stencil-weight generation (moment matrix, dual solve), level-6 material — while `form_optimizer` GOVERNS forms and is model discovery, \((Q,\bot)\to(Q,R)\), level-8 inference; the releveling commit should consider separating the two sectors | `class_form_pruner`, `class_fitted_balance` |
| `class_graph_linearization` | 1, stratification | **6** | \(\partial R/\partial Q\) is operator machinery: a discrete operator derived from a discrete operator | `class_graph_newton` |
| `class_graph_step` | 1, stratification | **6** | one application of a time-recurrence operator — a discrete operator, not a process | `class_graph_marcher` |
| `class_graph_marcher` | 2, stratification | **8** | forward inference in time: realizes Q from R, state by state | none (leaf) |
| `class_graph_mesh` | 1, stratification | **5** | metric measurements are fields over \(\Gamma\)'s carriers — geometric data, Q-like material | `class_conduction`, `class_advection`, `class_diffusion_statement`, `class_robin_condition`, `class_fitted_balance`, `class_mesh_builder` |
| `graph_forms` | 1 (bare) | **6** | the audit read the form's owned symbols as basis evaluation/derivative ACTIONS (OPERATOR, with the subset_set parentage secondary): shape machinery from which level-8 discovery assembles candidate R terms | `graph_fitting`, `class_polynomial_form`, `class_harmonic_form`, `class_form_pruner`, `class_diffusion_statement`, `class_fitted_balance` |
| `class_conduction`, `class_advection`, `class_robin_condition` | 3, stratification | **9** | constitution: material/model law | statements |
| `class_fitted_balance`, `class_diffusion_statement` | 3, stratification | **9 / 10** | the split between constitutive machinery (9) and the problem asked (10) needs the audit's per-symbol categories before the line is drawn | each other |
| `graph_grammar` | legacy compatibility contract | not a level | the four-root grammar is compatibility, retired by Phase 11, never renumbered | twenty-one legacy consumers (see audit) |
| `graph_calculus` | legacy compatibility | not a level | same | thirteen legacy consumers (see audit) |

## What the releveling commit must also do

1. Retire every "LEVEL n OF THE STRATIFICATION" banner in favour of one
   tower; no header may cite the old numbering except as history.
2. Restate `graph_minimization`'s header in the canonical direction
   vocabulary (`solve for data`); restate `graph_fitting`'s in
   `discover operator`. Neither is a rename of any public symbol.
3. Decide `graph_state`'s import surface before seating it at 7: the
   module may not import fitting, minimization or any inference
   machinery, ever — dependencies point upward into it, not out of it.
4. Leave `forward`/`reverse` out of every type name; directions are
   process vocabulary (COMPUTATIONAL-GRAPH.md 6).
5. Hold the certified tower suites (calculator, learning,
   derivative-action, adjoint, visualization, time-integration,
   partitioned) as unchanged oracles through the move.
