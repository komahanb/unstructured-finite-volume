# Roadmap status after the dissolution

The roadmap this file once carried classified the GTI prototype's
assumptions and ordered their retirements. The prototype is gone
(PR #40) and every retirement on that list either happened inside
the living modules or became unnecessary when its host was deleted.
This file records what happened to each item, which design rules
remain binding, and what is still open.

## Retirements, item by item

1. Hand-tabulated derivative patterns -> retired. The chain rule
   generates integer partitions of the degree recursively
   (`class_graph_chain_rule`); no degree is privileged, and the
   int64 bound on the multinomial count (degree > 20) is checked.

2. Affine-only design paths -> retired. Parameter paths arrive as
   `argument_path` values carrying full derivative sequences, both
   in `chain_rule % assemble` and in `march_directional`.

3. Uniform-only BDF rows -> retired. `set_bdf(k, steps)` prices the
   order-2 row on the steps actually taken and assigns the exact
   uniform constants when the steps are equal.

4. Missing adaptivity -> retired. `march_adaptive` plus the
   `step_policy` family (propose / judge / retry) choose steps at
   run time; the marcher measures, the policy decides.

5. Execution/solver seats -> retired. Linear solves live in the
   `graph_minimization` family (`class_graph_dense_direct`); no
   separate solver hierarchy exists.

6. Typed evaluation components, GEOM and time channels -> deleted
   with the prototype. State and parameters are input fields; a
   time-dependent statement takes time as one more parameter slot.

## Rules that remain binding

- No gti_* module may exist. Solvers of every kind join
  `graph_minimization`; a new linear-algebra hierarchy is not
  accepted.
- No `solve_transpose` method: a transposed system is solved by
  handing `transpose(a)` to the same solver.
- Every BDF coefficient is written in `set_bdf` and nowhere else.
- One Jacobian per edge serves every derivative order in
  `march_directional`; only the right-hand side changes with the
  order.
- Linearization dispatch happens only in `tangent_of`; callers do
  not inspect statement types themselves.
- A step policy reads scalars only (estimate, step, attempt count),
  never a state, a graph, or a statement.
- A coupled multi-stage (block implicit) scheme is one state of
  stacked stages solved by one `step_operator`; it is not a new
  marcher rule.
- Comments are operational: what is checked, what input is invalid,
  what happens on failure, why the check is needed.

## Open items

- BDF orders 3 and higher have no coefficient tables in `set_bdf`.
- The Adams corrector at nonuniform history (history-ratio
  coefficients) is not implemented; the current corrector reuses
  the backward step with a scaled hs.
- `march_adjoint` and `march_directional` assemble dense step
  Jacobians (`dense_matrix_of`); large states need a sparse or
  matrix-free route through the same minimizer family.
- The comment and naming rules above are enforced in the modules
  listed in `doc/GTI.md`; older modules elsewhere in `src/` have
  not been rewritten to them.
