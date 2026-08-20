# Time integration and differentiation

This document describes the time-integration and differentiation
modules as they exist in `src/` today. The GTI prototype this
document once described was dissolved into these modules (PR #40)
and its test directories were replaced (PR #41); a table at the end
maps the old module names to their current locations.

## The mathematical objects

Time is a chain graph: instants are vertices, steps are edges, and a
step size h is a number stored per edge. A physical statement is an
operation S on that graph; by convention S returns minus the
velocity, so the explicit update is q <- q - h S(q).

One step of a recurrence is the operation

    a0 q + a1 qold + a2 qolder + hs S(q)

whose coefficients (a0, a1, a2, hs) select the scheme. A derivative
of a composed quantity is assembled from partitions of the degree
(Faa di Bruno's formula). A linear operator is represented, when a
dense matrix is required, by evaluating it on the standard basis.

## Modules

### operation_action: max_degree and partial_action

Every operation reports, beyond `apply`, its exact partial actions:

    partial_action(input_graph, input_data, slots, directions, output)

is the statement differentiated once per entry of `slots(:)`,
contracted against the matching direction field, returned on the
statement's domain. `max_degree()` declares how deep this calculus
goes; the default is 0, and the default `partial_action` stops the
program. No derivative tensor is stored.

### operation_linearization

One derived operation, `linearization`: the tangent J v of a
statement at a frozen state, built by `tangent_of(action)`. It takes
the exact road, J v = D S(q)[v] as one partial action in input slot
1, when the statement's `max_degree` is at least one, and finite
differences of two residuals otherwise; its name begins `exact
derivative of` or `derivative of` accordingly. Every caller (newton,
the marcher) builds it through this one function.

### operation_chain_rule

`chain_rule % assemble(statement, input_graph, input_data, degree,
paths, output)` computes the total derivative d^n/ds^n of the
composition S(x_1(s), ..., x_m(s)).

An `argument_path` names one input slot and carries the derivative
sequence of that argument's path: `derivative(k)` is a
`path_derivative` holding x^(k) as a direction field, or unoccupied,
which reads as zero. Terms are indexed by integer partitions of the
degree, each carrying the multinomial count

    c(d) = n! / ( prod_i d_i! * prod_j multiplicity_j! ),

and every partition entry ranges over all paths providing that
derivative order, in ordered tuples. Partitions are generated, not
tabulated, so no degree is privileged.

Stops the program on: a negative degree, a path naming a slot the
statement does not take, two paths naming one slot, a partition
needing a partial past `max_degree`, and a multinomial count outside
int64 (degree 21 and beyond).

### operation_step

`scheme` is the step recurrence as an operation. Constructors
`backward_euler(action, h)`, `bdf(k, action, h)`, and
`bdf_variable(k, action, steps)` fill the coefficients; every
coefficient comes from `set_bdf(k, steps)`, which reconfigures a
standing operator between edges. With h0 = h_n and h1 = h_{n-1} the
order-2 row is

    a0 = (2 h0 + h1)/(h0 + h1)
    a1 = -(h0 + h1)/h1
    a2 = h0^2 / (h1 (h0 + h1)),

assigned as the exact constants (3/2, -2, 1/2) when h0 = h1. Only
orders 1 and 2 have tables. A diagonal RK stage and an Adams
corrector are this same backward step with a scaled hs against an
externally assembled base, so they need no constructors.
`dependencies()` returns the fan-in of reach + 1 instants onto the
newest, including the self-edge for the implicit unknown.

### operation_marching

The integrator. Rules: `MARCH_FORWARD` (explicit), `MARCH_BACKWARD`
(backward Euler, one solve per edge by the held minimizer `inner`),
`MARCH_BDF2` (order 2, started by one backward step).

- `march(action, on, q, nsteps [, steps] [, trajectory])` integrates
  forward; `steps(:)` gives one step per edge, `trajectory` records
  the state at every instant.
- `march_adjoint(transposed, on, lambda, nsteps [, steps] [, action]
  [, trajectory] [, seeds])` runs the chain in reverse. Explicit
  rule: the caller's transposed statement is applied edge by edge.
  Implicit rules: backward substitution; at each edge

      J_e^T lambda_e = seed_e,      J_e = tangent_of(scheme),

  the tangent of the step equation taken at the recorded state,
  compiled to a `stencil`, transposed, and solved by `dense_direct`;
  the couplings to earlier instants are the scheme's constant
  coefficients a1 and a2. Requires `action` and
  `trajectory`. `seeds(:, k)` is added when instant k is reached; on
  return `lambda` is the sensitivity at the first instant.
- `march_directional(action, on, nsteps, trajectory, order,
  sensitivities [, steps] [, parameters] [, paths])` computes
  forward directional derivatives of every order up to `order`. At
  each edge

      J_e q_e^(s) = -( a1 qprev^(s) + a2 qolder^(s)
                       + degree-s composition over the scheme with
                         the order-s state derivative set to zero ),
      J_e = tangent_of(scheme) at the solved state,

  so one tangent per edge, attached to `dense_direct`, serves every
  order. Parameter paths arrive on input slots 2 and higher as
  `argument_path` values; the state's path is computed, never
  supplied.
- `march_adaptive(action, on, q, duration, policy, max_attempts,
  steps_taken, completed)` chooses steps at run time: the policy
  proposes, the marcher computes a trial without modifying q,
  measures the distance from the extrapolating predictor, and the
  policy judges. A rejected trial changes nothing; a spent attempt
  budget returns `completed = .false.` with q at the last accepted
  state. What returns is `steps_taken`; trajectory and adjoint come
  from marching again with those steps.

### operation_step_policy

The decision half of adaptive marching: `propose` (first step at a
new edge), `judge` (accept or reject from the estimate alone),
`retry` (step after a rejection). A policy reads scalars only.
`halving_policy` is the reference member: start at `first_step`,
accept at or below `tolerance`, halve on rejection.

### operation_dense_direct

Gaussian elimination with partial pivoting as a concrete minimizer
in the `operation_minimization` family. The matrix is assembled by
applying the attached operation's matvec to each basis vector. A
pivot at or below `singular_tolerance` stops the program. It owns no
matrix representation: `stencil(action, on, width)` compiles any
operation onto a stencil by evaluation on the standard basis, and a
stencil's `transpose()` reverses its edges, so a dense array is
never built by hand.

There is no `solve_transpose`: the transposed stencil is attached
and solved like any other operation.

### operation_newton

The nonlinear governor. Its Jacobian action is `tangent_of(action)`,
so a differentiable statement is linearized exactly and any other
statement by differences, with no dispatch of newton's own.

### map_change_protocol, map_value, map_value_change

The reversible mutation stack. `run_change` owns the
lifecycle apply -> check -> keep | revert; a failure reported by
apply or check is reverted and returned; a step that returns without
marking its work stops the program. `value_map` stores, per graph
identity, a status (UNATTACHED / UNKNOWN / KNOWN) and a value field
on the graph's own domain, keyed on identity tokens copied at
attach, never on position. `value_change` is the concrete change
that updates a value map through the controller, with revert
restoring the prior row exactly.

## Where the old GTI modules went

| Old gti_* content                          | Current location                          |
|--------------------------------------------|-------------------------------------------|
| differentiable form contract               | `operation_action` (max_degree, partial_action) |
| form evaluator, evaluation point, buffers  | deleted; fields and domain checks cover it |
| chain-rule assembler                       | `operation_chain_rule`                   |
| motifs, motif builders, variable-step rows | `operation_step` (set_bdf)               |
| time graph, forward driver                 | `operation_marching` (march)              |
| local newton driver                        | `operation_newton` + `tangent_of`        |
| tangent driver                             | `operation_linearization`          |
| adjoint, reverse, seed drivers             | `march_adjoint` (backward substitution)    |
| degree-2 and degree-r directional drivers  | `march_directional`                        |
| adaptive growth driver, controller         | `march_adaptive` + `operation_step_policy`|
| linear solve helpers                       | `operation_dense_direct`                 |
| change protocol                            | `map_change_protocol`                    |
| attached value map, value change           | `map_value`, `map_value_change`    |

## Tests

- `test/graph-differentiation`: set_bdf rows, tangent_of dispatch,
  chain rule to degree 8 against analytic and Taylor-convolution
  values, march_directional to degree 8 in exact rationals,
  march_adjoint closed forms, march_adaptive closed forms, the
  nonuniform implicit march.
- `test/graph-change`: the controller lifecycle, the value map
  transitions and storage rules, value_change restorations.
- `test/graph-dense-direct`: the elimination, a stencil compiled
  from an operation, its transpose, and static checks on the source.
- `test/graph-marching`, `test/graph-minimization`: the marcher's
  integration rules and the minimizer family.
