# The graph time integrator

One program, `graph_time_integrator`, built from the single source
`module_graph_time_integrator.f90` over the library in `../src`.

## The problem it solves

Let q(t) be the state of the governing equation in descriptor form

    R( q^(N), q^(N-1), ..., q ; nu, t ) = 0        on  [0, T],

with N = `state_degree`, design parameter nu, and a time functional

    f(nu) = integral over [0, T] of  F( q, q', ..., q^(N-1) ; nu ) dt.

The program discretises [0, T], marches the equation with every scheme
the configuration names, and computes

    f,  df/dnu,  d^2 f/dnu^2,  ...,  d^m f/dnu^m,      m = max_derivative_degree,

printing one row per scheme. Every derivative is exact with respect to
the discrete problem: the residual and the functional are stated as
expressions over the state's components and the design, and their
partials of any order are obtained by evaluating those expressions
over an arithmetic that carries mixed derivatives - nothing is
differentiated by hand and nothing is differenced.

## Building

    cd ..            # the repository root
    ./build.sh       # the library, into lib/
    cd application
    ./build.sh       # the program

`PRECISION=quad ./build.sh`, at both levels, builds the tower over
real128 into `lib_quad/` and places the binary in `quad/`; the two
builds coexist.

## Running

    ./graph_time_integrator                          # config/homogeneous.cfg, the default
    ./graph_time_integrator --config=uniform         # config/uniform.cfg
    ./graph_time_integrator --config=uniform instants=41 "families=bdf dirk"

A configuration is named, not pathed. Every later `key=value` argument
overrides one setting; a value holding blanks is quoted whole. An
argument naming no setting stops the program: the input language has
no silently ignored word.

| configuration | what it states |
|---|---|
| `homogeneous` | every family alone, the default |
| `heterogeneous` | chains whose blocks change family along the horizon |
| `uniform` | the same schemes on the uniform partition |
| `linear` | the linear physics of a spatial field |
| `field` | a two-dimensional field: mesh, operator, march |
| `accounting` | the cost of each derivative order, measured |

## The objects of a run, and the keys that state them

**The equation.** `physics = vanderpol`, `state_degree = N`,
`design = nu`. The residual is one expression, stated once in the
module `gti_physics`:

    r = stated( derivative(q, N)
              - nu * (1 - derivative(q, 0)**2) * derivative(q, N-1)
              + derivative(q, 0),  N, 'van der pol residual' )

`derivative(q, d)` is the state's component of degree d - along the
instants the derivatives are unknowns the scheme relates, so the
symbol selects and does not differentiate. A new equation is a new
function beside this one, built from the components, the design,
constants, the four arithmetic operations, integer and real powers,
and sin, cos, exp, log, sqrt. The functionals F are expressions of the
same kind.

**The initial state.** `initial_state` lists q(0), q'(0), ...,
q^(N-1)(0) as blank-separated words, short lists padded with zeros.
The highest component q^(N)(0) is then solved from R itself, so the
initial state satisfies the equation rather than approximating it.

**The grid.** A grid is the discretisation of the coordinate [0, T]:
a finite measure, points t_k and weights h_k. `grid` is

| kind | the measure |
|---|---|
| `uniform` | h_k = T/(n-1), the instants equidistant |
| `random` | a reproducible drawn spacing from `seed`, each weight within [1/2, 3/2] of uniform |
| `adaptive` | the steps an error-controlled march discovers to `tolerance`, then frozen; `instants` is set by the result |

with n = `instants`. When `designs` names `grid`, the weights h_k join
nu as designs and the table reports df/dh beside df/dnu, together with
the identity sum over k of h_k df/dh_k = 0, the steps being
homogeneous of degree zero in their weights. The same map carries a
quadrature kind in the library: the Gauss rule's points and weights,
exact on polynomials of degree below 2n, demonstrated in
`assembled_tower`.

**The schemes.** `families = bdf adams dirk`, orders up to
`max_discretization_order`: backward differences and Adams-Moulton at
any order, diagonally implicit Runge-Kutta at orders two to four (the
implicit midpoint rule and the two Crouzeix tableaux).
`combinations = homogeneous | pairs | triples` builds, beside the
single-family rows, chains whose blocks change family along the
horizon, joined at their shared instants. A family whose constraint
reaches back over r > 1 instants is started by a stage block over
steps refined by `startup_refinement`, so every row integrates the
same initial-value problem.

**The functionals and their derivatives.**
`functionals = energy dissipation`, `designs = physics [grid]`,
`max_derivative_degree = m`. The derivatives are computed by a forward
expansion in the design; with several designs the tangent and adjoint
routes are chosen by counting - the forward route costs one solve per
design, the reverse one per functional - and either can be checked
against the other. `check` names a comparison against something known:

| check | the statement verified |
|---|---|
| `routes` | the tangent and adjoint routes agree over the whole derivative table |
| `sinks` | J_ii lambda_i = g_i on every unknown no row reads, at every degree |
| `ode` | a field at kappa = 0 equals one node's ordinary equation |
| `mode` | a rectangle at nu = 0 against the separated solution of the heat equation |
| `operator` | the fitted balance against kappa times the laplacian of the mode |

**The solvers.** Newton drives every block; `linear_solver` is
`direct` or `iterative` (GMRES), refined by `assembly`, `storage`,
`multigrid`, `krylov_restart`, `smoothing_sweeps`,
`max_linear_iterations`. `tolerance` with
`tolerance_criterion = relative | absolute` and
`iteration_criterion = by_rate | by_count` govern every stopping
question: each tolerance is relative or declared absolute, each budget
by rate or by count, and no threshold is a chosen number.

**A spatial field.** `spatial_counts = n1 n2` above zero lay a mesh
under the march: the unit square partitioned by the same grid
machinery as time, mapped by `spatial_geometry`
(`cartesian | circular | elliptical`), with the diffusion operator a
fitted polynomial balance of degree `spatial_order` and conductivity
`diffusion`. `sweep = space-time | time | space` solves each block
whole, instant by instant, or node by node to a fixed point.
`export = paraview` writes one `.vtu` per instant to `export_path`.

**Accounting.** `accounting = T` files what the run spends -
`measurements` among wall_time, primal_loops, tangent_loops,
adjoint_loops, newton_solves, linear_solves, factorisations - under
the level of the hierarchy it was spent in (expansion, horizon, block,
stage) and the derivative order, with a model-against-count table for
the sensitivity substitutions.

## The demonstrations

Each former standalone driver survives inside the program as a
demonstration. Each prints the quantity it checks and the measured
departure beside the floor it is held to; every floor is derived from
the arithmetic, none is chosen. Arguments after the demonstration's
name pass through to it.

    ./graph_time_integrator --list-demos

| demonstration | what it verifies | how to run |
|---|---|---|
| `adaptive_grid` | step-doubling grids at four tolerances; the step count scales as tol^(-1/(p+1)) and the two sensitivity routes agree on every grid | `./graph_time_integrator --demo=adaptive_grid` |
| `assembled_tower` | the expansion graph assembled and read back whole; the designed grid's partials; the Gauss rule exact on every power below 2n | `./graph_time_integrator --demo=assembled_tower` |
| `chained_horizon` | a chain across families: each derivative of the functional against a difference of the one below | `./graph_time_integrator --demo=chained_horizon` |
| `constraint_rows` | one block's rows; the physics partials against their closed form and a central difference | `./graph_time_integrator --demo=constraint_rows` |
| `coupling_relation` | the relations a scheme's coupling carries, and the weights on them | `./graph_time_integrator --demo=coupling_relation` |
| `expansion_check` | the expansion's derivatives to fourth order, each against a difference of the one below | `./graph_time_integrator --demo=expansion_check` |
| `family_coefficients` | multistep coefficients on non-uniform steps against the Lagrange functionals | `./graph_time_integrator --demo=family_coefficients` |
| `function_identities` | the exact arithmetic's elementary functions: fifty-five identities to five directions | `./graph_time_integrator --demo=function_identities` |
| `grid_design_check` | derivative tables in the step weights, checked three ways at every order | `./graph_time_integrator --demo=grid_design_check` |
| `jacobian_shape` | how far a block's rows reach, and how much of the square is empty | `./graph_time_integrator --demo=jacobian_shape` |
| `level_maps` | what each level of the tower carries, read back by one traversal | `./graph_time_integrator --demo=level_maps` |
| `level_shape` | two blocks built through the level storage and read back through the level view | `./graph_time_integrator --demo=level_shape` |
| `marched_block` | one block against the cosine of the harmonic oscillator; steps halved, the observed order against 2^p | `./graph_time_integrator --demo=marched_block` |
| `marched_horizon` | a horizon marched and its sensitivities, tangent against adjoint | `./graph_time_integrator --demo=marched_horizon` |
| `marched_stages` | a stage block against the cosine, order observed under halving | `./graph_time_integrator --demo=marched_stages` |
| `memory_shape` | what the representation costs, part by part | `./graph_time_integrator --demo=memory_shape` |
| `randomized_checks` | the splitting, route and order invariants over drawn parameters | `./graph_time_integrator --demo=randomized_checks 7 2` |
| `scheme_weights` | a row's weights against the polynomial the scheme reproduces | `./graph_time_integrator --demo=scheme_weights` |
| `sensitivity` | df/dnu by tangent, by adjoint, and by difference | `./graph_time_integrator --demo=sensitivity` |
| `solve_cost` | what one formation and one solve cost | `./graph_time_integrator --demo=solve_cost` |
| `tolerance_form` | where a march's tolerance floor lies: eps times the norms the solve carries | `./graph_time_integrator --demo=tolerance_form` |

`randomized_checks` takes a seed and a case count; run without them it
stops.
