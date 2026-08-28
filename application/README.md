# The graph time integrator

The main program, `graph_time_integrator`, is built from a single
source, `module_graph_time_integrator.f90`, over the library in
`../src`.

# Mathematics and Architecture

**Pass 1 — steady: no time, a point domain.**

For the design state nu, the physical field Q must satisfy the
governing law R without residue, while producing an output f read on
the design. For example, given the airfoil shape, the flow state
variables must satisfy the conservation laws without residue to
produce the lift that sustains the required maneuver. The
correspondence graph between f and nu is drawn below:

```mermaid
graph LR
  NU["nu — free: input"] --> R["R — governing law"]
  R -->|determines| Q["Q — state"]
  NU --> F["F — output relation"]
  Q --> F
  F --> f["f — sink"]
```

A control knob is turned to pitch the airfoil; higher-order
information about the change of lift with the pitch input is of
enormous value. A high-fidelity computation of the lift and its m-th
rate with respect to the design is non-trivial because of the
complexity of the framework that simulates the output: the larger the
codebase, the harder it becomes to extend the output's derivative
order - a code entropy barrier that limits the pace of technology
advancement.

The remedy is a separation of concerns into basis and coefficients:
the basis is the topology of the graph - which vertex reads which,
fixed once by the statement of the problem - and the coefficients are
the numbers its edges carry, the partials the arithmetic evaluates.
The rate df/dnu is then not code but a coefficient of the graph: it
captures how the output moves along the one free vertex, and the
m-th rate is reached by adding levels to the derivative tower below,
never by writing new code paths. That is what removes the barrier by
construction: the derivatives are properties of the graph, not code
to be written.

In this pass every vertex lives on a point domain: that is what
steady means in this reading, and it is why discretisation has
nothing to act on here and the square's vertical arrows are
identities. The model is
a relation between knob and state:

    R(Q; nu) = 0,        f = F(Q; nu).

Which vertex is the input is not a convention but an accounting: the
input is whatever the constraints leave free; everything else is
determined, and the determined vertices split into the state - read
again downstream - and the outputs, read by nothing. R is square
against Q and consumes exactly its freedom; nothing closes over nu;
so nu is the source, Q interior, f a sink. Where no model connects two
vertices, the identity relation Q = I(nu) is the default operator
between them. An edge is a read and carries the relation, never a
rate: the partials are not edges but the weights the linearised
traversal assigns to these same edges.

**The derivative tower.** The rates of change are not a chain hanging
off f; they are a second layer of vertices over the same topology,
with known connectivity. Writing V^(s) = d^s Q / dnu^s and
differentiating the one functional R(Q(nu); nu) = 0 repeatedly along
the design, every level obeys the same equation:

    J V^(s) = - R^(s),        s = 1, ..., m,

with J = R_Q the one operator, factored once, and R^(s) the s-th
design-derivative of the residual taken with the highest tangent
V^(s) withheld - not a new functional at any level, but R
differentiated. The first few, written out, show every entry to be a
partial of R contracted against the lower tower:

    R^(1) = R_nu
    R^(2) = R_QQ[V1, V1] + 2 R_Qnu[V1] + R_nunu
    R^(3) = R_QQQ[V1, V1, V1] + 3 R_QQ[V1, V2] + ...           (Faa di Bruno)

The structure is strictly lower triangular and dense below the
diagonal, not diagonal: V^(s) reads all of V^(1..s-1), never its
predecessor alone. The outputs are parallel sinks over the shared
tower - f^(m) reads V^(1..m) and nu, and no output reads another
output. The edges now carry the weights they always deserved, on the
vertices they belong to:

```mermaid
graph TD
  NU["nu"] -->|"R_nu"| V1
  subgraph T["the tangent tower — one J, factored once"]
    V1["V1 :  J V1 = − R^(1)"]
    V2["V2 :  J V2 = − R^(2)"]
    V3["V3 :  J V3 = − R^(3)"]
    V1 -->|"R_QQ[V1,·], R_Qnu"| V2
    V1 -->|"R_QQQ, R_QQnu, ..."| V3
    V2 -->|"R_QQ[·,V2]"| V3
  end
  V1 -->|"F_Q"| f1["df/dnu"]
  V1 -->|"F_QQ[·,V1], ..."| f2["d2f/dnu2"]
  V2 -->|"F_Q"| f2
  V1 --> f3["d3f/dnu3"]
  V2 --> f3
  V3 -->|"F_Q"| f3
  NU -->|"F_nu"| f1
  NU --> f2
  NU --> f3
```

**Three systems, one after another.** U is the primal state (Q of
section 2, written U here to sit beside its own derivatives), V the
tangent, W the adjoint. Only the first is nonlinear.

    nonlinear (state U, Newton):
    
    J^k dU  =  - R(U^k; nu),        U^(k+1) = U^k + dU,        J^k = R_U(U^k; nu)

One linearisation per iterate, solved and discarded, until R(U; nu)
falls to tolerance; J at the converged U is the one operator every
later system reuses, unrefactored.

    tangent (V, forward substitution):
    
    [ J              ] [ V^(1) ]      [ R^(1) ]
    [ L21  J         ] [ V^(2) ]  = - [ R^(2) ]
    [ L31  L32  J    ] [ V^(3) ]      [ R^(3) ]
    
    adjoint (W, back substitution - the transpose):
    
    [ J^T  L21^T  L31^T ] [ W^(1) ]     [ g^(1) ]
    [      J^T    L32^T ] [ W^(2) ]  =  [ g^(2) ]
    [             J^T   ] [ W^(3) ]     [ g^(3) ]

with the L's the partial-of-R couplings drawn above, and g^(s) the
seed the functional releases at order s - the F-partials read into
V^(s) in the tower diagram above. Back substitution reads the upper
system bottom row first: J^T W^(3) = g^(3), then
J^T W^(2) = g^(2) - L32^T W^(3), then
J^T W^(1) = g^(1) - L21^T W^(2) - L31^T W^(3) - the same J,
transposed, at every level, exactly as the forward system reused it
untransposed, and exactly the J the nonlinear solve above already
factored once.

The last two are the generalized sensitivity statement: one matrix,
its two triangles the two traversal directions. Each output closes
either way - as the contraction g . V along the tower (tangent), or
as W . (-R) along the seeded column (adjoint) - and the two are the
two evaluations of one bilinear form, W^T (block matrix) V. In the
code this is literal: one operator, one boolean naming the
orientation, and their agreement to round-off is what
`check = routes` measures.

Seen from linear algebra alone, this is an LU solve whose
factorisation was performed by causality: the time ordering of the
graph is an elimination ordering, so the space-time operator arrives
already triangular - the march is the forward substitution, the
adjoint the back substitution of the transpose, as in a Cholesky
solve where one factor serves both directions. And what the two
substitutions compute is an entry of the inverse: every sensitivity
is the bilinear form g . T^-1 r, evaluated either as a column of the
inverse (tangent, one solve per input) or as a row (adjoint, one
solve per output - the discrete Green's function of the functional).
The route gate chooses rows against columns of an inverse never
formed, and the routes check is, in classical language, a reciprocity
test: g . (T^-1 r) = (T^-T g) . r.

The inverse itself is never the thing to compute - it costs a
factorisation but applies with worse rounding and fills the sparsity
in - and what stands in for it is the factorisation of J, held under
a stamp and reused exactly across every level of the tower, the
adjoint included, since they share J to the last bit. Where J drifts,
the old factors become a preconditioner rather than an answer, and
each drift names an economy: frozen factors across Newton iterates,
factors carried from one step to the next until the inner iteration
count degrades, factors carried from one design to the neighbouring
one. And preconditioning the whole space-time matrix T is meaningful
in exactly one circumstance - when the causal order is deliberately
broken to solve across time in parallel, an approximate triangular
solve preconditioning the exact one; that is the mathematical content
of parallel-in-time iteration, and the sequential march is the limit
in which the preconditioner is exact and one sweep suffices.

The functional has the same triangle - F_Q repeated on its diagonal,
the higher F-partials below - but where the R-tower is solved, the
F-triangle is only ever applied: nothing in an output is implicit.
The framework therefore never assembles it; it forms its action by
evaluation - the tower is loaded into the arithmetic's subsets, the
functional's expression is evaluated over them, and the coefficient
that comes back is one block-row already contracted. In the joint
matrix of states and outputs the observation rows carry the identity
on their diagonal, which is the precise sense in which outputs are
sinks - their costate equation is trivial, the functional's own
multiplier is one - and why assembling the F-triangle would be
assembling rows whose solves are free.



```mermaid
flowchart TD
  F0["f^(0) = F(Q,nu)"]

  FQ["F_Q"]
  Fn["F_nu"]

  FQQ["F_QQ"]
  FQn["F_Qnu"]
  Fnn["F_nunu"]

  Fh["higher F_{Q^a nu^b}"]

  F0 --> FQ
  F0 --> Fn

  FQ --> FQQ
  FQ --> FQn
  Fn --> FQn
  Fn --> Fnn

  FQQ --> Fh
  FQn --> Fh
  Fnn --> Fh
```

**The functional as one more state.** The same fact, written
implicitly so that one structure serves both inputs: make the output
an unknown y with a residual in R's own form,

    R_aug(Q, y; nu) = [ R(Q; nu), y - F(Q; nu) ] = 0,
    
    J_aug = [ J     0 ]
            [ -F_Q  I ]

The diagonal of the y-rows is the identity - a zero diagonal would be
singular; the structural zero is the upper block, R never reading y,
which is "outputs are sinks" said in algebra and what keeps J_aug
triangular, so the augmentation is free: the tower recursion on J_aug
delivers V^(s) and f^(s) together, and the adjoint seeded with e_y
gives W_y = 1 by the identity block - the functional's
multiplier held at one is the trivial back-substitution. In time the
y-row is a running sum, y_k - y_(k-1) - h_k F_k = 0: a scheme row of
the most degenerate family, one-step reach, the quadrature weights
its coefficients. Both R and F already enter as expressions; what the
uniform treatment costs is exactly the several-unknowns extension,
so that (q, y) share a slice - after which the functional's separate
machinery collapses into the ordinary rows, the terminal seed on y
replaces g, and the sinks check certifies the functional rows with no
new code. The extension deletes a code path rather than adding one.

The continuous statement and its discrete image commute: discretising
the state and discretising the model are the two sides of one square,

    Q  ———————  R          (continuous)
    |           |
    Q̄  ———————  R̄          (discretised)

and the program lives on the bottom row. Every derivative it reports
is exact with respect to the bottom row — the discrete problem — not
an approximation of the top one.

**Pass 2 — unsteady: time, a line domain.**

Into the same graph we plug the time-dependent model, and apply the
one rule of a pass to *every* vertex: open it into its members.

The knob opens, and its members carry their own domains. What was the
point {nu} is now the set

    {nu}  <——  { mu, t },        mu on a point domain, t on the line [0, T],

mu the van der Pol parameter and t the independent variable, with R
stated over both:

    R(Q; mu, t) = 0        on  [0, T].

The domain of a member decides everything that follows. A point domain
carries one value and needs no grid, so mu passes to the discrete row
unchanged. A one-dimensional domain must be discretised, so t is the
member the grid acts on: it opens into the instants {t_1, ..., t_n},
equivalently the steps {h_k}, which is precisely how the domain joins
mu as a design (`designs = grid`) and why T is not the only parameter.
For van der Pol the edge from t into R carries zero — the equation is
autonomous — but the vertex is there systematically. The state q is a
field over t's line, so it opens twice: degree-wise here, instant-wise
under the grid. And f closes the ledger: the integral collapses the
line back to a point, so the output lives where the knob mu does —
which is what makes df/dmu a number.

The state opens. What was the point {Q} is now the set, ordered lowest
degree to highest,

    {Q}  <——  { q, q', ..., q^(N) },      N = state_degree,

with the edges carrying the scheme's coefficients that relate the
members: along the instants the derivatives are unknowns, and the
integration stencil of the configured family and order connects them.

The output opens. What was the value f is now the time functional and
its derivatives in the opened knob,

    f(mu) = integral over [0, T] of  F( q, q', ..., q^(N-1) ; mu ) dt,

and the residual and the functional are two trees over the same
leaves:

```mermaid
graph TD
  MU["mu — parameter"] --> R["R = q^(N) − mu (1 − q²) q^(N−1) + q"]
  T["t — time"] -.->|"zero for van der Pol: autonomous"| R
  Q0["q"]  --> R
  QN1["q^(N−1)"] --> R
  QN["q^(N)"] --> R
  Q0 --> F["F = (q² + q'²)/2"]
  Q1["q'"] --> F
  T -->|"discretised: instants t_k, steps h_k"| QT
  R -.->|"marched to zero at every instant"| QT["q(t_k), all instants"]
  QT --> INT["sum of h_k F_k"] --> f["f, df/dmu, ..., d^m f/dmu^m"]
  MU --> f
```

The program discretises [0, T], marches the equation with every scheme
the configuration names, and computes

    f,  df/dmu,  d^2 f/dmu^2,  ...,  d^m f/dmu^m,      m = max_derivative_degree,

and, when the steps are designs, df/dh_k beside them. Every derivative
is exact with respect to the discrete problem: the residual and the
functional are stated as expressions over the state's components and
the knob's — mu is a leaf of both trees, t a leaf of R — and their
partials of any order are obtained by evaluating those expressions
over an arithmetic that carries mixed derivatives. Nothing is
differentiated by hand and nothing is differenced. In the
configuration's flat vocabulary the one physics parameter mu is the
key `design`, written nu in the code.

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

Each key states one member of the problem graph, and the member's
domain decides its treatment: a point passes to the discrete row
unchanged, an extended domain is discretised by a grid, and an
operator is discretised by the scheme or stencil that binds it to that
grid's arithmetic.

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
initial state satisfies the equation rather than approximating it. In
the graph this is the degenerate junction: the boundary with the empty
past, treated by the same mechanism as the junction between two
blocks, not by a separate one.

**The grid.** The grid is the data-side vertical arrow of the square,
acting on the one extended domain the ODE problem has: t's line. It is
the discretisation of the coordinate [0, T] — a finite measure, points
t_k and weights h_k — and nothing else; the operator is discretised
separately, by the scheme. `grid` is

| kind | the measure |
|---|---|
| `uniform` | h_k = T/(n-1), the instants equidistant |
| `random` | a reproducible drawn spacing from `seed`, each weight within [1/2, 3/2] of uniform |
| `adaptive` | the steps an error-controlled march discovers to `tolerance`, then frozen; `instants` is set by the result |

with n = `instants`. When `designs` names `grid`, the weights h_k join
nu as designs and the table reports df/dh beside df/dnu, together with
the identity sum over k of h_k df/dh_k = 0, the steps being
homogeneous of degree zero in their weights. The same map carries a
quadrature kind in the library — the Gauss rule's points and weights,
exact on polynomials of degree below 2n, demonstrated in
`assembled_tower` — which is the measure a knob would take were its
point domain widened to an axis with a distribution: the third pass
this document does not yet write.

**The schemes.** The scheme is the operator-side vertical arrow of
the square: R bound to the grid's arithmetic. Its exactness class is
what the order p names — the square commutes identically on the
polynomials the scheme reproduces (`constraint_rows` shows the derived
rows vanish there) and to order h^p off them (`marched_block` shows
the defect fall as 2^p under halving). `families = bdf adams dirk`,
orders up to
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
routes are chosen by counting sources against sinks - the tangent
sweeps forward from the free vertices and costs one solve per design,
the adjoint sweeps backward from the outputs and costs one per
functional - and either can be checked against the other, the two
being the same linearised operator traversed with and against its
edges. `check` names a comparison against something known:

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

**A spatial field.** `spatial_counts = n1 n2` above zero widen the
state's domain a second time: beside t's line, two parametric lines
xi and eta, discretised by the very same grid machinery — extended
domains get grids, points pass through — and mapped into the plane by
`spatial_geometry` (`cartesian | circular | elliptical`), whose
mapping pushes the parametric measure forward into cell volumes. The
operator side follows its own arrow, as always: the diffusion operator
is a fitted polynomial balance of degree `spatial_order` with
conductivity `diffusion`, bound to that mesh.
`sweep = space-time | time | space` solves each block whole, instant
by instant, or node by node to a fixed point.
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
| `memory_shape` | the peak memory of one part of the representation: none, vertices, edges, field or block | `./graph_time_integrator --demo=memory_shape block 41` |
| `randomized_checks` | the splitting, route and order invariants over drawn parameters | `./graph_time_integrator --demo=randomized_checks 7 2` |
| `scheme_weights` | a row's weights against the polynomial the scheme reproduces | `./graph_time_integrator --demo=scheme_weights` |
| `sensitivity` | df/dnu by tangent, by adjoint, and by difference | `./graph_time_integrator --demo=sensitivity` |
| `solve_cost` | what one formation and one solve cost | `./graph_time_integrator --demo=solve_cost` |
| `tolerance_form` | where a march's tolerance floor lies: eps times the norms the solve carries | `./graph_time_integrator --demo=tolerance_form` |

`memory_shape` takes a part and an instant count, `block` and 41 unless
given; `randomized_checks` takes a seed and a case count, 7 and 2 unless
given.
