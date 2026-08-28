# Physics and functionals as graphs of typed operators

A governing equation and a functional are each stated once, at one
instant, as an expression in the state, the design and the elementary
operations. The expression is a graph: its vertices are typed
operators, its edges are the reads between them, its leaves are the
components of the state and the design, its root is the residual (or
the integrand). Nothing about the physics is declared beside the
expression; the degree, the partials in either argument, and the
spatial structure are derived from the graph.

    residual R   =   q'' - nu (1 - q^2) q' + q

              (-)                               root
             /   \
           (-)   q(0)                           leaf: state, degree 0
          /   \
       q(2)   (*)
             /   \
          (*)    q(1)
         /   \
       nu    (-)
            /   \
           1    (^2)
                  |
                q(0)

    energy F     =   (q^2 + q'^2) / 2          a second root over the same leaves

## What the framework already supplies

Two facts settle most of the design.

1. **The partials are already exact and automatic.** `nodal_integrand`
   (`physics/physics_integrand.f90`) evaluates a rule over
   `derivative_terms` (`src/util_derivative_terms.f90`), which carries
   the value and every mixed partial along the directions seeded on the
   state and the design. `gti_block` reads the tangent through
   `partial_action`, one direction per degree. So a Jacobian, a design
   partial, or a mixed partial of any degree costs one evaluation of the
   expression over `derivative_terms`. **No symbolic derivative is
   built**: differentiating an expression tree by rewriting (product
   rule, chain rule) grows the tree with every product and needs a
   simplifier; evaluating it over `derivative_terms` does not.

2. **A time derivative is a read, not an operator.** The state holds one
   component per degree per instant; the scheme (BDF, Adams, DIRK) is
   what relates the degrees. So the vertex `derivative(q, k)` selects the
   component `q(k)` at the instant. A spatial derivative is the stencil
   `gti_block` adds as its `spatial` rows, linear in `q(0)`.

What the framework lacks, and the draft was right to want:

- elementary functions (`sin`, `cos`, `exp`, `log`, `sqrt`, real powers)
  over `derivative_terms`: the rule today is restricted to `+ - * /` and
  integer powers;
- a statement of the physics as data rather than as a Fortran type per
  equation, so that a new equation is a line and the degree is derived;
- one type for a residual and for a functional integrand.

## The expression

One concrete type, parented to `operation` through the type that is
already the rule at one instant: `nodal_integrand`
(`physics/physics_integrand.f90`) declares the two arguments, the
domain (one value per instant), `max_degree` and `partial_action` over
`derivative_terms`, and leaves only `at_instant` deferred. The
expression is that type with the rule held as data in place of the
deferred binding, so the abstract and every physics type extending it
are replaced by one concrete `expression`, and the module moves to
`src/operation_expression.f90` following doc/coding-standards.md
(`operation_<subject>`, type in English order like `stencil`, `scheme`,
`fit`, `walk`).

A vertex has a kind, reads at most two earlier vertices, and holds one
number (a constant, or an integer exponent) and one integer (a leaf
component, or a function index). The vertices are stored in evaluation
order and the root is the last one. This is the same choice
`operation_walk` makes: a new operator kind costs a case, not a class.

    module operation_expression

      type, extends(operation) :: expression
         integer , allocatable, private :: kind(:)
         integer , allocatable, private :: first(:), second(:)   ! the vertices read; 0 if none
         integer , allocatable, private :: argument(:)           ! for a leaf: which argument is read
         integer , allocatable, private :: order(:)              ! leaf component (state degree), or function index
         real(dp), allocatable, private :: coefficient(:)        ! constant, or integer exponent
         integer               , private :: degree = 0           ! of the state it is bound to
      contains
         procedure :: name, domain, apply, max_degree, partial_action   ! the operation contract, as nodal_integrand has it
         procedure :: evaluated          ! the rule at one point, over derivative_terms
         procedure :: highest_degree     ! the largest state component read
         procedure :: reads_design       ! whether a design leaf is present
         procedure :: written            ! the expression as text, for reports
      end type expression

    Vertex kinds:   LEAF(argument j, component i)  CONSTANT(c)     the state's component is its degree
                    SUM  DIFFERENCE  PRODUCT  QUOTIENT  POWER(n)
                    FUNCTION(f)  with f in {sin, cos, exp, log, sqrt}
                    EXPANSION(fit)  PIECE(grid)                     see "Forms" below

Construction is by the intrinsic operators and the constructors below,
so the physics reads as it is written. A binary operator appends the
right operand's vertices after the left's, renumbers them, and adds one
vertex reading both roots.

    q  = unknown()                        ! the state
    nu = design()                         ! the design
    derivative(q, k)                      ! the component of degree k, along the instants
    1.0_dp - derivative(q, 0)**2          ! real * expression, expression ** integer
    sin(x), cos(x), exp(x), log(x), sqrt(x)

    operator(+), operator(-), operator(*), operator(/) :   expression x expression,
                                                           real x expression, expression x real
    operator(**)                                       :   expression x integer

Two consequences of the storage:

- The table *is* the reads relation R ⊆ V x V, so it can be handed to
  the framework's directed view when a walk over it is wanted
  (topological order, display). No walk is needed to evaluate it: the
  vertices are already in an order every read precedes.
- A subexpression used twice is stored twice. That is a tree, not a
  DAG. Sharing would need identity (tokens) on vertices; at the size of
  a physics rule it saves no evaluation, and it is noted here so that it is a
  decision and not an oversight.

### Evaluation

One loop over the vertices, over `derivative_terms`. Directions seeded
on `q` and `nu` propagate through every operation, so the root carries
the value and every requested mixed partial at once.

    pure function evaluated(this, arguments) result(r)
      type(derivative_terms), intent(in) :: arguments(:)   ! components of every argument, in order
      type(derivative_terms) :: v(size(this % kind))
      do i = 1, size(this % kind)
         select case (this % kind(i))
         case (LEAF);       v(i) = read(arguments, this % argument(i), this % order(i))   ! q(order) or nu
         case (CONSTANT);   v(i) = derivative_terms(this % coefficient(i), arguments(1))
         case (SUM);        v(i) = v(this % first(i)) + v(this % second(i))
         case (DIFFERENCE); v(i) = v(this % first(i)) - v(this % second(i))
         case (PRODUCT);    v(i) = v(this % first(i)) * v(this % second(i))
         case (QUOTIENT);   v(i) = v(this % first(i)) / v(this % second(i))
         case (POWER);      v(i) = integer_power(v(this % first(i)), nint(this % coefficient(i)))
         case (FUNCTION);   v(i) = composed(v(this % first(i)), this % order(i))
         end select
      end do
      r = v(size(this % kind))
    end function evaluated

Refused: a state degree read past the degree the statement is bound to;
`log`, `sqrt` or a real power at a nonpositive value (no derivative
there); `abs` is not offered, for the same reason.

## Elementary functions over derivative_terms

The one mathematical addition. For `g = f(a)` with `a` carrying
coefficients on subsets of the n directions, the coefficient of `g` on
a subset `m` is a sum over the set partitions of `m`:

    g_m  =  sum over partitions {B_1..B_p} of m  of  f^(p)(a_0) * a_{B_1} * ... * a_{B_p}

computed without enumerating partitions. Let `i` be the lowest element
of `m` and define `G(k, m)` as the same sum with `f^(p+k)` in place of
`f^(p)`. The block containing `i` is `S ∪ {i}` for some `S ⊆ m \ {i}`,
and the remaining blocks partition what is left, with one more
derivative taken:

    G(k, ∅)  =  f^(k)(a_0)
    G(k, m)  =  sum over S ⊆ m\{i}  of  a_{S ∪ {i}} * G(k+1, m \ {i} \ S)
    g_m      =  G(0, m)

`G` is filled in increasing `m`, for `k = 0 .. n - |m|`: a table of
`(n+1) x 2^n` numbers, each entry a subset sum, so the cost is the same
`3^n` the product already pays. The function supplies its derivatives
at the value `a_0` for `k = 0 .. n`:

    exp     e^x for every k
    sin     sin, cos, -sin, -cos, repeating
    cos     cos, -sin, -cos, sin, repeating
    log     log x, then (-1)^(k-1) (k-1)! / x^k
    x^p     p (p-1) ... (p-k+1) x^(p-k)         sqrt is p = 1/2

So `util_derivative_terms` gains one kernel, `composed(a, f)`, five
short derivative tables, and generic interfaces on the intrinsic names.
`integer_power` stays: it is exact at a zero base, where the real power
is not.

## Where it goes in the tower

    util_derivative_terms      +  composed, sin cos exp log sqrt, real power        (~80 lines)
    operation_expression  NEW     physics_integrand moved to src and made concrete:
                                  + vertices, operators, derivative(), evaluated    (349 → ~500 lines)
                                  - nodal_integrand (abstract), at_instant, zero_integrand
    physics_vanderpol          →  two functions returning an expression            (149 → ~25 lines)
    gti_block, gti_chain, gti_expansion, gti_field, gti_march, gti_stage, gti_taylor,
    chained_horizon, constraint_rows      one rename, class(nodal_integrand) → type(expression), 40 sites

Where `derivative_terms` sits. Differentiation is an operation, and in
the tower it already is one: `partial_action`, `linearization`,
`tangent_of`. `derivative_terms` is not the derivative but the number
the rule is evaluated in so that the derivative comes out - the role
`real(dp)` plays for a stencil. It stays beside `util_precision` as the
arithmetic, and the expression's `partial_action` is the operation that
uses it.

The same theorem is stated three times on three substrates:
`derivative_terms` does the product rule on subsets of directions,
`operation_chain_rule` does the total derivative of a composition with
integer partitions and multinomial counts, and `composed` does a
function of a `derivative_terms` value with set partitions. All three
are Faa di Bruno's formula; the set-partition form is the general one
and the integer-partition form is its symmetric case (`set_symmetric`).
`operation_chain_rule` could therefore be derived from
`derivative_terms` with symmetric seeding. That is a possible later
reduction, to be measured, not part of this plan.

The operation contract of `nodal_integrand` (two arguments, one value
per instant, `max_degree`, `partial_action`) is kept exactly; the only
change above the physics is the type name in signatures. An expression
is bound to the degree of the state it reads when it is stated:

    function stated(rule, degree, label) result(this)
       ! refused: rule % highest_degree() > degree - the state has no such component.
       ! a residual must read degree exactly, since that is the one degree no
       ! derived row determines; an integrand may read less.

The van der Pol file becomes the two rules and nothing else; the
comment block listing the partials by hand is deleted, since they are
no longer stated anywhere:

    function van_der_pol(n) result(r)
      q = unknown(); nu = design()
      r = stated( derivative(q, n) - nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, n-1) + derivative(q, 0), &
                & degree = n, label = 'van der pol residual')
    end function

    function van_der_pol_energy(n) result(f)
      q = unknown()
      f = stated( 0.5_dp * (derivative(q, 0)**2 + derivative(q, 1)**2), degree = n, label = 'van der pol energy')
    end function

`chained_horizon` and `constraint_rows` call these as they do today.
Once `gti_configuration` can pick a rule by name, the file can go
entirely and the rules live beside the other configured choices.

### The hierarchy this completes

The tower already evaluates each level as a loop over a reads graph:
the march is the loop over the blocks' reads graph forward and the
adjoint the same loop in reverse; a block is the loop over its
instants. The expression is the same structure one level down, inside
one instant:

    expansion  ⊃  sweep  ⊃  block  ⊃  slice (instant)  ⊃  expression vertex

and the derivative propagation is the same loop with `derivative_terms`
in place of numbers. Every constant vertex can be promoted to a design
leaf, which is how the network of `doc/NEURALNETWORK.md` extends to
the physics: a coefficient of the equation becomes a parameter, and its
gradient comes from the existing sensitivity sweeps with no new code.

## Coordinates: the continuum axis and its discretisation

Reviewed against `pspace/core.py` (branch `bsf`), which states the
probabilistic axis the way this framework states the temporal and
spatial ones. The two are the same construction with different
measures, and each has the half the other lacks.

A coordinate is the uniform axis xi on [0, 1], identified by token and
labelled. Every physical axis is its image under a mapping, and the
measure on the physical axis is the pushforward of the uniform measure
under that mapping - the Lebesgue measure on [0, T] is uniform up to
the scale T, which the mapping's derivative carries:

    time           t = T xi,  or a designed grid t(xi)        dt = t'(xi) d xi
    space          x = Phi(xi, eta), the geometry              |det D Phi| d xi d eta
    probability    y = F^-1(xi), the inverse distribution      rho(y) dy

So an integral over any axis is an expectation over xi of the
integrand times the mapping's Jacobian, and nothing else
distinguishes the three. pspace's `CoordinateType` enum
(PROBABILISTIC, SPATIAL, TEMPORAL) is not needed. What makes an axis
*temporal* is not a property of the axis but of the reads relation
laid along it - every scheme reads backward, so the relation is
acyclic and the level is marched; along a spatial axis the stencil
reads both ways and the level is solved. That is already how
`gti_block` derives its sweep order, and it is the right place for the
distinction.

Grading a mesh is choosing a density: a non-uniform partition of the
physical axis is the equiprobable partition of the pushforward
density, and designing the grid is designing that density through the
mapping's design leaves - the same object as designing a shape.

The basis and the quadrature are chosen for the pushforward density,
not mapped from xi. For an affine mapping the two coincide (Legendre,
Gauss-Legendre). For the normal, F^-1 has an unbounded derivative at
xi = 0 and 1 and mapped Gauss-Legendre points converge poorly; the
family orthogonal under the pushforward is Hermite, and pspace uses it
directly. The coordinate carries neither the measure nor the basis;
both follow from the mapping.

A discretisation of a coordinate is a finite measure on it: points and
weights, agreeing with the continuum measure on a class of functions.

    partition       cells and their measures     exact on piecewise constants     operation_grid + partitioned
    quadrature      Gauss points and weights     exact on polynomials to 2n-1     pspace getQuadraturePointsWeights

These are one map, coordinate -> (points, weights); `partitioned(grid,
n, dt, t)` and pspace's rule per degree are its two concretions. The
mesh is the product of per-axis discretisations followed by the
mapping: `gti_space` already does exactly this (`partitioned` along
xi and eta, then `mapped`, with each cell's measure read from its
mapped corners); pspace's `build_quadrature` is the same product with
the weights multiplied, and no mapping because its axes are
independent. The mapping belongs to the coordinate *system*, not the
coordinate: the polar mesh couples xi and eta, which pspace's per-axis
frames (physical y, standard z, quadrature x) and diagonal covariance
cannot express.

The basis is the third object and is kept apart from both. pspace's
psi_k are polynomials orthonormal under the axis measure (Hermite for
normal, Legendre for uniform, Laguerre for exponential); this
framework's `field_forms` holds the same objects - a table of functions
and their derivatives with an active membership - as `polynomial_form`
and `harmonic_form`, with no measure to be orthogonal under. The
orthogonal families are further concretions of `form`, one per
measure, and the fitted operator gains a well-conditioned basis
wherever the mapping determines the measure. pspace stores the basis degree on
the coordinate (`Coordinate.degree`, `CoordinateSystem.basis`); the
2026-08-21 ruling here keeps basis and coordinates as two branches of
one graph, queried through it, and that separation is kept.

Identity: pspace mints coordinate ids from a counter and preserves them
across `make_cs` so a `PolyFunction` remains valid when the degree
changes; the name is a sympy symbol beside the id. This is
`token_identity` (assigned once) and `map_label` (a name is not an
address), so a coordinate here is identified by token and labelled, and
the expression refers to axes by position in the coordinate system,
never by name.

    coordinate          the uniform axis, token + label            NEW: the input of operation_grid, in place of a bare span
    grid                coordinate -> (points, weights)            EXISTS: uniform, random, designed; NEW concretion: gauss
    coordinate system   product of coordinates + mapping           EXISTS in gti_space as (xi, eta, geometry); the mapping becomes an expression
    measure             pushforward of uniform under the mapping   DERIVED: the mapping's Jacobian; cell measures are its discrete form
    mesh                product of the grids, mapped               EXISTS: gti_space spatial_mesh, view_mesh
    form                basis orthonormal under the measure        EXISTS: polynomial, harmonic; NEW: legendre, hermite, laguerre, chosen by the measure
    expression          derivative(q, alpha), alpha by position    this document

The link to the physics statement: `derivative(q, alpha)` names axes of
the coordinate system the statement is laid on, by position. Along the
marched axis the derivative of degree k is a state component; along a
spatial axis it is compiled to the mesh's stencil; along a
probabilistic axis there is no derivative vertex, since the state is
sampled there rather than differentiated - the design `nu` *is* the
coordinate, and its discretisation is either the quadrature (pspace:
one trajectory per Gauss point, the functional's expectation a weighted
sum) or the Taylor coefficients this framework already carries (the
mixed partials of the functional in `nu` to any degree at one point).
Collocation and Taylor are then two grids on the same coordinate,
compared on the same functional, which is the comparison
`doc/NEURALNETWORK.md` §"Stochastic parameters" asks for.

### Geometry is the mapping

`gti_space` already states the geometry as a mapping of two parametric
coordinates xi, eta on [0, 1] into the domain, as a `select case` over
three formulas. With coordinates and the expression type, the mapping
is an expression in the coordinates with the extents as leaves:

    cartesian    x = a xi,               y = b eta
    circular     x = a xi cos 2 pi eta,  y = a xi sin 2 pi eta
    elliptical   x = a xi cos 2 pi eta,  y = b xi sin 2 pi eta

The polar coordinates are the parametric axes up to scale, r = a xi
and theta = 2 pi eta, so no coordinate is named for them. The measure
on physical space is the pushforward |det D Phi| d xi d eta; the mesh
reads each cell's measure from its mapped corners as a polygon, which
is the discrete form of the same measure, and the continuum form comes
from evaluating the mapping over `derivative_terms` - the metric,
exact, through the same `composed` kernel. `view_mesh_geometry` is
unchanged: the corners are the mapping evaluated at the corner
parameters.

A constant of the mapping promoted to a design leaf gives the
functional's partial in the shape by the same route as its partial in
`nu`. Caveat: the polygon geometry (areas, centroids, normals) is
written over `real(dp)`; a shape partial through the mesh needs it
evaluated over `derivative_terms` too. It is arithmetic on corners, so
it can be; it is not available today. The polar degeneracy (the first
ring collapsed to one cell) is det D Phi = 0 on the face xi = 0, a
property of the mapping and detectable from its derivatives; the
collapse itself stays as written.

So the expression's leaf is *component i of argument j*, not "state
degree" and "design" specifically: the physics reads (state, design),
the mapping reads (parametric point, shape parameters), and
`designed_grid` is the one-dimensional case - a mapping of [0, 1] onto
[0, T] with design leaves.

### Forms: generic function shapes

Two degrees must not be confused. `max_degree` on `operation` is the
order of exact partial action; for an expression it is
`max_subset_width()` and limits nothing in practice. The polynomial
degree is pspace's `Coordinate.degree`, and no coordinate here carries
one: the coordinate is the uniform axis. What shapes a mapping, a
physics coefficient or a grid may be built from is a question about
the *form*, and the tower has the type: `field_forms` holds a family
of functions of position with `values(x, at)` and `slopes(x, at, n)`,
`polynomial_form` and `harmonic_form` are its concretions, and
`operation_fitting`'s `fit` is form + coefficients - the object
G = (form basis, coefficients) of the fitted-expansion slice.

The generic function form is an expansion in a form with the
coefficients as leaves:

    f(xi)  =  sum over k of  c_k phi_k(xi)          c_k design leaves, phi_k the form's members

For derivative propagation every mixed partial of each phi_k is
needed, and `slopes` gives one first derivative. Rather than extend
each concretion by hand, each member of a form is itself an expression
in the position leaves and is evaluated in the same loop: a monomial is
PRODUCT and POWER vertices, a harmonic is sin or cos of a PRODUCT,
`values` and `slopes` are the zero- and one-direction evaluations, and
every higher derivative comes with them. `polynomial_form` and
`harmonic_form` shrink to two constructors of member lists. This is a
candidate merge of `field_forms` into `expression`, to be measured on
the fitting hot path before it is taken.

Two vertex kinds carry the expansion and the piecewise case:

    EXPANSION(fit)     sum of c_k phi_k over the form's active members (restrict = pspace's adaptive basis)
    PIECE(grid)        the cell of the coordinate value selects a sub-expression

A spline is a piecewise polynomial on a grid of the coordinate; a
B-spline member of degree p has support on p+1 cells, so each cell's
sub-expression holds only the members nonzero there, and evaluation is
local without a search. Continuity across breakpoints is a property of
the coefficients, or of the B-spline basis by construction; it is
verified at the breakpoints, not assumed. At a breakpoint the
derivative is one-sided, by the grid's half-open cell convention.

So a geometry may be a mapping whose components are PIECE vertices over
a grid of xi with polynomial sub-expressions and the control points as
design leaves - the CAD form - and its shape partials come by the same
route as the partial in `nu`. A designed time grid t(xi), a spatially
varying coefficient k(x) in the physics, and a pspace basis are the
same object with different leaves.

pspace's `Operation` (forward = decompose, inverse = reconstruct,
residual = the round trip's error) corresponds to `operation_fitting`:
the fit's coefficients are the decomposition, its evaluation the
reconstruction, and `tolerance_form` measures the residual. Nothing
is imported from it; the correspondence is recorded so that the two
codebases are read as one construction.

## Extensions, in order

1. **Space.** `derivative(q, [a, b])`: a vertex kind `SPATIAL(alpha)`, a derivative
   multi-index over the spatial axes of the coordinate system the
   statement is laid on (positions in the axis sequence; labels come
   from the coordinates, not from a fixed list of names). Compile: every `SPATIAL` vertex
   that enters the root linearly is moved into the stencil `gti_space`
   builds today, and the rest of the tree is the nodal rule. A
   `SPATIAL` vertex inside a product with the state (Burgers,
   advection) is not nodal; it needs the stencil's image at the point
   as a third input to the rule, seeded through the stencil's linear
   map. That is one more argument on `nodal_integrand` and one more
   gather in `gti_block`; it is not in the first slice.
2. **Several unknowns and several designs.** `unknown(i)`, `design(j)`;
   the state layout becomes (instant, unknown, degree) and the design
   (instant, design). `gti_block`'s gather and the tangent assembly
   index one more dimension. The expression does not change.
3. **Configured physics.** A rule chosen by name in `gti_configuration`,
   which deletes `physics_vanderpol.f90`.

## What is not built

- No symbolic differentiation, Jacobian or Hessian by tree rewriting:
  the partials come from evaluation over `derivative_terms`.
- No string parser: expressions are built by Fortran operators, so the
  compiler checks them.
- No fixed list of coordinate names (`t, x, y, z, r, theta, phi`): the
  time degree is an integer, the space axes are positions on the mesh.
- No `evaluate(values)` by name lookup: a leaf is a position in the
  state or the design, matched by the arguments `nodal_integrand`
  already declares.

## Verification

Each stays a percentage against theory, in the existing suites.

- `composed`: `exp`, `sin`, `log`, `sqrt` against their closed-form
  mixed partials for n = 1, 2, 3 seeded directions, and `sqrt(x)**2`
  against `x`.
- The stated van der Pol residual against the hand-written one: the
  value, `partial_action` in every degree, and the design partial, to
  the arithmetic's floor; then the hand-written type is deleted.
- `chained_horizon` and `constraint_rows` unchanged in output.

## Line count

Roughly +80 in `util_derivative_terms`, +150 net in the moved
`operation_expression` (the table and operators added, the abstract and
`zero_integrand` deleted), -125 in `physics_vanderpol`. Net about +100,
one type fewer, and new capability - every future equation costs a line
and the elementary functions become available - but the count is
stated here rather than assumed.

## Status, 2026-08-26

Done on tolerance-form-identification-via-graph: affabcb (composition
kernel and elementary functions on derivative_terms, 55 identities
checked to five directions), adcde6b (expression type, van der Pol as
three one-line rules), 8e581b5 (src/operation_expression.f90 as one
concrete operation, physics_integrand and its abstract deleted, 40
sites renamed, check_naming rule 5 recognises extended generics).
Outputs bitwise identical to the unmodified tree on chained_horizon,
randomized_checks and constraint_rows.

Measured line count, not the estimate above: +1227 / -515 across the
three commits, net +712, of which 156 is the identities driver;
operation_expression.f90 is 769 lines (the estimate said ~500) and the
derivative_terms additions 245 (the estimate said ~80). The
extensions - spatial vertex, several unknowns and designs, physics by
name in the configuration - are not started.

2026-08-27: the coordinates section's unification is implemented - the
grid and the quadrature are one map. GRID_GAUSS is a kind of the one
grid type beside uniform, random, designed and fixed: the same
partitioned machinery yields its weights, abscissae yields each kind's
points (a partition's instants, the rule's Legendre nodes), and the
assembled_tower demo holds the exactness law - the n-point rule
integrates every power below 2n over the span to the arithmetic's
floor. The probabilistic axis of the pspace review now has its
collocation half: a quadrature is an ordinary grid to everything that
takes one.
