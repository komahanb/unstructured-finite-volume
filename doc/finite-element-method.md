# The finite element method

Status: implemented, 2026-10-02: on a line the Hermite cubic, for the
beam's physics and its identification; on the cells of a mesh of two
and three coordinates the Lagrange elements of degree one, for the heat
equation. The note continues
`doc/unified-residuals.md`: section 5 keeps the derivative operators of
a law unexpanded so that the concretion chooses its form, and section 6
lists finite elements with the essential conditions at boundary nodes
as their open rows.

## The method

1. The finite element method is the Lagrangian restricted to a basis:
   the unknown and its multiplier are expanded in the same basis, and
   every integral is a sum over the quadrature points of the cells. An
   inner product is an integral, and an integral is a quadrature.
2. The law is written once, in flux form. Finite differences expand an
   outer derivative of a flux by the chain rule; finite elements move
   it onto the multiplier, so that each row is the law tested with one
   basis function.
3. The coefficients are components of the jet at the nodes: the value
   for a Lagrange element, the value and the slope for the Hermite
   cubic of the beam. The components the law reads above them are
   formed from the basis at the quadrature points, so the
   representation laws along space hold without rows.
4. An essential condition fixes a coefficient in place of that
   coefficient's row. A natural condition belongs to the boundary terms
   of the integration by parts; it is not yet stated (see Scope).
5. The transpose of the Galerkin system is the Galerkin system of the
   adjoint. Along time the adjoint remains the transpose of the family.

The architecture is that of the semi-intrusive stochastic Galerkin
method (dissertation, section 7.3), with space in place of probability:
the law at a point is the deterministic element, the basis and the
quadrature of the cells are the space, and the projection evaluates
the law at the quadrature points and tests it against the basis. A
point of a collocation is an element of one node, whose basis is the
identity on the stored jet.

On the beam, with m = EI w_xx,

    <lambda, rhoA w_tt + m_xx - f>
        = integral of lambda (rhoA w_tt - f) + lambda_xx m  dx
          + [lambda m_x - lambda_x m] at x = 0 and x = L.

The clamped ends fix w and w_x, the two coefficients of each end node,
so the bracket enters no row.

## Three choices the implementation made

- **The test function of a coefficient is its basis function divided
  by its moment of the coefficient's order about the node:** the
  integral of the basis function for the value, of the basis function
  times (x - x_i) for the slope. The row of the value is then the law's
  average about the node, and the row of the slope the law's slope
  there. The rows have the scale of the rows of finite differences, and
  the multiplier of the value's row is the adjoint times the node's
  volume, as under finite differences. Unnormalized, the stiffness
  (3e4) and the slope's mass (3e-5) put the smallest singular value of
  a block near 1e-13. dJ/ds by the adjoint and by the tangent then
  agreed to 1e-5 relative; normalized, they agree to 1e-12.
- **The initial data are projected onto the functions that satisfy the
  essential conditions.** At the first instant a clamped end takes the
  rows of the coefficients it fixes in place of the projection's. With
  the projection unconstrained, the damped beam converged at order 2.8;
  constrained, at order 4.
- **On a line, a block of the Gauss-Seidel preconditioner is the whole
  line at one moment,** as on the open line. The mass matrix couples
  every node of a moment, so blocks of one node leave Newton's
  iteration stopping with the slopes unconverged. With the line as the
  block, the solves of the default tolerances are accurate. The line is
  solved whole on every image.

## The beam programs

`physics/beam/euler_bernoulli_beam.f90` has two cases more: case 6 is
case 2 by finite elements, and case 7 is case 5, the Kelvin-Voigt
damping, by finite elements.

```diff
+    else if (by_elements) then
+       L_h = Lagrangian % discretize(omega_h, &
+            & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
+            &            from=[(k, k = 1, num_instants)]), space=finite_element(degree=3))
```

`physics/beam/euler_bernoulli_beam_optimization.f90` identifies the
stiffness by finite elements:

```diff
-  use operation_finite_difference, only : finite_difference
+  use operation_finite_element   , only : finite_element
 ...
-       &            from=[(k, k = 1, num_instants)]), space=finite_difference(order=6))
+       &            from=[(k, k = 1, num_instants)]), space=finite_element(degree=3))
```

Nothing else in either program changes.

## Measured

The error against the exact solution at 801 instants, over 9, 17 and
33 points:

| Case | Finite elements | Finite differences |
|---|---|---|
| 2, clamped | 1.37e-4, 8.47e-6, 4.70e-7 | 6.68e-2, 8.72e-4, 2.14e-5 |
| 5, damped | 4.52e-5, 3.49e-6, 3.57e-7 | 3.06e-2, 4.07e-4, 9.03e-6 |

Under finite elements the error falls by 16 per halving: order 4. At 65
points over 21 to 161 instants the error of case 6 falls as 1.69e-3,
2.50e-4, 4.14e-5, 8.15e-6, the order of the families in time.

The identification, from s = 0: the gradient falls as 1.20e-1, 2.20e-2,
1.57e-3, 3.74e-5, 7.13e-7, 2.00e-11, quadratically. The design found is
s1 = 0.51178039, 0.50183579, 0.50033454 on 21 x 17, 41 x 33 and 81 x
33, where finite differences found 0.50961827, 0.50205024, 0.50054911.
At the true design on 21 x 17, dJ/ds by the adjoint and by the tangent
agree to 4e-12 relative. The central differences of step 1e-4 agree to
1e-6, which is their own error.

The element's audit: on a cell of length 0.5 its mass and stiffness are
the Hermite cubic's, h/420 [156 22h 54 -13h; ...] and 1/h^3 [12 6h -12
6h; ...], to 2.8e-17 and 3.6e-15.

Every digit of the other cases is unchanged, serially and on 1, 2 and 4
images. The finite element cases give the same digits on 1, 2 and 4
images, but for dJ/ds at the design found, which differs at 1e-11, as
it does under finite differences.

## What changed in src/

- `operation_derivative_approximation.f90`: a space method answers
  three questions: the relations among the stored components
  (`derivative_rows`), whether it projects the law on the cells
  (`projects`), and the elements (`elements`), the cells' nodes,
  quadrature points, weights, the basis's jet and the test functions.
- `operation_finite_element.f90`, new: `finite_element(degree)`. On a
  line the cells are the segments between consecutive points; in two
  and three coordinates, the cells a mesh of corners keeps. The basis
  is dual to the coefficients, a column of the inverse of the matrix of
  the coefficients of the polynomials of `polynomial_form`, factorised
  by `dense_factorisation`. The quadrature is Gauss-Legendre,
  `gauss_node_weight` of `operation_grid.f90`, now public.
- `view_mesh_builder.f90` and `operation_manifold.f90`: `corners(m)`,
  the corners of a mesh as the points of a region; a face of a region
  read from a geometry file at a coordinate's value.
- `operation_expression.f90`: `flux_form`, the law's terms grouped by
  the derivative of the multiplier that tests them; `coefficient_form`,
  the layout of the coefficients a node stores.
- `operation_residual.f90`: the residual operator projects the law on
  the cells when the space method gives elements. The rules are
  evaluated at the quadrature points, on the jet the basis forms from
  the nodes' coefficients. The values are tested into the rows of the
  coefficients, and the Jacobian is the elements' sum. The coefficients
  couple through the elements, and the parts of several images read
  the nodes of their instances. The sweeps relate each coefficient's jet
  in time by the family. The conditions take the coefficients' rows. The
  objective is integrated at the quadrature points, and its
  differential is returned to the nodes by the transpose of the basis.
  The sources of the expansions and the design products are evaluated
  at the quadrature points. The Jacobian is summed element by element
  and written once per element.

## Scope

- One law per unknown, as elsewhere.
- A law whose unknown keeps more derivatives along space than the
  element's continuity allows is refused, as is a flux differentiated
  along two coordinates.
- A natural condition, a condition on a derivative above the element's
  order, is refused rather than satisfied through the boundary terms.
  The simply supported beam, whose w_xx = 0 is natural, is not yet
  stated by finite elements.
- A condition on the time factor (the integral of an unknown), an
  objective on a face, and the checkpoint and snapshot schedules are
  refused under a projection.
- The representation laws along space are refused: the basis states
  them, and `disagreement` reads their multipliers.

## Two and three coordinates

The points of the region are the corners of a mesh read from a file,
and its cells are the elements:

```fortran
omega_h = omega % discretize(time=points(n), space=corners(mesh('square_q8.msh')))
L_h     = L % discretize(omega_h, time=chain(...), space=finite_element(degree=1))
```

- `corners(m)` is the mesh whose cells are the corners some cell of m
  reads, each with a share of the volumes of the cells around it, its
  faces joining the corners of each cell. It keeps the cells of m as
  the elements of their corners. A periodic mesh is refused: its
  corners on the two sides of a periodic face are not identified.
- The elements are of degree one: triangles and tetrahedra with the
  polynomials of degree one, quadrangles and hexahedra with the
  bilinear and trilinear ones. The basis on the reference cell is dual
  to the values at its corners, by the same factorisation as on the
  line, and each cell is the image of the reference cell under its own
  basis. The quadrature is the tensor product of Gauss-Legendre rules
  of three points, collapsed onto the triangle and the tetrahedron.
- A face of a region read from a geometry file is where a coordinate
  takes a value; its points are the corners that lie there. A corner on
  several faces belongs to the condition stated first.
- Over several images the part of an image states the rows of its own
  corners and reads the corners of every cell around them.

`test/graph-scheme/heat_conduction_study.f90` states the heat equation
u_t - div(k grad u) = f, k = 1 + x/2, on the unit square and the unit
cube, with u = 0 on every face and f such that u = sin(pi x) sin(pi y)
cos(t) is the solution, times sin(pi z) on the cube. The error in L2
over [0, 1/2] x B:

| Cells | n = 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| quadrangles, 81 instants | 1.98e-2 | 4.94e-3 | 1.23e-3 | 3.08e-4 |
| triangles, 81 instants | 5.16e-2 | 1.37e-2 | 3.47e-3 | |

| Cells | n = 2 | 4 | 8 |
|---|---|---|---|
| hexahedra, 41 instants | 6.41e-2 | 1.54e-2 | 3.81e-3 |
| tetrahedra, 41 instants | 1.67e-1 | 6.58e-2 | 1.92e-2 |

The error falls by 4 per halving, order 2, the order of the elements of
degree one. On the tetrahedra, six to each cube of an extruded mesh, it
falls by 2.5 and then 3.4, still approaching that order at n = 8. The
digits are the same on 1, 2 and 4 images.

## Next

- The element interface of the dissertation as one type: the point of
  a collocation as the element of one node, so that the placement at
  the points and the projection are one code path.
- Natural conditions through the boundary terms, and elements of
  higher degree in two and three coordinates, which need the incidence
  of the cells with their edges and faces.
- Found in passing, read and not run: with no time chain,
  `quadrature_nodes` gives every node the weight one, which omits the
  cell volume from a steady objective under finite differences.
