# Unified residuals: physics, representation and conditions on one manifold

Status: design note, not implemented. Revised 2026-09-28 after the
owner's review. Written for two Euler-Bernoulli beam examples in
`test/graph-scheme`: an analysis, then an optimization of the same
Lagrangian.

The rank measurements quoted here are reproduced by
`doc/unified-residuals-rank.py [P] [N]`.

## Status (2026-09-29)

| Phase | State | Commit |
|---|---|---|
| 1. Positions in every rule | done; variable EI beam order 4 in x, 2 in t, dJ/ds by tangent and adjoint | 1c898b8 |
| (type-bound procedures of the residual) | done | 20b0c1f |
| 2. Geometry factors: `region(interval)`, `points(n)`, `boundary(c, value)` | done | ff210ce |
| 3. Jet component and derivative operator; `jet_form`; representation laws as residuals | done | 13f8cf7, 52ed4b4, 1a2728b |
| 4. Conditions on the faces of the region, with closed fitted differences | done | ddfa240 |
| 4. The open line: the trapezoid rule on the segments, the conditions of both ends in its open rows, closed differences at the explicit moments, Gauss-Seidel on whole lines | done; order 2, stable at 801 x 129 | 1866adc |
| 5. The beam analysis example | done | 102c263 |
| Multipliers of face conditions and of representation laws | done | 1dd35e0, 92f41dc |
| 6. Mixed components: a leaf of orders along the first coordinate and one later coordinate, each defined along space from the time component | done; the Kelvin-Voigt beam, order 4+ in x and 2 in t, dJ/dc by tangent and adjoint | (this commit) |
| 7. The beam optimization example (identification by Newton on the design) | done | bd608da |
| 8. The disagreement estimate `disagreement(solution, other)` | done; effectivity 1.0007 on the open line, 0.35 to 0.78 under closed differences | (this commit) |
| 8. The grid as a design: one grading coordinate of the beamline, Newton on the sum of squares of the estimate's contributions | done; error of J 2.34e-2 to -7.4e-4 at a fixed 17 points; gradient by central differences | (this commit) |

Every commit above was compared with the suites, the demos and a
two-image run against the state before it; the only differences are
the new examples' own output.

## 1. Purpose

Every law of a problem is stated the same way: a residual on its own
domain, paired with its own multiplier, summed into one Lagrangian.

    L = J + <lambda_r, r> + <mu_t, r_t> + <mu_x, r_x> + <lambda_ic, ic> + <lambda_0, bc0> + <lambda_L, bcL>

The kinds of residual are:

- `r`, the physics;
- `r_t` and `r_x`, the representation of the derivatives along each
  coordinate;
- `ic`, the initial conditions;
- `bc0` and `bcL`, the boundary conditions at the two ends.

Each residual is a `continuous_residual` on a manifold, and each has a
multiplier.

- The stationarity of L in the multipliers is the forward problem.
- The stationarity in the unknowns is the adjoint.
- The stationarity in a design field is the optimality condition of the
  design.

Time and space are coordinates of the same kind, and each has boundary
faces.

## 2. Geometry and its discretization

A continuous manifold (a geometry) and its discrete image are separate
objects. A `.geo` file is one instance of a geometry, and a `.msh` file
one instance of a discretized geometry.

| Geometry | Discrete image | Today | After this change |
|---|---|---|---|
| `interval(a, b)` of a coordinate | points on it, ends included | time only (`instants(n)`) | any coordinate, `points(n)` (`instants` renamed, one type for time and space) |
| `region('file.geo')`, or `region(interval(a, b))` | `mesh('file.msh')`, or `points(n)`: the vertex-centred line mesh | gmsh files | both |
| a product of factors, `beamline x timeline` | the product of the images | time x region only | any factors |
| curves, surfaces, volumes from the coordinate algebra | generated points and cells | none | later; the factor interface admits them |

What each factor provides:

- **Its coordinates.** These are numbered with the time first, when
  present, then the others in the order the factors are given.
- **Its boundary faces.** `omega % boundary(c, value)` is the face where
  the coordinate `c` takes `value`, an end of its interval. The face is a
  manifold of the other coordinates, and records its parent, the axis
  and the value.
  - `boundary(t, 0.0_dp)` replaces `boundary(time=0.0_dp)`, and every
    caller migrates.
  - `boundary(x, 0.0_dp)` and `boundary(x, L)` are the two ends of the
    beamline.
- **Its discretization.** Points on an interval are generated in `src/`,
  with no file. Cells come from a mesh.

## 3. Fields: unknown and given

A field on the manifold is one kind of object, whether it is unknown
or given.

| Field | Kind | Source of its values |
|---|---|---|
| `w` | unknown | the stationarity of L |
| `EI`, `rhoA` | given | a coordinate expression, data on the points, or a design field |
| `f` | given | a coordinate expression; a manufactured solution's forcing |
| a design field `d` | unknown | its own stationarity, section 14 |

**Measured gap: rules are evaluated without positions.** The rules on
the whole manifold are applied over the points without their positions,
so a law that reads a coordinate stops.

- `operation_expression`: `values_over` and `evaluated_over` stop with
  "the expression reads a coordinate, but no position was given".
- Conditions already receive a position (`face_relation`).
- Measured: the periodic beam with `EI = 1 + 0.5 cos(x)` and a
  manufactured forcing is accepted when stated, and stops on its first
  evaluation.
- Consequence: variable coefficients, a forcing `f(x, t)` and
  manufactured solutions all require phase 1 (section 12).

## 4. The jet: components as unknowns

Every component of every unknown's jet is an unknown of the discrete
problem.

**Multi-indices.** A component is named by a multi-index
`alpha = (alpha_t, alpha_x, ...)`, the order along each coordinate. The
jet of a field is a set of multi-indices, containing zero and every
multi-index the laws read.

- Today only pure multi-indices are stored: the orders along t, then the
  pure orders along each further axis (`order_at`, `components_per_field`).
  `field_derivative` refuses a mixed one with "a mixed derivative is not
  a component of the jet".
- After this change, mixed components such as `w_xt` are admitted. Only
  the system grows.

**Commutativity is assumed.** The component is indexed by the multi-index
alone, not by the order in which the derivatives are taken, so
`d_x d_t w` and `d_t d_x w` name the same component `w_xt`.

This is the equality of mixed partials for a function with continuous
second derivatives. It is also the reason the Taylor expansion is written

    w(z0 + h) = sum over alpha of  (d^alpha w)(z0) h^alpha / alpha!,
    alpha! = alpha_t! alpha_x! ...

The factor `alpha!` counts the orderings of one multi-index: the
`(|alpha|)! / alpha!` sequences of first derivatives that reach the same
component are one term. The design jets already use this convention.
`util_derivative_terms` stores one coefficient per subset of directions,
and `symmetric_terms` stores one per multiset.

**What commutativity requires of the discretization.** See section 9,
measured.

**A component and an operator are different things.**

- `w % jet_component(alpha)` is the stored unknown, a leaf of a law.
- `u % derivative(alpha)` is the operator `d^alpha` applied to any field,
  stored in the law as a vertex. `jet_form` replaces it by the jet
  components it denotes, once, in `discretize`; the concretion along
  the coordinate may choose another form (section 6).
- Today `w % derivative([t])` is a component: it selects a stored value
  and does not differentiate, so `d_t w - w_t` would read `0 = 0`. The
  representation laws need the operator.

## 5. Laws

**A point law** is algebraic in the components at one point. It may read
given fields and the coordinates.

| Law | Beam | Domain |
|---|---|---|
| Physics | `rhoA w_tt + d_xx(EI w_xx) - f` | Omega |
| Initial | the jet of `w` at t = 0 (section 8) | `beamline x {0}` |
| Boundary, clamped | `w`, `w_x` | each end `x = 0`, `x = L`, times the timeline |
| Boundary, simply supported | `w`, `w_xx` | the same |
| Gauge (Taylor-Green) | `integral of p over the region` | the time factor |

**A representation law** states what a stored component is: the
derivative of the root field along its multi-index.

    q_alpha - d^alpha w = 0,      one law per stored alpha other than zero

The law names no path. It does not say "the next component", nor
through which intermediate components the derivative is reached. The
discretization chooses the path, and it may relate any subset of the
components:

- finite differences give `w_xxxx` from the values of `w` directly;
- Newmark relates `w`, `w_t` and `w_tt` together;
- the time families step one order at a time.

What must hold is the count: one discrete relation per stored non-root
component per point, less the rows left open for the conditions
(section 7).

**The laws are stated as residuals on Omega, like the physics:**

    r   = continuous_residual(omega, [rhoA*w % jet_component([t, t]) + derivative(EI*w % jet_component([x, x]), [x, x]) - f])

    r_t = continuous_residual(omega, [derivative(w, [t])    - w % jet_component([t]),                    &
                                      derivative(w, [t, t]) - w % jet_component([t, t])])

    r_x = continuous_residual(omega, [derivative(w, [x])          - w % jet_component([x]),              &
                                      derivative(w, [x, x])       - w % jet_component([x, x]),           &
                                      derivative(w, [x, x, x])    - w % jet_component([x, x, x]),        &
                                      derivative(w, [x, x, x, x]) - w % jet_component([x, x, x, x])])

    ic  = continuous_residual(omega % boundary(t, 0.0_dp), [...])

    bc0 = continuous_residual(omega % boundary(x, 0.0_dp), [w, w % jet_component([x])])
    bcL = continuous_residual(omega % boundary(x, L),      [w, w % jet_component([x])])

    L   = energy + lambda_r*r + mu_t*r_t + mu_x*r_x + lambda_ic*ic + lambda_0*bc0 + lambda_L*bcL

A law that does not state its representation residuals receives them
generated from the components it reads. This is the present behaviour,
which becomes a convenience instead of a separate mechanism.

**Operators stay unexpanded.** `derivative(EI*w_xx, [x, x])` is kept as
an operator node in the law. Today the total derivative expands it by
the chain rule when the law is written (`total_derivative`,
`differentiated`): `EI w_xxxx + 2 EI' w_xxx + EI'' w_xx`. That fixes the
discretization before `discretize` runs. Kept unexpanded, the
concretion chooses:

- the chain-rule expansion, for finite differences;
- the flux form, `EI w_xx` at faces then differenced, for finite volumes
  and finite elements.

## 6. Discretization: the abstract law and its concretions

The continuous representation residual is the abstract statement. A
scheme is its discrete analog.

    L_h = L % discretize(omega_h, [discretization(r_t, chain([dirk(2), bdf(2), adams(2)], from=[1, 5, 8])), &
                                   discretization(r_x, finite_difference(order=6))])

The present `time=` and `space=` keywords are the case where `r_t` and
`r_x` are implicit.

**Each concretion owns its discrete domain, and the number of rows it
leaves open at each face of its axis.**

| Concretion | Discrete domain | Rows open at the faces of its axis |
|---|---|---|
| Families along t (Adams, BDF, DIRK, Newmark) | the steps between instants | at t = 0: one per component below the top order |
| A family along x (for example the trapezoid rule, or a Hermite relation) | the segments between points | at both ends, together: one per component below the top order |
| Fitted finite differences, closed | the points | none; open rows come from the physics' domain instead (section 7) |
| Finite volumes | the cells, with their faces | none per cell; conditions enter as face values |
| Finite elements | the elements | the essential conditions at boundary nodes |

The points and the steps, segments, cells and faces are the two sides
of the incidence form in `/home/komahan/git/AGENTS.md`. A concretion is
a residual on one side that reads the other.

## 7. Counting

The forward solve is `A w = b`, where A is the Jacobian of all the
constraints. The adjoint solve is `A^T lambda = -grad J`. Both need A
to be square and nonsingular.

The number of multipliers equals the number of constraints, one per
row. What must hold is

    number of constraints = number of unknowns in the jet stack

and A must have full rank.

**Along one coordinate with an open concretion.** Take `P` points and a
chain of degree `D`:

- `D (P - 1)` relations;
- the physics determines the top component at each point;
- this leaves `D` rows open, which the conditions on the faces of that
  coordinate fill.

**With a closed concretion.** The relations fill every point. The rows
for the conditions come from the physics instead, which acts on the
points away from the faces: two fewer at each end for the clamped beam.

**Measured.** Small beam grids, `w_tt + w_xxxx = 0`, clamped, backward
Euler along t:

| Variant | P, N | Unknowns | Rows | Rank | Smallest singular value |
|---|---|---|---|---|---|
| Family along x; initial `w, w_t`; boundary rows at every instant | 6, 4 | 168 | 168 | 165 | 1e-17 |
| The same, initial `w` dropped at the ends | 6, 4 | 168 | 166 | 165 | 2e-16 |
| Family along x; initial jet; x relations and boundary rows for t > 0 | 6, 4 | 168 | 168 | 168 | 1.2e-3 |
| Fitted differences, closed; physics on the points away from the ends for t > 0 | 6, 4 | 168 | 168 | 168 | 5.7e-5 |

The same pattern holds on `P, N = 5, 3` and `7, 6`.

## 8. The initial slice

**With an open concretion along x, the initial data is the whole jet.**

- At t = 0 the initial condition fixes `w` at every point, so the
  relations along x only differentiate known values.
- Each of the four levels (`w_x` from `w`, up to `w_xxxx` from
  `w_xxx`) is a relation on `P - 1` segments for `P` unknowns, and so
  needs one datum.
- The clamped conditions give two data at the first level (one
  surplus) and none at the other three. That leaves three components
  undetermined, as the rank deficiency of 3 in section 7 shows.
- The rule: the initial condition prescribes the jet of `w` along x at
  t = 0, and the x relations and boundary conditions act on
  `beamline x (0, T]`. The corners belong to the initial data.
- The jet comes from the coordinate algebra, differentiating `w0(x)`
  exactly, or from a closed concretion applied to its values, as the
  application's `consistent_states` computes.

**With a closed concretion along x, no rule is needed.** It determines
the x components from the values at every instant, t = 0 included. The
initial condition states `w` and `w_t` alone (section 7, last row).

**Compatibility.** The data must satisfy the boundary conditions at the
corners, for example `w0(0) = 0`.

## 9. Mixed components and commutativity (measured)

`w_xt` added to the jet of the family variant with the initial jet:

| Defining law of `w_xt` | P, N | Unknowns | Rows | Rank |
|---|---|---|---|---|
| Along t, from `w_x` (open at t = 0: `w_xt(x, 0)` joins the initial jet) | 6, 4 | 192 | 192 | 192 |
| Along x, from `w_t` (open at an end: one condition per instant) | 6, 4 | 192 | 192 | 192 |
| Both | 6, 4 | 192 | 212 | 192 |

- **Either path alone gives a square system of full rank.**
- **Both paths add 20 rows and no rank.** Every added row is a
  combination of the others. On a product of axes, the concretion along
  t acts on the index of the instant and the one along x on the index of
  the point. Their matrices are Kronecker factors, `D_t (x) I` and
  `I (x) D_x`, and so commute. The discrete jet is then commutative too,
  and the second path repeats the first.
- **Rule: exactly one defining law per stored component.** The path is
  the concretion's choice.
- On an unstructured region, the fitted operators along two axes do not
  commute. The second path would differ from the first by the truncation
  error and contradict it. A single law, the mixed fit `q_alpha` from
  the values directly, remains the rule.

**Implemented.** A leaf of the jet carries a mixed order beside its
order: the component of order `a` along the first coordinate and `o`
along one later coordinate. `w % derivative([x, x, t])` and
`w % jet_component([x, x, t])` name `w_xxt`; the entries commute.

- **Layout.** Each field stores, after its pure components, one block of
  the later coordinates' orders for each mixed order `a = 1 .. A`. A law
  with no mixed leaf has the layout it had.
- **The path.** The concretion along space defines the mixed component
  of orders `(a, o)` as the difference of order `o` of the time
  component `w_(a)`, the path along x of the table above. The families
  along time state the pure time components alone. The rows are those
  of the pure components with the source column moved from `w` to
  `w_(a)`.
- **Refused.** A representation law of a mixed component, the
  stationarity and the partials in a mixed component, a mixed component
  of two later coordinates, and the open line (the chain of `w_t` along
  the line would need conditions at the ends).
- **Measured** on the clamped Kelvin-Voigt beam
  `rhoA w_tt + d_xx(EI w_xx + c w_xxt) = f`, `c(x) = c_0 (1 + (x/L)^2)`,
  whose law reads `w_xxt`, `w_xxxt` and `w_xxxxt`, with fitted
  differences of order 6:

  | grid (instants x points) | error |
  |---|---|
  | 801 x 9, 17, 33 | 3.06e-2, 4.07e-4, 9.03e-6 |
  | 21, 41, 81, 161 x 65 | 8.91e-4, 1.51e-4, 3.13e-5, 8.33e-6 |

  One Newton step per block. With `c_0` the design, `dJ/dc_0` of the
  misfit at `c_0 = 0.5` on 21 x 17 is 0.1096767998 by the tangent,
  0.1096767982 by the adjoint, 0.1096764971 by the central difference
  of step 1e-3, whose own error `h^2/6 J'''` is -3e-7.

## 10. Multipliers and the adjoint

Each residual's multiplier is returned under its own name, at the
points of its domain.

- `lambda_r`: the adjoint of the physics, as now.
- `lambda_ic`: the sensitivity of J to the initial data, as now.
- `lambda_0` and `lambda_L`: the sensitivity of J to each boundary
  datum. This is the mechanical reaction (shear force, bending moment)
  only when L is the action of the beam. With J the energy integral that
  reading is not demonstrated, and the example prints both.
- `mu_t` and `mu_x`: the multipliers of the representation laws, returned
  under their names. Law k along time is the family's row of the k-th
  component it determines; law k along a space coordinate is the
  difference row of that component, or the row it moved to where a face
  condition took it. Measured on the simply supported beam:
  `mu_x(4) = -EI adjoint` at every point after the first instant, to
  zero difference, as the adjoint equation of the column `w_xxxx` states.

## 11. What must not change

The following must reproduce their present output, with timing lines
excluded:

- the time families and the chain, with the same rows and blocks;
- the periodic Taylor-Green problems (closed fitted differences, no
  boundary conditions);
- every `verify.sh` suite except accuracy-contract (its D07/D08 failure
  is independent), and the demos.

## 12. Phases, each ending in a comparison against the commit before it

1. **Positions in every rule.**
   - The whole-manifold rules receive the position of each point.
   - Check: every output unchanged, and the periodic beam with
     `EI = 1 + 0.5 cos(x)` and its manufactured forcing converging.
2. **Geometry factors.**
   - `interval` for any coordinate, `points(n)`, the generated
     one-dimensional mesh, and `boundary(c, value)`.
   - Migrate `boundary(time=...)`.
   - Check: every output unchanged.
3. **Component and operator.**
   - `w % jet_component(alpha)` and `derivative(u, alpha)`.
   - Operator nodes kept unexpanded.
   - Representation residuals as `continuous_residual`s, generated when
     not stated.
   - Along t, the rows the families build today.
   - Check: every output unchanged, and `mu_t` reported.
4. **Concretions owning their domain and open rows.**
   - A family along x on segments, the two-point rows within each time
     block, and the initial-slice rule.
   - Fitted differences with the physics on the points away from the
     faces.
   - Check: the rank on small grids, and the beam converging at the
     concretion's order.
5. **The beam analysis example** (section 13).
6. **Mixed components**, with one defining law each.
   - Check: the rank table of section 9 on the assembled system.
7. **The beam optimization example** (section 14).
8. **Later: disagreement estimates and the grid as a design**
   (sections 15 and 16).
   - The estimate `eta` checked against the known error of the beam,
     with its effectivity measured.
   - Then the grid moved by the gradient of `eta`.

## 13. The analysis example

The equation, with the material and geometric parameters named in the
program:

    rhoA(x) w_tt + d_xx( EI(x) w_xx ) = f      on [0, L] x [0, T]

- **Simply supported, uniform `EI` and `rhoA`, `f = 0`.** The exact
  solution is `w = sin(pi x / L) cos(omega t)`, with
  `omega = (pi / L)^2 sqrt(EI / rhoA)`.
- **Clamped, with `EI(x)` and `rhoA(x)` varying.** Checked against a
  manufactured solution that satisfies the conditions, with the forcing
  it implies.
- **Objective.** `J = integral over Omega of (rhoA w_t^2 + EI w_xx^2) / 2`.
- **Printed:**
  - the error;
  - J against its closed form, where one exists;
  - the boundary multipliers;
  - the observed orders along x and along t.

## 14. The optimization example

The same Lagrangian, with a stiffness that is a design field.

**The design as a field.** `EI(x)` does not vary in time. It is an
unknown field on Omega with two laws:

- the representation law `d_t EI = 0` on Omega;
- the initial condition `EI - EI0(x) = 0` on `beamline x {0}`.

Nothing new is needed. The multiplier of the initial condition is the
functional derivative `dJ / dEI0(x)`, the gradient of the objective in
the design, delivered by the reverse sweep that computes every other
multiplier.

**The stationarity in the design.** Setting `dL / dEI0 = 0` is the
optimality condition. It is solved in the example by an outer iteration
on `EI0`:

- each iteration is one forward and one reverse sweep;
- the step comes from the gradient.

The one-shot alternative, with the design stationarity among the rows of
one system, is left open.

**The objective must have a well-posed optimum.** Two candidates:

- minimum compliance under a fixed mass, with the mass constraint
  `integral of rhoA = m0` as one more residual with its own multiplier;
- the identification of `EI(x)` from a measured deflection, with a
  regularization term.

The choice is the owner's.

**Size.** The gradient is available once phases 1 to 5 exist. The
example then adds the design field's two laws, the objective and the
outer iteration. It is a short extension of the analysis example, not
of the code as it stands today.

## 15. Disagreement between representations as an error estimate

**The principle.** Two representations of the same law disagree by an
amount of the size of their truncation error. That disagreement is a
measured quantity, and it vanishes as the grid is refined. Instances:

| Two representations of one law | Disagreement |
|---|---|
| Two paths to a mixed component (section 9) | the discrete commutator `(D_x D_t - D_t D_x) w` |
| Two concretions of one derivative | for example `finite_difference` against `finite_volume` |
| One concretion at orders `p` and `p + 1` | the lower order's truncation error, to leading order |
| A Taylor truncation at orders `n` and `n + 1` | the first term omitted |

The code already uses one instance. The accuracy contract's
functional-error estimate (cases G22 and G23) estimates BDF-3 by BDF-4
and checks the effectivity ratio.

**In the Lagrangian.** The disagreement enters L as a resolution term.

- It can be a weighted soft term, `1/2 || q_alpha(one) - q_alpha(other) ||^2`.
- Or it can be a slack field with its own multiplier.

The stationary point then divides the disagreement between the
representations, instead of taking one path and discarding the other.
The reverse sweep delivers the multiplier of the resolution term, which
is the adjoint weight of the disagreement.

**The goal-oriented estimate.**

    eta = sum over the points, steps or cells of  (multiplier x local disagreement)

This is a dual-weighted residual estimate of the error in J. It is
localized to the discrete domain of the representations, which is what
goal-oriented refinement needs. It comes from the same stationarity as
every other multiplier, with no separate estimation machinery.

**Limits.**

- **It is blind where the two representations agree exactly.** On a
  product grid the commutator is zero, so it must be paired with a
  representation whose error it can see, such as `p` against `p + 1`.
- **It measures the difference of two errors, not either one.** The
  two can partly cancel. Its effectivity is measured, as the accuracy
  contract does for G22 and G23, and not assumed.

**Measured on the beam.** `L_h % disagreement(solution, other)`
returns one contribution per point: for each representation law along
space, minus its multiplier times the law evaluated on the jet that the
concretion `other` forms from the values of w. At a point the law is
`D^k w - w_(k)`. On the open line it is the trapezoid relation of the
segment, evaluated on that jet, at the segment's upper point. The
mode of the simply supported beam, where `J - J_h` is known:

| concretion | estimator | grid (instants x points) | eta / (J - J_h) |
|---|---|---|---|
| open line | differences of order 6 | 801 x 9, 17, 33, 65 | 1.0435, 1.0101, 1.0018, 1.0007 |
| open line | differences of order 8 | 41 x 17 | 1.009 |
| differences of order 4 | order 6 | 801 x 9, 17, 33 | 0.37, 0.47, 0.62 |
| differences of order 4 | order 8 | 801 x 17, 33, 65 | 0.46, 0.62, 0.76 |
| differences of order 6 | order 8 | 41 x 17; 801 x 17, 33 | 0.78; 0.41, 0.35 |

- **The open line is estimated to 1 %.** Every row of a law carries its
  multiplier; the conditions of the ends take the open rows.
- **Under closed differences the estimate is partial.** Next to each
  end a condition takes the row of a representation law.
  The mismatch of that law then has no multiplier. At 801 x 17 the
  points beside the ends contribute 1.7e-6, the next ones -3.8e-4,
  the interior at most 1.8e-5.
- **Time is not estimated.** The laws along time are those of the
  families; `r_t` has no second concretion here.

**The structural alternative.** Mimetic, or discrete-exterior-calculus,
concretions are built so that the discrete derivatives commute exactly.
There the paths agree by construction, and commutativity is a property
required of the concretion, not one resolved afterwards.

**For the demos.** The concretions choose one path (section 9).

## 16. The grid as a design

The positions of the points and of the instants can be design fields,
like `EI(x)` in section 14. Starting from a coarse grid in space and
time, the Lagrangian's stationarity in the positions moves them.

- Each position has the representation law `d_t x_i = 0` along the
  optimization, and an initial condition, as `EI` has.
- The gradient of the objective in each position comes from the
  reverse sweep.
- The update is the same outer iteration as in section 14.

**The objective cannot be J_h alone.** The application's
`grid_stationary_partition` (`adaptive_check = grid_stationarity`)
already measures the stationarity of the discrete functional in the step
sizes, and records the failure mode.

- A grid on which `F_h` is stationary need not be accurate.
- Implicit midpoint on the oscillator gives
  `F_h = sum_k h_k / (2 (1 + h_k^2 / 4))`. This is stationary on every
  uniform grid, while `F_h - F = -1/65` at `h = 1/4` (`test/gti-contract`).
- Minimizing `J_h` over the grid can therefore move the grid towards
  where the discrete functional is wrong in a favourable direction.

The objective is instead the estimate `eta` of section 15, the estimated
error in J, at a fixed number of points.

**Constraints on the grid.** Each is one more residual with its own
multiplier.

- **Order:** `x_(i+1) - x_i > 0` and `t_(k+1) - t_k > 0`, an inequality
  whose treatment (a barrier, or a slack with a sign constraint) is
  open.
- **The ends stay on the geometry:** `x_0 = 0`, `x_n = L`, `t_0 = 0`,
  `t_N = T`. These are conditions on the boundary faces, as for the
  solution.
- **The scale:** the direction that scales every step by one factor. The
  application's comment notes that `F_h` is homogeneous of degree zero
  in the step fractions, so the gradient along this direction vanishes
  identically. That direction is fixed, or projected away, as
  `grid_stationary_partition` does.

**Measured on the beam** (`euler_bernoulli_beam_grid`). The clamped
beam on the open line with the solution `sin(pi (x/L)^2)^2 cos(omega t)`,
whose curvature grows along the beamline, and the goal `J = integral of
w^2 / 2`. One design coordinate `s` grades the points,
`x_i = L (xi_i + s xi_i (1 - xi_i))`; the ends stay on the geometry and
the points ascend for `|s| < 1`, so neither constraint needs a
multiplier, and the scale direction does not arise.

- **The objective is `Phi = sum of eta_ik^2`** over points and instants.
  - `eta` itself changes sign: its zero, near `s = 0.41`, is a grid
    where contributions cancel.
  - `sum of |eta_ik|` has its minimum near the same grid but is not
    differentiable.
  - `sum over i of |sum over k of eta_ik|` does not follow the error:
    it is smaller at `s = -0.3`, where the error of w is twice the
    uniform grid's, than at `s = 0`.
- **Newton's iteration on `dPhi/ds = 0`** from the uniform grid, with
  steps bounded by 0.2, converges in four iterations to `s = 0.372` at
  101 x 17 (0.379 at 201 x 17):

  | s | Phi | J - J_h | error of w |
  |---|---|---|---|
  | 0 | 1.46e-5 | 2.34e-2 | 5.49e-2 |
  | 0.372 | 1.03e-5 | -7.4e-4 | 3.14e-2 |

- **The gradient is by central differences** of step 0.02, three solves
  per iteration. The gradient through the reverse sweep needs a second
  adjoint, since `eta` reads the multipliers, and the derivatives of the
  concretion's weights in the positions: the segment relations are
  linear in the lengths, the fitted differences are not.
- **The estimate fails on grids graded the other way:**
  `eta / (J - J_h) = 0.37` at `s = -0.6`.

**Relation to refinement.** Moving the positions (r-adaptation) keeps
the size of the system, while refinement (h-adaptation) changes it.
Both read the same localized `eta`: refinement splits where the
contribution is largest, and moving the grid transfers size between
cells by the gradient.

## 17. Open questions

1. **The family along x of order above two.** A two-point Hermite
   relation reading the higher components, or collocation on each
   segment.
2. **Several spatial coordinates on an unstructured region.** Open rows
   and conditions on its faces. The note covers products of intervals,
   and regions with closed concretions.
3. **A condition at `t = T`**, already a two-point problem in time,
   together with boundary conditions along x. Not measured.
4. **The objective and constraint of the optimization example**
   (section 14).
5. **Names.** `component`, `derivative`, `discretization`, `points` and
   `boundary(c, value)` are proposals, subject to the naming law and its
   hook.
6. **Stability of the fitted differences next to a face (measured).**
   Order-6 fitted differences, one-sided at the ends, give a
   semi-discrete operator with growing modes on fine grids:
   - the simply supported beam at 801 instants and 129 points diverges
     (error 6.6e77) at order 6 and converges (5.1e-5) at order 4;
   - the clamped case diverges at the same grid from the Adams-Moulton
     blocks on, with no damping in the family to hide the growth;
   - the periodic problems have no end and do not show it.

   Candidates: a closure of lower order at the ends, a closure that
   satisfies summation by parts, or the open family along x of phase 4.
   The open line (`segment_trapezoid`) is stable on the same grids,
   8.64e-5 at 801 x 129, of order two; a relation of higher order on
   the segments (Hermite) is the next step for its accuracy.
7. **Which of two disagreeing paths holds.** There are three answers:
   - choose one (the demos);
   - resolve both through a term in L (section 15);
   - require concretions that commute by construction.

   Which one is the default for regions where the concretions do not
   commute is open.
8. **The resolution term's form and weight.** A soft penalty or a slack
   with a multiplier, and the scale of its weight relative to J.
9. **The inequality constraints of the grid as a design** (section 16),
   and the gradient of `eta` in the positions through a second adjoint.
10. **Mixed components of two later coordinates.** A leaf carries one
    mixed order, along the first coordinate, so `w_xy` is stored on a
    region without time and not with it. A multi-index of one order per
    coordinate is the general form.
