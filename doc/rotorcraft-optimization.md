# The rotorcraft optimization in the Lagrangian framework

Status: plan, 2026-10-02; the rigid rotor built, 2026-10-03, in
`physics/rotorcraft-dynamics/` (section "The rigid rotor"); inequalities
on the design in `minimize`, 2026-10-03 (stage 4). The source
is komahanb/rotorcraft-optimization, `deterministic-optimization`:
Python scripts that build the case in TACS (`tacs_builder.py`,
`scripts/multiscenario.py` and its three scenarios). The facts below
are read from those scripts, the TACS sources the clone vendors, and
its run log (`logs/multiscenario.log`).

## The rigid rotor

Built 2026-10-03, with the decisions of the author: each blade at the
azimuth of its cap; the blade model chosen at setup, rigid first, then
beam, plate and solid; the in-house minimizer in place of ALGENCAN; the
time history drawn for paraview.

**The statement.** `multibody_dynamics.f90` states a rigid body and its
joints as fields of time. A body's action is its kinetic energy less
its potential energy, with the unit norm of its Euler parameters paired
with a multiplier; a joint is a list of conditions, each paired with a
reaction. The Lagrangian density of the mechanism is the sum of the
actions and of the products of the conditions with their reactions,
and its stationarity in each unknown is that unknown's law: 355
unknowns of time, 26 bodies of eight and 147 reactions.
`rotorcraft_dynamics.f90` states the 26 bodies (the table
`rotor_bodies.txt`, read from the case's `data/inp`) and the 34 joints
of `collective.py`. A blade is a plate of span 1.56 from the radius
0.44, chord 0.121, density 2500 and the thickness of the setup, along
the azimuth of its cap and linked to it.

**Two rotations agree** where the vector part of the rotation between
them vanishes, eta_a eps_b - eta_b eps_a - eps_a x eps_b. The rigid
link and the hub's driver are stated so: three conditions independent
at every orientation. Equal parameters eps_a = eps_b leave eta to the
unit norm, whose derivative 2 eta vanishes at half a turn; the block of
that instant is singular. TACS adds terms of a factor delta to the unit
norms for the same reason.

**The first instant.** The initial data fix the positions and the
rates, so a joint's condition reads nothing the first instant leaves
free. The library states such a law there by its derivative along time
of the least order that reads the accelerations: the condition hidden
in its holding at every instant, which determines the reactions at the
first instant.

**The march** is bdf(1) over 360 steps of the revolution, as TACS
marches it, each step from the solution of the step before
(`minimize(..., continued=.true.)`).

**The setup**, `rotorcraft.nml`: the blade model, the scenario
(collective, longitudinal, lateral), the thickness, the revolutions,
the steps, the directory of the bodies' meshes and the stride of the
instants drawn.

**Verified.**
- The triple pendulum of the dissertation (`triple_pendulum.f90`):
  the joints' conditions below 1e-15 at every instant; the change of
  the energy at t = 2 is -6.74e-2, -1.16e-2 and -2.20e-3 at 401, 801
  and 1601 instants under bdf(2), approaching order two.
- The rotor: the initial energy 4.99144494e5 J, TACS's 4.99144e5; the
  hub's angle to 1e-16 and every condition of the joints to 2e-14 at
  every instant of the revolution.
- Against TACS with the same rigid blades (36 steps, BDF of order
  one): the first step agrees to 5.7e-10. At every step of the
  revolution the motion here satisfies the Euler-Lagrange equations of
  the bodies, projected on the conditions, to 3e-10 relative,
  evaluated independently; TACS's motion misses them by up to 4.5e-3.
  TACS's rigid link applies the moment of its position reaction on
  body A with the wrong sign (`TACSRigidLink::addResidual` adds
  `addEMatTransProduct(-1.0, t, ...)` where the gradient of
  `uB - uA - t + C_A^T t` gives +1). Every other element of TACS
  reacts along the gradient of its condition, to 1e-10. With that sign
  reversed in the equations, TACS's states satisfy them to 3e-10 at
  every step and these miss them by up to 5e-3: the sign is the whole
  difference. The upper swashplate's turn about the shaft drifts from
  it: the energies after a revolution are 4.739e5 J here and 4.719e5
  in TACS. TACS's runs of the case have the error in the caps' links,
  to the pitch links and to the blades.
- The sensitivities (`triple_pendulum_sensitivity.f90`): the pendulum
  with a design in a mass and one in a joint's condition, which the
  first instant states by its second derivative. dJ/ds by the adjoint
  and by the expansion agree to 2.9e-13 and 2.3e-14, 7.6e-15 and
  2e-16, where J itself varies by 6e-15 between designs 1e-9 apart,
  and central differences of step 1e-4 to 1.5e-9 and 8.3e-10, their
  own error; the same digits on 1, 2 and 4 images.
- 360 steps: the energy falls to 4.94255e5 J over the revolution, where
  36 steps lose 2.5e4: the dissipation of bdf(1). 168 seconds, 232 MB.

**The drawing.** Every rigid body is drawn with its surface mesh (the
case's `data/mesh`, the corners of its CQUAD9 cells), each blade with
the box of its plate, moved with the body to each instant drawn: x =
r0 + u + C^T (X - r0). One file of 168326 points per instant,
`rotor_<instant>.vtu`, and `rotor.pvd`, the time series paraview reads.

**What changed in src/.**
- `operation_residual.f90`: the laws at the first instant (above);
  `minimize(..., continued)`; the rule formed in one pass.
- `operation_expression.f90`: the stationarity of a Lagrangian that is
  a sum of multipliers times laws is the law the multiplier multiplies,
  read on the Lagrangian's tuple, where before each row evaluated the
  whole Lagrangian; the symbolic derivative skips a subtree that does
  not read its variable; the offsets of the fields are stored; the
  rule is evaluated on one array of coefficients.
- `util_derivative_terms.f90`: each rule on arrays of coefficients,
  the operators of the type those rules on its own.
- `operation_field.f90`: `at` of a field that is not an integral is
  its value at every point of the solution.
- `operation_dense_direct.f90`: a pivot is singular at or below the
  tolerance times the largest entry of its column; the reactions of
  the index-three system have pivots near the square of the step.
- `view_paraview_writer.f90`: a writer of cells alone, the file in
  binary, the time series; `view_bdf_loader.f90`: the surface of a
  Nastran bulk data file.

Every value of the 92 cases is unchanged, serially and on 1, 2 and 4
images. The count of nonzeros in Newton's report falls, the explicit
zeros of the rows no longer listed. The laws of the rotor are formed
in 1.1 s, where they took 225 s, and a step takes 0.46 s.

**Decision before the beam.** A blade that turns a full revolution is
not a linear beam, plate or solid in fixed axes. Either it is
geometrically exact in fixed axes, as in the dissertation and TACS,
linked to its cap at the root; or it is linear in the axes of its cap,
its kinetic energy read in those turning axes and its stiffening by
the axial stress of the rotation added to its potential energy.
Recommended: the second. It keeps the beam, the plate and the solid
linear and alike, one model replacing another at setup, and the cap
turns the revolution. Decided by the author, 2026-10-03: the second.

**The library piece, the same for both.** Either way the bodies are
unknowns of time alone and the blade is a field of time and of the
blade, on one manifold: the cap's coordinates are read by the blade's
law (in the cap's axes) or by its root condition (in fixed axes), and
the cap's law reads an integral over the blade. Built 2026-10-03:
- An unknown declared with time alone among its arguments, on a
  manifold of time and space, `omega % unknown('theta', [t])`, is
  stored at every point; its equation is stated once per instant, at
  the first point, its terms of time once and its integral over the
  region at the quadrature points of the elements; at every other point
  the copy equals the one before, so the basis forms the unknown
  itself.
- A field may contain an integral over the region,
  `I*theta_tt + integral(g, over=omega % space()) - tau`, entering
  linearly, differentiated along time inside the integral.
- Newton ends converged where no halving of its step decreases the
  residual and the step is below half the digits of the iterate: an
  iterate started from the instant before, whose residual a relative
  tolerance would place below the rounding of the residual's own
  evaluation.
- Verified on `euler_bernoulli_beam.f90`, case 8: the clamped beam of
  case 6 on a turning hub, the hub's angle an unknown of time alone
  whose law reads the integral of the beam's moment of momentum. At
  801 instants the error of the deflection falls as 1.51e-4, 9.47e-6,
  5.40e-7 over 9, 17 and 33 points, order 4, and of the angle as
  4.18e-5, 2.66e-6, 5.10e-7; with the angle's law reading no integral
  the deflection is case 6's to every printed digit. dJ/ds through the
  hub's inertia by the adjoint and by the expansion agree to 5e-13,
  central differences to their own error, on 1, 2 and 4 images.
- The stationarity of a density along several coordinates,
  `density % stationarity(w, [t, x])`, keeps each derivative an
  operator, so that the elements read the outer derivatives along space
  as fluxes: case 8's laws formed from its action give its errors.
- At the first instant an equation of time alone that reads no
  component of top degree, a joint's condition, is stated by its
  derivative along time under the projection too; a march begun from
  the instant before is solved to the rounding of its residual, and
  its first instant begins from the estimate's derivatives along time.
  Before, it began from zero rates, which laws quadratic in the rates
  reached by a long damped iteration: the rigid rotor's first instant
  now takes two Newton steps, its revolution the same digits.

**The blade as a beam in the cap's axes** (`multibody_dynamics.f90`,
`beam`), the frame the author chose. The blade is a field of time and of its span: the
displacement along its span, chord and thickness and its twist, the
section's centre at root + (x + u) e_s + v e_c + w e_f in the cap's
reference axes. Its density over the span is the kinetic energy of the
section's centre and of its turning (the cap's angular velocity and the
twist's rate, with the section's moments of inertia), less the strain
energy with the axial strain to second order in the slopes, which
stiffens the turning blade, and less the potential in gravity. The
bodies' coordinates and the reactions become unknowns of time alone;
the cap's laws gain the integral over the span of the stationarity of
the blade's density. `blade_on_hub.f90` verifies it on one blade on a
driven hub: as E grows the driver's reactions approach the rigid
plate's (-225.9008 and -225.91163 at E = 7e15 and 7e17 against
-225.91174), and the one about z falls with the step for both, the
dissipation of bdf(1); the same digits on 1, 2 and 4 images.

**The rotor with beam blades** (`blade = 'beam'` in `rotorcraft.nml`,
the span's points `span_points`). Every point of the span stores
every field, the 355 unknowns of time alone among them, so the
solver eliminates the rows that tie their copies - with the history,
the family's relations and the zero coefficients - and solves the
rest directly: about 355 unknowns of time alone and 32 per point of
the span at each moment. The element matrices are formed field by
field, an expression's degrees are derived in one pass over its
leaves, and the rules are evaluated on the real values the basis
forms. The setup takes 5.5 s and a step about 2 s at 5 points of the
span; a revolution of 360 steps takes 22 minutes, the tip stretched
by 7.5e-4 by the turning, the hub's angle to 1e-16. At the tenth step
the tip of the blade at zero degrees has moved along the chord by
-7.10262e-4, -7.07935e-4, -7.07820e-4 and -7.07813e-4 at 3, 5, 9 and
17 points of the span, order four; the same digits on 1, 2 and 4
images. A march of 36 steps a revolution does not resolve the blades'
motion: its tip motion is ten times that of 360 steps. Drawn, each
point of the box of a blade moves with its section of the beam, and
`rotor_<instant>.vtu` states the beam's deflection as a vector at
every point, for paraview's warp by vector.

A law of time alone that reads no integral over the span - of a body
the blades are not fixed to, or of a joint - is tested at one site of
each moment with one, so that it is evaluated there alone; every
rule is evaluated only where its test functions are other than zero,
an expression forms the place of each leaf in its tuple once with its
domain, and the reverse sweep's sources are formed in the components
a rule reads alone. At 36 steps and 5 points of the span the march
takes 2 min 20 s where it took 5 min, the motion the same to every
printed digit and the smooth maximum to 7e-15. Still to come: the
blade's energy over the span in the printed energy.

Each block of the march begins from the instant before with the
unknowns of time alone - the bodies' coordinates - extrapolated along
time to the block's last instant by the taylor polynomial of their
derivatives below the top one; the hub turns by 10 degrees in a step
of 36, and Newton began every block that far from its solution. The
top derivative is not read: at the end of a step of a tableau the law
alone determines it, 2e2 where the acceleration is 8 on the beam of
case 8. Nor is the beam's field extrapolated: its stiff modes' rates
times the step begin Newton at six times the residual, its test
relative to that residual then ends at the edge of the rounding, and
case 8's misfit varied by 5e-15 between designs 1e-9 apart where it
varies by 8e-16. The march of 36 steps takes 1 min 7 s where it took
2 min 17 s, of 12 steps 2 min 2 s where it took 8 min 24 s, the motion
the same to every printed digit and the smooth maximum to 1.2e-14 and
4.4e-14; the rigid rotor's revolution 2 min 9 s where it took 2 min
38 s, the same energies; `blade_on_hub` as a beam 7.6 s where it took
15.7 s. The direct solve of a block with unknowns of time alone scales
its rows and then its columns to a largest entry of one before it
factors them: along the eight strips the derivatives of the smooth
maximum by the reverse sweep and by the expansion agree to 4e-14 to
2.5e-13 where they agreed to 1.6e-11 to 1.1e-9, and `blade_on_hub` at
E = 7e15, whose factorisation found a singular pivot, runs again,
-225.9008 and 3531.5 as its header records.

**The span's points and the cost of a march.** The 48 strips of the
case need a point of the span at the edge of each, 49 points, where
the rest of this note takes 5. At 36 steps, 48 strips at 2 cm, the
march takes:

| points of the span | time | memory | smooth maximum |
|---|---|---|---|
| 13 | 8 min 55 s | 0.66 GB | 1.877 |
| 25 | 19 min 59 s | 0.95 GB | 1.920 |
| 49 | 36 min 30 s | 1.41 GB | 1.935 |

against 1.676 at 5 points: the designs found on 5 points rest on a
stress some 15 per cent below the one the span resolves. The time
grows as the points, not as the cube of the block's 355 + 32 np
unknowns: at 25 points, of 18 min 53 s, the residual's 1394
evaluations take 607 s, the 372 jacobians 424 s, and the eliminated
solve, its dense factorisation within it, 56 s. With the thickness
uniform the same march takes 12 min 52 s and with eight strips 13 min
49 s, so that the 47 unit steps of the 48 strips cost a third of it.
Newton takes about ten iterations a block, and at 5 points 5.9, of
which some two stand at the rounding of the residual, each halving its
step up to eight times before the iteration ends. An exit at the
rounding by the course of the iteration - a full step that does not
decrease the residual after one of full length that cut it a
hundredfold, below 1e-6 of the first - was tried and is not taken: at
25 points it takes a quarter of the time, 14 min 11 s, 298 jacobians
and 945 residuals, and the cases of the suite are the same to every
digit, but case 8's misfit then varies by 3.6e-15 between designs
1e-9 apart where it varies by 8e-16, the iterations it ends polishing
the hub's angle. A revolution of 360 steps at 49 points and its
reverse sweep would take, by these marches of 36 steps, hours for each
scenario, and an optimization of the three days. At 49 points the
eliminated solve met 1e-6 to 2e-8 of the residual it began at, so that
Newton's steps were ruled by its error; it refines its solution now:
of the 227 solves of a march of 10 steps, 198 meet their tolerance,
most to 1e-14, and 29 end between 1.5e-8 and 1.6e-7, against 1.5e-8 -
16 in a block that converges in 29 iterations, 12 in a step of 36
degrees that does not converge.

**The thickness, the mass and the stress** (`strips`, `thickness`,
`sensitivities` in `rotorcraft.nml`). With strips, the thickness of
each strip of the span is a design coordinate, the same for the four
blades: h = h_1 plus, at the inner edge x_j of each strip after the
first, (h_j - h_(j-1)) H(x - x_j), H the unit step, a function of the
coordinates alone whose derivative along a coordinate it reads is
refused. The elements read the thickness inside the fluxes alone, so
nothing differentiates it along the span. The functions:

- the blades' mass over the revolution, the integral of 4 rho c h:
  2.1737794328871449 kg s at 2 cm, TACS's 2.173779432887122;
- the smooth maximum (`operation_field`, `maximum`) of the von Mises
  function (sigma^2 + 3 tau^2) / sigma_y^2 at the corners of each
  section and the middles of its long sides, over the four blades and
  the revolution: log( integral of sum exp(rho f) ) / rho, evaluated
  as m + log( integral of sum exp(rho (f - m)) ) / rho, m the largest
  f at the quadrature nodes. The square of the ratio, which has no
  root to differentiate where the stress vanishes, at the first
  instant.

With eight strips at 2 cm the motion is the same to every printed
digit as with the uniform thickness. At 36 steps, 5 points of the
span and rho = 100, the derivative of the smooth maximum along the
thickness of the first strip by the expansion is -166.45582952328755;
central differences of steps 1e-5 and 2e-5 miss it by 1.318e-4 and
5.271e-4, ratio 4.00, their own error, and their extrapolation is
-166.4558295191, within 2.5e-11. Along the root strip the stress falls
and along the outer strips it rises, the centrifugal load they add.

**The stress of the beam blades and the march.** The largest von Mises
ratio of the blades at 2 cm, collective, 5 points of the span:

| march | largest ratio | its parts: axial, lag, flap |
|---|---|---|
| 36 steps, bdf(1) | 1.29 | 0.16, 0.85, 0.28 |
| 72 steps, bdf(1) | 0.83 | 0.16, 0.51, 0.16 |
| 36 steps, bdf(2) | 0.52 | 0.16, 0.34, 0.02 |
| 72 steps, bdf(2) | 0.29 | 0.22, 0.06, 0.01, at the start |
| 360 steps, bdf(1) | 0.33 | 0.16, 0.14, 0.03 |

The lag bending is the error of the march: bdf(1) forms the
centripetal acceleration of the turning blade a step late, a load
along the chord of about w dt times the centrifugal load, which bends
the blade at the root; a rigid blade takes it as a reaction of its
link. The axial stress, 0.16 of the yield, is that of the turning,
and neither it nor the lag bending depends on a uniform thickness: at
36 steps and 9 points of the span the ratio is 1.34 at 1 cm and 1.36
at 2 cm. Only the
thickness's distribution changes them, a blade thinner outboard
loading its root less. KS <= 1/4 is not met by a uniform blade at
any of these marches; TACS's run ends where it is not met either, KS
0.37 at about 1 cm.

The smooth maximum here is the logarithm of an integral over the
revolution and the span, of the square of the ratio, at `ks_weight` =
100: it differs from the largest value by the logarithm of the measure
near it, summed over the blades and the points of the sections, over
the weight. At 360 steps it is 0.0538 for the uniform
blade at 2 cm, where the square of the largest ratio is 0.109, and
0.0185 for the thickness falling from 2 cm at the root to 1 cm at the
tip in eight strips; at 36 steps, where the lag bending of bdf(1)
lasts the revolution, 1.676 against 1.66. TACS bounds 4 KS - 1 with
rho_KS = 1e5 on the ratio itself, the largest value to 1e-4: the case
as TACS states it takes `ks_weight` near 1e5.

## The case

- **Bodies.** The case has 26 bodies:
  - 22 rigid: the hub; the lower and upper swashplates and the sphere;
    four lower and four upper pitch links; four blade caps; three push
    rods; the base plate; two push horns.
  - 4 flexible blades: flat plates 1.56 m by 0.121 m at z = 0.9, each
    192 geometrically nonlinear MITC9 shell elements (48 spanwise by 4
    chordwise). The material is isotropic first-order shear (E = 70
    GPa, rho = 2500, yield 350 MPa).
- **Joints.** 34 joints are stated as 66 multiplier elements, with no
  boundary condition anywhere:
  - two drivers: the hub turns at 109.12 rad/s, and each push rod is
    driven along z by a sin(27.28 t);
  - revolute, spherical and cylindrical joints;
  - rigid links from each cap to the 9 root nodes of its blade;
  - one fixed joint at the base plate.
- **Unknowns.** At each node: the displacement, the Euler parameters,
  and the multiplier of their unit norm. In all, 3580 nodes give 28640
  unknowns.
- **Time.** One revolution, 360 steps of BDF1, with Newton's iteration
  on a direct factorization. The run takes two or three Newton steps
  per step.
- **Design.** 48 thicknesses, one per spanwise strip, shared by the
  four blades, between 1 and 2 cm.
- **Objective.** The blades' mass integrated over the revolution: 2.1738
  kg s at 2 cm.
- **Constraints.** Per scenario, 4 KS - 1 <= 0: a KS of the von Mises
  ratio over the blades and the instants, rho_KS = 1e5, so that the
  stress stays below a quarter of the yield stress. Also 94 linear
  smoothness constraints between adjacent strips.
- **Scenarios.** Three: collective, longitudinal cyclic and lateral
  cyclic. They differ only in the push-rod amplitudes. The objective is
  the collective's.
- **Optimizer.** ALGENCAN, an augmented Lagrangian method, on gradients
  by the discrete adjoint.

## The case in this framework

| TACS | Here |
|---|---|
| an element's residual R(t, q, q', q'') | a law on its manifold |
| a kinematic constraint and its multiplier element | a condition paired with its multiplier, the joint's reaction |
| a driver | a condition stating the motion |
| a function, sum over steps of w_k f_k | an objective, an integral over the manifold |
| the BDF integrator | a chain of families |
| the adjoint and the total derivative | the reverse sweep and the design's multipliers |
| the assembler | the residual operator |
| three scenarios on one design | three Lagrangians on one design |

The thesis derives the equations of motion from the stationarity of
the action (section 4.1). Here they come from that action itself: the
kinetic and potential energies of the bodies, and each joint's
condition paired with its reaction. `stationarity` forms the
Euler-Lagrange equation in each coordinate, as it does for any
Lagrangian density. The joints' conditions hold at every instant, so
they belong to the law on the manifold of time, and their reactions
are unknowns of that manifold.

One scenario's Lagrangian:

```
L = J(s) + <adjoint, r> + <lambda, ic> + <kappa, s - s_0>
```

- J is the blades' mass integrated over the revolution.
- r collects the Euler-Lagrange equations of the action and the
  joints' conditions; ic the initial state.
- Removing the design condition makes the design free, as in the
  beam's identification (`doc/one-lagrangian-four-problems.md`).
- The constraints KS <= 1/4 are inequalities on the design factor,
  each with a multiplier that is zero where it is not active.
- The three scenarios are three such Lagrangians on one design, the
  non-intrusive form of optimization under uncertainty with three nodes
  of weight one.

## Stages

Each stage ends with the same digits serially and on several images,
and with dJ/ds by the adjoint against the tangent and central
differences.

1. **The rigid mechanism** (`collective_rigid.py`: the same model with
   the blades rigid). It needs bodies as unknowns of time alone, the
   action of each body, and the joints and the two drivers as
   conditions with their reactions. The result is an index-3
   differential-algebraic system marched by BDF1, as TACS marches it.
   Verified against TACS: the initial kinetic energy, the hub's driving
   torque, and the joints' reactions over the revolution.
2. **The blades as rotating beams** along the span by the Hermite
   elements: centrifugal stiffening, the root fixed to its cap's motion,
   and the 48 strips' thicknesses as design coordinates in EI and rhoA.
   The beam is a model of the shell, so it is verified by its own
   convergence and compared with TACS's tip motion and stresses, not to
   the digit.
3. **The functions.** The mass and the KS of the bending stress as
   integrals over the manifold; dJ/ds and dKS/ds by the adjoint,
   verified against central differences.
4. **Inequality constraints in `minimize`.** Each constraint KS <= 1/4
   and each smoothness constraint is a condition with a nonnegative
   multiplier, solved by an augmented Lagrangian as ALGENCAN does, or
   by an active-set Newton on the reduced space. Verified first on the
   beam's identification with a bound on the stiffness.
   Built 2026-10-03 for inequalities on the design factor,
   `continuous_residual(omega % design(), [c], inequality=.true.)`
   paired with its multiplier: each Newton step on the design minimizes
   the quadratic model of J under the inequalities linearized at the
   iterate, its active set settled within the step, and the multipliers
   returned are those of the design returned. The beam's identification
   with s1 <= 0.45 finds (0.45, -0.22396067) in six iterations, the
   multiplier 2.4205141386e-5, -dJ/ds1 at that design by the adjoint to
   every digit; a lower bound, two bounds of which one is active, and
   s1^2 + s2^2 <= 1/4 converge in six iterations each, the multipliers
   satisfying the stationarity of L in the design. The smoothness
   constraints are of this kind. KS <= 1/4 is not: it is a functional of
   the state, and needs its own reverse sweep.
   Built 2026-10-03 for functionals of the state: an inequality on an
   integral or a smooth maximum over the manifold has its value and its
   gradient from the expansion of the state, and enters the reverse
   sweep times its multiplier of the step before, so that the model's
   hessian is the Lagrangian's. Each iterate is tested by the merit J +
   penalty * violations; a step the merit refuses is taken again from
   the design before it with the hessian shifted and each violated
   inequality asked for a quarter of the reduction, and a hessian that
   is not positive definite - the mass's, zero - is shifted from the
   start. The steps of the bounds above are unchanged. The beam's least
   mean stiffness with its misfit at most m (`misfit=m`): 13
   iterations at m = 1e-4 and 32 at 1e-6, the multiplier the same along
   s1 and s2 to 3e-12 and 6e-11.
   The hessian by the expansion costs a reverse sweep's expansion along
   every design coordinate: on the rotor at 36 steps and eight strips
   251 s forward and 1998 s backward an iteration, growing with the
   strips. With the design discretized by the expansion of order zero
   the gradients are those of the reverse sweeps alone, J's and each
   functional's of its own on the same forward sweep, and the hessian
   their damped BFGS update, its cost apart from the number of strips:
   on the same rotor 69 s forward and 20 s backward. The beam's bound
   and misfit cases converge to the same designs (`bfgs`, a fourth
   argument), in 19, 15 and 45 iterations where the expansion took 6,
   13 and 32. The rotor optimizes by this update.
   Each inequality has its own penalty in the merit, at least twice its
   own multiplier. One penalty, twice the largest multiplier, weighed
   the rotor's smooth maximum of the stress - of no units and of
   multiplier 0.04 - by the 34.6 of the bounds on the thickness, so
   that a step that lowered the mass by 3.4e-3 and exceeded the
   stress's bound by 3.4e-3 was refused, and taken again unchanged, its
   length fixed by its active inequalities; a refused step is now taken
   again at most half as long. The beam's cases, of one inequality
   each, are the same to every digit. The rotor's least mass at 36
   steps of bdf(1), 5 points of the span and eight strips from 2 cm,
   the bound of the von Mises ratio 1.25 - the march's own stress is
   1.29 at 2 cm -, takes six iterations where it took 21, eleven
   minutes: the root strip 1.0203272556788471 cm, the seven others at
   their bound of 1 cm, the mass 1.0896514020886789 kg s, the smooth
   maximum on its bound to 5e-15, and its multiplier 0.037891294249,
   the 21 iterations' to 1.6e-12.
   Each step is found by the dual method of Goldfarb and Idnani: from
   the minimizer of the model without inequalities, the most violated
   inequality joins the active set, the normals of the active ones
   remaining independent - one that depends on them moves the multipliers
   alone, and one leaves where its multiplier vanishes - and the step
   and the multipliers are those of one solve of the active set's
   equalities at every move. The primal active set it replaces took
   dependent normals in: the rotor's three scenarios at 2 cm, its
   eight upper bounds active and three smooth maxima violated, made its
   system singular. Where no step meets the inequalities linearized,
   each violated one is required to remove a quarter of its fraction,
   three times, and then only not to grow, which the null step meets. Of 2000
   random problems, every step meets its conditions to 1.3e-11; each of
   the 230 without one has the weights of inequalities whose
   normals cancel and whose constants sum above zero; and where the
   primal active set settled, 1738, the step and the multipliers are
   its own to every digit. The beam's cases and the rotor's are the
   same to every digit.
5. **Three scenarios on one design.** The sum of three discrete
   Lagrangians, each with its own state and adjoint, and one design
   iteration.
   Built 2026-10-03: `minimize_together(residuals, estimates,
   solutions)`, residuals on manifolds of the same design coordinates,
   each with its forward and reverse sweeps at every iterate; the
   gradient and the hessian the sums of theirs, the inequalities those
   of all, the merit the sum of the objectives, any residual having
   none. The beam's identification as two scenarios at two frequencies
   converges in six iterations, the derivatives of the two misfits
   along the design summing to 5e-11; the least mean stiffness of one
   scenario with the misfit of the other at most 1e-4 in thirteen, the
   multiplier the same along s1 and s2 to 2e-12.
   The rotor's scenarios (`scenario` lists several in
   `rotorcraft.nml`), each stated on a manifold of its own and marched
   in turn; to optimize, one design of all, the first stating the mass
   and the inequalities on the design alone and every one its bound on
   the stress. At 36 steps, 5 points of the span and eight strips from
   2 cm, their smooth maxima are 1.676, 2.106 and 2.991 at 2 cm and
   1.638, 1.811 and 2.323 at 1 cm, so that a bound of 1.8 leaves every
   strip at 1 cm, in two iterations. With ks_bound = 1.5 the least mass
   takes six iterations, 33 minutes, none refused: the root strip
   1.0308636977683841 cm and the seven others at 1 cm, the mass
   1.0910828959080567 kg s, the lateral's smooth maximum on its bound
   to 7e-15 and its multiplier 0.0578856747211, the collective's and
   the longitudinal's 1.525 and 1.739 and zero. The lateral's derivative
   along the root strip by central differences of step 1e-6,
   -234.706109775, gives that multiplier as the mass's over it to
   5e-9.
6. **The blades as plates** (later): first-order shear plates on the
   span-by-chord rectangle, by the Lagrange elements of two
   coordinates, nearer TACS's shells.

## What exists, what is missing

- **Exists:** the manifold of time with design coordinates; laws in
  sin, cos, exp and products; `stationarity`; conditions with their
  multipliers; objectives as integrals; the reverse sweep and the
  expansion along the design; a free design found by Newton; the finite
  element method on lines and meshes; runs on several images;
  inequalities on the design factor; unknowns of time alone on a
  manifold of time and space, their laws reading integrals over the
  region; inequalities on functionals of the state, each with its own
  reverse sweep; several Lagrangians on one design, the rotor's three
  scenarios among them; the smooth maximum with its largest value
  subtracted; the blades as beams, their section given by its
  properties per length, as TACS's TimoshenkoStiffness takes them, or
  implied by a uniform solid rectangle, which alone states its stress;
  ten objectives of the blades as beams - the mass, the smooth maximum
  of the stress, the strain energy, the mean squares of the flap, lag
  and twist and of the pitch-link, push-rod and hub loads - with their
  derivatives along every strip by the expansion and one's by the
  reverse sweep, which agree to 1e-13.
- **Missing:**
  - a way to write a body's 7 coordinates and its action once, not as
    7 scalar laws by hand, without moving solver logic into `physics/`;
  - an initial state consistent with the joints (TACS starts the push
    rods at rest against their drivers); the families on the index-3
    system are measured in "A step of a DIRK tableau on the joints";
  - the case as TACS states it: 48 strips, three scenarios, 360 steps,
    KS <= 1/4 with `ks_weight` near 1e5;
  - the blades as plates and as solids;
  - the fan plot of each scenario, with TACS's uniform HART-II section.
    TACS's own (examples/rotating_beam) holds no input: one blade on a
    driven hub at 2.5 degrees of precone, the frequencies of its laws
    linearised about the state one revolution reaches, at 0.1 to 1.2
    times 109.12 rad/s. A collective held is a steady state, its
    frequencies the eigenvalues of the laws linearised about it; a
    cyclic held, a periodic one, its frequencies Floquet's, of the map
    over a revolution, each up to a multiple of the rate. Both need
    the state with the inputs held, and the eigenvalues.

## Objectives, constraints and their reductions

Every objective and constraint of the rotor is a reduction of a field
over the manifold - the revolutions and the span - and over the field's
components, the blades', the links', the rods' or the hub's
(`operation_field`): the mass an integral, the stress a smooth maximum,
the others means. The reductions are one operation of several
generators, phi^-1 of the integral of phi(f): the integral's phi the
identity, the mean's the identity over the measure made one, the smooth
maximum's the exponential, of negative weight the smooth minimum. The
smooth maximum is so not an integral of an integral but one integral
in the semiring whose sum is log(exp(rho a) + exp(rho b)) / rho, the
maximum as rho grows and, of the mean's normalized measure, the mean as
rho falls to zero; an integral of an integral is a reduction over a
product of factors, which of one generator is the reduction over each in
turn (the smooth maximum over the revolutions of the smooth maximum over
the span is the smooth maximum over both, exactly), and of two
generators a composition whose order matters, not yet stated. The
derivative of each reduction is a measure over the points: uniform for
the mean, the weights exp(rho f) / integral for the smooth maximum. A
constraint is the same reduction of the field less its limit, and an
optimization demonstration differs from another only in its
`objective`, `constraints` and `limits`.

## The published optimization

The dissertation (section 5.4.3, eq. 5.1) optimizes this rotor, the
TACS case's mechanism: the least mass of the blades over the 48 panels'
thicknesses, 10 to 20 mm, adjacent panels within 1 mm, the smooth
maximum of the von Mises ratio (TACS's KSFailure, the ratio itself, of
weight 1e5) at most 0.25 in each of the three scenarios, from 2 cm, by
SLSQP in 73 iterations; all three stress constraints active, the mass
some 0.56 of the initial, the thickness some 17 mm at the root falling
by about 1 mm a panel to 10 mm at 40 per cent of the span and on to the
tip (Figs. 5.15 and 5.16, read off).

Stated here - `scenario = 'collective', 'longitudinal', 'lateral',
steps = 72, order = 2, strips = 8, thickness = 8*0.02, optimize =
.true., ks_bound = 0.25, ks_weight = 1.0e4, bounds = 0.01, 0.02,
smoothness = 6.0e-3` - the blades are beams of eight strips on five
points of the span, adjacent strips within 6 mm, the dissertation's 1
mm over the six panels of a strip, and the smooth maximum is of the
squared ratio of weight 1e4, about 5e3 on the ratio itself. In 68 min
41 s it finds 11.79 mm in the root strip and the lower bound outboard,
the mass 0.511 of the initial; the longitudinal bound active, mu_ks
1.31, the collective's 0.12 per cent and the lateral's 1.8 per cent
below theirs. This is not the published optimum, and what decides it
is not the blade's stress:

- **The largest stress is the start's.** The blades begin unstretched
  at the full rate, and the axial wave this sets off, near 850 Hz and
  49 periods a revolution, gives the largest ratio at the second or
  third instant: at the design found 0.250, 0.207, 0.227 and 0.264 at
  72, 144, 288 and 576 steps, not converging. The wave lasts as long
  as the march damps it: at 576 steps the ratio is 0.21 to 0.25 up to
  60 degrees and 0.17 to 0.18 up to 120, where at 288 steps it is 0.16.
  The blade itself is undamped. Past 150 degrees the ratio is 0.152 to
  0.159 at 288 and 576 steps, near the centrifugal stress at the root
  of the second strip, 0.153. A uniform blade's centrifugal stress at
  the root is rho Omega^2 (R^2 - r0^2) / (2 sigma_y), 0.162 of the
  yield at any thickness: the 0.293 to 0.299 measured at 72 steps,
  read before as that stress, is the start's. The bound 0.25 and the
  limits set from those measurements bound the march's error.
- **The reactions' mean squares include the start's impulse.** The
  push rods' doubles as the step halves, 6.7e8, 1.3e9, 2.7e9 and 5.3e9
  N^2: it grows as 1/dt. The flap's and the lag's converge, the
  flap's at the second order: 9.2e-7, 3.5e-7, 2.7e-7 and 2.5e-7 m^2.
- **Five points of the span resolve the stress roughly** where the
  thickness changes from strip to strip: one cubic element spans two
  strips, its strain continuous. At 288 steps, past the start, the
  ratio is 0.152 to 0.159 at five points and 0.160 to 0.166 at nine,
  one element a strip; the start's peak is 0.227 and 0.234. At 72
  steps the design found exceeds its bound at nine points in every
  scenario by 12 per cent of the squared ratio (0.0703, 0.0704,
  0.0690): the start's peak, which the span's resolution changes too.
- **TACS's shells state a stress at the root that the beam does not.**
  Their smooth maximum at 2 cm over the revolution, their own start
  included, is 0.438 to 0.441 at 360 steps ("Targets from the TACS
  run"); the beam's over the revolution at the design found is 0.227
  and 0.264 at 288 and 576 steps, and 0.16 once the start has passed.
  Their optimum keeps 17 mm at the root.

At the design found, rounded to nine digits, the reverse sweep and the
expansion give the smooth maximum's derivatives in the longitudinal
scenario to 2.2e-14; the root strip's stationarity, dm/dh_1 + mu_ks
dKS/dh_1 = 0, gives mu_ks to 1.0e-8 of the one returned, and the seven
lower bounds' multipliers, dm/dh_k + mu_ks dKS/dh_k, are 14.1 to 19.0,
all positive: the discrete problem is solved.
The comparison with the dissertation waits on initial data consistent
with the rotation ("Decisions for the author", 5) and a march
converged in the step. The least mean square of the flap under the
mass and the stress (`objective = 'flap', constraints = 'mass',
'stress', limits = 1.6303345746653587, 0.0849`, 72 steps), whose
report first stopped the program (23b8294), ends with both constraints
inactive, its limits also set from the march's error; at 72 steps its
objective is mostly that error, the flap at the published design 2.6
times its value at 144 steps. Its design is not reported.

## Other published optimizations

Chattopadhyay and Walsh (NASA TM-100569, 1988; the journal paper,
Journal of the American Helicopter Society 34(4), 1989) minimize the
weight of the blades of a four-bladed articulated rotor in vacuum, 22
ft in radius at 293 rpm, its blade a box beam inside the airfoil in ten
segments. The design is each segment's wall thicknesses and
nonstructural weight and, for the tapered blade, the box's height at
the root and its taper ratio. The constraints are the first three
lead-lag and the first two flap frequencies within 1 per cent of the
reference blade's; the autorotational inertia, the sum of each
segment's weight times the square of its radius, at least the
reference's; and in each segment the centrifugal stress, the weight
outboard times Omega^2 r over its area, at most the allowable over a
factor of safety of 3. The reference blade weighs 98.27 lb, the optimum
95.62 lb with 30 design variables, 93.61 with 40 and 92.16 tapered with
42, the first lead-lag frequency at its upper bound and the inertia at
its lower in every case. The paper tabulates the reference blade
segment by segment and the optima (its Tables 1 to 3 and 5 to 7): a
validation for the blades as beams. Its mass and its stress are
reductions the rotor states already, and the inertia is another, the
integral of the mass per length times r^2. Missing are the frequencies,
the eigenvalues of the laws linearised about the steady rotation; a
bound from below, the inertia at least its reference, where a limit
bounds from above; the box beam's section; the hinges at 0.833 ft; and
the nonstructural masses.

Peters, Ko, Korn and Rossow place the frequencies of a blade by its
mass and stiffness, also in vacuum (39th Forum of the American
Helicopter Society, 1983; NASA grant NAG-1-250, 1982 to 1985), which
also needs the eigenvalues. The optimizations of Ganguli and Chopra, of
Friedmann and of the HART II blades reduce the vibratory loads of the
hub in forward flight and need the air.

## HART II

The blades of HART II are the TACS case's in radius (2 m), chord (0.121
m) and cutout (0.44 m), and add a precone of 2.5 degrees, a linear twist
of -8 degrees over the radius zero at 1.5 m, and the NACA 23012 with a
trailing-edge tab (the geometry's `--preset hart2`). The rotor takes
them at setup: `precone` cones each blade about its chord at the root,
its beam's axes, its rigid plate's centre and inertia and its drawing;
`twist` and `twist_zero` pitch the beam's section along the span, its
stiffnesses taking the curvatures along the turned axes and its inertia
and stress points turning with them, so that a constant pitch of 90
degrees marches as the section with its flap and lag properties
exchanged, the tip's displacements to every printed digit and the
driver's reactions to 7.7e-12; `section_kind = 'given'` takes TACS's
uniform HART-II section; `blade_meshes` draws the blades with their
meshes. With precone and twist zero every output is the same to every
byte.

## A step of a DIRK tableau on the joints

At the end of a step of a DIRK tableau the quadrature of the stages
fixes the positions and the velocities, as the initial data fix them
at the first instant. A joint's condition there read nothing the step
left free, its reaction had no row, and every tableau stopped the
pendulum and the rotor at its first step with a singular pivot. The
block now states the laws at such an explicit moment by the rule of
the first instant, each condition by its second derivative along
time; projected on the elements, an instance's columns include those
that rule reads. The 92 cases and the heat cases are the same to
every digit, serially and on 1, 2 and 4 images, and the rotor's
marches, derivatives and least mass to every printed digit. Through
the new rows the pendulum's dJ/ds by the adjoint and by the expansion
agree to 2.2e-15 by dirk(3) and to 1.5e-13 by dirk(4,5), and the
rotor's along two strips to 5e-14 from a start by dirk(3). A march by
a tableau gives on 1 and 2 images the digits it gives serially, but
for the last two of the rotor's mass on 2 images, which the sum over
the images takes in another order, as it did by bdf(1) before.

The tableaux on conditions of position, a system of index three, the
rotor's blades as beams at 2 cm, 36 steps and 5 points of the span:

| march | the pendulum, 201 to 801 instants | the rotor's largest von Mises ratio |
|---|---|---|
| dirk(2) x12, bdf(2) x12, adams(2) x12 | diverges | 2.46, the energy swinging threefold |
| dirk(3) x12, bdf(2) x12, adams(2) x12 | energy of order 2.6, conditions off by O(h^2) at the instants | 0.462 |
| dirk(4,5) x12, bdf(2) x12, adams(2) x12 | energy of order one, conditions exact | 0.394 |
| bdf(1) x1, bdf(2) x23, adams(2) x12 | | 0.568 |
| bdf(2) | | 0.520 |
| bdf(1) | | 1.290 |

The implicit midpoint, dirk(2), and adams(2), the trapezoidal rule,
multiply a stiff error by -1 at every step and damp none: the
midpoint's step ends at twice its stage less its start, off the
conditions, and its stages' return to them grows; adams(2) leaves the
rigid rotor's energy oscillating by 1 per cent. The usual remedy is
the conditions stated at the velocity, a choice of the statement left
open. The SDIRK of Hairer and Wanner weighs two of its stages
negatively, -49/48 and -85/12, and the smooth maximum, the logarithm
of a quadrature over the stages, is then not defined: NaN.

## Decisions for the author

Decided 2026-10-03: each blade at its cap; the blade model chosen at
setup, rigid, beam, plate and solid in that order; the in-house
minimizer; the blade linear in its cap's axes; the optimization with
the beam, the plate and the solid, on a reduced case first (about 8
strips, one scenario, 36 steps) and then on the case as TACS states
it; TACS compared as it is, its rigid-link sign noted.

Open, with what the work has measured:

1. **The march of the reduced case.** At 36 steps of bdf(1) the stress
   of the blades is 1.29 times the yield at 2 cm, its lag bending the
   centripetal acceleration a step late, so that KS <= 1/4 is met by no
   uniform blade; 72 steps of bdf(2) give 0.29, 360 of bdf(1) 0.33. The
   reduced optimizations above take ks_bound = 1.25 and 1.5 in its
   place. Extrapolated in the step, bdf(1) from 72 and 360 steps at
   the first order and bdf(2) from 36 and 72 at the second, the ratio
   converges near 0.21, 0.205 and 0.213. The marches of three families
   are measured in "A step of a DIRK tableau on the joints"; a march
   chosen by its measured convergence is proposed.
2. **The plate theory**, before the plate and the solid: Kirchhoff
   plates by bicubic Hermite elements, which need the mixed
   derivatives along space in the elements; first-order shear plates
   by the Lagrange elements of degree two, Q9, which serve the solid
   too and need the cells' edges numbered as members of the mesh, not
   a mesh of the second order, proposed; or by Q1, which locks.
3. **The case as TACS states it**: its 48 strips need 49 points of the
   span, where a revolution of 360 steps and its reverse sweep would
   take, by the marches of 36 steps measured, hours for each scenario
   and the optimization of the three days; or
   a coarser case - 48 strips on 25 points, or fewer strips - which is
   not TACS's.
4. **The weight of the smooth maximum** of that case: near 1e5, so
   that it bounds the largest stress as TACS's rho_KS = 1e5 does, the
   weight of 100 here lying a factor of two below it at 360 steps.
5. **The initial data.** The rotor starts at the full rate with its
   blades unstretched, and its push rods at rest while their drivers
   move them, a sin(Omega t / 4) at the rate a Omega / 4, as TACS's
   case does. The start then decides the largest stress and the
   reactions' mean squares, and neither converges in the step ("The
   published optimization"): every optimization stated so far bounds
   the march's error. Proposed: initial data consistent with the
   rotation and the drivers - the blades in the steady rotation,
   stretched by their centrifugal load, and the rods moving at their
   drivers' velocity. The steady rotation is also the state about which
   the fan plot's frequencies are the eigenvalues, and its stress is the
   published minimum-weight blades' (TM-100569). The alternative, the reductions over a window of time
   after the start - a part of the time factor, as a face is a part of
   the boundary - converges only as fast as the march damps the start's
   wave, which the blade does not.

## Targets from the TACS run

At 2 cm thickness the initial kinetic energy is 4.99144e5 J, the
objective is 2.173779432887122 kg s, and KS is 0.437839, 0.440713 and
0.434083 for the collective, longitudinal and lateral scenarios. These
are targets for the shell model; the rigid stage reproduces the
kinematics and the beam stage approximates the rest.
