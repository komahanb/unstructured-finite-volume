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
  and by the expansion agree to 1.1e-14 and 1.1e-14, and central
  differences of step 1e-4 to 3e-10 and 7e-10, their own error; the
  same digits on 1, 2 and 4 images.
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
5. **Three scenarios on one design.** The sum of three discrete
   Lagrangians, each with its own state and adjoint, and one design
   iteration.
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
  region.
- **Missing:**
  - inequalities on functionals of the state, KS <= 1/4, each with its
    own reverse sweep (those on the design factor exist);
  - several Lagrangians on one design;
  - a way to write a body's 7 coordinates and its action once, not as
    7 scalar laws by hand, without moving solver logic into `physics/`;
  - a check that the families march an index-3 system, and an initial
    state consistent with the joints (TACS starts the push rods at rest
    against their drivers);
  - the KS aggregation with its maximum shift;
  - the rotating beam of the blades.

## Decisions for the author

Decided 2026-10-03: each blade at its cap; the blade model chosen at
setup, rigid, beam, plate and solid in that order; the in-house
minimizer; the blade linear in its cap's axes; the optimization with
the beam, the plate and the solid, on a reduced case first (about 8
strips, one scenario, 36 steps) and then on the case as TACS states
it; TACS compared as it is, its rigid-link sign noted.

1. **The blades at 90 and 270 degrees** are attached across the hub:
   `blade90.inp` reads `rotor2_192dv.bdf`, which lies at y from -2.0 to
   -0.44, but it is linked to the cap at (0, 0.429, 0.9). `blade270.inp`
   is the mirror case. Collective input is symmetric and unaffected;
   under cyclic input the two blades are 180 degrees out of phase with
   their names. Reproduce this, or place each blade at its cap?
2. **The blade model of the first optimization:** the rotating beam of
   stage 2, or the plate of stage 6.
3. **The optimizer:** an augmented Lagrangian in `minimize`, keeping
   the design iteration in the library.

## Targets from the TACS run

At 2 cm thickness the initial kinetic energy is 4.99144e5 J, the
objective is 2.173779432887122 kg s, and KS is 0.437839, 0.440713 and
0.434083 for the collective, longitudinal and lateral scenarios. These
are targets for the shell model; the rigid stage reproduces the
kinematics and the beam stage approximates the rest.
