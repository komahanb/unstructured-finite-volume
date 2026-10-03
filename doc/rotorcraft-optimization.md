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
  in TACS. TACS's runs of the case carry the error in the caps' links,
  to the pitch links and to the blades.
- The sensitivities (`triple_pendulum_sensitivity.f90`): the pendulum
  with a design in a mass and one in a joint's condition, which the
  first instant states by its second derivative. dJ/ds by the adjoint
  and by the expansion agree to 8e-14 and 9e-15, and central
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
carries the revolution.

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
  element method on lines and meshes; runs on several images.
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
minimizer.

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
