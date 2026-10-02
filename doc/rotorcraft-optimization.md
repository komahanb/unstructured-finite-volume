# The rotorcraft optimization in the Lagrangian framework

Status: plan, 2026-10-02. The source is komahanb/rotorcraft-optimization,
`deterministic-optimization`: Python scripts that build the case in TACS
(`tacs_builder.py`, `scripts/multiscenario.py` and its three scenarios).
The facts below are read from those scripts, the TACS sources the clone
vendors, and its run log (`logs/multiscenario.log`).

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
  - inequality conditions on the design, and an optimizer for them;
  - several Lagrangians on one design;
  - a way to write a body's 7 coordinates and its action once, not as
    7 scalar laws by hand, without moving solver logic into `physics/`;
  - a check that the families march an index-3 system, and an initial
    state consistent with the joints (TACS starts the push rods at rest
    against their drivers);
  - the KS aggregation with its maximum shift;
  - the rotating beam of the blades.

## Decisions for the author

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
