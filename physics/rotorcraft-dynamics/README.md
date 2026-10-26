# The dynamics of a rotorcraft

Consider choosing the material thickness along a rotor blade. Reducing
its mass also changes its resistance to deformation and its capacity
to carry load. We need its motion and the stresses it experiences
whilst in motion before we can choose thicknesses of least mass under
a stress bound. The program computes these quantities, their rates of
change with the strip thicknesses, and a design for one scenario or
for three together.

We use the rotor of the TACS case of `komahanb/rotorcraft-optimization`
(dissertation, section 5.4). With rigid blades it has 26 bodies and
35 joints: TACS's 34 and the rotating scissors, which turns the upper
swashplate with the hub. The hub turns at 109.12 rad/s, and three push
rods drive the swashplates. The bodies begin at the velocities the
drivers require. At setup we choose the model of each of the four
blades: a rigid body with the plate's mass and inertia, or a beam
fixed at its root to its cap. The choices may differ between blades.

With any beam present, its motion varies with time and position along
the span, so the manifold is time x span, with design coordinates
added for the strip thicknesses. The rigid bodies, including any rigid
blade, retain unknowns of time alone on that manifold. `rigid_link`
joins each rigid blade to its cap through `cap_blade_<deg>`, as in the
all-rigid rotor. Only the beams contribute span fields and root
conditions. The reported `mass` integrates the beams' mass alone over
the revolutions, in kg s; a rigid blade's constant mass is not counted.

| program | what it does |
|---|---|
| `rotorcraft_dynamics` | the rotor; its setup read from `rotorcraft.nml` |
| `triple_pendulum` | the triple pendulum of the dissertation, section 5.1 |
| `triple_pendulum_sensitivity` | the pendulum's dJ/ds by the adjoint and by the expansion |
| `blade_on_hub` | one blade on a driven hub, as a beam or as a rigid plate |

`rotor_bodies.txt` is the table of the rigid bodies, read by
`rotorcraft_dynamics`; the header of each program states what it
computes and what it was measured to give.

## Building and running

From this directory,

    make                              # the library, if its sources changed, and the four programs
    ./rotorcraft_dynamics             # reads rotorcraft.nml in the directory it runs in
    ./triple_pendulum 201             # the instant count (default 201)
    ./triple_pendulum_sensitivity 1 1 1   # the two design values and the order of the expansion
    ./blade_on_hub beam 36 5 70e9     # beam [steps] [span points] [E], or: rigid [steps]

`./run.sh` builds the library from the repository's root and runs the
four programs with their defaults. On several images, build with
`make COARRAY=lib` and run with `mpirun -np 2 ./rotorcraft_dynamics`;
one image needs `OMPI_MCA_osc=pt2pt mpirun -np 1 ...`. The rotor's march
gives the same digits on 1 and 2 images, the pendulum's on 1, 2 and 4.

To run a setup of your own, copy `rotor_bodies.txt` and a
`rotorcraft.nml` into a directory and run the program there: each run
writes its drawing into the directory it runs in.

Each program prints on image 1 the convergence of every solve, as
`van_der_pol` does: for each block of the march its Newton steps - the
residual, the step and its factor, the jacobian, the linear solve -
and where there is an adjoint, its blocks in reverse; to optimize,
each design iteration's objective, conditions and stationarity. The
rotor names each scenario before its march, and the scenarios before
the design's iteration, whose every iterate marches them in turn; the
results follow the iterations. The title names the models in azimuth
order when they differ, for example `rigid, beam, beam, beam blades`.
The motion report gives the rigid bodies' energy, the hub's angle and,
with beams, the tip of the first beam with its azimuth printed.

## The setup, `rotorcraft.nml`

    &setup
      blade = 'beam', scenario = 'collective', steps = 144, span_points = 5
    /

| name | meaning | default |
|---|---|---|
| `blade` | `'rigid'` or `'beam'` for all four blades when `blades` is omitted | `'rigid'` |
| `blades` | four models, each `'rigid'` or `'beam'`, at azimuths 0, 90, 180 and 270 degrees; `'plate'` and `'solid'` are refused by error stop | every blade takes `blade` |
| `scenario` | `'collective'`, `'longitudinal'` or `'lateral'`, or several of them | `'collective'` |
| `thickness` | `thickness(1)` for rigid blades and beams without strips; with strips, `thickness(k)` initializes strip k of every beam | 2 cm |
| `revolutions` | of the hub | 1 |
| `steps` | per revolution | 360 |
| `order` | of the backward differences of the march, 1 or 2 | 2 |
| `span_points` | of a blade as a beam | 9 |
| `strips` | of each beam's span, the thickness of each a design coordinate; 0, uniform | 0 |
| `each_blade` | with strips: `.false.`, shared thicknesses; `.true.`, each beam's own thicknesses and reported `blade_ks`, with `optimize` bounding each `blade_ks` separately | `.false.` |
| `sensitivities` | with strips: the derivatives along each thickness | `.true.` |
| `optimize` | with strips: the least mass whose stress is bounded | `.false.` |
| `ks_bound` | the bound of the von Mises ratio, its square bounding the smooth maximum | 0.25 |
| `ks_weight` | the weight of the smooth maximum | 100 |
| `bounds` | the least and the greatest thickness | 1 cm, 2 cm |
| `smoothness` | the largest difference of adjacent strips of the same blade | 1 mm |
| `yield` | the yield stress | 350e6 |
| `meshes` | the directory of the bodies' surface meshes; none, no drawing | none |
| `every` | the stride of the instants drawn | 10 |
| `lagrangian` | L of each scenario written as a formula, one equation a line, to `rotor_lagrangian.txt` | `.false.` |

With `strips`, we can vary the beam thicknesses together or separately.
By default, `each_blade = .false.` gives every beam the same coordinates
`h1` through `h<strips>`. The smooth maximum `ks` of the von Mises
function is taken over all the beams, their section corners and the
middles of the long sides, and the span and revolutions. With
`optimize`, this one maximum must satisfy `ks <= ks_bound^2`.

Set `each_blade = .true.` to give each beam its own strips. The
function `strip` names the coordinates `h<deg>_<k>`: for example,
`h90_1` through `h270_8` for eight strips on the beams at 90, 180 and
270 degrees. Each beam's strip k begins at `thickness(k)`. Its own
smooth maximum, `blade_ks`, is reported; with `optimize`, each
`blade_ks <= ks_bound^2` has a component of `mu_ks`. The thickness
`bounds` apply to every strip, and `smoothness` compares adjacent
strips within the same blade. The optimized report prints every
thickness and every component of `mu_ks`, ordered by beam azimuth.

For sensitivities, we still use `ks` over all beams as the one
objective: a sum of smooth maxima is not a functional the library
forms. The expansion and the reverse sweep give its derivatives along
every design coordinate. With `each_blade`, the derivative table names
that coordinate in the first column, then gives the mass derivative,
the two derivatives of `ks`, and the expansion derivative of each
`blade_ks` in beam azimuth order. This lets us read both the change in
the combined stress measure and the change in each beam's measure.

## Setups, and how long they take

The following measurements use `blades` and `each_blade` omitted,
on one image of this workstation, four setups running at once:

- **The rigid rotor over a revolution**, with `blade = 'rigid'`,
  360 steps of bdf(2): about 4 min. Printed: the energy and the hub's
  angle every 30 degrees.
- **The blades as beams**, `blade = 'beam', steps = 144, span_points =
  5`: about 3.5 min. Printed also: the motion of the tip of the blade
  at zero degrees, the blades' mass and the smooth maximum of their
  stress.
- **The derivatives along the thickness of 8 strips**, adding
  `strips = 8, thickness = 8*0.02`: about 3.7 hours, the expansion
  along every strip. Printed for each strip: the mass's derivative, and the
  smooth maximum's by the expansion and by the reverse sweep.
- **The least mass, one scenario**, adding `optimize = .true.`: 2
  design iterations, about 15 min.
- **The least mass, three scenarios**, `scenario = 'collective',
  'longitudinal', 'lateral'`: 2 design iterations, about 46 min.
  Printed for each scenario, at the design found: its motion, mass and
  smooth maximum, the thickness of each strip and the multiplier of its
  bound.

At 144 steps of bdf(2) the blades' stress is the centrifugal stress of
the turning blade, rho w^2 (R^2 - r^2)/2, 0.162 of the yield at the
root at uniform thickness: the smooth maximum is 2.515e-2 at 2 cm and
2.513e-2 at 1 cm, under the 0.0625 of `ks_bound = 0.25`, so that the
least mass is every strip at its least thickness. A thinner root strip
under the outer strips' load raises it: its derivative along the root
strip's thickness is -1.2065, and the derivatives along the eight
strips sum to 1.1e-3. Five points of the span resolve the smooth maximum
to 2.2e-4 of its value at 17; 72 steps overstate it by 58 per cent,
144 and 288 agree to 0.46 per cent. At 36 steps of bdf(1) without the
rotating scissors the upper swashplate lagged the hub by 64 degrees
and the stress was 1.29 times the yield; the costs of the case as TACS
states it in `doc/rotorcraft-optimization.md` were measured on that
march.

To check the model and thickness choices together, consider a rigid
blade at zero degrees and beams at the other three azimuths, each with
its own strips:

    &setup
      blades = 'rigid', 'beam', 'beam', 'beam', each_blade = .true.,
      scenario = 'collective', steps = 12, order = 2, span_points = 5,
      strips = 8, thickness = 8*0.02, sensitivities = .true.
    /

These measurements use 12 steps of bdf(2) from one step of bdf(1).
The march is coarse: the smooth maximum is near 5 where 144 steps
give 2.5e-2. They verify the statement and its derivatives, not the
physics. At this same resolution, `blade = 'beam'` with `blades` and
`each_blade` omitted prints every line as before, with sensitivities
and with `optimize` at `ks_bound = 2.24`; only the times differ.

- **Three beams with separate strips:** 24 design coordinates, with
  the tip at 90 degrees reported. The march, its 24 expansions and
  its adjoint took 4077 s and 633 MB. The beams' mass is
  1.6303345746653470 kg s, three quarters of the four beams'
  2.1737794328871334 at 12 steps. The derivatives of `ks` by expansion
  and reverse sweep agree to 1.4e-13 relative along all 24 strips.
- **The same rotor with `optimize = .true.`:** at `ks_bound = 2.24`,
  two iterations put every strip at its bound of 1 cm. The beams'
  mass is 8.1516728733267363e-1 kg s; each stress bound is inactive
  and every component of `mu_ks` is zero. At `ks_bound = 2.2`, all
  three stress bounds are active. Five iterations take 4 min 26 s,
  giving each `blade_ks` at 4.84 to 1.4e-14 and mass
  8.1559522908988735e-1 kg s. Only the root strips remain above
  1 cm; the program's header records their thicknesses and multipliers.

## Drawing the motion in ParaView

Set `meshes` to the TACS case's surface meshes,

    meshes = '/home/komahan/_git/rotorcraft-optimization/deterministic-optimization/data/mesh'

and `every` to the stride of the instants drawn. The table's rigid
bodies are drawn with their meshes, moved with the bodies to each
instant. A rigid blade is drawn as a box at `thickness(1)`. Each beam
is drawn as a box with the thickness of each strip, its points moved
with their beam sections. Initially these are the setup thicknesses;
at the design found, each beam uses its own strips with `each_blade`,
or the shared strips otherwise. The run writes:

| setup | files |
|---|---|
| one scenario | `rotor_<instant>.vtu`, `rotor.pvd` |
| several scenarios | `rotor_<scenario>_<instant>.vtu`, `rotor_<scenario>.pvd` |
| to optimize | the initial design marched first: `rotor_initial.pvd` and `rotor_optimized.pvd`, or `rotor_<scenario>_initial.pvd` and `rotor_<scenario>_optimized.pvd` |

The least mass of one scenario, drawn at the initial design and at
the design found, every twelfth of the 145 instants, about 19 min:

    &setup
      blade = 'beam', scenario = 'collective', steps = 144, span_points = 5,
      strips = 8, thickness = 8*0.02, optimize = .true., every = 12,
      meshes = '/home/komahan/_git/rotorcraft-optimization/deterministic-optimization/data/mesh'
    /

writes `rotor_initial.pvd` and `rotor_optimized.pvd`, 13 instants each:
the blades 2 cm thick in the first and 1 cm in the second.

The mixed run with `meshes`, `each_blade` and `optimize` also writes
`rotor_initial_*` and `rotor_optimized_*`: each contains the rigid
blade's box and the three beam boxes, with each beam drawn at its own
strip thicknesses.

Each `.vtu` has `body`, the number of the body of each cell, and
with beams `deflection`, the beam's deflection at each point in the
fixed axes, zero on the rigid bodies. In ParaView:

1. File, Open, the `.pvd`; Apply. The time slider runs through the
   instants drawn; Play animates the revolution.
2. Colour by `body` to tell the bodies apart.
3. With beams, Filters, Warp By Vector, the vector `deflection`, a
   scale factor of 100 or more to see the blades bend, their tips by
   some 2 mm; colour by `deflection` to see its magnitude.
4. To compare the initial and the optimized design, open both `.pvd`
   files, split the view, show one in each and link the cameras
   (right-click in a view, Link Camera); the time slider moves both.

Without the window, `pvpython` draws an instant into a picture; on
this workstation it needs the system's Python path,

    PYTHONPATH=/usr/lib/python3/dist-packages pvpython --force-offscreen-rendering draw.py

with `draw.py`:

    from paraview.simple import *
    rotor = PVDReader(FileName='rotor_optimized.pvd')
    warped = WarpByVector(Input=rotor, Vectors=['POINTS', 'deflection'], ScaleFactor=50.0)
    view = GetActiveViewOrCreate('RenderView')
    display = Show(warped, view)
    ColorBy(display, ('CELLS', 'body'))
    view.ViewTime = rotor.TimestepValues[-1]
    view.ResetCamera()
    SaveScreenshot('rotor_optimized.png', view, ImageResolution=[1200, 900])
