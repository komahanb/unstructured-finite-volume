# The dynamics of a rotorcraft

The rotor of the TACS case of `komahanb/rotorcraft-optimization`
(dissertation, section 5.4): 26 bodies and 35 joints - TACS's 34 and
the rotating scissors, which turns the upper swashplate with the hub -
the hub turning at 109.12 rad/s, three push rods driving the
swashplates, and four blades, each a rigid plate or a beam, chosen at
setup. The bodies begin at the velocities the drivers require. Its motion over
the revolutions, the blades' mass and the smooth maximum of their
stress, the derivatives of these along the thickness of each strip of
the span, and the thickness of least mass whose stress is bounded in
one scenario or in three.

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
results follow the iterations.

## The setup, `rotorcraft.nml`

    &setup
      blade = 'beam', scenario = 'collective', steps = 144, span_points = 5
    /

| name | meaning | default |
|---|---|---|
| `blade` | `'rigid'` or `'beam'` (plate and solid to come) | `'rigid'` |
| `scenario` | `'collective'`, `'longitudinal'` or `'lateral'`, or several of them | `'collective'` |
| `thickness` | of the blades; with strips, of each strip | 2 cm |
| `revolutions` | of the hub | 1 |
| `steps` | per revolution | 360 |
| `order` | of the backward differences of the march, 1 or 2 | 2 |
| `span_points` | of a blade as a beam | 9 |
| `strips` | of the span, the thickness of each a design coordinate; 0, uniform | 0 |
| `sensitivities` | with strips: the derivatives along each thickness | `.true.` |
| `optimize` | with strips: the least mass whose stress is bounded | `.false.` |
| `ks_bound` | the bound of the von Mises ratio, its square bounding the smooth maximum | 0.25 |
| `ks_weight` | the weight of the smooth maximum | 100 |
| `bounds` | the least and the greatest thickness | 1 cm, 2 cm |
| `smoothness` | the largest difference of adjacent strips | 1 mm |
| `yield` | the yield stress | 350e6 |
| `meshes` | the directory of the bodies' surface meshes; none, no drawing | none |
| `every` | the stride of the instants drawn | 10 |

## Setups, and how long they take

On one image of this workstation, four setups running at once:

- **The rigid rotor over a revolution**, with `rotorcraft.nml` as it is,
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
root whatever the thickness: the smooth maximum is 2.515e-2 at 2 cm and
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

## Drawing the motion in ParaView

Set `meshes` to the TACS case's surface meshes,

    meshes = '/home/komahan/_git/rotorcraft-optimization/deterministic-optimization/data/mesh'

and `every` to the stride of the instants drawn. Every rigid body is
drawn with its mesh, moved with the body to each instant; each blade
with the box of its plate, its thickness that of each strip - the
setup's, and for the design found the thicknesses found - and as a
beam each point of the box moved with its section of the beam. The run
writes:

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
