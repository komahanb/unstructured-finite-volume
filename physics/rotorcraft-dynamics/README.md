# The dynamics of a rotorcraft

The rotor of the TACS case of `komahanb/rotorcraft-optimization`
(dissertation, section 5.4): 26 bodies and 34 joints, the hub turning
at 109.12 rad/s, three push rods driving the swashplates, and four
blades, each a rigid plate or a beam, chosen at setup. Its motion over
the revolutions, the blades' mass and the smooth maximum of their
stress, the derivatives of these along the thickness of each strip of
the span, and the thickness of least mass whose stress is bounded in
one scenario or in three.

| program | what it does |
|---|---|
| `rotorcraft_dynamics` | the rotor; its setup read from `rotorcraft.nml` |
| `triple_pendulum` | the triple pendulum of the dissertation, section 5.1, marched by bdf(2) or by a DIRK tableau |
| `triple_pendulum_sensitivity` | the pendulum's dJ/ds by the adjoint and by the expansion, by either march |
| `blade_on_hub` | one blade on a driven hub, as a beam or as a rigid plate |
| `flexible_mechanism` | the template of a problem of a flexible mechanism - its assembly, its design, its scenarios, the march and the audit - compiled with the case `hart2_case`: the HART II rotor from rest on the ground, spun up by its controls, its eight strips' least mass under a bound of the stress |

`rotor_bodies.txt` is the table of the rigid bodies, read by
`rotorcraft_dynamics`; the header of each program states what it
computes and what it was measured to give.

## Building and running

From this directory,

    make                              # the library, if its sources changed, and the four programs
    ./rotorcraft_dynamics             # reads rotorcraft.nml in the directory it runs in
    ./triple_pendulum 201             # the instant count (default 201)
    ./triple_pendulum 201 'dirk(3)'   # and the march: bdf(2) (default), dirk(3) or dirk(4,5)
    ./triple_pendulum_sensitivity 1 1 1   # the two design values and the order of the expansion
    ./triple_pendulum_sensitivity 1 1 1 'dirk(4,5)'   # and the march
    ./blade_on_hub beam 36 5 70e9     # beam [steps] [span points] [E], or: rigid [steps]

`./run.sh` builds the library from the repository's root and runs the
four programs with their defaults. On several images, build with
`make COARRAY=lib` and run with `mpirun -np 2 ./rotorcraft_dynamics`;
one image needs `OMPI_MCA_osc=pt2pt mpirun -np 1 ...`. The rotor's march
gives the same digits on 1 and 2 images, the pendulum's on 1, 2 and 4.

To run a setup of your own, copy `rotor_bodies.txt` and a
`rotorcraft.nml` into a directory and run the program there: each run
writes its drawing into the directory it runs in.

## The setup, `rotorcraft.nml`

    &setup
      blade = 'beam', scenario = 'collective', steps = 36, span_points = 5
    /

| name | meaning | default |
|---|---|---|
| `blade` | `'rigid'` or `'beam'` (plate and solid to come) | `'rigid'` |
| `scenario` | `'collective'`, `'longitudinal'` or `'lateral'`, or several of them | `'collective'` |
| `thickness` | of the blades; with strips, of each strip | 2 cm |
| `revolutions` | of the hub | 1 |
| `steps` | per revolution | 360 |
| `order` | of the backward differences of the march, 1 or 2 | 1 |
| `span_points` | of a blade as a beam | 9 |
| `precone` | of the blades, degrees, each span inclined upward about the chord at the root; HART II's is 2.5 | 0 |
| `twist`, `twist_zero` | of the blades as beams: degrees over the blade's radius, linear, nose up positive, zero at the radius `twist_zero`; HART II's -8 and 1.5 m | 0, 1.5 |
| `section_kind` | of a blade as a beam: `'solid'`, the plate of the setup's thickness; `'given'`, by the properties below, constant along the span | `'solid'` |
| `section_mass` | given: the mass per length, kg/m | 0.95 |
| `section_moments` | given: the moments of inertia per length about the span, the chord and the thickness, kg m | 7.1828e-4, 1.7e-5, 7.0128e-4 |
| `section_axial` | given: EA, N | 1.17e7 |
| `section_flap` | given: EI_f, in bending along the thickness, N m^2 | 250 |
| `section_lag` | given: EI_c, in bending along the chord, N m^2 | 5200 |
| `section_torsion` | given: GJ, N m^2 | 160 |
| `strips` | of the span, the thickness of each a design coordinate; 0, uniform | 0 |
| `sensitivities` | with strips: the derivatives along each thickness | `.true.` |
| `optimize` | with strips: the least mass whose stress is bounded | `.false.` |
| `ks_bound` | the bound of the von Mises ratio, its square bounding the smooth maximum | 0.25 |
| `ks_weight` | the weight of the smooth maximum | 100 |
| `objective` | the functional of the Lagrangian: to optimize, the one minimized; with strips and sensitivities, the one the reverse sweep differentiates too | `'mass'` to optimize, `'stress'` otherwise |
| `constraints` | to optimize: the functionals bounded in every scenario | `'stress'` |
| `limits` | their limits, each functional at most its limit | the stress's `ks_bound`^2 |
| `bounds` | the least and the greatest thickness | 1 cm, 2 cm |
| `smoothness` | the largest difference of adjacent strips | 1 mm |
| `yield` | the yield stress | 350e6 |
| `meshes` | the directory of the bodies' surface meshes; none, no drawing | none |
| `blade_meshes` | the directory of the blades' surface meshes, `blade0.bdf` .. `blade270.bdf` (the geometry's `--preset hart2`), drawn in place of the plates' boxes | none |
| `every` | the stride of the instants drawn | 10 |

## The objectives

Functionals of the motion of the blades as beams, printed after each
march, each a reduction (`operation_field`) of a field over the
manifold - the revolutions and the span - and over the field's
components: the blades', the links', the rods' or the hub's. With
strips and `sensitivities`, the derivatives of every one along the
thickness of each strip by the expansion of the motion, and those of
the one named `objective` by the reverse sweep too. A constraint is the
same reduction of the field less its limit, at most zero, with its
multiplier in the Lagrangian; an optimization demonstration differs
from another only by its `objective`, `constraints` and `limits`.

| name | the reduction | its field |
|---|---|---|
| `mass` | the integral, kg s | the blades' mass per length |
| `stress` | the smooth maximum of weight `ks_weight`; the solid section alone states it | the squared von Mises ratio at the section's corners and long sides |
| `strain_energy` | the mean, J/m | the blades' strain energy per length |
| `flap`, `lag`, `twist` | the mean, m^2 and rad^2 | the squares of the flap and lag displacements and of the elastic twist |
| `pitch_link` | the mean, N^2 | the squared forces of the four pitch links on the upper swashplate |
| `push_rod` | the mean, N^2 | the squared forces of the three push rods' drivers |
| `hub_force` | the mean, N^2 | the squared force of the hub's driver |
| `hub_moment` | the mean | the squared reactions of the hub driver's conditions on its tilt about x and y |

Not yet: the driver's torque and power, in vacuum mostly the
dissipation of the march; the placement of the frequencies on the
fan plot, which needs the eigenvalues; anything of the air.

## Setups, and how long they take

On one image of this workstation:

- **The rigid rotor over a revolution**, with `rotorcraft.nml` as it is:
  about 2 min. Printed: the energy and the hub's angle every 30
  degrees.
- **The blades as beams**, `blade = 'beam', steps = 36, span_points =
  5`: about 1 min. Printed also: the motion of the tip of the blade at
  zero degrees, the blades' mass and the smooth maximum of their
  stress.
- **The blades as beams of a given section**, adding `section_kind =
  'given'`: about 1 min. The defaults are the properties of TACS's
  uniform HART-II case (its `examples/rotating_beam`). Printed as for
  the beams, without the smooth maximum: a section given by its
  properties states no stress, and the strips, which design the
  plate's thickness, are refused with it.
- **The derivatives along the thickness of 8 strips**, adding
  `strips = 8, thickness = 8*0.02`: about half an hour, the expansion
  along every strip. Printed for each strip: the mass's derivative, and the
  smooth maximum's by the expansion and by the reverse sweep.
- **The least mass, one scenario**, adding `optimize = .true.,
  ks_bound = 1.25`: 6 design iterations, about 11 min.
- **The least mass, three scenarios**, `scenario = 'collective',
  'longitudinal', 'lateral'` with `ks_bound = 1.5`: 6 design
  iterations, about 33 min. Printed for each scenario, at the design
  found: its motion, mass and smooth maximum, the thickness of each
  strip and the multiplier of its bound.

At 36 steps of bdf(1) the blades' stress is 1.29 times the yield at
2 cm and no thickness meets `ks_bound = 0.25`; the bounds 1.25 and 1.5
above stand in for it until the march of the reduced case is decided.
That stress, as the 0.29 at 72 steps of bdf(2), is the start's: the
blades begin unstretched at the full rate, and the axial wave this
sets off gives the largest stress, which does not converge in the
step; once it has passed, the ratio is the centrifugal stress, some
0.16 (`doc/rotorcraft-optimization.md`, "The published optimization").
Five points of the span resolve the stress at the root to some 15 per
cent; the case as TACS states it - 48 strips, 49 points of the span,
360 steps, three scenarios - would take days (see
`doc/rotorcraft-optimization.md`).

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
the design found, every third of the 37 instants, about 12 min:

    &setup
      blade = 'beam', scenario = 'collective', steps = 36, span_points = 5,
      strips = 8, thickness = 8*0.02, optimize = .true., ks_bound = 1.25, every = 3,
      meshes = '/home/komahan/_git/rotorcraft-optimization/deterministic-optimization/data/mesh'
    /

writes `rotor_initial.pvd` and `rotor_optimized.pvd`, 13 instants each:
the blades 2 cm thick in the first, and in the second 1.02 cm at the
root strip and 1 cm outboard.

Each `.vtu` has `body`, the number of the body of each cell, and
with beams `deflection`, the beam's deflection at each point in the
fixed axes, zero on the rigid bodies. In ParaView:

1. File, Open, the `.pvd`; Apply. The time slider runs through the
   instants drawn; Play animates the revolution.
2. Colour by `body` to tell the bodies apart.
3. With beams, Filters, Warp By Vector, the vector `deflection`, a
   scale factor of 10 to 100 to see the blades bend - at 36 steps of
   bdf(1) mostly in lag, the centripetal acceleration a step late;
   colour by `deflection` to see its magnitude.
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
