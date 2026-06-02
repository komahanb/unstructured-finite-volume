# Towards a General Purpose Finite Volume Solver

Right now the assembler *is* the Laplace equation with a particular box's boundary
conditions baked in. That has to stop. The assembler should only *discretize*; the
physics, the boundaries and the run should all come from outside.

## Principles/Goals

1. Truly Fortran -std=f2018 :sunglasses:
2. Object Oriented Design for separation of **geometry**, **physics**, **boundaries**, **solution** and **run**
3. Nothing problem-specific in the library - the *driver* and a *config file* decide the problem
4. Coarrays for distributed memory parallelism
5. *pure*, *elemental*, *do concurrent* for shared memory parallelism

## What is wrong today (`src/class_assembler.f90`)

- boundary conditions are literals keyed to tag numbers - `phifront=5`, `tag 1 -> 5`, `tag 3 -> 10`, ... only dirichlet, only this one box
- the operator is a unit-diffusivity laplacian - `farea*(qN-qP)/fdelta`, no `gamma`, no convection, no reaction
- `evaluate_source = 0.0` - source is dead code
- one variable per cell - `num_variables = 1`, `num_state_vars = num_cells`
- `write_solution` writes 2d tecplot titled `"FVM-Laplace"` with one variable `"T"`

The mesh already gives us names: `tag_info`, `tag_numbers`, `tag_physical_dimensions`,
`face_tags`, and `face_group`/`create_face_groups` (currently commented out at
`class_mesh.f90:492`). Per-face geometry is all there (`face_areas`, `face_deltas`,
`cell_face_normals`, `cell_face_tangents`, `face_cells`, `lvec`, `cell_volumes`).
So a boundary should be addressed by *name*, never by a magic integer.

We also already have the right abstractions sitting unused (not in `src/OBJECTS`):
`interface_physics`, `interface_assembler`, `class_finite_volume_assembler`,
`class_physics_unsteady_heat`, and the time/nonlinear stack
(`integrator_interface`, `backward_differences_class`, `nonlinear`). The last three are
broken - they `use dynamic_physics_interface, only : dynamics` and that module was never
written. (KB: `dynamics` vs `physics` naming was never resolved - fix with a shim.)

## Design - who owns what

- **mesh** owns geometry and named tags. Add `find_tag_by_name(name)` and switch
  `create_face_groups` on. (nothing else changes here)
- **graph** owns connectivity and dofs. a graph of `vertex`es and `edge`s - the cells are
  vertices, a shared interior face is an edge. it hands the assembler the dof map
  (`(cell,variable) -> global`) and the sparsity (who couples to whom), and it is the thing
  we partition. `module_mesh_utils`'s `reverse_map` / `transpose_connectivities` /
  `sparse_transpose_matmul` are already graph operations - fold them in here. plenty of
  other uses: dof reordering for bandwidth, coloring for the smoothers, halo lists.
  (KB: later hand the graph to parmetis and exchange halos over coarrays - goal 4)
- **equation** owns the physics. an abstract `equation` that hands the assembler, per cell,
  the diffusion tensor `K(x)`, the source, (later) a convecting velocity and reaction, and
  the number of variables. `class_diffusion` is the first concrete one.
- **boundary_condition** owns one named boundary - dirichlet / neumann / robin and its data.
  the assembler holds a table of these, one per face tag.
- **assembler** owns nothing physical. it walks faces and cells and *assembles*:
  matrix-free `J*x`, residual and rhs - asking the equation for fluxes/sources and the bc
  table for the boundary contribution. it also implements `dynamics` so the integrators and
  newton can drive it.
- **config** owns the run. a plain text input file says which mesh, which equation and
  coefficients, the bc on each named boundary, the solver, the time scheme and the output.
- **solvers / integrator / writer** stay as they are - they already work on flat vectors and
  3d multi-field vtu.

a boundary condition, in our usual style:

```fortran
type :: boundary_condition
   type(string) :: name              ! physical group name from the mesh
   integer      :: tag               ! resolved tag number
   integer      :: kind              ! BC_DIRICHLET / BC_NEUMANN / BC_ROBIN
   real(dp)     :: a, b, c           ! robin: a*phi + b*dphi/dn = c ; dirichlet/neumann fall out
 contains
   procedure :: apply                ! contribution to (lhs diagonal, rhs)
end type boundary_condition
```

and the physics behind a flux, deferred so any pde can plug in:

```fortran
type, abstract :: equation
   integer :: num_variables
 contains
   procedure(diffusion_tensor_interface), deferred :: diffusion_tensor   ! K(x) per cell
   procedure(source_interface)          , deferred :: source             ! per cell
end type equation
```

dofs and connectivity as a graph - vertices and edges, in the usual style:

```fortran
type :: vertex
   integer :: number                 ! global dof / cell number
   integer :: part                   ! owning image after partitioning
end type vertex

type :: edge
   integer :: tail, head             ! coupled vertices (cells sharing a face)
   integer :: face                   ! the interior face that made this edge
end type edge

type :: graph
   type(vertex), allocatable :: vertices(:)
   type(edge)  , allocatable :: edges(:)
 contains
   procedure :: dof                  ! (cell, variable) -> global dof
   procedure :: neighbours           ! adjacency of a vertex
   procedure :: partition            ! later: parmetis
end type graph
```

## File Groups

### New
- `src/class_graph.f90` - `vertex`/`edge`/`graph`; owns dof numbering, sparsity and (later) partitioning
- `src/class_boundary_condition.f90` - dirichlet/neumann/robin, addressed by name
- `src/interface_equation.f90` - abstract pde (K tensor, source, #variables, later convection/reaction)
- `src/class_diffusion.f90` - first concrete equation (anisotropic K, per-cell source)
- `src/class_config.f90` - read the run from a text file (reuse `class_file` + `class_string % tokenize`)
- `src/dynamic_physics_interface.f90` - the missing `dynamics` shim so the time/nonlinear stack builds
- `examples/solver/main.f90` + `*.cfg` - one generic driver, problems live in config

### Generalize
- `src/class_assembler.f90` - operator from the equation, boundaries from the table, multi-field, 3d output
- `src/class_mesh.f90` - turn on `create_face_groups`, add `find_tag_by_name`

### Reuse as-is
- `src/class_paraview_writer.f90` (already 3d, multi-field)
- `src/interface_linear_solver.f90` + `class_conjugate_gradient/sor/gauss_seidel/gauss_jacobi`
- `src/class_file.f90`, `src/class_string.f90`, `src/scalar.fpp`

### Repair and put in `src/OBJECTS`
- `src/interface_physics.f90`, `src/integrator_interface.f90`,
  `src/backward_differences_class.f90`, `src/nonlinear.f90`

## Milestones

Each milestone builds clean and keeps box-3 / sphere / cylinder working - they are the
regression. (KB: the first three already kill the hardcoding - that is the important bit.)

0. name lookup + switch on `create_face_groups`; snapshot box-3/sphere/cylinder as golden
1. named, typed boundary conditions - delete every `phifront`/`tag .eq. 1` branch from
   `get_source` and the boundary branch of `get_jacobian_vector_product`
2. equation layer - anisotropic `K` and per-cell source; flux `area*(n^T K n)/delta` and the
   skew correction projected through `K` (generalize the minimum-correction code already in
   `get_skew_source`)
3. graph + multi-field - introduce `class_graph` (cells -> vertices, interior faces -> edges)
   and get the dof map from it - `dof = graph % dof(icell, ivar)`; loop fields in jvp/source/bc;
   solvers untouched. fold `reverse_map`/`transpose_connectivities` onto the graph.
4. config file - generic `examples/solver/main.f90`; box-3/sphere/cylinder become `*.cfg`
5. time + nonlinear - write `dynamic_physics_interface`, make the assembler a `dynamics`, wire
   bdf (transient) and newton (nonlinear, reuse `add_jacobian_fd`)
6. output + cleanup - `write_solution` 3d and n named fields through the paraview writer; drop
   the 2d tecplot/`"FVM-Laplace"`/`"T"` block; tidy drivers, examples, readme
7. distributed - partition the graph with parmetis, one image per part, exchange halo dofs over
   coarrays (goal 4). (KB: the long game - the graph is what makes this tractable)

## Verification

- build the usual way - `make -C src F90=gfortran && make -C src install` (now `-std=f2018`),
  then the drivers (`-fcoarray=single` is already in `Makefile.in`)
- regression every milestone: box-3 (dirichlet), sphere (constant -> exact 5.0 ~1e-13),
  cylinder-coarse (skewed -> converges, no NaN) - unchanged when the config reproduces the
  old problem
- operator sanity: `A = L + U + D` (assembly test) still holds
- new capability: mms error-vs-refinement for anisotropic `K` (box-36 -> box-fine); a
  zero-flux neumann symmetry check; a 2-field reaction-diffusion mms; transient decay vs
  analytic (bdf); a nonlinear mms (newton). everything exports 3d multi-field `.vtu`

## Notes

- big job, done in independent milestones - 0..2 remove the hardcoding, 3..5 add the breadth
- multi-field only touches the dof map in the matrix-free jvp; the linear solvers are flat-vector and stay put
- all of this is on `mesh_changes_3d`. the python `ufvm/` port on `master` is out of scope here
