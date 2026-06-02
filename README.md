![airfoil](https://github.com/komahanb/unstructured-finite-volume/blob/master/misc/airfoil.png)

# Unstructured Finite Volume Solver

An object-oriented, unstructured finite-volume solver for partial differential
equations. The original library is written in modern Fortran (`src/`); the
`ufvm/` package is a faithful 1:1 Python + NumPy port of its working core, with
the same classes, methods, and argument names.

## Design

The framework cleanly separates **geometry** (mesh), **assembly** (the discrete
operator and right-hand side), and **solution** (linear solvers):

```
GmshLoader ──> Mesh ──> Assembler ──> {ConjugateGradient, Sor, GaussSeidel, GaussJacobi}
```

## Install

```bash
pip install numpy
# run from the repository root, or add it to PYTHONPATH
export PYTHONPATH=/path/to/unstructured-finite-volume
```

## Quick start — solve Poisson on an unstructured mesh

```python
import numpy as np
from ufvm import GmshLoader, Mesh, Assembler, ConjugateGradient

# 1. Load a GMSH mesh and build the topology + geometry (centers,
#    volumes, face areas/normals, interpolation weights, ...)
grid = Mesh(GmshLoader("examples/poisson/square-40.msh"))

# 2. Assemble the finite-volume system over the mesh
fvm = Assembler(grid)

# 3. Solve  A x = b  iteratively
cg = ConjugateGradient(FVAssembler=fvm, max_it=100,
                       max_tol=100.0 * np.finfo(np.float64).eps,
                       print_level=-1)
phi = cg.solve()                       # cell-centered solution

# 4. Post-process and export for Tecplot
fvm.write_solution("poisson.dat", phi)
```

## Swappable linear solvers

All solvers share the same interface (`LinearSolver.solve`), so they are
interchangeable:

```python
from ufvm import ConjugateGradient, Sor, GaussSeidel, GaussJacobi

tol = 100.0 * np.finfo(np.float64).eps
solvers = [
    ConjugateGradient(fvm, max_it=100, max_tol=tol, print_level=-1),
    Sor(fvm, omega=1.8545, max_it=100, max_tol=tol, print_level=-1),
    GaussSeidel(fvm, max_it=100, max_tol=tol, print_level=-1),
    GaussJacobi(fvm, max_it=100, max_tol=tol, print_level=-1),
]
for s in solvers:
    phi = s.solve()                    # all converge to the same field
```

## Matrix-free assembly & operator splitting

The Jacobian is matrix-free — `Assembler` exposes `J·x` products and filtered
assembly (full, diagonal, lower, upper) used by the iterative solvers:

```python
Ax = np.zeros(fvm.num_state_vars)
fvm.get_jacobian_vector_product(Ax, phi)        # action of J on a vector

A = fvm.get_jacobian()                          # dense J  (= L + D + U)
L = fvm.get_jacobian(filter=fvm.LOWER_TRIANGLE)
U = fvm.get_jacobian(filter=fvm.UPPER_TRIANGLE)
D = fvm.get_jacobian(filter=fvm.DIAGONAL)
assert np.allclose(A, L + U + D)
```

## Flux interpolation & mesh queries

```python
phi_v = fvm.evaluate_vertex_flux(phi)           # cell -> vertex values
phi_f = fvm.evaluate_face_flux(phi)             # cell -> face   values

grid.to_string()                                # print topology + geometry

# Geometry/topology arrays live on the mesh:
grid.num_cells, grid.num_vertices, grid.num_faces
grid.cell_centers      # (3, num_cells)
grid.cell_volumes      # (num_cells,)
grid.face_centers      # (3, num_faces)
grid.face_areas        # (num_faces,)
grid.cell_face_normals # (3, max_cell_faces, num_cells)
```

## Package layout

| Component        | Classes / module |
|------------------|------------------|
| Mesh & IO        | `GmshLoader`, `MeshLoader`, `Mesh`, `module_mesh_utils` |
| Assembly         | `Assembler` |
| Linear solvers   | `LinearSolver`, `ConjugateGradient`, `Sor`, `GaussSeidel`, `GaussJacobi` |
| Utilities        | `String`, `File`, `Set`, `List` |

## Examples & tests

Runnable scripts mirror the Fortran drivers:

```bash
python examples/poisson/test.py     # Poisson on a square, 4 solvers, RMSE vs exact
python test/assembly/test.py        # Jacobian filters and A = L + U + D check
python test/mesh/test.py            # GMSH loading across several meshes
python test/solve/test.py           # CG solve on an airfoil mesh
python test/unsteady/test.py        # CG solve on a rectangle mesh
```

The Python port is validated against the original Fortran (built with
`gfortran`): solutions and RMSE match to machine precision (~1e-16).

## Original Fortran library

The reference implementation lives in `src/` and builds to `libufvm.a`:

```bash
cd src && make            # requires gfortran + LAPACK
```

Design goals of the Fortran core: modern `-std=f2008/f2018`, object-oriented
separation of geometry/physics/solution, coarrays for distributed-memory
parallelism, and `pure`/`elemental`/`do concurrent` for shared-memory parallelism.
