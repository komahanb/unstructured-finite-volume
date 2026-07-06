![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/airfoil.png)
![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/cylinder.png)
![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/poisson.png)

# Unstructured Finite Volume Solver

Unstructured Finite Volume Solver for Partial Differential Equations

## Principles/Goals

1. Truly Fortran -std=f2018
2. Object Oriented Design for separation of geometry, physics and solution
3. Coarrays for distributed memory parallelism
4. *pure*, *elemental*, *do concurrent* for shared memory parallelism

## Architecture

every concrete class extends one abstract interface file (`interface_*.f90`):

| interface | concrete classes |
|---|---|
| `interface_graph` | `class_graph` (mesh), `class_chain` (rule-generated) |
| `interface_marcher` | ancestor of the three solver families below |
| `interface_assembler` | `class_assembler`, `class_partitioned_assembler` |
| `interface_flux` | `class_diffusion_flux`, `class_advection_flux` |
| `interface_linear_solver` | `class_conjugate_gradient`, `class_gmres`, `class_normal_cg`, `class_sor`, `class_gauss_seidel`, `class_gauss_jacobi` |
| `interface_nonlinear_solver` | `class_newton_solver` |
| `interface_integrator` | `class_bdf` |
| `interface_multigrid` | `class_algebraic_multigrid` |
| `interface_physics` | user-supplied physics |

**marcher contract**: linear solver, newton and the integrator are one
concept — advance a deferred one-step map until termination — and share
one ancestor (`interface_marcher`: tolerance, budget, verbosity). each
family binds its context name through a generic (`solve => march`,
`integrate => march`); the provided march is the residual-minimization
iteration (linear: `converge` around the deferred `iterate`; time:
`integrate` around the deferred `step`).

**system contract**: the assembler is the system. a solver asks it three
things only — `get_residual(r, x)`, the unified jacobian-vector product
`get_jacobian_residual_product(w, v, mode, part)` (mode = forward/
transpose, part = whole/diagonal/triangles), and `inner_product(a, b)`.
discretization words (source, skew) never reach the solver; two analytic
consistency checks (`verify_transpose_consistency`,
`verify_parts_consistency`) audit the product at machine precision.
solvers are stateless w.r.t. the system; `integrator % integrate(mode)`
owns both primal and adjoint trajectories. flux projections
(`normal_diffusivity`, `normal_speed`) live in `interface_flux`.

**graph contract**: the graph owns the retained adjacency, visit
orders and partitioning (`graph % partition()`, `partition_rcb`);
neighbour queries are its compiler-enforced deferred contract. the
DIRECTED graph (`digraph`) adds out/in queries, the topological
`dependency_order` (with a cycle refusal and a pure `is_acyclic`), the
discrete adjoint's structure (`accumulate_adjoint`, direction read from
structure, objective vertex explicit) and its type-bound witness. the
mesh graph is an undirected subclass — the directed protocol does not
compile on it, deliberately; the chain (iterate/step sequences) is
directed. test/graph exercises both standalone.

## Build

```
./build.sh           # builds lib/libufvm.{a,so}; auto-detects gfortran
bash build_parallel.sh   # coarray build into lib_par/ (caf, -fcoarray=lib)
```

## Meshes

meshes are generated on the fly with the gmsh python api (4.15+, msh 4.1) - none
are committed. the loader reads 4.1 only; `meshgen/generate.py` reproduces every
geometry with the physical-group names the configs/tests use:

```
python meshgen/generate.py box-3 test/box-3.msh   # box / sphere / cylinder / square
```

## Run

the run-scripts generate the mesh, build, and solve:

```
bash examples/solver/run.sh box.cfg    # also sphere.cfg / cylinder.cfg / box-transient.cfg
bash test/regression/run.sh            # analytic-oracle suites (3d + 2d)
```

## Solvers

the steady linear solve picks a method by config string: `cg` (matrix-free
conjugate gradient, the default), `sor` / `gs` / `gj` (stationary), or
`pcg-amg` - smoothed-aggregation algebraic multigrid as a CG preconditioner,
h-independent iteration counts (src/class_algebraic_multigrid.f90, built on an
assembled csr).

distributed-memory solves (goal 3) live on the SYSTEM side
(src/class_partitioned_assembler.f90): the mesh is partitioned by
recursive coordinate bisection (`graph % partition_rcb`), each image
computes only its owned rows of the jacobian-vector product and its
owned share of every inner product (reduced with co_sum) — and the
solver is the ordinary `conjugate_gradient`, textually identical in
serial and parallel; one image reduces exactly to the serial CG. a
per-image AMG on the owned block (`block_preconditioner`) is an
additive-schwarz preconditioner that plugs into the usual
`pre_conditioner` slot.
run it with:

```
bash test/parallel/run.sh              # serial + cafrun -np 2,4
```

nonsymmetric operators (e.g. advection) need a non-CG krylov method:
src/class_gmres.f90 (restarted GMRES(m), transpose-free, optional right
preconditioner so the AMG plugs in) and src/class_normal_cg.f90 (cgnr / cgne on
the normal equations - a cheap robust baseline, but they converge on the
*squared* condition number, so GMRES wins decisively once advection dominates;
see test/krylov/run.sh for the comparison).

advection is the first such operator: src/class_advection_flux.f90
(`advection_flux` F = v q, `advection_diffusion_flux` F = v q - K grad q) plugs
into the same flux seam - the assembler picks up the normal advection speed
v.n (central differencing -> a skew-symmetric, non-symmetric operator). solve
it with GMRES on the assembled csr; test/advection/run.sh checks 2nd-order
accuracy against the exact 1d boundary-layer solution.

## Post-processing

each solve writes the field as both a paraview `.vtu` and a gmsh `.msh`:

```
paraview examples/solver/box.vtu
gmsh     examples/solver/box.msh
```