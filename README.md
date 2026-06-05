![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/airfoil.png)
![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/cylinder.png)
![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/poisson.png)

# Unstructured Finite Volume Solver

Unstructured Finite Volume Solver for Partial Differential Equations

## Principles/Goals

1. Truly Fortran -std=[f2008,f2018] :sunglasses:
2. Object Oriented Design for separation of geometry, physics and solution
3. Coarrays for distributed memory parallelism
4. *pure*, *elemental*, *do concurrent* for shared memory parallelism

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
h-independent iteration counts (src/class_amg.f90, built on an assembled csr).

distributed-memory CG (goal 3) lives in src/class_distributed_cg.f90: the mesh
is partitioned by recursive coordinate bisection (src/class_graph.f90
`partition_rcb`, organised by src/class_partition.f90), each image owns its
rows, the matvec exchanges only the halo over coarrays and the dot products
reduce with co_sum; one image reduces exactly to the serial CG. a per-image AMG
on the owned block is a restricted-additive-schwarz preconditioner. run it with:

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