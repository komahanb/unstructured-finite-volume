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
./build.sh        # builds lib/libufvm.{a,so}; auto-detects gfortran
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

## Post-processing

each solve writes the field as both a paraview `.vtu` and a gmsh `.msh`:

```
paraview examples/solver/box.vtu
gmsh     examples/solver/box.msh
```