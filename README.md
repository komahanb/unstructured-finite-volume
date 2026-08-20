![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/airfoil.png)
![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/cylinder.png)
![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/poisson.png)

# Unstructured Finite Volume Solver

Unstructured finite volume solver for partial differential equations,
built as a tower of graph mathematics.

## Principles

1. Fortran, latest released standard the compiler speaks (f2023 on gcc 13+)
2. A small set of prime types generates every higher type by composition
3. `pure` / `elemental` where the mathematics is pure
4. Coarrays for distributed memory (serial build is `-fcoarray=single`)

## Architecture

`doc/ARCHITECTURE.md` is the map, generated from the code: eight prime
types (graph, token, relation, field, operation, transform, map, view)
and the composition arithmetic that produces everything else — stencils,
time steps, linearizations, minimizers, the marcher, the adjoint, and
the mesh itself. `AGENTS.md` is the constitution behind it;
`doc/coding-standards.md` holds the writing rules.

The mesh is derived, not stored: a mesh file supplies member sets, the
cell-to-vertex relation, vertex coordinates, and boundary tags;
everything else (faces, adjacencies, centres, areas, normals, volumes,
weights) is computed by relation algebra and geometry and seated on the
one directed graph type.

## Build

```
./build.sh               # builds lib/libufvm.{a,so}; auto-detects gfortran
bash build_parallel.sh   # coarray build into lib_par/ (caf, -fcoarray=lib)
```

## Test

Each directory under `test/` with a `run.sh` is a suite; the tower
suites (`*-tower/`) carry per-level import gates. Run one suite:

```
(cd test/graph-ordinary && bash run.sh)
```

## Meshes

Sample meshes are generated with the gmsh python api (4.15+, msh 4.1):

```
python meshgen/generate.py box-3 test/box-3.msh   # box / sphere / cylinder / square
```
