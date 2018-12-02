![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/doc/poisson.png)

# Unstructured Finite Volume Solver

Unstructured Finite Volume Solver for Partial Differential Equations

## Principles/Goals

1. Truly Fortran -std=[f2008,f2018] :sunglasses:
2. Object Oriented Design for separation of geometry, physics and solution
3. Coarrays for distributed memory parallelism
4. *pure*, *elemental*, *do concurrent* for shared memory parallelism