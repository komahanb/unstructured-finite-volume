![alt text](https://github.com/komahanb/unstructured-finite-volume/blob/master/misc/airfoil.png)

# Unstructured Finite Volume Solver

Unstructured Finite Volume Solver for Partial Differential Equations

## Principles/Goals

1. Truly Fortran -std=[f2008,f2018] :sunglasses:
2. Object Oriented Design for separation of geometry, physics and solution
3. Coarrays for distributed memory parallelism
4. *pure*, *elemental*, *do concurrent* for shared memory parallelism

# File Groups

##  Assembler

- unstructured-finite-volume/src/class_assembler.f90 (KB: What is the difference between assembler and finite volume assembler?)

## Physics

- unstructured-finite-volume/src/interface_physics.f90
- unstructured-finite-volume/src/class_physics_unsteady_heat.f90

## Time integrators

- unstructured-finite-volume/src/integrator_interface.f90
- unstructured-finite-volume/src/backward_differences_class.f90

## Nonlinear System Solvers

- unstructured-finite-volume/src/interface_nonlinear_solver.f90
- unstructured-finite-volume/src/class_nonlinear_solver.f90

## Linear System Solvers

- unstructured-finite-volume/src/interface_linear_solver.f90
- unstructured-finite-volume/src/class_sor.f90
- unstructured-finite-volume/src/class_conjugate_gradient.f90
- unstructured-finite-volume/src/class_gauss_jacobi.f90
- unstructured-finite-volume/src/class_gauss_seidel.f90

## Mesh, Readers and Writers

- unstructured-finite-volume/src/interface_mesh_loader.f90
- unstructured-finite-volume/src/mesh/class_gmsh_loader.f90
- unstructured-finite-volume/src/module_mesh_utils.f90
- unstructured-finite-volume/src/class_mesh.f90

## Utilities

- unstructured-finite-volume/src/class_string.f90
- unstructured-finite-volume/src/class_file.f90
- unstructured-finite-volume/src/class_set.f90
- unstructured-finite-volume/src/class_list.f90

## Examples

- unstructured-finite-volume/examples/unsteady-heat/test.f90
- unstructured-finite-volume/examples/poisson/test.f90

## Unit Tests

### Test Assembly of Linear System
- unstructured-finite-volume/test_sparse.f90
- unstructured-finite-volume/test/assembly/test.f90
- unstructured-finite-volume/test/unsteady/test.f90

### Unit Test Mesh IO
- unstructured-finite-volume/test/mesh/test.f90
- unstructured-finite-volume/test/mesh/class_test_mesh_loader.f90

### Unit Test Solve
- unstructured-finite-volume/test/solve/test.f90
