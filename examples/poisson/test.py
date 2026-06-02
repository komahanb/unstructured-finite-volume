#=====================================================================#
# Test loading of mesh and mesh pre-processing
#
# Port of examples/poisson/test.f90 (program test_mesh).
#=====================================================================#

import os
import sys

import numpy as np

# Make the ufvm package importable when run from this directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from ufvm import (
    GmshLoader, Mesh, Assembler,
    ConjugateGradient, GaussJacobi, GaussSeidel, Sor,
)


def get_exact_solution(x):
    # Arguments: x(:,:) coordinates -> returns fexact(npts)
    PI = 4.0 * np.arctan(1.0)
    alpha = 16.0 / (PI ** 4.0)

    npts = x.shape[1]
    fexact = np.zeros(npts, dtype=np.float64)

    # Exact solution as a function of coordinate
    for i in range(1, npts + 1):
        # Evalute f(x) as a summation
        for ii in range(0, 100):
            for jj in range(0, 100):
                mm = 2 * ii + 1
                nn = 2 * jj + 1
                fexact[i - 1] = fexact[i - 1] + \
                    alpha * np.sin(float(mm) * x[0, i - 1] * PI) * np.sin(float(nn) * x[1, i - 1] * PI) \
                    / (float(mm * nn) * float(mm ** 2 + nn ** 2))

    return fexact


def main():
    filename = "square-40.msh"

    # meshing : Create a mesh object
    gmsh = GmshLoader(filename)
    grid = Mesh(gmsh)
    del gmsh

    # assembly : Create an assembler object for the linear system
    FVMAssembler = Assembler(grid)

    # cg_solver
    max_tol = 100.0 * np.finfo(np.float64).eps
    max_it = 100
    print_level = -1

    solver = ConjugateGradient(
        FVAssembler=FVMAssembler, max_tol=max_tol, max_it=max_it, print_level=print_level)

    # Solve using solver method
    fhatc = solver.solve()
    print(' cg solution = ')
    for i in range(1, min(10, len(fhatc)) + 1):
        print(i, fhatc[i - 1])

    # Writes the mesh for tecplot
    FVMAssembler.write_solution("poission-cg-40.dat", fhatc)
    fhatv = FVMAssembler.evaluate_vertex_flux(fhatc)
    fhatf = FVMAssembler.evaluate_face_flux(fhatc)

    # Get exact cell center values of solution
    ff = get_exact_solution(FVMAssembler.grid.face_centers[0:2, :])
    fc = get_exact_solution(FVMAssembler.grid.cell_centers[0:2, :])
    fv = get_exact_solution(FVMAssembler.grid.vertices[0:2, :])

    # Write tecplot output of error
    ec = FVMAssembler.create_vector()
    ec = np.abs(fc - fhatc)
    FVMAssembler.write_solution("poission-ec-40.dat", ec)

    print(" num_faces      ", FVMAssembler.grid.num_faces)
    print(" num_cells      ", FVMAssembler.grid.num_cells)
    print(" num_vertices   ", FVMAssembler.grid.num_vertices)
    print(" num_state_vars ", FVMAssembler.num_state_vars)

    npts = (FVMAssembler.grid.num_cells
            + FVMAssembler.grid.num_vertices
            + FVMAssembler.grid.num_faces)

    rmse = np.sum(np.abs(fv - fhatv) ** 2)
    rmse = np.sqrt(rmse / float(FVMAssembler.grid.num_vertices))

    print(" rmse = ", rmse)
    print(" cell volume = ",
          np.sum(FVMAssembler.grid.cell_volumes) / float(FVMAssembler.grid.num_cells),
          FVMAssembler.grid.cell_volumes.min(),
          FVMAssembler.grid.cell_volumes.max())

    # The Fortran program issues `stop` here; the sor/seidel/jacobi
    # blocks below it are unreachable. They are reproduced as helper
    # functions so the parity with the source is explicit.
    return


def sor_solver(FVMAssembler):
    max_tol = 100.0 * np.finfo(np.float64).eps
    solver = Sor(FVAssembler=FVMAssembler, omega=1.8545, max_tol=max_tol, max_it=100, print_level=-1)
    x = solver.solve()
    print(' sor solution = ')
    for i in range(1, min(10, len(x)) + 1):
        print(i, x[i - 1])
    FVMAssembler.write_solution("poission-sor-40.dat", x)


def seidel_solver(FVMAssembler):
    max_tol = 100.0 * np.finfo(np.float64).eps
    solver = GaussSeidel(FVAssembler=FVMAssembler, max_tol=max_tol, max_it=100, print_level=-1)
    x = solver.solve()
    print(' seidel solution = ')
    for i in range(1, min(10, len(x)) + 1):
        print(i, x[i - 1])
    FVMAssembler.write_solution("poission-seidel-40.dat", x)


def jacobi_solver(FVMAssembler):
    max_tol = 100.0 * np.finfo(np.float64).eps
    solver = GaussJacobi(FVAssembler=FVMAssembler, max_tol=max_tol, max_it=100, print_level=-1)
    x = solver.solve()
    print(' jacobi solution = ')
    for i in range(1, min(10, len(x)) + 1):
        print(i, x[i - 1])
    FVMAssembler.write_solution("poission-jacobi-40.dat", x)


if __name__ == '__main__':
    main()
