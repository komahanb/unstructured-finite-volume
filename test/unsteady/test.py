#=====================================================================#
# Test loading of mesh and mesh pre-processing
#
# Port of test/unsteady/test.f90 (program test_mesh).
#=====================================================================#

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from ufvm import GmshLoader, Mesh, Assembler, ConjugateGradient


def main():
    filename = "../rectangle.msh"

    max_tol = 100.0 * np.finfo(np.float64).eps
    max_it = 100
    print_level = 1

    # Create a mesh object
    grid = Mesh(GmshLoader(filename))

    # Assembler Object coordinating geometry and physics
    FVMAssembler = Assembler(grid)

    # Linear Solution
    linear = ConjugateGradient(
        FVAssembler=FVMAssembler, max_tol=max_tol, max_it=max_it, print_level=print_level)

    # Solve using solver method
    x = linear.solve()
    print('cg solution = ')
    for i in range(1, min(10, len(x)) + 1):
        print(i, x[i - 1])

    # Writes the mesh for tecplot
    FVMAssembler.write_solution("mesh-cg.dat", x)


if __name__ == '__main__':
    main()
