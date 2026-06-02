#=====================================================================#
# Test loading of mesh and mesh pre-processing
#
# Port of test/assembly/test.f90 (program test_mesh).
#=====================================================================#

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from ufvm import GmshLoader, Mesh, Assembler


def print_matrix(matrix):
    m = matrix.shape[0]
    n = matrix.shape[1]
    for i in range(1, min(10, m) + 1):
        print("".join("%15.1f" % matrix[i - 1, j - 1] for j in range(1, min(10, n) + 1)))


def main():
    filename = "../triangle.msh"

    # meshing : Create a mesh object
    gmsh = GmshLoader(filename)
    grid = Mesh(gmsh)
    del gmsh

    # assembly
    FVMAssembler = Assembler(grid)

    # Assemble parts of the jacobian and test
    print('getting upper triangle')
    U = FVMAssembler.get_jacobian(filter=FVMAssembler.UPPER_TRIANGLE)
    print_matrix(U)

    print('getting lower triangle')
    L = FVMAssembler.get_jacobian(filter=FVMAssembler.LOWER_TRIANGLE)
    print_matrix(L)

    print('getting diagonal matrix')
    D = FVMAssembler.get_jacobian(filter=FVMAssembler.DIAGONAL)
    print_matrix(D)

    print('getting full jacobian')
    A = FVMAssembler.get_jacobian()
    print_matrix(A)

    # Check consistency of matrix assembly
    if np.abs(A - L - U - D).max() > np.finfo(np.float64).tiny:
        raise SystemExit("error in assembly")
    else:
        print('passed assembly test')

    print('performing symmetry test')
    AT = FVMAssembler.get_transpose_jacobian()
    print_matrix(AT)
    print("asymmetry", (np.abs(A - AT).max() > np.finfo(np.float64).tiny))


if __name__ == '__main__':
    main()
